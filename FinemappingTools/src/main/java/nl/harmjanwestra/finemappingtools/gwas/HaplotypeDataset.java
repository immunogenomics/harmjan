package nl.harmjanwestra.finemappingtools.gwas;

import cern.colt.matrix.tbit.BitVector;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import nl.harmjanwestra.finemappingtools.gwas.tasks.LRTestHaploTask;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.console.ProgressBar;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Triple;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;

/**
 * Created by hwestra on 11/2/16.
 */
public class HaplotypeDataset {

	private final ArrayList<VCFVariant> variants;
	private final BitVector[][] haplotypePairs;
	private final ArrayList<BitVector> availableHaplotypesList;
	private final DiseaseStatus[][] finalDiseaseStatus;
	private final DoubleMatrix2D finalCovariates;
	private final ArrayList<Integer> originalIndsIncluded;
	private double[] haplotypeFrequenciesCases;
	private double[] haplotypeFrequenciesControls;
	private BitVector referenceHaplotype;
	private int nrHaplotypes;
	private DenseDoubleMatrix2D haplotypeMatrix;
	private HashMap<BitVector, Integer> hapToInt;
	private int minorHapId = -1;
	private double minAlleleFreq;


	public HaplotypeDataset(BitVector[][] haplotypePairs,
							ArrayList<BitVector> availableHaplotypesList,
							DiseaseStatus[][] finalDiseaseStatus,
							DoubleMatrix2D finalCovariates,
							ArrayList<Integer> originalIndIdsIncluded,
							ArrayList<VCFVariant> variants) {
		this.finalDiseaseStatus = finalDiseaseStatus;
		this.finalCovariates = finalCovariates;
		this.haplotypePairs = haplotypePairs;
		this.availableHaplotypesList = availableHaplotypesList;
		this.variants = variants;
		this.originalIndsIncluded = originalIndIdsIncluded;
		index();
		calcHapFreq();
		populateMatrix();
	}


	private void index() {
		hapToInt = new HashMap<>();
		for (int i = 0; i < getNrHaplotypes(); i++) {
			hapToInt.put(availableHaplotypesList.get(i), i);
		}
	}

	private void populateMatrix() {
		haplotypeMatrix = new DenseDoubleMatrix2D(getNrIndividuals(), getNrHaplotypes());
		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			BitVector[] haps = haplotypePairs[i];

			Integer hapId1 = hapToInt.get(haps[0]);
			Integer hapId2 = hapToInt.get(haps[1]);

			if (hapId1.equals(hapId2)) {
				haplotypeMatrix.set(i, hapId1, 2);
			} else {
				haplotypeMatrix.set(i, hapId1, 1);
				haplotypeMatrix.set(i, hapId2, 1);
			}
		}
	}

	public ArrayList<VCFVariant> getVariants() {
		return variants;
	}

	public BitVector[][] getHaplotypePairs() {
		return haplotypePairs;
	}

	public ArrayList<BitVector> getAvailableHaplotypesList() {
		return availableHaplotypesList;
	}

	public DiseaseStatus[][] getFinalDiseaseStatus() {
		return finalDiseaseStatus;
	}

	public DoubleMatrix2D getFinalCovariates() {
		return finalCovariates;
	}

	public double[] getHaplotypeFrequenciesCases() {
		return haplotypeFrequenciesCases;
	}

	public double[] getHaplotypeFrequenciesControls() {
		return haplotypeFrequenciesControls;
	}

	public BitVector getReferenceHaplotype() {
		return referenceHaplotype;
	}

	public static HaplotypeDataset create(ArrayList<VCFVariant> variants,
										  DiseaseStatus[][] finalDiseaseStatus,
										  DoubleMatrix2D finalCovariates,
										  double genotypeprobthreshold,
										  boolean[] includeTheseIndividuals,
										  ExecutorService exService) {
// get a list of available haplotypes
		HashSet<BitVector> availableHaplotypeHash = new HashSet<BitVector>();
		ArrayList<BitVector> availableHaplotypesList = new ArrayList<>();

		CompletionService<Triple<BitVector[], Integer, Boolean>> jobHandler = new ExecutorCompletionService<Triple<BitVector[], Integer, Boolean>>(exService);

		int submit = 0;
		int[] iToIndex = new int[finalDiseaseStatus.length];
		ArrayList<Integer> originalIndIdsIncluded = new ArrayList<>();
		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			if (includeTheseIndividuals == null || includeTheseIndividuals[i]) {
				jobHandler.submit(new LRTestHaploTask(i, variants, genotypeprobthreshold));
				iToIndex[i] = submit;
				originalIndIdsIncluded.add(i);
				submit++;
			} else {
				iToIndex[i] = -1;
			}
		}

		int returned = 0;
		BitVector[][] haplotypePairs = new BitVector[submit][];
		ProgressBar pb = new ProgressBar(submit);
		while (returned < submit) {

			try {
				Triple<BitVector[], Integer, Boolean> fut = jobHandler.take().get();
				BitVector[] haps = fut.getLeft();

				int i = fut.getMiddle();
				int index = iToIndex[i];

				if (!fut.getRight()) {

					haplotypePairs[index] = haps;

					if (!availableHaplotypeHash.contains(haps[0])) {
						availableHaplotypesList.add(haps[0]);
						availableHaplotypeHash.add(haps[0]);
					}
					if (!availableHaplotypeHash.contains(haps[1])) {
						availableHaplotypesList.add(haps[1]);
						availableHaplotypeHash.add(haps[1]);
					}

				} else {
					haplotypePairs[index] = null;
				}

				pb.iterate();
				returned++;
				if (availableHaplotypeHash.size() > 1E5) {
					System.out.println("Can't work with this number of haplotypes..");
					System.exit(-1);
				}

//				System.out.print("Individual: " + returned + "\t" + availableHaplotypeHash.size() + " available haplotypes so far.\r");
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}
		pb.close();

		if (submit < finalDiseaseStatus.length) {

			// prune the disease status and covariates
			DiseaseStatus[][] statuses = new DiseaseStatus[submit][finalDiseaseStatus[0].length];
			DoubleMatrix2D covariates = new DenseDoubleMatrix2D(submit, finalCovariates.columns());
			for (int i = 0; i < finalDiseaseStatus.length; i++) {
				int index = iToIndex[i];
				if (index != -1) {
					for (int j = 0; j < finalDiseaseStatus[0].length; j++) {
						statuses[index][j] = finalDiseaseStatus[i][j];
					}

					for (int j = 0; j < finalCovariates.columns(); j++) {
						covariates.setQuick(index, j, finalCovariates.getQuick(i, j));
					}
				}
			}

			return new HaplotypeDataset(haplotypePairs,
					availableHaplotypesList,
					statuses,
					covariates,
					originalIndIdsIncluded,
					variants);
		} else {
			return new HaplotypeDataset(haplotypePairs,
					availableHaplotypesList,
					finalDiseaseStatus,
					finalCovariates,
					originalIndIdsIncluded,
					variants);
		}


	}


	public HaplotypeDataset prune(double frequencyThreshold) {
		System.out.println("Variant allele frequencies: ");
		for (int i = 0; i < variants.size(); i++) {
			VCFVariant var = variants.get(i);
			System.out.println(var.getId() + "\t" + var.getAlleleFrequencies()[0] + "\t" + var.getAlleleFrequencies()[1]);
		}


//		get a list of haplotypes below the threshold
		HashSet<Integer> keepTheseIndividuals = new HashSet<>();
		ArrayList<BitVector> haplotypesAboveThreshold = new ArrayList<>();
		for (int b = 0; b < availableHaplotypesList.size(); b++) {
			if (haplotypeFrequenciesControls[b] > frequencyThreshold) {
				haplotypesAboveThreshold.add(availableHaplotypesList.get(b));
			}
		}

		HashSet<BitVector> allowedHaplotypes = new HashSet<>();
		allowedHaplotypes.addAll(haplotypesAboveThreshold);
		ArrayList<Integer> origIds = new ArrayList<Integer>();
		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			BitVector[] haps = haplotypePairs[i];
			if (haps != null) {
				if (allowedHaplotypes.contains(haps[0]) && allowedHaplotypes.contains(haps[1])) {
					keepTheseIndividuals.add(i);
					origIds.add(originalIndsIncluded.get(i));
				}
			}
		}

		System.out.println(keepTheseIndividuals.size() + " individuals will remain after pruning.");
		if (keepTheseIndividuals.isEmpty()) {
			return null;
		}

		// prune haplotypes, covariates and disease status
		DiseaseStatus[][] remainingdiseaseStatus = new DiseaseStatus[keepTheseIndividuals.size()][1];
		DoubleMatrix2D remainingCovariates = new DenseDoubleMatrix2D(keepTheseIndividuals.size(), finalCovariates.columns());

		BitVector[][] remainingHaplotypes = new BitVector[keepTheseIndividuals.size()][2];

		int ictr = 0;
		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			BitVector[] haps = haplotypePairs[i];
			if (keepTheseIndividuals.contains(i)) {
				remainingHaplotypes[ictr] = haps;

				remainingdiseaseStatus[ictr][0] = finalDiseaseStatus[i][0];
				for (int c = 0; c < finalCovariates.columns(); c++) {
					remainingCovariates.setQuick(ictr, c, finalCovariates.getQuick(i, c));
				}
				ictr++;
			}
		}
		return new HaplotypeDataset(remainingHaplotypes, haplotypesAboveThreshold, remainingdiseaseStatus, remainingCovariates, origIds, variants);
	}

	public void calcHapFreq() {
		// index haplotypes

		HashMap<BitVector, Integer> index = new HashMap<BitVector, Integer>();
		for (int b = 0; b < availableHaplotypesList.size(); b++) {
			index.put(availableHaplotypesList.get(b), b);
		}

		// determine their frequencies
		this.haplotypeFrequenciesCases = new double[availableHaplotypesList.size()];
		this.haplotypeFrequenciesControls = new double[availableHaplotypesList.size()];

		int totalChromosomesCases = 0;
		int totalChromosomesControls = 0;

		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			BitVector[] haps = haplotypePairs[i];
			DiseaseStatus diseaseStatus = finalDiseaseStatus[i][0];
			if (haps != null) {

				Integer index1 = index.get(haps[0]);
				Integer index2 = index.get(haps[1]);

				if (diseaseStatus.equals(DiseaseStatus.CASE)) {
					haplotypeFrequenciesCases[index1]++;
					haplotypeFrequenciesCases[index2]++;
					totalChromosomesCases += 2;
				} else {
					haplotypeFrequenciesControls[index1]++;
					haplotypeFrequenciesControls[index2]++;
					totalChromosomesControls += 2;
				}
			}
		}

		BitVector reference = null;
		double maxfreq = -1;
		System.out.println("Haplotype frequencies:");
		System.out.println("Haplotype\tCases\tControls");
		minAlleleFreq = 1;
		for (int i = 0; i < availableHaplotypesList.size(); i++) {
			haplotypeFrequenciesCases[i] /= totalChromosomesCases;
			haplotypeFrequenciesControls[i] /= totalChromosomesControls;
			if (haplotypeFrequenciesControls[i] > maxfreq) {
				reference = availableHaplotypesList.get(i);
				maxfreq = haplotypeFrequenciesControls[i];
			}
			if (haplotypeFrequenciesControls[i] < minAlleleFreq) {
				minorHapId = i;
				minAlleleFreq = haplotypeFrequenciesControls[i];
			}
			System.out.println(getHaplotypeDesc(i) + "\t" + haplotypeFrequenciesCases[i] + "\t" + haplotypeFrequenciesControls[i]);
		}
		referenceHaplotype = reference;
	}

	public String getHaplotypeDesc(int i) {

		BitVector bitVector = availableHaplotypesList.get(i);
		return getHaplotypeDesc(bitVector);
	}

	public SNPFeature asFeature(int hap, Integer referenceHapId) {
		Chromosome chr = variants.get(0).getChrObj();
		int start = variants.get(0).getPos();
		int stop = variants.get(variants.size() - 1).getPos();
		SNPFeature snp = new SNPFeature(chr, start, stop);
		snp.setAFCases(haplotypeFrequenciesCases[hap]);
		snp.setAFControls(haplotypeFrequenciesControls[hap]);
		snp.setName(getVariantsAsString());
		if (referenceHapId == null) {
			snp.setAlleles(new String[]{"All", getHaplotypeDesc(hap)});
			if (!availableHaplotypesList.get(hap).equals(referenceHaplotype)) {
				snp.setMinorAllele(getHaplotypeDesc(hap));
			} else {
				snp.setMinorAllele(getHaplotypeDesc(hapToInt.get(referenceHaplotype)));
			}
		} else {
			String[] alleles = new String[availableHaplotypesList.size()];
			for (int h = 0; h < availableHaplotypesList.size(); h++) {
				alleles[h] = getHaplotypeDesc(h);
			}
			snp.setAlleles(alleles);
			snp.setMinorAllele(getHaplotypeDesc(hap));
		}

		return snp;
	}

	public SNPFeature asFeature(ArrayList<BitVector> haps, Integer referenceHapId) {
		Chromosome chr = variants.get(0).getChrObj();
		int start = variants.get(0).getPos();
		int stop = variants.get(variants.size() - 1).getPos();

		double[] afsControls = new double[haps.size()];
		double[] afsCases = new double[haps.size()];
		String[] hapDescs = new String[haps.size()];
		for (int b = 0; b < haps.size(); b++) {
			Integer id = hapToInt.get(haps.get(b));
			afsControls[b] = haplotypeFrequenciesControls[id];
			afsCases[b] = haplotypeFrequenciesCases[id];
			hapDescs[b] = getHaplotypeDesc(id);
		}

		SNPFeature snp = new SNPFeature(chr, start, stop);
		snp.setAFCases(afsCases);
		snp.setAFControls(afsControls);
		snp.setName(getVariantsAsString());


		snp.setAlleles(new String[]{getHaplotypeDesc(referenceHapId), Strings.concat(hapDescs, Strings.comma)});
		snp.setMinorAllele(getHaplotypeDesc(minorHapId));


		return snp;
	}

	public int getNrIndividuals() {
		return finalCovariates.rows();
	}

	public DoubleMatrix2D getIntercept() {
		DoubleMatrix2D intercept = DoubleFactory2D.dense.make(getNrIndividuals(), 1);
		intercept.assign(1);
		return intercept;
	}

	public double[] getDiseaseStatus() {
		double[] diseaseStatus = new double[getNrIndividuals()];
		for (int d = 0; d < finalDiseaseStatus.length; d++) {
			diseaseStatus[d] = finalDiseaseStatus[d][0].getNumber();
		}
		return diseaseStatus;
	}

	public String getVariantsAsString() {
		String variantStr = "";
		for (VCFVariant variant : variants) {
			if (variantStr.length() == 0) {
				variantStr += variant.getId();
			} else {
				variantStr += "," + variant.getId();
			}
		}

		return variantStr;
	}

	public int getNrHaplotypes() {
		return availableHaplotypesList.size();
	}

	public DoubleMatrix2D getHaplotype(int hap) {
		DoubleMatrix2D haps = haplotypeMatrix.viewColumn(hap).reshape(getNrIndividuals(), 1);
		return haps;
	}

	public Integer getHaplotypeId(BitVector referenceHaplotype) {
		return hapToInt.get(referenceHaplotype);
	}

	public ArrayList<Integer> getOriginalIndsIncluded() {
		return originalIndsIncluded;
	}

	public HaplotypeDataset selectHaplotypes(ArrayList<Integer> hapsToSelect) {

//		get a list of haplotypes below the threshold
		HashSet<Integer> keepTheseIndividuals = new HashSet<>();
		ArrayList<BitVector> haplotypesAboveThreshold = new ArrayList<>();
		for (int b = 0; b < hapsToSelect.size(); b++) {
			haplotypesAboveThreshold.add(availableHaplotypesList.get(hapsToSelect.get(b)));
		}

		HashSet<BitVector> allowedHaplotypes = new HashSet<>();
		allowedHaplotypes.addAll(haplotypesAboveThreshold);
		ArrayList<Integer> origIds = new ArrayList<>();
		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			BitVector[] haps = haplotypePairs[i];
			if (haps != null) {
				if (allowedHaplotypes.contains(haps[0]) && allowedHaplotypes.contains(haps[1])) {
					keepTheseIndividuals.add(i);
					origIds.add(originalIndsIncluded.get(i));

				}
			}
		}

		System.out.println(keepTheseIndividuals.size() + " individuals will remain after pruning.");
		if (keepTheseIndividuals.isEmpty()) {
			return null;
		}

		// prune haplotypes, covariates and disease status
		DiseaseStatus[][] remainingdiseaseStatus = new DiseaseStatus[keepTheseIndividuals.size()][1];
		DoubleMatrix2D remainingCovariates = new DenseDoubleMatrix2D(keepTheseIndividuals.size(), finalCovariates.columns());

		BitVector[][] remainingHaplotypes = new BitVector[keepTheseIndividuals.size()][2];

		int ictr = 0;
		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			BitVector[] haps = haplotypePairs[i];
			if (keepTheseIndividuals.contains(i)) {
				remainingHaplotypes[ictr] = haps;

				remainingdiseaseStatus[ictr][0] = finalDiseaseStatus[i][0];
				for (int c = 0; c < finalCovariates.columns(); c++) {
					remainingCovariates.setQuick(ictr, c, finalCovariates.getQuick(i, c));
				}
				ictr++;
			}
		}
		return new HaplotypeDataset(remainingHaplotypes, haplotypesAboveThreshold, remainingdiseaseStatus, remainingCovariates, origIds, variants);
	}

	public HaplotypeGroup collapse(Integer haplotype, boolean[] varCombo, int ctr) throws IOException {

		HashSet<BitVector> haplotypesToCollapseOn = new HashSet<>();
		if (haplotype != null) {
			// collapse on a single haplotype
			haplotypesToCollapseOn.add(availableHaplotypesList.get(haplotype));
		} else {
			// use set to determine on which alleles to collapse
			for (int a = 0; a < availableHaplotypesList.size(); a++) {
				BitVector v = availableHaplotypesList.get(a);
				ArrayList<Integer> positions = new ArrayList<>();
				for (int i = 0; i < varCombo.length; i++) {
					if (varCombo[i]) {
						positions.add(i);
					}
				}

				boolean include = true;
				for (int pos : positions) {
					if (!v.get(pos)) {
						include = false;
					}
				}
				if (include) {
					haplotypesToCollapseOn.add(v);
				}
			}
		}

		if (haplotypesToCollapseOn.isEmpty()) {
			return null;
		}

		// now collapse on selected haplotypes
		DoubleMatrix2D haplotypeMatrixCollapsed = new DenseDoubleMatrix2D(getNrIndividuals(), 1);
		ArrayList<BitVector> group0 = new ArrayList<>();
		ArrayList<BitVector> group1 = new ArrayList<>();


		System.out.println(haplotypesToCollapseOn.size() + " haps to collapse on.");
		for (int a = 0; a < availableHaplotypesList.size(); a++) {
			BitVector hap = getAvailableHaplotypesList().get(a);
			if (haplotypesToCollapseOn.contains(hap)) {
				group1.add(hap);
				System.out.println("including:\t" + getHaplotypeDesc(a));
			} else {
				group0.add(hap);
			}
		}

		TextFile outf = new TextFile("/Data/tmp/tnfaip3/haplogroupout-" + ctr + ".txt", TextFile.W);

		String header = "ind";
		for (int a = 0; a < availableHaplotypesList.size(); a++) {
			header += "\t" + getHaplotypeDesc(a);
		}

		HaplotypeGroup groupaaaa = new HaplotypeGroup(haplotypeMatrixCollapsed, group0, group1, variants, finalDiseaseStatus);
		String desc = groupaaaa.getHaplotypeDesc(group1);
		header += "\thap1\thap2\tcollapsed (" + desc + ")";

		outf.writeln(header);


		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			BitVector[] haps = haplotypePairs[i];

			boolean hap1Selected = haplotypesToCollapseOn.contains(haps[0]);
			boolean hap2Selected = haplotypesToCollapseOn.contains(haps[1]);

			String hapStr = "";
			for (int a = 0; a < availableHaplotypesList.size(); a++) {
				double haplo = getHaplotype(a).get(i, 0);
				if (hapStr.length() == 0) {
					hapStr += "" + haplo;
				} else {
					hapStr += "\t" + haplo;
				}
			}


			if (hap1Selected && hap2Selected) {
				haplotypeMatrixCollapsed.set(i, 0, 2);
			} else if (hap1Selected || hap2Selected) {
				haplotypeMatrixCollapsed.set(i, 0, 1);
			}

			outf.writeln(i + "\t" + hapStr + "\t" + getHaplotypeDesc(haps[0]) + "\t" + getHaplotypeDesc(haps[1]) + "\t" + haplotypeMatrixCollapsed.get(i, 0));
		}
		outf.close();

		HaplotypeGroup group = new HaplotypeGroup(haplotypeMatrixCollapsed, group0, group1, variants, finalDiseaseStatus);


		return group;
	}

	private String getHaplotypeDesc(BitVector hap) {
		String hapdesc = "";
		for (int b = 0; b < hap.size(); b++) {
			String all = "";
			if (hap.get(b)) {
				all = variants.get(b).getAlleles()[1];
			} else {
				all = variants.get(b).getAlleles()[0];
			}
			if (hapdesc.length() == 0) {
				hapdesc += all;
			} else {
				hapdesc += "|" + all;
			}
		}
		return hapdesc;

	}

	public Pair<DoubleMatrix2D, ArrayList<BitVector>> getHaplotypesExcludeReference() {

		Integer id = hapToInt.get(referenceHaplotype);
		DoubleMatrix2D matrix = new DenseDoubleMatrix2D(haplotypeMatrix.rows(), haplotypeMatrix.columns() - 1);
		for (int i = 0; i < haplotypeMatrix.rows(); i++) {
			int ctr = 0;
			for (int j = 0; j < haplotypeMatrix.columns(); j++) {
				if (j != id) {
					matrix.setQuick(i, ctr, haplotypeMatrix.getQuick(i, j));
					ctr++;
				}
			}
		}

		ArrayList<BitVector> remainingHaps = new ArrayList<>();
		for (int i = 0; i < availableHaplotypesList.size(); i++) {
			if (i != id) {
				remainingHaps.add(availableHaplotypesList.get(i));
			}
		}

		return new Pair<DoubleMatrix2D, ArrayList<BitVector>>(matrix, remainingHaps);
	}

	public Integer getReferenceHaplotypeId() {
		return hapToInt.get(referenceHaplotype);
	}

	public void createGroups(ArrayList<Integer> groupVars) {

		ArrayList<ArrayList<BitVector>> hapspervar = new ArrayList<>();
		for (int i = 0; i < groupVars.size(); i++) {
			Integer v = groupVars.get(i);
			ArrayList<BitVector> haps = new ArrayList<>();
			for (int b = 0; b < availableHaplotypesList.size(); b++) {
				BitVector vector = availableHaplotypesList.get(b);
				if (vector.get(v)) {
					haps.add(vector);
				}
			}
		}

		// create a new haplogroup is any of the variants
		// overlap haplotypes.. i.e. when a hap is present for > 1 variants
		int[] hapsctr = new int[availableHaplotypesList.size()];
		for (int i = 0; i < hapspervar.size(); i++) {
			ArrayList<BitVector> haps = hapspervar.get(i);
			for (int j = 0; j < haps.size(); j++) {
				Integer id = hapToInt.get(haps.get(j));
				hapsctr[id]++;
			}
		}

		HashSet<Integer> ctrHash = new HashSet<Integer>();
		for (int i = 0; i < hapsctr.length; i++) {
			ctrHash.add(hapsctr[i]);
		}

		ArrayList<ArrayList<BitVector>> groups = new ArrayList<>();
		for (Integer q : ctrHash) {
			ArrayList<BitVector> group = new ArrayList<>();
			for (int z = 0; z < hapsctr.length; z++) {
				if (q.equals(hapsctr[z])) {
					group.add(availableHaplotypesList.get(z));
				}
			}
			groups.add(group);
		}

		DoubleMatrix2D matrix = new DenseDoubleMatrix2D(haplotypeMatrix.rows(), groups.size());
		for (int g = 0; g < groups.size(); g++) {
			ArrayList<BitVector> group = groups.get(g);

			for (int i = 0; i < haplotypeMatrix.rows(); i++) {
				for (int z = 0; z < group.size(); z++) {
					BitVector v = group.get(z);
					Integer id = hapToInt.get(v);
					double val = haplotypeMatrix.get(i, id);
					double val2 = matrix.getQuick(i, g);
					matrix.setQuick(i, g, val + val2);
				}
			}
		}
	}
}
