package nl.harmjanwestra.gwas;

import cern.colt.matrix.tbit.BitVector;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import nl.harmjanwestra.gwas.tasks.LRTestHaploTask;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Triple;

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
	private double[] haplotypeFrequenciesCases;
	private double[] haplotypeFrequenciesControls;
	private BitVector referenceHaplotype;
	private int nrHaplotypes;
	private DenseDoubleMatrix2D haplotypeMatrix;
	private HashMap<BitVector, Integer> hapToInt;


	public HaplotypeDataset(BitVector[][] haplotypePairs,
							ArrayList<BitVector> availableHaplotypesList,
							DiseaseStatus[][] finalDiseaseStatus,
							DoubleMatrix2D finalCovariates,
							ArrayList<VCFVariant> variants) {
		this.finalDiseaseStatus = finalDiseaseStatus;
		this.finalCovariates = finalCovariates;
		this.haplotypePairs = haplotypePairs;
		this.availableHaplotypesList = availableHaplotypesList;
		this.variants = variants;
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
										  ExecutorService exService) {
// get a list of available haplotypes
		HashSet<BitVector> availableHaplotypeHash = new HashSet<BitVector>();
		ArrayList<BitVector> availableHaplotypesList = new ArrayList<>();

		CompletionService<Triple<BitVector[], Integer, Boolean>> jobHandler = new ExecutorCompletionService<Triple<BitVector[], Integer, Boolean>>(exService);
		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			jobHandler.submit(new LRTestHaploTask(i, variants, genotypeprobthreshold));
		}

		int returned = 0;
		BitVector[][] haplotypePairs = new BitVector[finalDiseaseStatus.length][];
		ProgressBar pb = new ProgressBar(finalDiseaseStatus.length);
		while (returned < finalDiseaseStatus.length) {

			try {
				Triple<BitVector[], Integer, Boolean> fut = jobHandler.take().get();
				BitVector[] haps = fut.getLeft();

				if (!fut.getRight()) {
					haplotypePairs[fut.getMiddle()] = haps;


					if (!availableHaplotypeHash.contains(haps[0])) {
						availableHaplotypesList.add(haps[0]);
						availableHaplotypeHash.add(haps[0]);
					}
					if (!availableHaplotypeHash.contains(haps[1])) {
						availableHaplotypesList.add(haps[1]);
						availableHaplotypeHash.add(haps[1]);
					}

				} else {
					haplotypePairs[fut.getMiddle()] = null;
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

		return new HaplotypeDataset(haplotypePairs,
				availableHaplotypesList,
				finalDiseaseStatus,
				finalCovariates,
				variants);
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
		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			BitVector[] haps = haplotypePairs[i];
			if (haps != null) {
				if (allowedHaplotypes.contains(haps[0]) && allowedHaplotypes.contains(haps[1])) {
					keepTheseIndividuals.add(i);
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
		return new HaplotypeDataset(remainingHaplotypes, haplotypesAboveThreshold, remainingdiseaseStatus, remainingCovariates, variants);
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
		for (int i = 0; i < availableHaplotypesList.size(); i++) {
			haplotypeFrequenciesCases[i] /= totalChromosomesCases;
			haplotypeFrequenciesControls[i] /= totalChromosomesControls;
			if (haplotypeFrequenciesControls[i] > maxfreq) {
				reference = availableHaplotypesList.get(i);
				maxfreq = haplotypeFrequenciesControls[i];
			}
			System.out.println(getHaplotypeDesc(i) + "\t" + haplotypeFrequenciesCases[i] + "\t" + haplotypeFrequenciesControls[i]);
		}
		referenceHaplotype = reference;
	}

	public String getHaplotypeDesc(int i) {

		BitVector bitVector = availableHaplotypesList.get(i);
		String hapdesc = "";
		for (int b = 0; b < bitVector.size(); b++) {
			String all = "";
			if (bitVector.get(b)) {
				all = variants.get(b).getAlleles()[1];
			} else {
				all = variants.get(b).getAlleles()[0];
			}
			if (hapdesc.length() == 0) {
				hapdesc += all;
			} else {
				hapdesc += ";" + all;
			}
		}
		return hapdesc;
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
			if(!availableHaplotypesList.get(hap).equals(referenceHaplotype)){
				snp.setMinorAllele(getHaplotypeDesc(hap));
			} else {
				snp.setMinorAllele(getHaplotypeDesc(hapToInt.get(referenceHaplotype)));
			}
		} else {
			String[] alleles = new String[availableHaplotypesList.size()];
			for(int h=0;h<availableHaplotypesList.size();h++){
				alleles[h] = getHaplotypeDesc(h);
			}
			snp.setAlleles(alleles);
			snp.setMinorAllele(getHaplotypeDesc(hap));
		}

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
				variantStr += ";" + variant.getId();
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

	public HaplotypeDataset selectHaplotypes(ArrayList<Integer> hapsToSelect) {

//		get a list of haplotypes below the threshold
		HashSet<Integer> keepTheseIndividuals = new HashSet<>();
		ArrayList<BitVector> haplotypesAboveThreshold = new ArrayList<>();
		for (int b = 0; b < hapsToSelect.size(); b++) {
			haplotypesAboveThreshold.add(availableHaplotypesList.get(hapsToSelect.get(b)));
		}

		HashSet<BitVector> allowedHaplotypes = new HashSet<>();
		allowedHaplotypes.addAll(haplotypesAboveThreshold);
		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			BitVector[] haps = haplotypePairs[i];
			if (haps != null) {
				if (allowedHaplotypes.contains(haps[0]) && allowedHaplotypes.contains(haps[1])) {
					keepTheseIndividuals.add(i);
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
		return new HaplotypeDataset(remainingHaplotypes, haplotypesAboveThreshold, remainingdiseaseStatus, remainingCovariates, variants);
	}


//	public String getHaplotypeComboDescription(ArrayList<BitVector> comboToTest, ArrayList<VCFVariant> variants) {
//
//		String[] haps = new String[comboToTest.size()];
//		for (int i = 0; i < comboToTest.size(); i++) {
//			haps[i] = getHaplotypeDesc(comboToTest.get(i), variants);
//		}
//
//		return Strings.concat(haps, Strings.semicolon);
//	}
}
