package nl.harmjanwestra.gwas;

import cern.colt.matrix.tbit.BitVector;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.gwas.tasks.LRTestHaploTestTask;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.math.LogisticRegressionOptimized;
import nl.harmjanwestra.utilities.math.LogisticRegressionResult;
import nl.harmjanwestra.utilities.sets.Bits;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import org.apache.commons.math3.util.CombinatoricsUtils;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ChiSquare;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.concurrent.Executors;

/**
 * Created by hwestra on 7/5/16.
 */
public class LRTestHaplotype extends LRTest {

	private final DiseaseStatus[][] finalDiseaseStatus;
	private final DoubleMatrix2D finalCovariates;
	LRTestHaploTestTask lrht = new LRTestHaploTestTask();
	LogisticRegressionOptimized reg = new LogisticRegressionOptimized();


	public LRTestHaplotype(LRTestOptions options) throws IOException {
		super(options);
		exService = Executors.newFixedThreadPool(options.getNrThreads());
		this.finalDiseaseStatus = sampleAnnotation.getSampleDiseaseStatus();
		this.finalCovariates = sampleAnnotation.getCovariates();
		testHaplotype();

		exService.shutdown();
	}


	public void testHaplotype() throws IOException {


		System.out.println("Testing haplotypes");
		options.setSplitMultiAllelic(true);

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(options.getBedfile());
		ArrayList<VCFVariant> variants = readVariants(options.getOutputdir() + "variantlog.txt", regions);
		if (variants.isEmpty()) {
			System.out.println("Sorry no variants found.");
			System.exit(-1);
		}

		// perform univariate test
		System.out.println();

		HaplotypeDataset ds = HaplotypeDataset.create(variants, finalDiseaseStatus, finalCovariates, options.getGenotypeProbThreshold(), exService);
		int nrHapsPrev = ds.getAvailableHaplotypesList().size();

		HaplotypeDataset dspruned = null;
		int delta = 1;
		int pruningiter = 0;
		while (delta > 0) {
			System.out.println();
			System.out.println("Pruning iteration " + pruningiter);
			dspruned = ds.prune(options.getHaplotypeFrequencyThreshold());
			delta = nrHapsPrev - dspruned.getAvailableHaplotypesList().size();
			nrHapsPrev = dspruned.getAvailableHaplotypesList().size();
			pruningiter++;
		}

//		univariate(dspruned);

//		exhaustive(variants);

		collapse(variants);

		System.exit(-1);
	}

	private void collapse(ArrayList<VCFVariant> variants) throws IOException {
		System.out.println("Collapsing tests / defining haplogroups");
		String outloc = options.getOutputdir() + "haptest-collapsed.txt";
		TextFile out = new TextFile(outloc, TextFile.W);

		AssociationFile f = new AssociationFile();
		out.writeln(f.getHeader());
		int haptr = 0;
		ArrayList<Pair<HaplotypeGroup, AssociationResult>> allResults = new ArrayList<Pair<HaplotypeGroup, AssociationResult>>();
		int nrMaxVariants = variants.size() + 1;
		for (int nrVariants = 2; nrVariants < nrMaxVariants; nrVariants++) {
			System.out.println("Nr Variants: " + nrVariants);
			Iterator<int[]> combos = CombinatoricsUtils.combinationsIterator(variants.size(), nrVariants);
			while (combos.hasNext()) {
				int[] combo = combos.next();

				ArrayList<VCFVariant> selectedVars = new ArrayList<>();
				for (int q = 0; q < combo.length; q++) {
					selectedVars.add(variants.get(combo[q]));
				}

				HaplotypeDataset ds = HaplotypeDataset.create(selectedVars, finalDiseaseStatus, finalCovariates, options.getGenotypeProbThreshold(), exService);

				int nrHapsPrev = ds.getAvailableHaplotypesList().size();
				HaplotypeDataset dspruned = null;
				int delta = 1;
				int pruningiter = 0;
				while (delta > 0) {
					System.out.println();
					System.out.println("Pruning iteration " + pruningiter);
					dspruned = ds.prune(options.getHaplotypeFrequencyThreshold());
					delta = nrHapsPrev - dspruned.getAvailableHaplotypesList().size();
					nrHapsPrev = dspruned.getAvailableHaplotypesList().size();
					pruningiter++;
				}

				// make combinations of haplotypes
				int[] values = Bits.valuesForBitCombinations(nrVariants);
				for (int v = 0; v < values.length; v++) {
					boolean[] set = Bits.convertAsBoolean(values[v], nrVariants);

					boolean allfalse = true;
					boolean alltrue = true;

					String setStr = "";
					for (int s = 0; s < set.length; s++) {
						if (set[s]) {
							allfalse = false;
							setStr += "1";
						} else {
							alltrue = false;
							setStr += "0";
						}
					}

					if (!alltrue && !allfalse) {
						HaplotypeGroup group = dspruned.collapse(null, set, haptr);
						if (group != null && group.getFreqControls() != 1d) {
							System.out.println(setStr + "\t" + group.toString());
							AssociationResult r = testHaploGroup(dspruned, group);
							addHaploGroupResult(new Pair<HaplotypeGroup, AssociationResult>(group, r), allResults);
						}
						haptr++;
					}
				}

				// test individual haplotypes
				for (int v = 0; v < nrVariants; v++) {
					HaplotypeGroup group = dspruned.collapse(v, null, haptr);
					if (group != null && group.getFreqControls() != 1d) {
						AssociationResult r = testHaploGroup(dspruned, group);
						addHaploGroupResult(new Pair<HaplotypeGroup, AssociationResult>(group, r), allResults);
					}
					haptr++;
				}

			}
		}
		out.close();

		HaplotypeGroupPlot plot = new HaplotypeGroupPlot();
		try {
			plot.run(allResults, variants, "/Data/tmp/tnfaip3/plot.pdf");
		} catch (DocumentException e) {
			e.printStackTrace();
		}

	}

	private void addHaploGroupResult(Pair<HaplotypeGroup, AssociationResult> haplotypeGroupAssociationResultPair, ArrayList<Pair<HaplotypeGroup, AssociationResult>> allResults) {
		if (allResults.isEmpty()) {
			allResults.add(haplotypeGroupAssociationResultPair);
		} else {

			ArrayList<VCFVariant> variants2 = haplotypeGroupAssociationResultPair.getLeft().getVariants();
			ArrayList<BitVector> group0 = haplotypeGroupAssociationResultPair.getLeft().getGroup0();
			ArrayList<BitVector> group1 = haplotypeGroupAssociationResultPair.getLeft().getGroup1();

			HashSet<VCFVariant> variantSet = new HashSet<VCFVariant>();
			variantSet.addAll(variants2);

			HashSet<BitVector> group0Set = new HashSet<BitVector>();
			group0Set.addAll(group0);
			HashSet<BitVector> group1Set = new HashSet<BitVector>();
			group1Set.addAll(group1);

			boolean include = true;
			for (Pair<HaplotypeGroup, AssociationResult> p : allResults) {
				boolean sharedVariants = false;
				boolean sharedHaplogroups = false;

				ArrayList<VCFVariant> variants1 = p.getLeft().getVariants();
				ArrayList<BitVector> groupP0 = p.getLeft().getGroup0();
				ArrayList<BitVector> groupP1 = p.getLeft().getGroup1();


				if (variants1.size() == variants2.size()) {
					// check whether result shares same variants
					int nrShared = 0;
					for (VCFVariant v : variants1) {
						if (variantSet.contains(v)) {
							nrShared++;
						}
					}

					if (nrShared == variantSet.size()) {
						sharedVariants = true;
						// now check whether the haplogroups are the same
						if (group0.size() == groupP0.size() && group1.size() == groupP1.size()) {
							int nrGroup0Shared = 0;
							int nrGroup1Shared = 0;

							for (BitVector v : groupP0) {
								if (group0Set.contains(v)) {
									nrGroup0Shared++;
								}
							}
							for (BitVector v : groupP1) {
								if (group1Set.contains(v)) {
									nrGroup1Shared++;
								}
							}

							if (nrGroup0Shared == groupP0.size() && nrGroup1Shared == groupP1.size()) {
								sharedHaplogroups = true;
							}
						}
					}
				}

				if (sharedHaplogroups && sharedVariants) {
					include = false;
				}
			}

			if (include) {
				allResults.add(haplotypeGroupAssociationResultPair);
			}
		}
	}

	private AssociationResult testHaploGroup(HaplotypeDataset dspruned, HaplotypeGroup group) {

//							System.exit(-1);
		DoubleMatrix2D intercept = dspruned.getIntercept();
		DoubleMatrix2D xprime = DoubleFactory2D.dense.appendColumns(intercept, dspruned.getFinalCovariates());
		double[] y = dspruned.getDiseaseStatus();
		LogisticRegressionResult resultCovars = reg.univariate(y, xprime);

		DoubleMatrix2D interceptWHaps = DoubleFactory2D.dense.appendColumns(intercept, group.getHaplotypeMatrixCollapsed());
		DoubleMatrix2D x = DoubleFactory2D.dense.appendColumns(interceptWHaps, dspruned.getFinalCovariates());
		LogisticRegressionResult resultX = reg.univariate(y, x);

		int nrRemaining = 1;
		double devx = resultX.getDeviance();
		double devnull = resultCovars.getDeviance();
		double[] betasmlelr = new double[nrRemaining];
		double[] stderrsmlelr = new double[nrRemaining];
		double[] or = new double[nrRemaining];
		double[] orhi = new double[nrRemaining];
		double[] orlo = new double[nrRemaining];
		int ctr = 0;

		double beta = -resultX.getBeta()[1];
		double se = resultX.getStderrs()[1];
		betasmlelr[ctr] = beta;
		stderrsmlelr[ctr] = se;

		double OR = Math.exp(beta);
		double orLow = Math.exp(beta - 1.96 * se);
		double orHigh = Math.exp(beta + 1.96 * se);
		or[ctr] = OR;
		orhi[ctr] = orHigh;
		orlo[ctr] = orLow;

		double deltaDeviance = devnull - devx;
		int df = 1;
		double p = ChiSquare.getP(df, deltaDeviance);

		AssociationResult result = new AssociationResult();
		result.setDevianceNull(devnull);
		result.setDevianceGeno(devx);
		result.setDf(df);
		result.setDfalt(x.columns());
		result.setDfnull(xprime.columns());

		result.setBeta(betasmlelr);
		result.setSe(stderrsmlelr);
		result.setPval(p);

		SNPFeature snp = group.asFeature();
		result.setSnp(snp);
		result.setN(x.rows());

		snp.setImputationQualityScore(1d);
		return result;
	}

	private void exhaustive(ArrayList<VCFVariant> variants) throws IOException {
		System.out.println("Conditional tests");

		String outloc = options.getOutputdir() + "haptest-exhaustive.txt";
		String outloc2 = options.getOutputdir() + "haptest-exhaustive2.txt";

		System.out.println("Output is here: " + outloc);
		TextFile out = new TextFile(outloc, TextFile.W);
		TextFile out2 = new TextFile(outloc2, TextFile.W);
		AssociationFile f = new AssociationFile();
		out.writeln(f.getHeader());
		out2.writeln(f.getHeader());


		for (int nrVariants = 2; nrVariants < variants.size(); nrVariants++) {
			System.out.println("Nr Variants: " + nrVariants);
			// create all possibilities of 2 variants
			Iterator<int[]> combos = CombinatoricsUtils.combinationsIterator(variants.size(), nrVariants);
			while (combos.hasNext()) {
				int[] combo = combos.next();

				ArrayList<VCFVariant> selectedVars = new ArrayList<>();
				for (int q = 0; q < combo.length; q++) {
					selectedVars.add(variants.get(combo[q]));
				}

				HaplotypeDataset ds = HaplotypeDataset.create(selectedVars, finalDiseaseStatus, finalCovariates, options.getGenotypeProbThreshold(), exService);

				int nrHapsPrev = ds.getAvailableHaplotypesList().size();
				HaplotypeDataset dspruned = null;
				int delta = 1;
				int pruningiter = 0;
				while (delta > 0) {
					System.out.println();
					System.out.println("Pruning iteration " + pruningiter);
					dspruned = ds.prune(options.getHaplotypeFrequencyThreshold());
					delta = nrHapsPrev - dspruned.getAvailableHaplotypesList().size();
					nrHapsPrev = dspruned.getAvailableHaplotypesList().size();
					pruningiter++;
				}

				for (int hap = 0; hap < dspruned.getNrHaplotypes(); hap++) {
					DoubleMatrix2D intercept = dspruned.getIntercept();
					DoubleMatrix2D xprime = DoubleFactory2D.dense.appendColumns(intercept, dspruned.getFinalCovariates());
					double[] y = dspruned.getDiseaseStatus();
					LogisticRegressionResult resultCovars = reg.univariate(y, xprime);

					DoubleMatrix2D haps = dspruned.getHaplotype(hap);
					DoubleMatrix2D interceptWHaps = DoubleFactory2D.dense.appendColumns(intercept, haps);
					DoubleMatrix2D x = DoubleFactory2D.dense.appendColumns(interceptWHaps, dspruned.getFinalCovariates());
					LogisticRegressionResult resultX = reg.univariate(y, x);

					int nrRemaining = 1;
					double devx = resultX.getDeviance();
					double devnull = resultCovars.getDeviance();
					double[] betasmlelr = new double[nrRemaining];
					double[] stderrsmlelr = new double[nrRemaining];
					double[] or = new double[nrRemaining];
					double[] orhi = new double[nrRemaining];
					double[] orlo = new double[nrRemaining];
					int ctr = 0;

					double beta = -resultX.getBeta()[1];
					double se = resultX.getStderrs()[1];
					betasmlelr[ctr] = beta;
					stderrsmlelr[ctr] = se;

					double OR = Math.exp(beta);
					double orLow = Math.exp(beta - 1.96 * se);
					double orHigh = Math.exp(beta + 1.96 * se);
					or[ctr] = OR;
					orhi[ctr] = orHigh;
					orlo[ctr] = orLow;

					double deltaDeviance = devnull - devx;
					int df = 1;
					double p = ChiSquare.getP(df, deltaDeviance);

					AssociationResult result = new AssociationResult();
					result.setDevianceNull(devnull);
					result.setDevianceGeno(devx);
					result.setDf(df);
					result.setDfalt(x.columns());
					result.setDfnull(xprime.columns());

					result.setBeta(betasmlelr);
					result.setSe(stderrsmlelr);
					result.setPval(p);

					SNPFeature snp = dspruned.asFeature(hap, null);
					result.setSnp(snp);
					result.setN(x.rows());

					snp.setImputationQualityScore(1d);
					out.writeln(result.toString());
				}


				Integer referenceIndex = dspruned.getHaplotypeId(dspruned.getReferenceHaplotype());
				for (int hap = 0; hap < dspruned.getNrHaplotypes(); hap++) {
					if (!referenceIndex.equals(hap)) {
						BitVector hapToRun = dspruned.getAvailableHaplotypesList().get(hap);
						ArrayList<Integer> hapsToSelect = new ArrayList<Integer>();

						hapsToSelect.add(referenceIndex);
						hapsToSelect.add(hap);
						System.out.println();
						System.out.println(referenceIndex + "\tvs\t" + hap);
						HaplotypeDataset refDs = dspruned.selectHaplotypes(hapsToSelect);

						DoubleMatrix2D intercept = refDs.getIntercept();
						double[] y = refDs.getDiseaseStatus();
						DoubleMatrix2D xprime = DoubleFactory2D.dense.appendColumns(intercept, refDs.getFinalCovariates());
						LogisticRegressionResult resultCovars = reg.univariate(y, xprime);

						Integer newHapIndex = refDs.getHaplotypeId(hapToRun);
						Integer newRefHapIndex = refDs.getHaplotypeId(dspruned.getReferenceHaplotype());
						DoubleMatrix2D haps = refDs.getHaplotype(newHapIndex);
						DoubleMatrix2D interceptWHaps = DoubleFactory2D.dense.appendColumns(intercept, haps);
						DoubleMatrix2D x = DoubleFactory2D.dense.appendColumns(interceptWHaps, refDs.getFinalCovariates());
						LogisticRegressionResult resultX = reg.univariate(y, x);

						int nrRemaining = 1;
						double devx = resultX.getDeviance();
						double devnull = resultCovars.getDeviance();
						double[] betasmlelr = new double[nrRemaining];
						double[] stderrsmlelr = new double[nrRemaining];
						double[] or = new double[nrRemaining];
						double[] orhi = new double[nrRemaining];
						double[] orlo = new double[nrRemaining];
						int ctr = 0;

						double beta = -resultX.getBeta()[1];
						double se = resultX.getStderrs()[1];
						betasmlelr[ctr] = beta;
						stderrsmlelr[ctr] = se;

						double OR = Math.exp(beta);
						double orLow = Math.exp(beta - 1.96 * se);
						double orHigh = Math.exp(beta + 1.96 * se);
						or[ctr] = OR;
						orhi[ctr] = orHigh;
						orlo[ctr] = orLow;

						double deltaDeviance = devnull - devx;
						int df = 1;
						double p = ChiSquare.getP(df, deltaDeviance);

						AssociationResult result = new AssociationResult();
						result.setDevianceNull(devnull);
						result.setDevianceGeno(devx);
						result.setDf(df);
						result.setDfalt(x.columns());
						result.setDfnull(xprime.columns());

						result.setBeta(betasmlelr);
						result.setSe(stderrsmlelr);
						result.setPval(p);

						SNPFeature snp = refDs.asFeature(newHapIndex, newRefHapIndex);
						result.setSnp(snp);
						result.setN(x.rows());

						snp.setImputationQualityScore(1d);
						out2.writeln(result.toString());
					}
				}


			}
		}

		out2.close();
		out.close();


	}


	private void univariate(HaplotypeDataset dspruned) throws IOException {

		System.out.println("Univariate tests");

		String outloc = options.getOutputdir() + "haptest-univariate.txt";
		System.out.println("Output is here: " + outloc);
		TextFile out = new TextFile(outloc, TextFile.W);
		AssociationFile f = new AssociationFile();
		out.writeln(f.getHeader());

		for (int hap = 0; hap < dspruned.getNrHaplotypes(); hap++) {
			DoubleMatrix2D intercept = dspruned.getIntercept();
			DoubleMatrix2D xprime = DoubleFactory2D.dense.appendColumns(intercept, dspruned.getFinalCovariates());
			double[] y = dspruned.getDiseaseStatus();
			LogisticRegressionResult resultCovars = reg.univariate(y, xprime);

			DoubleMatrix2D haps = dspruned.getHaplotype(hap);
			DoubleMatrix2D interceptWHaps = DoubleFactory2D.dense.appendColumns(intercept, haps);
			DoubleMatrix2D x = DoubleFactory2D.dense.appendColumns(interceptWHaps, dspruned.getFinalCovariates());
			LogisticRegressionResult resultX = reg.univariate(y, x);

			int nrRemaining = 1;
			double devx = resultX.getDeviance();
			double devnull = resultCovars.getDeviance();
			double[] betasmlelr = new double[nrRemaining];
			double[] stderrsmlelr = new double[nrRemaining];
			double[] or = new double[nrRemaining];
			double[] orhi = new double[nrRemaining];
			double[] orlo = new double[nrRemaining];
			int ctr = 0;

			double beta = -resultX.getBeta()[1];
			double se = resultX.getStderrs()[1];
			betasmlelr[ctr] = beta;
			stderrsmlelr[ctr] = se;

			double OR = Math.exp(beta);
			double orLow = Math.exp(beta - 1.96 * se);
			double orHigh = Math.exp(beta + 1.96 * se);
			or[ctr] = OR;
			orhi[ctr] = orHigh;
			orlo[ctr] = orLow;

			double deltaDeviance = devnull - devx;
			int df = 1;
			double p = ChiSquare.getP(df, deltaDeviance);

			AssociationResult result = new AssociationResult();
			result.setDevianceNull(devnull);
			result.setDevianceGeno(devx);
			result.setDf(df);
			result.setDfalt(x.columns());
			result.setDfnull(xprime.columns());

			result.setBeta(betasmlelr);
			result.setSe(stderrsmlelr);
			result.setPval(p);

			SNPFeature snp = dspruned.asFeature(hap, null);
			result.setSnp(snp);
			result.setN(x.rows());

			snp.setImputationQualityScore(1d);
			out.writeln(result.toString());
		}

		out.close();

		outloc = options.getOutputdir() + "haptest-univariate-relative.txt";
		System.out.println("Output is here for relative univariate: " + outloc);
		out = new TextFile(outloc, TextFile.W);
		out.writeln(f.getHeader());

		// now rerun assuming a reference haplotype
		Integer referenceIndex = dspruned.getHaplotypeId(dspruned.getReferenceHaplotype());
		for (int hap = 0; hap < dspruned.getNrHaplotypes(); hap++) {
			if (!referenceIndex.equals(hap)) {
				BitVector hapToRun = dspruned.getAvailableHaplotypesList().get(hap);
				ArrayList<Integer> hapsToSelect = new ArrayList<Integer>();

				hapsToSelect.add(referenceIndex);
				hapsToSelect.add(hap);
				System.out.println();
				System.out.println(referenceIndex + "\tvs\t" + hap);
				HaplotypeDataset refDs = dspruned.selectHaplotypes(hapsToSelect);

				DoubleMatrix2D intercept = refDs.getIntercept();
				double[] y = refDs.getDiseaseStatus();
				DoubleMatrix2D xprime = DoubleFactory2D.dense.appendColumns(intercept, refDs.getFinalCovariates());
				LogisticRegressionResult resultCovars = reg.univariate(y, xprime);

				Integer newHapIndex = refDs.getHaplotypeId(hapToRun);
				Integer newRefHapIndex = refDs.getHaplotypeId(dspruned.getReferenceHaplotype());
				DoubleMatrix2D haps = refDs.getHaplotype(newHapIndex);
				DoubleMatrix2D interceptWHaps = DoubleFactory2D.dense.appendColumns(intercept, haps);
				DoubleMatrix2D x = DoubleFactory2D.dense.appendColumns(interceptWHaps, refDs.getFinalCovariates());
				LogisticRegressionResult resultX = reg.univariate(y, x);

				int nrRemaining = 1;
				double devx = resultX.getDeviance();
				double devnull = resultCovars.getDeviance();
				double[] betasmlelr = new double[nrRemaining];
				double[] stderrsmlelr = new double[nrRemaining];
				double[] or = new double[nrRemaining];
				double[] orhi = new double[nrRemaining];
				double[] orlo = new double[nrRemaining];
				int ctr = 0;

				double beta = -resultX.getBeta()[1];
				double se = resultX.getStderrs()[1];
				betasmlelr[ctr] = beta;
				stderrsmlelr[ctr] = se;

				double OR = Math.exp(beta);
				double orLow = Math.exp(beta - 1.96 * se);
				double orHigh = Math.exp(beta + 1.96 * se);
				or[ctr] = OR;
				orhi[ctr] = orHigh;
				orlo[ctr] = orLow;

				double deltaDeviance = devnull - devx;
				int df = 1;
				double p = ChiSquare.getP(df, deltaDeviance);

				AssociationResult result = new AssociationResult();
				result.setDevianceNull(devnull);
				result.setDevianceGeno(devx);
				result.setDf(df);
				result.setDfalt(x.columns());
				result.setDfnull(xprime.columns());

				result.setBeta(betasmlelr);
				result.setSe(stderrsmlelr);
				result.setPval(p);

				SNPFeature snp = refDs.asFeature(newHapIndex, newRefHapIndex);
				result.setSnp(snp);
				result.setN(x.rows());

				snp.setImputationQualityScore(1d);
				out.writeln(result.toString());
			}
		}
		out.close();
	}

}
