package nl.harmjanwestra.gwas;

import cern.colt.matrix.tbit.BitVector;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import nl.harmjanwestra.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.gwas.tasks.LRTestHaploTask;
import nl.harmjanwestra.gwas.tasks.LRTestHaploTestTask;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.math.LogisticRegressionOptimized;
import nl.harmjanwestra.utilities.math.LogisticRegressionResult;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ChiSquare;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.Executors;

/**
 * Created by hwestra on 7/5/16.
 */
public class LRTestHaplotype extends LRTest {

	LRTestHaploTestTask lrht = new LRTestHaploTestTask();

	private DiseaseStatus[][] finalDiseaseStatus;
	private DoubleMatrix2D finalCovariates;


	public LRTestHaplotype(LRTestOptions options) throws IOException {
		super(options);
		exService = Executors.newFixedThreadPool(options.getNrThreads());
		this.finalDiseaseStatus = sampleAnnotation.getSampleDiseaseStatus();
		this.finalCovariates = sampleAnnotation.getCovariates();
		testHaplotype();
		exService.shutdown();
	}


	public void testHaplotype() throws IOException {


		System.out.println("Testing haplotype thingies");
		options.setSplitMultiAllelic(true);
		ArrayList<VCFVariant> variants = readVariants(options.getOutputdir() + "variantlog.txt");
		Triple<DoubleMatrix2D, ArrayList<BitVector>, BitVector[][]> haplotypeData = createHaplotypes(variants);


		DoubleMatrix2D hapdata = haplotypeData.getLeft();
		for (int i = 0; i < hapdata.rows(); i++) {
			int sum = 0;
			for (int j = 0; j < hapdata.columns(); j++) {
				sum += hapdata.get(i, j);
			}
			if (sum != 2) {
				System.err.println("SUM!=2: " + sum + " individual " + i);
			}
		}

		BitVector[][] haplotypePairs = haplotypeData.getRight();
		ArrayList<BitVector> haplotypeWithFrequencyAboveThreshold = haplotypeData.getMiddle();
		Triple<ArrayList<BitVector>, BitVector, double[]> frequencyData = haplotypeFrequency(haplotypeWithFrequencyAboveThreshold, haplotypePairs, variants);

		BitVector referenceHaplotype = frequencyData.getMiddle();
		double[] haplotypeFrequencies = frequencyData.getRight();

		int start = variants.get(0).getPos();
		int stop = variants.get(variants.size() - 1).getPos();
		Chromosome chr = variants.get(0).getChrObj();

		// perform univariate test
		System.out.println();
		System.out.println("Univariate tests");
		System.out.println("Output is here: " + options.getOutputdir() + "haptest.txt");
		TextFile out = new TextFile(options.getOutputdir() + "haptest.txt", TextFile.W);
		AssociationFile f = new AssociationFile();
		out.writeln(f.getHeader());

		univariate(out, hapdata, haplotypeWithFrequencyAboveThreshold, variants);

		out.writeln("multivariate");
		int refhapcol = 0;
		double maxfreq = 0;
		for (int i = 0; i < haplotypeFrequencies.length; i++) {
			if (haplotypeFrequencies[i] > maxfreq) {
				refhapcol = i;
				maxfreq = haplotypeFrequencies[i];
			}
		}
		multivariate(out, hapdata, haplotypeWithFrequencyAboveThreshold, variants, refhapcol);

		// perform multivariate test


		Triple<String, AssociationResult, VCFVariant> result = lrht.calc(null,
				referenceHaplotype,
				haplotypeWithFrequencyAboveThreshold,
				null,
				haplotypePairs,
				finalDiseaseStatus,
				finalCovariates,
				start,
				stop,
				chr,
				variants);
		out.writeln(result.getMiddle().toString());
		System.out.println();
		System.out.println("multivariate:\t" + result.getMiddle().getLog10Pval());
		System.out.println("Reference haplotype: " + lrht.getHaplotypeDesc(referenceHaplotype, variants));
		out.writeln();
		System.out.println("Pairwise");

		// try other combinations of two haplotypes
		ArrayList<BitVector> availableHaplotypesForConditionalAnalysis = new ArrayList<>();
		for (BitVector b : haplotypeWithFrequencyAboveThreshold) {
			if (!b.equals(referenceHaplotype)) {
				availableHaplotypesForConditionalAnalysis.add(b);
			}
		}
		System.out.println(availableHaplotypesForConditionalAnalysis.size() + " haps for conditional analysis");

		System.out.println("Pairwise tests...");
		for (int i = 0; i < availableHaplotypesForConditionalAnalysis.size(); i++) {
			BitVector hap1 = availableHaplotypesForConditionalAnalysis.get(i);
			for (int j = i + 1; j < availableHaplotypesForConditionalAnalysis.size(); j++) {
				BitVector hap2 = availableHaplotypesForConditionalAnalysis.get(j);
				ArrayList<BitVector> toTest = new ArrayList<>();
				toTest.add(hap1);
				toTest.add(hap2);

				String haplotypeComboName = lrht.getHaplotypeComboDescription(toTest, variants);
				result = lrht.calc(null,
						referenceHaplotype,
						toTest,
						null,
						haplotypePairs,
						finalDiseaseStatus,
						finalCovariates,
						start,
						stop,
						chr,
						variants);
				System.out.println(haplotypeComboName + "\t" + result.getMiddle().getLog10Pval());
				out.writeln(result.getMiddle().toString());
			}
		}

		System.out.println();


		out.writeln("conditional");
		System.out.println("Conditional analysis");
		// perform conditional analyses
		// make all possible combinations of haplotypes


		HaplotypeCombiner combiner = new HaplotypeCombiner(availableHaplotypesForConditionalAnalysis);
		combiner.combine();
		ArrayList<List<BitVector>> haplotypeCombinations = combiner.getResults();


		for (int c = 0; c < haplotypeCombinations.size(); c++) {
			ArrayList<BitVector> comboToTest = new ArrayList<>();
			comboToTest.addAll(haplotypeCombinations.get(c));

			HashSet<BitVector> combosToTestHash = new HashSet<BitVector>();
			combosToTestHash.addAll(comboToTest);

			// make a list of available haps to condition on
			ArrayList<BitVector> possibleConditionalHaplotypes = new ArrayList<>();
			for (int j = 0; j < availableHaplotypesForConditionalAnalysis.size(); j++) {
				BitVector conditionalhaplotype = availableHaplotypesForConditionalAnalysis.get(j);

				boolean conditionOnThisHaplotype = true;

				if (combosToTestHash.contains(conditionalhaplotype)) {
					conditionOnThisHaplotype = false;
				}

				if (conditionOnThisHaplotype) {
					possibleConditionalHaplotypes.add(conditionalhaplotype);
				}
			}


			if (!possibleConditionalHaplotypes.isEmpty()) {
				String haplotypeComboName = lrht.getHaplotypeComboDescription(comboToTest, variants);

				combiner = new HaplotypeCombiner(possibleConditionalHaplotypes);
				combiner.combine();
				ArrayList<List<BitVector>> combinations = combiner.getResults();
				for (List<BitVector> combination : combinations) {
					ArrayList<BitVector> conditonalCombo = new ArrayList<>();
					conditonalCombo.addAll(combination);
					String conditionalHaplotypeComboName = lrht.getHaplotypeComboDescription(conditonalCombo, variants);
					result = lrht.calc(null,
							referenceHaplotype,
							comboToTest,
							conditonalCombo,
							haplotypePairs,
							finalDiseaseStatus,
							finalCovariates,
							start,
							stop,
							chr,
							variants);
					System.out.println(haplotypeComboName + "\t" + conditionalHaplotypeComboName + "\t" + result.getMiddle().getLog10Pval());
					out.writeln(result.getMiddle().toString());
				}
			}
		}
		out.close();
	}

	private void univariate(TextFile out, DoubleMatrix2D hapdata, ArrayList<BitVector> availableHaplotypes, ArrayList<VCFVariant> variants) throws IOException {

		// make the column for the intercept
		DoubleMatrix2D intercept = DoubleFactory2D.dense.make(hapdata.rows(), 1);
		intercept.assign(1);

//		DoubleMatrix2D diseaseStatus = new DenseDoubleMatrix2D(hapdata.rows(), 1);

		if (intercept.rows() != finalDiseaseStatus.length) {
			System.out.println("????");
			System.exit(-1);
		}

		double[] y = new double[finalDiseaseStatus.length];
		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			y[i] = finalDiseaseStatus[i][0].getNumber();
		}
		DoubleMatrix2D xprime = DoubleFactory2D.dense.appendColumns(intercept, finalCovariates);

		for (int hap = 0; hap < hapdata.columns(); hap++) {
			DoubleMatrix2D x = DoubleFactory2D.dense.appendColumns(intercept, hapdata.viewColumn(hap).reshape(hapdata.rows(), 1));
			x = DoubleFactory2D.dense.appendColumns(x, finalCovariates);


			System.out.println(x.rows() + "\t" + xprime.rows() + "\t" + y.length);
			System.out.println(x.columns() + "\t" + xprime.columns() + "\t" + 1);
			System.out.println();

			LogisticRegressionOptimized reg = new LogisticRegressionOptimized();
			LogisticRegressionResult resultCovars = reg.univariate(y, xprime);
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

			SNPFeature snp = new SNPFeature(Chromosome.ONE, 1, 1);


			result.setSnp(snp);
			result.setN(x.rows());
			snp.setMaf(0d);
			snp.setHwep(0d);

			Double imputationqualityscore = 1d;
			snp.setImputationQualityScore(imputationqualityscore);

			LRTestHaploTestTask lrht = new LRTestHaploTestTask();


			snp.setAlleles(new String[]{"", ""});
			snp.setName(lrht.getHaplotypeDesc(availableHaplotypes.get(hap), variants));
			snp.setMinorAllele(lrht.getHaplotypeDesc(availableHaplotypes.get(hap), variants));
//			} else {
//				ArrayList<String> haplotypeNames = new ArrayList<>();
//				haplotypeNames.add(getHaplotypeDesc(refHaplotype, variants));
//				for (int q = 0; q < haplotypesToTest.size(); q++) {
//					if (!haplotypesToTest.get(q).equals(refHaplotype)) {
//						haplotypeNames.add(getHaplotypeDesc(haplotypesToTest.get(q), variants));
//					}
//				}
//
//				snp.setAlleles(haplotypeNames.toArray(new String[0]));
//				snp.setMinorAllele("");
//			}
			out.writeln(result.toString());
			System.out.println(lrht.getHaplotypeComboDescription(availableHaplotypes, variants) + "\t" + Strings.concat(betasmlelr, Strings.semicolon) + "\t" + Strings.concat(or, Strings.semicolon) + "\t" + p);
		}
	}

	public void multivariate(TextFile out, DoubleMatrix2D hapdata, ArrayList<BitVector> availableHaplotypes, ArrayList<VCFVariant> variants, int referenceHaplotypeCol) throws IOException {
		// make the column for the intercept
		DoubleMatrix2D intercept = DoubleFactory2D.dense.make(hapdata.rows(), 1);
		intercept.assign(1);

//		DoubleMatrix2D diseaseStatus = new DenseDoubleMatrix2D(hapdata.rows(), 1);

		if (intercept.rows() != finalDiseaseStatus.length) {
			System.out.println("????");
			System.exit(-1);
		}

		double[] y = new double[finalDiseaseStatus.length];
		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			y[i] = finalDiseaseStatus[i][0].getNumber();
		}

		DoubleMatrix2D xprime = DoubleFactory2D.dense.appendColumns(intercept, finalCovariates);
		DenseDoubleAlgebra dda = new DenseDoubleAlgebra();
		int[] colindexes = new int[hapdata.columns() - 1];
		int ctr = 0;
		for (int i = 0; i < hapdata.columns(); i++) {
			if (i != referenceHaplotypeCol) {
				colindexes[ctr] = i;
				ctr++;
			}
		}
		DoubleMatrix2D submatX = dda.subMatrix(hapdata, 0, hapdata.rows() - 1, colindexes);

		DoubleMatrix2D x = DoubleFactory2D.dense.appendColumns(intercept, submatX);
		x = DoubleFactory2D.dense.appendColumns(x, finalCovariates);

		System.out.println(x.rows() + "\t" + xprime.rows() + "\t" + y.length);
		System.out.println(x.columns() + "\t" + xprime.columns() + "\t" + 1);
		System.out.println();

		LogisticRegressionOptimized reg = new LogisticRegressionOptimized();
		LogisticRegressionResult resultCovars = reg.univariate(y, xprime);
		LogisticRegressionResult resultX = reg.univariate(y, x);

		int nrRemaining = hapdata.columns() - 1;
		double devx = resultX.getDeviance();
		double devnull = resultCovars.getDeviance();
		double[] betasmlelr = new double[nrRemaining];
		double[] stderrsmlelr = new double[nrRemaining];
		double[] or = new double[nrRemaining];
		double[] orhi = new double[nrRemaining];
		double[] orlo = new double[nrRemaining];

		for (int i = 0; i < nrRemaining; i++) {
			int mlecol = i + 1;


			betasmlelr[i] = resultX.getBeta()[mlecol];
			stderrsmlelr[i] = resultX.getStderrs()[mlecol];
			double beta = betasmlelr[i];
			double se = resultX.getStderrs()[mlecol];
			double OR = Math.exp(beta);
			double orLow = Math.exp(beta - 1.96 * se);
			double orHigh = Math.exp(beta + 1.96 * se);
			or[i] = OR;
			orhi[i] = orHigh;
			orlo[i] = orLow;

		}

		double deltaDeviance = devnull - devx;
		int df = x.columns() - xprime.columns();
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

		SNPFeature snp = new SNPFeature(Chromosome.ONE, 1, 1);


		result.setSnp(snp);
		result.setN(x.rows());
		snp.setMaf(0d);
		snp.setHwep(0d);

		Double imputationqualityscore = 1d;
		snp.setImputationQualityScore(imputationqualityscore);

		LRTestHaploTestTask lrht = new LRTestHaploTestTask();


		snp.setAlleles(new String[]{"", ""});

		snp.setName(lrht.getHaplotypeComboDescription(availableHaplotypes, variants));
		snp.setMinorAllele("");
		out.writeln(result.toString());
		System.out.println(Math.exp(0.2244776321462565));
		System.out.println(lrht.getHaplotypeComboDescription(availableHaplotypes, variants) + "\t" + Strings.concat(betasmlelr, Strings.semicolon) + "\t" + Strings.concat(or, Strings.semicolon) + "\t" + p + "\t" + -Math.log10(p));
	}

	private Triple<DoubleMatrix2D, ArrayList<BitVector>, BitVector[][]> createHaplotypes(ArrayList<VCFVariant> variants) {

		System.out.println("Generating haplotype pairs from: " + variants.size() + " variants.");
		// get a list of available haplotypes

		HashSet<BitVector> availableHaplotypeHash = new HashSet<BitVector>();
		ArrayList<BitVector> availableHaplotypesList = new ArrayList<>();

		CompletionService<Triple<BitVector[], Integer, Boolean>> jobHandler = new ExecutorCompletionService<Triple<BitVector[], Integer, Boolean>>(exService);
		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			jobHandler.submit(new LRTestHaploTask(i, variants, options.getGenotypeProbThreshold()));
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

				pb.set(returned);
				returned++;
				if (availableHaplotypeHash.size() > 1E5) {
					System.out.println("Can't work with this number of haplotypes..");
					System.exit(-1);
				}


				System.out.print("Individual: " + returned + "\t" + availableHaplotypeHash.size() + " available haplotypes so far.\r");


			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}
		pb.close();
		System.out.println("Individual: " + returned + "\t" + availableHaplotypeHash.size() + " available haplotypes so far.");

		// now we have a collection of all haplotypes available in the data..


		// prune the data
		return pruneMissingSamplesAndReformat(haplotypePairs, availableHaplotypesList, variants);

	}

	private Triple<DoubleMatrix2D, ArrayList<BitVector>, BitVector[][]> pruneMissingSamplesAndReformat(BitVector[][] haplotypePairs, ArrayList<BitVector> availableHaplotypesList, ArrayList<VCFVariant> variants) {

		System.out.println("Variant allele frequencies: ");
		for (int i = 0; i < variants.size(); i++) {
			VCFVariant var = variants.get(i);
			System.out.println(var.getId() + "\t" + var.getAlleleFrequencies()[0] + "\t" + var.getAlleleFrequencies()[1]);
		}

		// determine frequencies
		Triple<ArrayList<BitVector>, BitVector, double[]> haplotypeFrequencyData = haplotypeFrequency(availableHaplotypesList, haplotypePairs, variants);

		if (haplotypeFrequencyData.getLeft().isEmpty()) {
			return null;
		}

		System.out.println();

		double[] haplotypeFreqs = haplotypeFrequencyData.getRight();

		System.out.println("Reference Haplotype: " + lrht.getHaplotypeDesc(haplotypeFrequencyData.getMiddle(), variants));

		// now we have selected the haplotypes to test, prune the individuals that don't have any of those haplotypes (and hets too)
		// make a list of individuals to remove
		ArrayList<Integer> keepTheseIndividuals = new ArrayList<>();
		ArrayList<BitVector> haplotypesAboveThreshold = haplotypeFrequencyData.getLeft();

		HashSet<BitVector> allowedHaplotypes = new HashSet<>();
		allowedHaplotypes.addAll(haplotypesAboveThreshold);
		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			BitVector[] haps = haplotypePairs[i];
			if (haps != null) {
				if (allowedHaplotypes.contains(haps[0]) && allowedHaplotypes.contains(haps[1])) {
					keepTheseIndividuals.add(i);
				} else {
					haplotypePairs[i] = null;
				}
			}
		}

		System.out.println(keepTheseIndividuals.size() + " individuals remain after pruning.");

		if (keepTheseIndividuals.isEmpty()) {
			return null;
		}
		// recalculate haplotype frequencies?
		haplotypeFrequencyData = haplotypeFrequency(availableHaplotypesList, haplotypePairs, variants);
		if (haplotypeFrequencyData.getLeft().isEmpty()) {
			return null;
		}

		System.out.println();
		haplotypeFreqs = haplotypeFrequencyData.getRight();
		System.out.println("Reference Haplotype: " + lrht.getHaplotypeDesc(haplotypeFrequencyData.getMiddle(), variants));

		// prune haplotypes, covariates and disease status
		DiseaseStatus[][] remainingdiseaseStatus = new DiseaseStatus[keepTheseIndividuals.size()][1];
		DoubleMatrix2D remainingCovariates = new DenseDoubleMatrix2D(keepTheseIndividuals.size(), finalCovariates.columns());
		DoubleMatrix2D haplotypeMatrix = new DenseDoubleMatrix2D(keepTheseIndividuals.size(), haplotypesAboveThreshold.size());
		BitVector[][] remainingHaplotypes = new BitVector[keepTheseIndividuals.size()][2];

		HashMap<BitVector, Integer> hapToInt = new HashMap<>();
		for (int i = 0; i < haplotypesAboveThreshold.size(); i++) {
			hapToInt.put(haplotypesAboveThreshold.get(i), i);
		}

		int ictr = 0;
		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			BitVector[] haps = haplotypePairs[i];
			if (haps != null) {
				Integer hapId1 = hapToInt.get(haps[0]);
				Integer hapId2 = hapToInt.get(haps[1]);

				if (hapId1.equals(hapId2)) {
					haplotypeMatrix.set(ictr, hapId1, 2);
				} else {
					haplotypeMatrix.set(ictr, hapId1, 1);
					haplotypeMatrix.set(ictr, hapId2, 1);
				}

				remainingHaplotypes[ictr] = haps;

				remainingdiseaseStatus[ictr][0] = finalDiseaseStatus[i][0];
				for (int c = 0; c < finalCovariates.columns(); c++) {
					remainingCovariates.setQuick(ictr, c, finalCovariates.getQuick(i, c));
				}
				ictr++;
			}
		}

		finalDiseaseStatus = remainingdiseaseStatus;
		finalCovariates = remainingCovariates;


		return new Triple<>(haplotypeMatrix, haplotypesAboveThreshold, remainingHaplotypes);
	}


	public Triple<ArrayList<BitVector>, BitVector, double[]> haplotypeFrequency(ArrayList<BitVector> availableHaplotypesList, BitVector[][] haplotypePairs, ArrayList<VCFVariant> variants) {

		// index haplotypes
		System.out.println();
		System.out.println("Haplotype frequencies: ");
		HashMap<BitVector, Integer> index = new HashMap<BitVector, Integer>();
		for (int b = 0; b < availableHaplotypesList.size(); b++) {
			index.put(availableHaplotypesList.get(b), b);
		}

		// determine their frequencies
		double[] haplotypeFrequencies = new double[availableHaplotypesList.size()];
		for (int i = 0; i < finalDiseaseStatus.length; i++) {
			BitVector[] haps = haplotypePairs[i];
			if (haps != null) {
				Integer index1 = index.get(haps[0]);
				Integer index2 = index.get(haps[1]);
				haplotypeFrequencies[index1]++;
				haplotypeFrequencies[index2]++;
			}
		}


		int totalChromosomes = finalDiseaseStatus.length * 2;
		double frequencythreshold = options.getHaplotypeFrequencyThreshold();
		int nrAboveThreshold = 0;
		Double maxFreq = 0d;
		BitVector referenceHaplotype = null;

		ArrayList<BitVector> haplotypeWithFrequencyAboveThreshold = new ArrayList<>();
		for (int i = 0; i < availableHaplotypesList.size(); i++) {
			double ctr = haplotypeFrequencies[i];

			haplotypeFrequencies[i] /= totalChromosomes;
			BitVector hap = availableHaplotypesList.get(i);
			int[] hapint = new int[variants.size()];
			for (int q = 0; q < variants.size(); q++) {
				if (hap.get(q)) {
					hapint[q] = 1;
				}
			}

			System.out.println(lrht.getHaplotypeDesc(availableHaplotypesList.get(i), variants) + "\t" + haplotypeFrequencies[i] + "\t" + ctr);
			if (haplotypeFrequencies[i] > frequencythreshold) {
				nrAboveThreshold++;
				haplotypeWithFrequencyAboveThreshold.add(availableHaplotypesList.get(i));
				if (haplotypeFrequencies[i] > maxFreq) {
					referenceHaplotype = availableHaplotypesList.get(i);
					maxFreq = haplotypeFrequencies[i];
				}
			}
		}
		System.out.println(nrAboveThreshold + " haplotypes with frequency > " + frequencythreshold);
		return new Triple<>(haplotypeWithFrequencyAboveThreshold, referenceHaplotype, haplotypeFrequencies);
	}

	public class HaplotypeCombiner {
		private final ArrayList<BitVector> inputHaplotypes;
		ArrayList<List<BitVector>> results = new ArrayList<List<BitVector>>();
		private List<BitVector> output = new ArrayList<BitVector>();

		public HaplotypeCombiner(final ArrayList<BitVector> haps) {
			inputHaplotypes = haps;
		}

		public void combine() {
			combine(0);
		}

		private void combine(int start) {
			for (int i = start; i < inputHaplotypes.size(); ++i) {
				output.add(inputHaplotypes.get(i));
				List<BitVector> tmpOutput = new ArrayList<BitVector>();
				tmpOutput.addAll(output);
				results.add(tmpOutput);
				if (i < inputHaplotypes.size()) {
					combine(i + 1);
				}
				output = output.subList(0, output.size() - 1);
			}
		}

		public ArrayList<List<BitVector>> getResults() {
			return results;
		}
	}

}
