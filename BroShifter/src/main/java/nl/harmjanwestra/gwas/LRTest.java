package nl.harmjanwestra.gwas;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import nl.harmjanwestra.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.gwas.tasks.LRTestExhaustiveTask;
import nl.harmjanwestra.gwas.tasks.LRTestTask;
import nl.harmjanwestra.gwas.tasks.LRTestVariantQCTask;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationFilePairwise;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.AssociationResultPairwise;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.individuals.Individual;
import nl.harmjanwestra.utilities.math.LogisticRegressionOptimized;
import nl.harmjanwestra.utilities.math.LogisticRegressionResult;
import nl.harmjanwestra.utilities.plink.PlinkFamFile;
import nl.harmjanwestra.utilities.vcf.SampleAnnotation;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.vcf.VCFVariantComparator;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.Console;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;

//import java.util.concurrent.Executors;

/**
 * Created by hwestra on 11/7/15.
 */
public class LRTest {

	protected ExecutorService exService;
	protected SampleAnnotation sampleAnnotation;
	LRTestOptions options;
	private int submitted;
	private int returned;
	private double highestLog10p;
	private HashMap<String, Integer> sampleToIntGenotypes;
	private HashSet<String> snpLimit;
	private ArrayList<Feature> bedRegions = null;
	private boolean[] genotypeSamplesWithCovariatesAndDiseaseStatus;
	private ProgressBar progressBar;

	public LRTest(LRTestOptions options) throws IOException {
		this.options = options;


		if (!initialize()) {
			System.err.println("Something went wrong during initialization.. Please check the input.");
			System.exit(-1);
		}

		switch (options.getAnalysisType()) {
			case CONDITIONAL:
				System.out.println("Will perform conditional logistic regression");
				System.out.println("Setting up threadpool with: " + options.getNrThreads() + " threads..");
				exService = Executors.newFixedThreadPool(options.getNrThreads());

				testConditional();
				exService.shutdown();
				break;
			case EXHAUSTIVE:
				System.out.println("Will perform exhaustive pairwise logistic regression");
				System.out.println("Setting up threadpool with: " + options.getNrThreads() + " threads..");
				exService = Executors.newFixedThreadPool(options.getNrThreads());
				testExhaustivePairwise();
				exService.shutdown();
				break;
			case HAPLOTYPE:
				System.out.println("Will perform haplotype logistic regression");
				break;
			case NORMAL:
				System.out.println("Will perform normal logistic regression");
				System.out.println("Setting up threadpool with: " + options.getNrThreads() + " threads..");
				exService = Executors.newFixedThreadPool(options.getNrThreads());
				testNormal();
				exService.shutdown();
				break;
		}

		// System.exit(0);
		System.out.println("Done.");
	}

	public static void main(String[] args) {


		String[] args4 = new String[]{
				"--gwas",
				"--conditional",
				"-c", "/Data/tmp/2016-06-24/2016-03-11-T1D-covarmerged.txtmergedCovariates-withPseudos.txt",
				"-d", "/Data/tmp/2016-06-24/2016-03-11-T1D-diseaseStatusWithPseudos.txt",
				"-f", "/Data/tmp/2016-06-29/T1D-recode-maf0005-ICRegions-samplenamefix-pseudo.vcf.gz-filtered-merged.fam",
//				"-f", "/Data/tmp/2016-06-24/T1D-recode-maf0005-ICRegions-samplenamefix-pseudo.vcf.gz.fam",
//				"-f", "/Data/tmp/2016-06-26/Allele+AA.NomissingNonRare0.0005.Include.fam", // xinli fam path
				"-i", "/Data/tmp/2016-06-24/T1DVCF.vcf",
				"-e", "/Data/tmp/2016-06-24/T1D-recode-regionsfiltered-allelesfiltered-samplenamefix-pseudo.vcf.gz-parents.txt",
				"-r", "/Data/tmp/2016-06-24/il2ra.bed",
				"--snplimit", "/Data/tmp/2016-06-24/limit.txt",
				"-t", "4",
				"-q", "0.3",
				"--maxiter", "5",
				"--splitmultiallelic",
				"--limittosamplesinfam",
				"-o", "/Data/tmp/2016-06-29/testoutgwas-hjfixed.txt"
		};
		LRTestOptions options0 = new LRTestOptions(args4);

//		try {
//			new LRTest(options0);
//
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

		args4 = new String[]{
				"--gwas",
				"--haplotype",
				"-c", "/Data/tmp/2016-06-24/2016-03-11-T1D-covarmerged.txtmergedCovariates-withPseudos.txt",
				"-d", "/Data/tmp/2016-06-24/2016-03-11-T1D-diseaseStatusWithPseudos.txt",
				"-f", "/Data/tmp/2016-06-29/T1D-recode-maf0005-ICRegions-samplenamefix-pseudo.vcf.gz-filtered-merged.fam",
				"-i", "/Data/tmp/2016-06-24/T1DVCF.vcf",
				"-e", "/Data/tmp/2016-06-24/T1D-recode-regionsfiltered-allelesfiltered-samplenamefix-pseudo.vcf.gz-parents.txt",
				"-r", "/Data/tmp/2016-06-24/il2ra.bed",
				"--snplimit", "/Data/tmp/2016-06-24/limit.txt",
				"-t", "4",
				"-q", "0.3",
				"--splitmultiallelic",
				"--limittosamplesinfam",
				"-o", "/Data/tmp/2016-06-24/testoutnormal.txt"
		};

		LRTestOptions options = new LRTestOptions(args4);

		try {
			new LRTestHaplotype(options);

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void waitForEnter(String message, Object... args) {
		Console c = System.console();
		if (c != null) {
			// printf-like arguments
			if (message != null)
				c.format(message, args);
			c.format("\nPress ENTER to proceed.\n");
			c.readLine();
		}
	}


	public void testNormal() throws IOException {

		ArrayList<VCFVariant> variants = readVariants(options.getOutputdir() + "variantlog.txt");


		// boot up threadpool
		String condition = options.getConditional();
		HashMap<Feature, Integer> conditionsPerIter = null;
		if (condition != null) {
			conditionsPerIter = new HashMap<Feature, Integer>();
			String[] conditionalSNPs = condition.split(",");
			for (int c = 0; c < conditionalSNPs.length; c++) {
				String[] snpelems = conditionalSNPs[c].split("-");
				Feature f = new Feature();
				f.setChromosome(Chromosome.parseChr(snpelems[0]));
				f.setStart(Integer.parseInt(snpelems[1]));
				f.setStop(Integer.parseInt(snpelems[1]));
				f.setName(snpelems[2]);
				conditionsPerIter.put(f, c + 1);
				System.out.println("Will condition on " + f.toString() + "-" + f.getName() + " in iteration " + (c + 1));
			}
			int maxNrIter = conditionsPerIter.size() + 1;
			options.setMaxIter(maxNrIter);
			System.out.println("Setting max iter to: " + maxNrIter);
		}


		// keep a list of genotypes to condition on
		int iter = 0;
		ArrayList<Pair<VCFVariant, Triple<int[], boolean[], Integer>>> conditional = new ArrayList<>();

		ArrayList<String> conditionalVariantIds = new ArrayList<String>();
		AssociationFile associationFile = new AssociationFile();
		String header = associationFile.getHeader();

		submitted = 0;
		returned = 0;

		int alleleOffsetGenotypes = 0;
		highestLog10p = 0;
		LogisticRegressionResult nullmodelresult = null;
		Integer nrNullColumns = null;
		if (options.assumeNoMissingData) {
			Pair<LogisticRegressionResult, Integer> nullmodelresultpair = getNullModel(variants.get(0), conditional, 1, 1 + (variants.get(0).getNrAlleles() - 1));
			nrNullColumns = nullmodelresultpair.getRight();
			nullmodelresult = nullmodelresultpair.getLeft();
		}


		CompletionService<Triple<String, AssociationResult, VCFVariant>> jobHandler = new ExecutorCompletionService<Triple<String, AssociationResult, VCFVariant>>(exService);
		for (int i = 0; i < variants.size(); i++) {
			VCFVariant variant = variants.get(i);
			// throw into a thread
			// TODO: conditional on dosages is separate from conditional on genotypes.. ?
			LRTestTask task = new LRTestTask(variant,
					iter,
					conditional,
					alleleOffsetGenotypes,
					sampleAnnotation,
					options);

			if (nullmodelresult != null) {
				task.setResultNullmodel(nullmodelresult, nrNullColumns);
			}

			jobHandler.submit(task);

			submitted++;
		}

		TextFile pvalout = new TextFile(options.getOutputdir() + "gwas-" + iter + ".txt", TextFile.W);
		System.out.println("Output will be written here: " + options.getOutputdir() + "gwas-" + iter + ".txt");
		pvalout.writeln(header);
		progressBar = new ProgressBar(submitted);
		clearQueue(null, pvalout, iter, variants, jobHandler, null);
		progressBar.close();
		pvalout.close();

	}

	protected boolean initialize() throws IOException {
		System.out.println("Assoc: " + options.getVcf());
		System.out.println("Covar: " + options.getCovariateFile());
		System.out.println("Disease: " + options.getDiseaseStatusFile());
		System.out.println("Out: " + options.getOutputdir());

		// index covariate samples to vcf samples
		HashSet<String> excludeTheseSamples = new HashSet<String>();
		if (options.getSamplesToExclude() != null) {
			TextFile excl = new TextFile(options.getSamplesToExclude(), TextFile.R);
			excludeTheseSamples = (HashSet<String>) excl.readAsSet(0, TextFile.tab);
			excl.close();
			System.out.println("Loaded: " + excludeTheseSamples.size() + " samples to exclude from " + options.getSamplesToExclude());
		}
		HashMap<String, DiseaseStatus[]> diseaseStatus = new HashMap<String, DiseaseStatus[]>();

		// load disease status
		TextFile tf = new TextFile(options.getDiseaseStatusFile(), TextFile.R);
		String[] elems = tf.readLineElems(Strings.tab);
		while (elems != null) {
			String sample = elems[0];
			if (excludeTheseSamples == null || !excludeTheseSamples.contains(sample)) {
				DiseaseStatus[] sampleDiseaseStatus;
				if (elems.length > 2) {
					sampleDiseaseStatus = new DiseaseStatus[elems.length - 1];
					for (int i = 1; i < elems.length; i++) {
						sampleDiseaseStatus[i - 1] = DiseaseStatus.parseStatus(elems[i]);
					}
				} else {
					String statusStr = elems[1];
					sampleDiseaseStatus = new DiseaseStatus[1];
					sampleDiseaseStatus[0] = DiseaseStatus.parseStatus(statusStr);
				}
				diseaseStatus.put(sample, sampleDiseaseStatus);
			}
			elems = tf.readLineElems(Strings.tab);
		}
		tf.close();

		System.out.println(diseaseStatus.size() + " disease status samples loaded");

		// load samples from vcf path
		VCFGenotypeData data = new VCFGenotypeData(options.getVcf());
		ArrayList<String> vcfSamples = data.getSamples();
		HashSet<String> vcfSamplesWithDiseaseStatus = new HashSet<String>();
		for (String sample : vcfSamples) {
			if (diseaseStatus.containsKey(sample)) {
				if (excludeTheseSamples == null || !excludeTheseSamples.contains(sample)) {
					vcfSamplesWithDiseaseStatus.add(sample);
				}
			}
		}

		// remove related samples based on fam path if any
		if (options.getFamfile() != null) {
			if (excludeTheseSamples == null) {
				excludeTheseSamples = new HashSet<String>();
			}

			if (options.isTestallsamplesinfam()) {
				excludeSamplesNotInFAMFile(excludeTheseSamples, vcfSamplesWithDiseaseStatus);
			} else {
				determineSamplesToRemoveFromFam(excludeTheseSamples, vcfSamplesWithDiseaseStatus);
			}

			HashSet<String> tmpVCFSamples = new HashSet<String>();
			for (String s : vcfSamplesWithDiseaseStatus) {
				if (!excludeTheseSamples.contains(s)) {
					tmpVCFSamples.add(s);
				}
			}
			vcfSamplesWithDiseaseStatus = tmpVCFSamples;
			System.out.println(vcfSamplesWithDiseaseStatus.size() + " samples left after filtering relateds");
		}

		// assume rows are samples and columns are covariates
		// this loads only samples that also have a disease status loaded..
		HashSet<String> covariatesToInclude = null;
		if (options.getCovariatesToInclude() != null) {
			// load hashet
			TextFile excl = new TextFile(options.getCovariatesToInclude(), TextFile.R);
			covariatesToInclude = (HashSet<String>) excl.readAsSet(0, TextFile.tab);
			excl.close();
			System.out.println(covariatesToInclude + " covariates selected from: " + options.getCovariatesToInclude());
		}

		// load covariates
		DoubleMatrixDataset<String, String> covariates = DoubleMatrixDataset.loadSubsetOfTextDoubleData(options.getCovariateFile(), "\t", vcfSamplesWithDiseaseStatus, covariatesToInclude);
		System.out.println("Covariate matrix: " + covariates.rows() + " samples " + covariates.columns() + " covariates");

		LinkedHashMap<String, Integer> covariateSampleHash = covariates.getHashRows();
		sampleToIntGenotypes = new HashMap<String, Integer>();
		int ctr = 0;
		genotypeSamplesWithCovariatesAndDiseaseStatus = new boolean[vcfSamples.size()];
		int vcfSampleCtr = 0;

		// index covariate samples
		HashSet<String> alternatecovariateSamples = new HashSet<String>();
		Set<String> keys = covariateSampleHash.keySet();
		HashMap<String, String> altToSample = new HashMap<String, String>();
		for (String sample : keys) {
			alternatecovariateSamples.add(sample + "_" + sample); // sometimes covariates have a weird samplename
			altToSample.put(sample, sample + "_" + sample);
		}

		// index the final list of samples
		ArrayList<String> samplesIntersect = new ArrayList<>();
		for (String sample : vcfSamples) {
			if (covariateSampleHash.containsKey(sample) || alternatecovariateSamples.contains(sample)) {
				sampleToIntGenotypes.put(sample, ctr);
				genotypeSamplesWithCovariatesAndDiseaseStatus[vcfSampleCtr] = true;
				samplesIntersect.add(sample);
				ctr++;
			}
			vcfSampleCtr++;
		}

		// load bed regions, if any

		if (options.getBedfile() != null) {
			BedFileReader reader = new BedFileReader();
			bedRegions = reader.readAsList(options.getBedfile());
			System.out.println(bedRegions.size() + " regions loaded from: " + options.getBedfile());
			Collections.sort(bedRegions, new FeatureComparator(false));
		}

		// load list of variants to limit upon if any
		if (options.getSnpLimitFile() != null) {
			TextFile excl = new TextFile(options.getSnpLimitFile(), TextFile.R);
			snpLimit = (HashSet<String>) excl.readAsSet(0, TextFile.tab);
			excl.close();
		}


		System.out.println(sampleToIntGenotypes.size() + " samples with disease status, covariates and genotypes");
		if (sampleToIntGenotypes.size() == 0) {
			System.out.println("Problem with matching samples...");
			return false;
		} else {
			DiseaseStatus[][] finalDiseaseStatus = new DiseaseStatus[sampleToIntGenotypes.size()][0];
			DoubleMatrix2D finalCovariates = new DenseDoubleMatrix2D(sampleToIntGenotypes.size(), covariates.columns());

			TextFile sampleListOut = new TextFile(options.getOutputdir() + "samplelist.txt", TextFile.W);
			System.out.println(options.getOutputdir() + "samplelist.txt");
			for (int i = 0; i < finalCovariates.rows(); i++) {
				sampleListOut.writeln(samplesIntersect.get(i));
			}
			sampleListOut.close();

			System.out.println("Final covariate array size: " + finalCovariates.rows());
			System.out.println("Final disease status array size: " + finalDiseaseStatus.length);

			// reorder the covariates to match the genotyped samples
			ArrayList<String> samplesFromCovariates = covariates.getRowObjects();
			int minstatus = Integer.MAX_VALUE;
			int maxstatus = -Integer.MAX_VALUE;
			int nrCases = 0;
			int nrControls = 0;
			int nrUnknown = 0;
			int nrTotal = 0;

			// make final list of covariates
			// remap disease status to 0 and 1
			for (int sid = 0; sid < samplesFromCovariates.size(); sid++) {
				String sample = samplesFromCovariates.get(sid);

				Integer id = sampleToIntGenotypes.get(sample);
				if (id == null) {
					id = sampleToIntGenotypes.get(altToSample.get(sample));
				}
				if (id != null) {


					finalDiseaseStatus[id][0] = diseaseStatus.get(sample)[0];
					if (finalDiseaseStatus[id][0].equals(DiseaseStatus.CONTROL)) {
						nrControls++;
					} else if (finalDiseaseStatus[id][0].equals(DiseaseStatus.CASE)) {
						nrCases++;
					} else {
						nrUnknown++;
					}
					for (int col = 0; col < covariates.columns(); col++) {
						finalCovariates.setQuick(id, col, covariates.getElement(sid, col));
					}
					nrTotal++;
				}
			}

//			for (int id = 0; id < finalDiseaseStatus.length; id++) {
//				System.out.println(id + "\t" + finalDiseaseStatus[id] + "\t" + finalDiseaseStatus[id].getNumber());
//			}
//			System.exit(-1);
			System.out.println(nrCases + " cases ");
			System.out.println(nrControls + " controls ");
			System.out.println(nrUnknown + " unknown ");
			System.out.println(nrTotal + " with covariates");

			sampleAnnotation = new SampleAnnotation();
//			sampleAnnotation.setIndividualGender(individualGender);
			sampleAnnotation.setCovariates(finalCovariates);
			sampleAnnotation.setSampleDiseaseStatus(finalDiseaseStatus);

			return true;
		}
	}

	private void excludeSamplesNotInFAMFile(HashSet<String> excludeTheseSamples, HashSet<String> vcfSamplesWithDiseaseStatus) throws IOException {

		PlinkFamFile famfile = new PlinkFamFile(options.getFamfile());
		ArrayList<Individual> samples = famfile.getSamples();
		ArrayList<String> samplesHash = new ArrayList<>();
		for (int i = 0; i < samples.size(); i++) {
			samplesHash.add(samples.get(i).getName());
		}

		for (String sample : vcfSamplesWithDiseaseStatus) {
			if (!samplesHash.contains(sample)) {
				excludeTheseSamples.add(sample);
			}
		}

		for (String sample : samplesHash) {
			if (!vcfSamplesWithDiseaseStatus.contains(sample)) {
				System.out.println(sample + "\t not present in my data?");
			}
		}
	}

	public void testConditional() throws IOException {


		if (options.getConditional() != null) {
			System.out.println("Using " + options.getConditional() + " for conditional analysis..");
		}

		AssociationFile associationFile = new AssociationFile();
		String header = associationFile.getHeader();

		// boot up threadpool
		ArrayList<HashMap<Feature, String>> variantsToConditionOn = null;
		if (options.getConditional() != null) {
			System.out.println("Parsing: " + options.getConditional());
			variantsToConditionOn = new ArrayList<>();

			// you can provide one file for each iteration.
			String file = options.getConditional();

			// format: Chr2_162960873-163361685        0       2_163110536_rs2111485   8.910715565551376

			// determine number of iterations in file
			TextFile tf = new TextFile(file, TextFile.R);
			String[] head = tf.readLineElems(TextFile.tab);
			if (head.length > 3) {
				if (!head[0].equals("Region") || !head[1].equals("Iter") || !head[2].equals("Variant")) {
					System.err.println("Error: " + file + " does not have expected format (Region\tIter\tVariant)");
					System.exit(-1);
				}
			}
			String[] elems = tf.readLineElems(TextFile.tab);
			int maxiter = 0;

			while (elems != null) {
				Integer iter = Integer.parseInt(elems[1]);
				if (iter > maxiter) {
					maxiter = iter;

				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

			// now load the data
			for (int q = 0; q < maxiter + 1; q++) {
				variantsToConditionOn.add(new HashMap<>());
			}

			System.out.println("Max iter in file: " + maxiter);

			tf.open();
			tf.readLine();
			elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length >= 3) {
					String region = elems[0];

					String[] regionElems = region.split("_");
					Chromosome chr = Chromosome.parseChr(regionElems[0]);
					String[] posElems = regionElems[1].split("-");
					Integer start = Integer.parseInt(posElems[0]);
					Integer stop = Integer.parseInt(posElems[1]);
					Feature reg = new Feature(chr, start, stop);
					Integer iter = Integer.parseInt(elems[1]);

					String varStr = elems[2];
					String[] varStrElems = varStr.split("_");
					varStr = Chromosome.parseChr(varStrElems[0]).toString() + "_" + varStrElems[1] + "_" + varStrElems[2];

					HashMap<Feature, String> toAdd = variantsToConditionOn.get(iter);
					toAdd.put(reg, varStr);
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

		}

		ArrayList<VCFVariant> allVariants = readVariants(options.getOutputdir() + "variantlog.txt");


		CompletionService<Triple<String, AssociationResult, VCFVariant>> jobHandler = new ExecutorCompletionService<Triple<String, AssociationResult, VCFVariant>>(exService);


		ArrayList<AssociationResult> assocResults = new ArrayList<>();


		System.out.println("Iteration " + 0 + " starting. Model: y ~ SNP + covar.");
		System.out.println("Output will be written here: " + options.getOutputdir() + "gwas-" + 0 + ".txt");
		TextFile pvalout = new TextFile(options.getOutputdir() + "gwas-" + 0 + ".txt", TextFile.W);
		TextFile logout = new TextFile(options.getOutputdir() + "gwas-" + 0 + "-log.txt.gz", TextFile.W);
		pvalout.writeln(header);

		int variantCtr = 0;
		submitted = 0;
		returned = 0;
		highestLog10p = 0;


		for (VCFVariant variant : allVariants) {
			ArrayList<Pair<VCFVariant, Triple<int[], boolean[], Integer>>> conditional = new ArrayList<>();

			LRTestTask task = new LRTestTask(variant,
					0,
					conditional,
					0,
					sampleAnnotation,
					options);

			jobHandler.submit(task);

			submitted++;
			if (submitted % 1000 == 0) {
				clearQueue(logout, pvalout, 0, null, jobHandler, assocResults);
			}
		}

		// clean up the queue
		clearQueue(logout, pvalout, 0, null, jobHandler, assocResults);
		pvalout.close();
		logout.close();

		System.out.println("Initial mapping done. Now performing conditional analysis");
		// run conditional per region..
		// get iter0 results
		// get best associated variant in region
		// testNormal conditional on best variant


		HashSet<Feature> visitedRegions = new HashSet<Feature>();
		for (VCFVariant v : allVariants) {
			for (int f = 0; f < bedRegions.size(); f++) {
				if (v.asFeature().overlaps(bedRegions.get(f))) {
					visitedRegions.add(bedRegions.get(f));
				}
			}
		}

		System.out.println(visitedRegions.size() + " regions are present in association output.");
		if (visitedRegions.isEmpty()) {
			System.out.println("No more work to do");

		} else {
			ArrayList<Feature> remainingRegions = new ArrayList<>();

			int nrMaxIter = 0;
			if (variantsToConditionOn != null) {

				HashSet<Feature> regionsToCondition = new HashSet<Feature>();
				for (int q = 0; q < variantsToConditionOn.size(); q++) {
					regionsToCondition.addAll(variantsToConditionOn.get(q).keySet());
				}
				System.out.println(regionsToCondition.size() + " total regions to run conditional analysis.");

				HashSet<Feature> regionsOnChr = new HashSet<Feature>();
				for (VCFVariant variant : allVariants) {
					for (Feature f : regionsToCondition) {
						if (variant.asFeature().overlaps(f)) {
							regionsOnChr.add(f);
						}
					}
				}
				remainingRegions.addAll(regionsOnChr);
				System.out.println(remainingRegions.size() + " regions remaining after filtering for variants");

				nrMaxIter = variantsToConditionOn.size();
			} else {
				remainingRegions.addAll(visitedRegions);
				nrMaxIter = options.getMaxIter();
			}

			System.out.println("Total number of iterations to run: " + nrMaxIter);


			ArrayList<ArrayList<VCFVariant>> conditionalVariants = new ArrayList<>();
			for (int i = 0; i < remainingRegions.size(); i++) {
				conditionalVariants.add(new ArrayList<>());
			}

			TextFile modelsout = new TextFile(options.getOutputdir() + "gwas-conditional-models.txt", TextFile.W);
			modelsout.writeln("Region\tIter\tModel");
			TextFile topvariantsout = new TextFile(options.getOutputdir() + "gwas-topvariants.txt", TextFile.W);
			topvariantsout.writeln("Region\tIter\tVariant\tPval");

			for (int iteration = 1; iteration < nrMaxIter; iteration++) {

				System.out.println("Output will be written here: " + options.getOutputdir() + "gwas-" + iteration + ".txt");
				pvalout = new TextFile(options.getOutputdir() + "gwas-" + iteration + ".txt", TextFile.W);

				pvalout.writeln(header);

				// start new result array
				ArrayList<AssociationResult> assocResultsIter = new ArrayList<>();
				submitted = 0;
				for (int regionId = 0; regionId < remainingRegions.size(); regionId++) {
					// get variants in region
					ArrayList<VCFVariant> variantsInRegion = filterVariantsByRegion(allVariants, remainingRegions.get(regionId));

					VCFVariant bestVariantLastIter = null;
					Pair<VCFVariant, AssociationResult> bestAssocLastIter = null;
					if (variantsToConditionOn != null) {

						HashMap<Feature, String> regionVars = variantsToConditionOn.get(iteration - 1);
						String variantToConditionOn = regionVars.get(remainingRegions.get(regionId));
						System.out.println("Looking for variant: " + variantToConditionOn + " in region " + remainingRegions.get(regionId).toString());
						for (VCFVariant var : variantsInRegion) {
							String varstr = Chromosome.parseChr(var.getChr()).toString() + "_" + var.getPos() + "_" + var.getId();
							if (varstr.equals(variantToConditionOn)) {
								bestVariantLastIter = var;
							}
						}

						bestAssocLastIter = getBestAssocForRegion(assocResults, remainingRegions.get(regionId), variantsInRegion, variantToConditionOn);
						if (bestVariantLastIter == null) {
							System.err.println("Error: could not find variant: " + variantToConditionOn);
							System.exit(-1);
						}
					} else {
						System.exit(-1);
						bestAssocLastIter = getBestAssocForRegion(assocResults, remainingRegions.get(regionId), variantsInRegion, null);
						bestVariantLastIter = bestAssocLastIter.getLeft();
					}
					topvariantsout.writeln(remainingRegions.get(regionId).toString() + "\t" + (iteration - 1) + "\t" + bestVariantLastIter.getChr() + "_" + bestVariantLastIter.getPos() + "_" + bestVariantLastIter.getId() + "\t" + bestAssocLastIter.getRight().getLog10Pval());

					// get conditional variants for this region
					ArrayList<VCFVariant> conditionalVariantsForRegion = conditionalVariants.get(regionId);
					if (conditionalVariantsForRegion == null) {
						conditionalVariantsForRegion = new ArrayList<>();
					}
					conditionalVariantsForRegion.add(bestVariantLastIter);

					// now put the data in the correct data structure.. recode and whatever
					ArrayList<Pair<VCFVariant, Triple<int[], boolean[], Integer>>> conditional = new ArrayList<>();


//						String model = "";
					LRTestVariantQCTask lqr = new LRTestVariantQCTask();
					int alleleOffsetGenotypes = 0;
					ArrayList<String> modelvariants = new ArrayList<>();
					for (int c = 0; c < conditionalVariantsForRegion.size(); c++) {
						VCFVariant variant = conditionalVariantsForRegion.get(c);
						modelvariants.add(variant.toString());
//							model += " + " + variant.getId();
						conditional.add(new Pair<>(variant,
								lqr.determineMissingGenotypes(variant.getGenotypeAllelesAsMatrix2D(),
										sampleAnnotation.getCovariates().rows())));
						alleleOffsetGenotypes += variant.getAlleles().length - 1;
					}
					modelsout.writeln(remainingRegions.get(regionId).toString() + "\t" + iteration + "\t" + Strings.concat(modelvariants, Strings.semicolon));

					int alleleOffsetDosages = 0;

					for (int v = 0; v < variantsInRegion.size(); v++) {
						// throw into a threads
						VCFVariant variant = variantsInRegion.get(v);
						LRTestTask task = new LRTestTask(variant,
								iteration,
								conditional,
								alleleOffsetGenotypes,
								sampleAnnotation,
								options);

						jobHandler.submit(task);
						submitted++;


						if (submitted % 1000 == 0) {
							clearQueue(null, pvalout, iteration, null, jobHandler, assocResultsIter);
						}
					}
				}

				clearQueue(null, pvalout, iteration, null, jobHandler, assocResultsIter);
				assocResults = assocResultsIter;
				pvalout.close();
				System.out.println("Done");
			}
			for (int regionId = 0; regionId < remainingRegions.size(); regionId++) {
				ArrayList<VCFVariant> variantsInRegion = filterVariantsByRegion(allVariants, remainingRegions.get(regionId));
				// get variants in region
				Pair<VCFVariant, AssociationResult> bestAssocLastIter = null;


				if (variantsToConditionOn != null) {
					VCFVariant bestVariantLastIter = null;
					int iteration = nrMaxIter;
					HashMap<Feature, String> regionVars = variantsToConditionOn.get(iteration - 1);
					String variantToConditionOn = regionVars.get(remainingRegions.get(regionId));
					System.out.println("Looking for variant: " + variantToConditionOn + " in region " + remainingRegions.get(regionId).toString());
					for (VCFVariant var : variantsInRegion) {
						String varstr = Chromosome.parseChr(var.getChr()).toString() + "_" + var.getPos() + "_" + var.getId();
						if (varstr.equals(variantToConditionOn)) {
							bestVariantLastIter = var;
						}
					}

//					bestAssocLastIter = getBestAssocForRegion(assocResults, remainingRegions.get(regionId), variantsInRegion, variantToConditionOn);
					if (bestVariantLastIter == null) {
						System.err.println("Error: could not find variant: " + variantToConditionOn);
						System.exit(-1);
					}

					String variant = variantsToConditionOn.get(nrMaxIter - 1).get(remainingRegions.get(regionId));
					bestAssocLastIter = getBestAssocForRegion(assocResults, remainingRegions.get(regionId), variantsInRegion, variant);
					if (bestVariantLastIter == null) {
						System.err.println("Error: could not find variant: " + variantToConditionOn);
						System.exit(-1);
					}
				} else {
					bestAssocLastIter = getBestAssocForRegion(assocResults, remainingRegions.get(regionId), variantsInRegion, null);
				}
				if (bestAssocLastIter != null) {
					VCFVariant bestVariantLastIter = bestAssocLastIter.getLeft();
					topvariantsout.writeln(remainingRegions.get(regionId).toString() + "\t" + (options.getMaxIter() - 1) + "\t" + bestVariantLastIter.getChr() + "_" + bestVariantLastIter.getPos() + "_" + bestVariantLastIter.getId() + "\t" + bestAssocLastIter.getRight().getLog10Pval());
				}
			}
			modelsout.close();
			topvariantsout.close();
		}

	}

	public ArrayList<VCFVariant> readVariants(String logfile) throws IOException {

		TextFile log = new TextFile(logfile, TextFile.W);
		System.out.println("Log will be written here: " + log.getFileName());
		log.writeln("Chr\tPos\tId\tMAF\tHWEP\tImpQual\tOverlap\tPassQC");

		int maxQueueSize = options.getNrThreads() * 10;
		CompletionService<Pair<VCFVariant, String>> jobHandler = new ExecutorCompletionService<Pair<VCFVariant, String>>(exService);

		TextFile vcfIn = new TextFile(options.getVcf(), TextFile.R);
		String vcfLn = vcfIn.readLine();
		while (vcfLn != null && vcfLn.startsWith("#")) {
			vcfLn = vcfIn.readLine();
		}
		int linessubmitted = 0;
		System.out.println("Parsing: " + options.getVcf());


		int totalsubmitted = 0;
		ArrayList<VCFVariant> variants = new ArrayList<>();
		while (vcfLn != null) {
			if (linessubmitted == maxQueueSize) {
				// get some of the output
				int returned = 0;
				while (returned < linessubmitted) {
					try {
						Future<Pair<VCFVariant, String>> future = jobHandler.take();
						if (future != null) {
							Pair<VCFVariant, String> result = future.get();

							if (log != null) {
								String logStr = result.getRight();
								log.writeln(logStr);
							}
							if (result.getLeft() != null) {
								VCFVariant variant = result.getLeft();
								int nrAlleles = variant.getNrAlleles();
								if (options.splitMultiAllelicIntoMultipleVariants && nrAlleles > 2) {
									for (int a = 1; a < nrAlleles; a++) {
										variants.add(variant.variantFromAllele(a));
									}
								} else {
									variants.add(variant);
								}
							}
						}
						returned++;

					} catch (InterruptedException e) {
						e.printStackTrace();
					} catch (ExecutionException e) {
						e.printStackTrace();
					}

				}
				linessubmitted = 0;

				System.out.print(totalsubmitted + " submitted. " + variants.size() + " variants in memory " + ManagementFactory.getThreadMXBean().getThreadCount() + " threads\r");
			}


			LRTestVariantQCTask task = new LRTestVariantQCTask(vcfLn,
					bedRegions,
					options,
					genotypeSamplesWithCovariatesAndDiseaseStatus,
					sampleAnnotation,
					snpLimit);
			jobHandler.submit(task);
			linessubmitted++;
			totalsubmitted++;
			vcfLn = vcfIn.readLine();
		}
		System.out.println();
		System.out.println(totalsubmitted + " total lines read.");
		vcfIn.close();

		int returned = 0;
		while (returned < linessubmitted) {
			try {
				Future<Pair<VCFVariant, String>> future = jobHandler.take();
				if (future != null) {
					Pair<VCFVariant, String> result = future.get();

					if (log != null) {
						String logStr = result.getRight();
						log.writeln(logStr);
					}
					if (result.getLeft() != null) {
						VCFVariant variant = result.getLeft();
						int nrAlleles = variant.getNrAlleles();
						if (options.splitMultiAllelicIntoMultipleVariants && nrAlleles > 2) {
							System.out.println("Splitting alleles");

							for (int a = 1; a < nrAlleles; a++) {
								variants.add(variant.variantFromAllele(a));
							}
						} else {
							variants.add(variant);
						}
					}
				}
				returned++;

			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}
		System.out.println(variants.size() + "  variants finally included.");

		log.close();

		// why not sort :)
		Collections.sort(variants, new VCFVariantComparator());
		return variants;

	}

	public void testExhaustivePairwise() throws IOException {
		System.out.println("Will perform exhaustive pairwise analysis");

		// get a list of variants for the region
		HashMap<Feature, ArrayList<VCFVariant>> variantsInRegions = new HashMap<Feature, ArrayList<VCFVariant>>();
		HashSet<Feature> availableRegions = new HashSet<Feature>();


		ArrayList<VCFVariant> allVariants = readVariants(options.getOutputdir() + "variantlog.txt");

		if (allVariants.isEmpty()) {
			System.out.println("Sorry. No work.");
		} else {
			for (int f = 0; f < bedRegions.size(); f++) {
				Feature region = bedRegions.get(f);
				for (VCFVariant variant : allVariants) {
					if (variant.asFeature().overlaps(region)) {
						ArrayList<VCFVariant> variants = variantsInRegions.get(region);
						if (variants == null) {
							variants = new ArrayList<>();
						}
						variants.add(variant);
						variantsInRegions.put(region, variants);
						availableRegions.add(region);
					}
				}
			}

			// combinatorial madness!
			if (availableRegions.isEmpty() || variantsInRegions.isEmpty()) {
				// print some fancy error message;
				System.out.println("Sorry. No work.");
			} else {
				System.out.println(availableRegions.size() + " regions to query.");
				CompletionService<AssociationResultPairwise> jobHandler = new ExecutorCompletionService<AssociationResultPairwise>(exService);

				LogisticRegressionResult nullmodelresult = null;
				Integer nrNullColumns = null;
				if (options.assumeNoMissingData) {
					System.out.println("Generating null model.");
					Feature[] regionsStr = availableRegions.toArray(new Feature[0]);
					ArrayList<VCFVariant> variants = variantsInRegions.get(regionsStr[0]);
					Pair<LogisticRegressionResult, Integer> nullmodelresultpair = getNullModel(variants.get(0), null, 1, 1 + (variants.get(0).getNrAlleles() - 1));
					nrNullColumns = nullmodelresultpair.getRight();
					nullmodelresult = nullmodelresultpair.getLeft();
				}

				System.out.println("Submitting jobs");
				int submitted = 0;
				for (Feature region : availableRegions) {
					ArrayList<VCFVariant> variants = variantsInRegions.get(region);
					for (int i = 0; i < variants.size(); i++) {
						for (int j = i + 1; j < variants.size(); j++) {
							// submit job to queue;
							LRTestExhaustiveTask lrtet = new LRTestExhaustiveTask(variants,
									i,
									j,
									genotypeSamplesWithCovariatesAndDiseaseStatus,
									sampleAnnotation,
									options);
							lrtet.setResultNullmodel(nullmodelresult, nrNullColumns);
							jobHandler.submit(lrtet);
							submitted++;

						}
					}
				}

				System.out.println("Submitted a total of " + submitted + " jobs");
				ProgressBar pb = new ProgressBar(submitted);
				int returned = 0;
				TextFile pvalOut = new TextFile(options.getOutputdir() + "pairwise.txt.gz", TextFile.W);
				AssociationFilePairwise assocFile = new AssociationFilePairwise();

				pvalOut.writeln(assocFile.getHeader());
				while (returned < submitted) {
					try {
						Future<AssociationResultPairwise> future = jobHandler.take();
						if (future != null) {
							AssociationResult result = future.get();
							if (result != null) {
								// write to disk
								pvalOut.writeln(result.toString());
							}
						}
						returned++;
						pb.iterate();
					} catch (InterruptedException e) {
						e.printStackTrace();
					} catch (ExecutionException e) {
						e.printStackTrace();
					}
				}
				pvalOut.close();
				pb.close();
			}
		}

	}

	private Pair<VCFVariant, AssociationResult> getBestAssocForRegion(ArrayList<AssociationResult> assocResults,
																	  Feature region,
																	  ArrayList<VCFVariant> variantsInRegion,
																	  String variantQuery) {

		AssociationResult topResult = null;
		if (variantQuery != null) {
			for (AssociationResult r : assocResults) {
				if (r.getSnp().overlaps(region)) {
					String id = r.getSnp().getChromosome().toString() + "_" + r.getSnp().getStart() + "_" + r.getSnp().getName();
					if (id.equals(variantQuery)) {
						topResult = r;
					}
				}
			}
		} else {
			for (AssociationResult r : assocResults) {
				if (r.getSnp().overlaps(region)) {

					if (topResult == null) {
						topResult = r;
					} else {
						if (topResult.getLog10Pval() < r.getLog10Pval()) {
							topResult = r;
						}
					}
				}
			}
		}


		// get the variant belonging to this assoc result
		if (topResult != null) {
			for (VCFVariant v : variantsInRegion) {
				if (v.asFeature().overlaps(topResult.getSnp())) {
					String alleles1 = Strings.concat(v.getAlleles(), Strings.comma);
					String alleles2 = Strings.concat(topResult.getSnp().getAlleles(), Strings.comma);

					if (v.getId().equals(topResult.getSnp().getName()) && alleles1.equals(alleles2)) {
						return new Pair<VCFVariant, AssociationResult>(v, topResult);
					}
				}
			}
		}

		return null;
	}

	private ArrayList<VCFVariant> filterVariantsByRegion(ArrayList<VCFVariant> variants, Feature region) {
		ArrayList<VCFVariant> output = variants.stream().filter(v -> v.asFeature().overlaps(region)).collect(Collectors.toCollection(ArrayList::new));
		return output;
	}

	private void clearQueue(TextFile logout, TextFile pvalout,
							int iter, ArrayList<VCFVariant> variants,
							CompletionService<Triple<String, AssociationResult, VCFVariant>> jobHandler,
							ArrayList<AssociationResult> associationResults) throws IOException {
//		System.out.println(submitted + " results to process.");
		while (returned < submitted) {
			try {
				Future<Triple<String, AssociationResult, VCFVariant>> future = jobHandler.take();
				if (future != null) {
					Triple<String, AssociationResult, VCFVariant> result = future.get();
					if (result != null) {
						AssociationResult assoc = result.getMiddle();
						VCFVariant variant = result.getRight();
						String logStr = result.getLeft();
						if (logStr != null && logout != null) {
							logout.writeln(logStr);
						}
						if (assoc != null) {
							double p = assoc.getLog10Pval();
							if (p > highestLog10p) {
								highestLog10p = p;
							}

							if (options.testMultiAllelicVariantsIndependently && variant.getNrAlleles() > 2) {
								AssociationResult[] subresults = assoc.getSubresults();
								for (AssociationResult r : subresults) {
									pvalout.writeln(r.toString());
								}
							} else {
								pvalout.writeln(assoc.toString());
							}


							if (associationResults != null) {
								associationResults.add(assoc);
							}
							if (variant != null && variants != null && iter == 0 && options.getMaxIter() > 1) {
								variants.add(variant);
							}
						}

					}
				}
				returned++;
				if (progressBar != null) {
					progressBar.iterate();
				}
//				System.out.print(returned + "/" + submitted);

			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}
//		System.out.println();
		returned = 0;
		submitted = 0;
	}

	private ArrayList<Pair<String, Triple<String, String, String>>> getTrios(String famfile) throws IOException {
		System.out.println("Loading trios from FAM path: " + famfile);
		ArrayList<Pair<String, Triple<String, String, String>>> output = new ArrayList<>();
		TextFile tf = new TextFile(famfile, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {

			if (elems.length >= 4) {
				String family = elems[0];

				String kid = elems[1];
				String dad = elems[2];
				String mom = elems[3];

				if (!dad.equals("0") && !mom.equals("0")) {
					Integer control = Integer.parseInt(elems[5]);
					if (control.equals(2)) {
						Triple<String, String, String> t = new Triple<String, String, String>(kid, mom, dad);
						Pair<String, Triple<String, String, String>> p = new Pair<String, Triple<String, String, String>>(family, t);
						output.add(p);
					}
				}
			}
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();
		System.out.println(output.size() + " trios found in FAM path");
		return output;
	}

	public void determineSamplesToRemoveFromFam(HashSet<String> samplesToRemove, HashSet<String> samplesWithGenotypeAndDiseaseStatus) throws IOException {

		// get the trios (only trios where kid is a case)
		ArrayList<Pair<String, Triple<String, String, String>>> famData = getTrios(options.getFamfile()); // format: familyname<kid, mom, dad>

		HashSet<String> visitedFamilies = new HashSet<String>();
		if (famData.size() > 0) {

			for (Pair<String, Triple<String, String, String>> p : famData) {
				Triple<String, String, String> trio = p.getRight();

				String fam = p.getLeft();
				String kid = trio.getLeft();
				if (visitedFamilies.contains(fam)) {
					String altkid = kid + "_" + kid;
					samplesToRemove.add(kid);
					samplesToRemove.add(altkid);
				} else {
					visitedFamilies.add(fam);
				}
				String mom = trio.getRight();
				String dad = trio.getMiddle();
				if (samplesWithGenotypeAndDiseaseStatus.contains(kid)) {
					samplesToRemove.add(mom);
					samplesToRemove.add(dad);
				} else {
					String altkid = kid + "_" + kid;
					if (samplesWithGenotypeAndDiseaseStatus.contains(altkid)) {
						String altDad = dad + "_" + dad;
						String altMom = mom + "_" + mom;
						samplesToRemove.add(altDad);
						samplesToRemove.add(dad);

						samplesToRemove.add(mom);
						samplesToRemove.add(altMom);
					}
				}
			}

			System.out.println(samplesToRemove.size() + " samples to remove using FAM path");
		}

	}

	public Pair<LogisticRegressionResult, Integer> getNullModel(VCFVariant variant,
																ArrayList<Pair<VCFVariant, Triple<int[], boolean[], Integer>>> conditional,
																int firstColumnToRemove,
																int lastColumnToRemove) {
		// get a random variant
		// prepare the matrix

		// generate pseudocontrol genotypes
		LRTestTask lrt = new LRTestTask(sampleAnnotation);
		Pair<DoubleMatrix2D, double[]> xandy = lrt.prepareMatrices(
				variant,
				conditional
		);

		// remove collinear variables and prune
		Pair<DoubleMatrix2D, boolean[]> pruned = lrt.removeCollinearVariables(xandy.getLeft());


		DoubleMatrix2D x = pruned.getLeft(); // x is now probably shorter than original X
		DenseDoubleAlgebra dda = new DenseDoubleAlgebra();
		boolean[] notaliased = pruned.getRight(); // length of original X

		// check if the alleles we put in are aliased
		int firstAllele = firstColumnToRemove; // intercept + allleles we conditioned on (original X indexing)
		int lastAllele = lastColumnToRemove;   // intercept + allleles we conditioned on + alleles for this variant (original X indexing)


		int nrRemaining = 0;
		int[] alleleIndex = new int[lastColumnToRemove - firstColumnToRemove];
		for (int i = 0; i < alleleIndex.length; i++) {
			alleleIndex[i] = -1;
		}

		int newXIndex = 0;
		ArrayList<Integer> colIndexArr = new ArrayList<>(x.columns());
		for (int i = 0; i < notaliased.length; i++) {
			if (notaliased[i]) {
				if (i >= firstAllele && i < lastAllele) {
					alleleIndex[i - firstAllele] = newXIndex;
					nrRemaining++;
				} else {
					colIndexArr.add(newXIndex);
				}
				newXIndex++;
			}
		}

		DoubleMatrix2D xprime = dda.subMatrix(x, 0, x.rows() - 1, Primitives.toPrimitiveArr(colIndexArr.toArray(new Integer[0])));
		LogisticRegressionOptimized reg = new LogisticRegressionOptimized();
		return new Pair<>(reg.univariate(xandy.getRight(), xprime), xprime.columns());
	}


}
