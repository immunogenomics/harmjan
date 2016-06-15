package nl.harmjanwestra.gwas;

import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.text.Strings;

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

	LRTestOptions options;
	private int submitted;
	private int returned;
	private double highestLog10p;
	private HashMap<String, Integer> sampleToIntGenotypes;
	private HashSet<String> snpLimit;
	private ArrayList<Feature> bedRegions = null;
	private boolean[] genotypeSamplesWithCovariatesAndDiseaseStatus;
	private DoubleMatrix2D finalCovariates;
	private DiseaseStatus[] finalDiseaseStatus;
	private ProgressBar progressBar;
	private final ExecutorService exService;

	public LRTest(LRTestOptions options) throws IOException {
		this.options = options;
		System.out.println("Setting up threadpool with: " + options.getNrThreads() + " threads..");
		exService = Executors.newFixedThreadPool(options.getNrThreads());
		if (options.isConditionalAnalysis()) {
			testConditional();
		} else if (options.isExhaustivePairwiseAnalysis()) {
			testExhaustivePairwise();
		} else {
			testNormal();
		}
		exService.shutdown();
	}

	public static void main(String[] args) {

		String[] args2 = new String[]{
				"-c", "/Data/tmp/2016-05-27/2016-03-11-T1D-covarmerged.txtmergedCovariates-withPseudos.txt",
				"-d", "/Data/tmp/2016-05-27/2016-03-11-T1D-diseaseStatusWithPseudos.txt",
				"-f", "/Data/tmp/2016-05-27/T1D-recode-maf0005-ICRegions-samplenamefix-pseudo.vcf.gz.fam",
				"-i", "/Data/tmp/2016-05-27/test.vcf",
				"-e", "/Data/tmp/2016-05-27/T1D-recode-regionsfiltered-allelesfiltered-samplenamefix-pseudo.vcf.gz-parents.txt",
				"-r", "/Data/tmp/2016-05-27/AllICLoci.bed",
				"-t", "1",
				"-q", "0.3",
				"-o", "/Data/tmp/2016-05-27/testout.txt"
		};

		String[] args3 = new String[]{
				"-c", "/Data/tmp/2016-06-10/covarmerged.txtmergedCovariates.txt",
				"-d", "/Data/tmp/2016-06-10/covarmerged.txtmergeddisease.txt",
				"-f", "/Data/tmp/2016-06-10/covarmerged.txtmergedfam.fam",
//				"-i", "/Data/tmp/2016-06-10/RA-Beagle1kg-regionfiltered-COSMO-ImpQualsReplaced-chr22.vcf.gz",
				"-i", "/Data/tmp/2016-06-10/test.vcf",
				"-r", "/Data/tmp/2016-06-10/loci.bed",
//				"--exhaustive",
//				"--maxiter", "5",
				"-t", "1",
				"-q", "0.3",
				"-o", "/Data/tmp/2016-06-10/testout.txt"
		};

		LRTestOptions options = new LRTestOptions(args3);

		try {
			new LRTest(options);

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
		System.out.println("Will perform logistic regression analysis");

		boolean initialized = initialize();

		// waitForEnter("Waiting for a response...");

		if (initialized) {

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
			ArrayList<Pair<VCFVariant, Triple<int[], boolean[], Triple<Integer, Double, Double>>>> conditional = new ArrayList<>();

			ArrayList<String> conditionalVariantIds = new ArrayList<String>();
			AssociationFile associationFile = new AssociationFile();
			String header = associationFile.getHeader();

			submitted = 0;
			returned = 0;

			int alleleOffsetGenotypes = 0;
			highestLog10p = 0;


			CompletionService<Triple<String, AssociationResult, VCFVariant>> jobHandler = new ExecutorCompletionService<Triple<String, AssociationResult, VCFVariant>>(exService);
			for (int i = 0; i < variants.size(); i++) {
				VCFVariant variant = variants.get(i);
				// throw into a thread
				// TODO: conditional on dosages is separate from conditional on genotypes.. ?
				LRTestTask task = new LRTestTask(variant,
						iter,
						finalDiseaseStatus,
						finalCovariates,
						conditional,
						alleleOffsetGenotypes,
						options);
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
	}

	private boolean initialize() throws IOException {
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
		HashMap<String, DiseaseStatus> diseaseStatus = new HashMap<String, DiseaseStatus>();

		// load disease status
		TextFile tf = new TextFile(options.getDiseaseStatusFile(), TextFile.R);
		String[] elems = tf.readLineElems(Strings.tab);
		while (elems != null) {
			String sample = elems[0];
			String statusStr = elems[1];
			if (excludeTheseSamples == null || !excludeTheseSamples.contains(sample)) {
				diseaseStatus.put(sample, DiseaseStatus.parseStatus(statusStr));
			}
			elems = tf.readLineElems(Strings.tab);
		}
		tf.close();

		System.out.println(diseaseStatus.size() + " disease status samples loaded");

		// load samples from vcf file
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

		// remove related samples based on fam file if any
		if (options.getFamfile() != null) {
			if (excludeTheseSamples == null) {
				excludeTheseSamples = new HashSet<String>();
			}
			determineSamplesToRemoveFromFam(excludeTheseSamples, vcfSamplesWithDiseaseStatus);

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
			finalDiseaseStatus = new DiseaseStatus[sampleToIntGenotypes.size()];
			finalCovariates = new DenseDoubleMatrix2D(sampleToIntGenotypes.size(), covariates.columns());

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
					finalDiseaseStatus[id] = diseaseStatus.get(sample);
					if (finalDiseaseStatus[id].equals(DiseaseStatus.CONTROL)) {
						nrControls++;
					} else if (finalDiseaseStatus[id].equals(DiseaseStatus.CASE)) {
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

			return true;
		}
	}


	public void testConditional() throws IOException {

		System.out.println("Will perform conditional analysis");

		boolean initialized = initialize();

		if (initialized) {

			AssociationFile associationFile = new AssociationFile();
			String header = associationFile.getHeader();

			// boot up threadpool

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
				ArrayList<Pair<VCFVariant, Triple<int[], boolean[], Triple<Integer, Double, Double>>>> conditional = new ArrayList<>();

				LRTestTask task = new LRTestTask(variant,
						0,
						finalDiseaseStatus,
						finalCovariates,
						conditional,
						0,
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
				remainingRegions.addAll(visitedRegions);


				System.out.println("Total number of iterations to run: " + options.getMaxIter());
				ArrayList<ArrayList<VCFVariant>> conditionalVariants = new ArrayList<>();
				for (int i = 0; i < remainingRegions.size(); i++) {
					conditionalVariants.add(new ArrayList<>());
				}

				LRTestVariantQCTask tasktmp = new LRTestVariantQCTask();
				TextFile modelsout = new TextFile(options.getOutputdir() + "gwas-conditional-models.txt", TextFile.W);
				modelsout.writeln("Region\tIter\tModel");
				TextFile topvariantsout = new TextFile(options.getOutputdir() + "gwas-topvariants.txt", TextFile.W);
				topvariantsout.writeln("Region\tIter\tVariant\tPval");

				for (int i = 1; i < options.getMaxIter(); i++) {

					System.out.println("Output will be written here: " + options.getOutputdir() + "gwas-" + i + ".txt");
					pvalout = new TextFile(options.getOutputdir() + "gwas-" + i + ".txt", TextFile.W);

					pvalout.writeln(header);

					// start new result array
					ArrayList<AssociationResult> assocResultsIter = new ArrayList<>();
					submitted = 0;
					for (int regionId = 0; regionId < remainingRegions.size(); regionId++) {
						// get variants in region
						ArrayList<VCFVariant> variantsInRegion = filterVariantsByRegion(allVariants, remainingRegions.get(regionId));
						Pair<VCFVariant, AssociationResult> bestAssocLastIter = getBestAssocForRegion(assocResults, remainingRegions.get(regionId), variantsInRegion);

						VCFVariant bestVariantLastIter = bestAssocLastIter.getLeft();
						topvariantsout.writeln(remainingRegions.get(regionId).toString() + "\t" + (i - 1) + "\t" + bestVariantLastIter.getChr() + "_" + bestVariantLastIter.getPos() + "_" + bestVariantLastIter.getId() + "\t" + bestAssocLastIter.getRight().getLog10Pval());
						// get conditional variants for this region
						ArrayList<VCFVariant> conditionalVariantsForRegion = conditionalVariants.get(regionId);
						if (conditionalVariantsForRegion == null) {
							conditionalVariantsForRegion = new ArrayList<>();
						}
						conditionalVariantsForRegion.add(bestVariantLastIter);

						// now put the data in the correct data structure.. recode and whatever
						ArrayList<Pair<VCFVariant, Triple<int[], boolean[], Triple<Integer, Double, Double>>>> conditional = new ArrayList<>();


						String model = "y ~ intercept ";
						LRTestVariantQCTask lqr = new LRTestVariantQCTask();
						int alleleOffsetGenotypes = 0;
						for (int c = 0; c < conditionalVariantsForRegion.size(); c++) {
							VCFVariant variant = conditionalVariantsForRegion.get(c);

							model += " + " + variant.getId();
							conditional.add(new Pair<>(variant,
									lqr.filterAndRecodeGenotypes(variant.getGenotypeAllelesAsMatrix2D(),
											finalDiseaseStatus,
											variant.getAlleles().length,
											finalCovariates.rows())));
							alleleOffsetGenotypes += variant.getAlleles().length - 1;
						}
						modelsout.writeln(remainingRegions.get(regionId).toString() + "\t" + i + "\t" + model);

						int alleleOffsetDosages = 0;

						for (int v = 0; v < variantsInRegion.size(); v++) {
							// throw into a thread
							// TODO: conditional on dosages is separate from conditional on genotypes.. ?
							VCFVariant variant = variantsInRegion.get(v);
							LRTestTask task = new LRTestTask(variant,
									i,
									finalDiseaseStatus,
									finalCovariates,
									conditional,
									alleleOffsetGenotypes,
									options);

							jobHandler.submit(task);
							submitted++;


							if (submitted % 1000 == 0) {
								clearQueue(null, pvalout, i, null, jobHandler, assocResultsIter);
							}
						}
					}

					clearQueue(null, pvalout, i, null, jobHandler, assocResultsIter);
					assocResults = assocResultsIter;
					pvalout.close();
				}
				for (int regionId = 0; regionId < remainingRegions.size(); regionId++) {
					// get variants in region
					ArrayList<VCFVariant> variantsInRegion = filterVariantsByRegion(allVariants, remainingRegions.get(regionId));
					Pair<VCFVariant, AssociationResult> bestAssocLastIter = getBestAssocForRegion(assocResults, remainingRegions.get(regionId), variantsInRegion);
					VCFVariant bestVariantLastIter = bestAssocLastIter.getLeft();
					topvariantsout.writeln(remainingRegions.get(regionId).toString() + "\t" + (options.getMaxIter() - 2) + "\t" + bestVariantLastIter.getChr() + "_" + bestVariantLastIter.getPos() + "_" + bestVariantLastIter.getId() + "\t" + bestAssocLastIter.getRight().getLog10Pval());
				}
				modelsout.close();
				topvariantsout.close();
			}
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
								variants.add(result.getLeft());
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

				System.out.println(totalsubmitted + " submitted. clearing queue. " + variants.size() + " variants in memory " + ManagementFactory.getThreadMXBean().getThreadCount());
			}

			LRTestVariantQCTask task = new LRTestVariantQCTask(vcfLn,
					bedRegions,
					options,
					genotypeSamplesWithCovariatesAndDiseaseStatus,
					finalDiseaseStatus,
					finalCovariates);
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
						variants.add(result.getLeft());
					}
				}
				returned++;

			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}
		System.out.println(variants.size() + " variants included.");

		log.close();
		return variants;

	}

	public void testExhaustivePairwise() throws IOException {
		System.out.println("Will perform exhaustive pairwise analysis");
		boolean initialized = initialize();
		if (initialized) {
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
					CompletionService<AssociationResult> jobHandler = new ExecutorCompletionService<AssociationResult>(exService);

					int submitted = 0;
					for (Feature region : availableRegions) {
						ArrayList<VCFVariant> variants = variantsInRegions.get(region);
						for (int i = 0; i < variants.size(); i++) {
							for (int j = i + 1; j < variants.size(); j++) {
								// submit job to queue;
								LRTestExhaustiveTask lrtet = new LRTestExhaustiveTask(variants, i, j, genotypeSamplesWithCovariatesAndDiseaseStatus, finalDiseaseStatus, finalCovariates, options);

								jobHandler.submit(lrtet);
								submitted++;

							}
						}
					}

					System.out.println("Submitted a total of " + submitted + " jobs");
					ProgressBar pb = new ProgressBar(submitted);
					int returned = 0;
					TextFile pvalOut = new TextFile(options.getOutputdir() + "-pairwise.txt.gz", TextFile.W);
					AssociationFile assocFile = new AssociationFile();
					assocFile.setPairWise(true);
					pvalOut.writeln(assocFile.getHeader());
					while (returned < submitted) {
						try {
							Future<AssociationResult> future = jobHandler.take();
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
	}

	private Pair<VCFVariant, AssociationResult> getBestAssocForRegion(ArrayList<AssociationResult> assocResults,
																	  Feature region,
																	  ArrayList<VCFVariant> variantsInRegion) {

		AssociationResult topResult = null;
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

		// get the variant belonging to this assoc result
		for (VCFVariant v : variantsInRegion) {
			if (v.asFeature().overlaps(topResult.getSnp())) {
				String alleles1 = Strings.concat(v.getAlleles(), Strings.comma);
				String alleles2 = Strings.concat(topResult.getSnp().getAlleles(), Strings.comma);

				if (v.getId().equals(topResult.getSnp().getName()) && alleles1.equals(alleles2)) {
					return new Pair<VCFVariant, AssociationResult>(v, topResult);
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
							pvalout.writeln(assoc.toString());
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
		System.out.println("Loading trios from FAM file: " + famfile);
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
		System.out.println(output.size() + " trios found in FAM file");
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

			System.out.println(samplesToRemove.size() + " samples to remove using FAM file");
		}

	}


}
