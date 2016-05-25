package nl.harmjanwestra.gwas;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
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

import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;

/**
 * Created by hwestra on 11/7/15.
 */
public class LRTest {

//	boolean useR = false;

	LRTestOptions options;
	private int submitted;
	private int returned;
	private double highestLog10p;
	private VCFVariant currentLowestVariant;
	private double highestLog10PProbs;
	private int nrTested;
	private int maxNrOfResultsInBuffer = 10000;
	private HashMap<String, Integer> sampleToIntGenotypes;
	private HashSet<String> snpLimit;
	private ArrayList<Feature> bedRegions = null;
	private boolean[] genotypeSamplesWithCovariatesAndDiseaseStatus;
	private DenseDoubleMatrix2D finalCovariates;
	private double[] finalDiseaseStatus;

	public LRTest(LRTestOptions options) throws IOException {
		this.options = options;
		if (options.isConditionalAnalysis()) {
			testConditional();
		} else if (options.isExhaustivePairwiseAnalysis()) {
			testExhaustivePairwise();
		} else {
			testNormal();
		}

	}

	public void testNormal() throws IOException {
		boolean initialized = initialize();

		if (initialized) {

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

			// TODO: other VCF info scores...

			// keep a list of genotypes to condition on
			int iter = 0;
			ArrayList<Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>>> conditional = new ArrayList<>();
			ArrayList<DoubleMatrix2D> conditionalDosages = new ArrayList<>();
			ArrayList<String> conditionalVariantIds = new ArrayList<String>();
			AssociationFile associationFile = new AssociationFile();
			String header = associationFile.getHeader();
			ArrayList<VCFVariant> variants = new ArrayList<>();

			// boot up threadpool
			System.out.println("Setting up threadpool with: " + options.getNrThreads() + " threads..");
			ExecutorService threadPool = Executors.newFixedThreadPool(options.getNrThreads());
			CompletionService<Triple<String, AssociationResult, VCFVariant>> jobHandler = new ExecutorCompletionService<Triple<String, AssociationResult, VCFVariant>>(threadPool);

			while (iter < options.getMaxIter()) {
				TextFile logout = null;
				if (iter == 0) {
					System.out.println("Iteration " + iter + " starting. Model: y ~ SNP + covar.");
					logout = new TextFile(options.getOutputdir() + "log-iter" + iter + ".txt", TextFile.W);
					String logoutheader = "SNP\tChr\tPos\tAlleles\tMinorAllele\tImputationQual\tMAF\tHWEP\tTested";
					logout.writeln(logoutheader);
					System.out.println("Log will be written here: " + options.getOutputdir() + "log-iter" + iter + ".txt");
				} else {
					System.out.println("Iteration " + iter + " starting. Model: y ~ " + Strings.concat(conditionalVariantIds, Strings.semicolon)
							+ " + SNP + covar.");
					System.out.println(variants.size() + " variants in total.");
				}

				TextFile pvalout = new TextFile(options.getOutputdir() + "gwas-" + iter + ".txt", TextFile.W);
				System.out.println("Output will be written here: " + options.getOutputdir() + "gwas-" + iter + ".txt");
				pvalout.writeln(header);

				currentLowestVariant = null;

				int variantCtr = 0;
				nrTested = 0;
				submitted = 0;
				returned = 0;

				int alleleOffsetGenotypes = 0;
				int alleleOffsetDosages = 0;
				highestLog10p = 0;
				highestLog10PProbs = 0;

				for (Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>> c : conditional) {
					alleleOffsetGenotypes += c.getLeft().rows();
				}
				for (DoubleMatrix2D c : conditionalDosages) {
					alleleOffsetDosages += c.columns();
				}

				String vcfLn = null;
				TextFile vcfIn = null;
				if (iter == 0) {
					vcfIn = new TextFile(options.getVcf(), TextFile.R);
					vcfLn = vcfIn.readLine();
					while (vcfLn != null && vcfLn.startsWith("#")) {
						vcfLn = vcfIn.readLine();
					}
				}

				while ((iter == 0 && vcfLn != null) || (iter > 0 && variantCtr < variants.size())) {

					VCFVariant variant = null;
					if (iter == 0) {
						variant = new VCFVariant(vcfLn, VCFVariant.PARSE.HEADER);
					} else {
						variant = variants.get(variantCtr);
					}

					String name = variant.getId();
					if (snpLimit == null || snpLimit.contains(name)) {

						Double imputationqualityscore = variant.getImputationQualityScore();
						if (imputationqualityscore == null) {
							imputationqualityscore = Double.NaN;
						}


						boolean impqualscoreOK = false;

						if ((Double.isNaN(imputationqualityscore) && options.isTestVariantsWithoutImputationQuality())
								|| (!Double.isNaN(imputationqualityscore) && imputationqualityscore >= options.getImputationqualitythreshold())) {
							impqualscoreOK = true;
						} else if (Double.isNaN(imputationqualityscore)) {
							System.err.println("No imputaton quality score for variant: " + variant.getChr() + "-" + variant.getPos() + "-" + variant.getId());
							System.err.println("In file: " + options.getVcf());
							if (iter == 0) {
//								SNP	Chr	Pos	ImputationQual	MAF	Region	OverlapOK	MAFOk	ImpQualOK
								logout.writeln(variant.getId()
										+ "\t" + variant.getChr()
										+ "\t" + variant.getPos()
										+ "\t" + Strings.concat(variant.getAlleles(), Strings.comma)
										+ "\t" + variant.getMinorAllele()
										+ "\t" + imputationqualityscore
										+ "\t" + null
										+ "\t" + null
										+ "\t" + null
										+ "\t" + false);
							}
						}

						boolean overlap = false;
						if (bedRegions == null) {
							overlap = true;
						} else if (bedRegions != null && impqualscoreOK) {
							Feature variantFeature = variant.asFeature();
//
//							if (variant.getId().equals("rs12022522")) {
//								System.out.println("Got it");
//							}

							for (Feature r : bedRegions) {
								if (r.overlaps(variantFeature)) {
									overlap = true;
									break;
								}
							}
//							if (variant.getId().equals("rs12022522")) {
//								System.out.println(overlap);
//								for (Feature r : bedRegions) {
//									System.out.println(r.toString() + "\t" + variantFeature.toString() + "\t" + r.overlaps(variantFeature));
//								}
//								System.exit(-1);
//							}
						}

						if (!impqualscoreOK || !overlap) {
							if (iter == 0) {
//								SNP	Chr	Pos	ImputationQual	MAF	Region	OverlapOK	MAFOk	ImpQualOK
								logout.writeln(variant.getId()
										+ "\t" + variant.getChr()
										+ "\t" + variant.getPos()
										+ "\t" + Strings.concat(variant.getAlleles(), Strings.comma)
										+ "\t" + variant.getMinorAllele()
										+ "\t" + imputationqualityscore
										+ "\t" + null
										+ "\t" + overlap
										+ "\t" + null
										+ "\t" + false);
							}
						} else {

							// throw into a thread
							// TODO: conditional on dosages is separate from conditional on genotypes.. ?
							LRTestTask task = new LRTestTask(vcfLn,
									variant,
									iter,
									genotypeSamplesWithCovariatesAndDiseaseStatus,
									finalDiseaseStatus,
									finalCovariates,
									conditional,
									conditionalDosages,
									alleleOffsetGenotypes,
									alleleOffsetDosages,
									options);
							jobHandler.submit(task);
							submitted++;

							if (submitted % 1000 == 0) {
								clearQueue(logout, pvalout, iter, variants, jobHandler, null);
							}
						}
					} // end snplimit.contains(snp)
					variantCtr++;

					if (variantCtr % 1000 == 0) {
						String currentLowestDosagePValSNPId = null;
						if (currentLowestVariant != null) {
							currentLowestDosagePValSNPId = currentLowestVariant.getId();
						}
						System.out.println("Iteration: " + iter + "\tNr variants: " + variantCtr + " loaded.\tSubmitted to queue: " + submitted + "\tNr Tested: " + nrTested + "\tHighest P-val: " + highestLog10p + "\tSNP: " + currentLowestDosagePValSNPId);
					}
					if (iter == 0) {
						vcfLn = vcfIn.readLine();
					}
				} // end data.hasnext

				if (iter == 0) {
					vcfIn.close();
				}

				if (submitted > 0) {
					clearQueue(logout, pvalout, iter, variants, jobHandler, null);
				}

				String currentLowestDosagePValSNPId = null;
				if (currentLowestVariant != null) {
					currentLowestDosagePValSNPId = currentLowestVariant.getId();
				}
				System.out.println("Iteration: " + iter + " is done. \tNr variants: " + variantCtr + "\tNr Tested: " + nrTested + "\tHighest P-val: " + highestLog10p + "\tBy variant: " + currentLowestDosagePValSNPId);

				// TODO: this code does not handle the absence of imputation dosages very well...
				LRTestTask tasktmp = new LRTestTask();
				if (conditionsPerIter != null) {
					System.out.println("Iterating " + variants.size() + " variants");

					// iterate the variants
					int nextIter = iter + 1;
					for (VCFVariant variant : variants) {
						Feature variantFeature = variant.asFeature();

						if (conditionsPerIter.containsKey(variantFeature)) {
							Integer iterForVariant = conditionsPerIter.get(variantFeature);
							System.out.println("Variant found: " + variantFeature.getName() + " iter: " + iterForVariant + " next: " + nextIter);

							if (iterForVariant != null && iterForVariant <= nextIter) {
								Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>> unfilteredGenotypeData = tasktmp.filterAndRecodeGenotypes(
										genotypeSamplesWithCovariatesAndDiseaseStatus,
										variant.getGenotypeAllelesAsMatrix2D(),
										finalDiseaseStatus,
										variant.getAlleles().length,
										finalCovariates.rows());
								conditional.add(unfilteredGenotypeData);

								DoubleMatrix2D dosages = variant.getImputedDosagesAsMatrix2D();
								conditionalDosages.add(dosages);
								conditionalVariantIds.add(variant.getId());
							}
						}
					}
				} else if (currentLowestVariant != null) {
					VCFVariant variant = currentLowestVariant;
					Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>> unfilteredGenotypeData = tasktmp.filterAndRecodeGenotypes(
							genotypeSamplesWithCovariatesAndDiseaseStatus,
							variant.getGenotypeAllelesAsMatrix2D(),
							finalDiseaseStatus,
							variant.getAlleles().length,
							finalCovariates.rows());
					conditional.add(unfilteredGenotypeData);


					DoubleMatrix2D dosages = variant.getImputedDosagesAsMatrix2D();
					conditional.add(unfilteredGenotypeData);
					conditionalDosages.add(dosages);
					conditionalVariantIds.add(variant.getId());
					System.out.println("Iteration: " + iter
							+ "\tNr variants: " + variantCtr
							+ "\tHighest P-val: " + highestLog10PProbs
							+ "\tBy variant: " + variant.getId());

				}

				if (iter == 0) {
					logout.close();
				}

				pvalout.close();
				iter++;
			}
			threadPool.shutdown();
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
		HashMap<String, Integer> diseaseStatus = new HashMap<String, Integer>();

		// load disease status
		TextFile tf = new TextFile(options.getDiseaseStatusFile(), TextFile.R);
		String[] elems = tf.readLineElems(Strings.tab);
		while (elems != null) {
			String sample = elems[0];
			String statusStr = elems[1];
			Integer status = Integer.parseInt(statusStr);
			if (excludeTheseSamples == null || !excludeTheseSamples.contains(sample)) {
				diseaseStatus.put(sample, status);
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
			finalDiseaseStatus = new double[sampleToIntGenotypes.size()];
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
					Integer status = diseaseStatus.get(sample);
					int remapstatus = (status - 1);
					finalDiseaseStatus[id] = remapstatus;

					if (status == 2) {
						nrCases++;
					} else if (status == 1) {
						nrControls++;
					} else {
						nrUnknown++;
					}

					if (remapstatus > maxstatus) {
						maxstatus = remapstatus;

					}
					if (remapstatus < minstatus) {
						minstatus = remapstatus;

					}
					for (int col = 0; col < covariates.columns(); col++) {
						finalCovariates.setQuick(id, col, covariates.getElement(sid, col));
					}
					nrTotal++;
				}
			}

			System.out.println(nrCases + " cases ");
			System.out.println(nrControls + " controls ");
			System.out.println(nrUnknown + " unknown ");
			System.out.println(nrTotal + " with covariates");

			System.out.println("max/min status: " + maxstatus + "/" + minstatus);
			return true;
		}
	}


	public void testConditional() throws IOException {

		boolean initialized = initialize();

		if (initialized) {

			// TODO: other VCF info scores...

			// keep a list of genotypes to condition on


			AssociationFile associationFile = new AssociationFile();
			String header = associationFile.getHeader();


			// boot up threadpool
			System.out.println("Setting up threadpool with: " + options.getNrThreads() + " threads..");
			ExecutorService threadPool = Executors.newFixedThreadPool(options.getNrThreads());
			CompletionService<Triple<String, AssociationResult, VCFVariant>> jobHandler = new ExecutorCompletionService<Triple<String, AssociationResult, VCFVariant>>(threadPool);


			boolean initialMappingPerformed = false;
			ArrayList<VCFVariant> variants = new ArrayList<>();
			ArrayList<AssociationResult> assocResults = new ArrayList<>();
			if (!initialMappingPerformed) {
				int iter = 0;
				System.out.println("Iteration " + iter + " starting. Model: y ~ SNP + covar.");
				TextFile logout = new TextFile(options.getOutputdir() + "log-iter" + iter + ".txt", TextFile.W);
				String logoutheader = "SNP\tChr\tPos\tAlleles\tMinorAllele\tImputationQual\tMAF\tHWEP\tTested";
				logout.writeln(logoutheader);
				System.out.println("Log will be written here: " + options.getOutputdir() + "log-iter" + iter + ".txt");

				System.out.println("Output will be written here: " + options.getOutputdir() + "gwas-" + iter + ".txt");
				TextFile pvalout = new TextFile(options.getOutputdir() + "gwas-" + iter + ".txt", TextFile.W);
				pvalout.writeln(header);

				currentLowestVariant = null;

				int variantCtr = 0;
				nrTested = 0;
				submitted = 0;
				returned = 0;


				highestLog10p = 0;
				highestLog10PProbs = 0;


				String vcfLn = null;
				TextFile vcfIn = null;

				vcfIn = new TextFile(options.getVcf(), TextFile.R);
				vcfLn = vcfIn.readLine();
				while (vcfLn != null && vcfLn.startsWith("#")) {
					vcfLn = vcfIn.readLine();
				}

				VCFVariant variant = null;

				while (vcfLn != null) {
					variant = new VCFVariant(vcfLn, VCFVariant.PARSE.HEADER);
					String name = variant.getId();
					if (snpLimit == null || snpLimit.contains(name)) {
						Double imputationqualityscore = variant.getImputationQualityScore();
						if (imputationqualityscore == null) {
							imputationqualityscore = Double.NaN;
						}

						boolean impqualscoreOK = false;
						if ((Double.isNaN(imputationqualityscore) && options.isTestVariantsWithoutImputationQuality())
								|| (!Double.isNaN(imputationqualityscore) && imputationqualityscore >= options.getImputationqualitythreshold())) {
							impqualscoreOK = true;
						} else if (Double.isNaN(imputationqualityscore)) {
							System.err.println("No imputaton quality score for variant: " + variant.getChr() + "-" + variant.getPos() + "-" + variant.getId());
							System.err.println("In file: " + options.getVcf());
							if (iter == 0) {
//								SNP	Chr	Pos	ImputationQual	MAF	Region	OverlapOK	MAFOk	ImpQualOK
								logout.writeln(variant.getId()
										+ "\t" + variant.getChr()
										+ "\t" + variant.getPos()
										+ "\t" + Strings.concat(variant.getAlleles(), Strings.comma)
										+ "\t" + variant.getMinorAllele()
										+ "\t" + imputationqualityscore
										+ "\t" + null
										+ "\t" + null
										+ "\t" + null
										+ "\t" + false);
							}
						}
						boolean overlap = false;
						if (impqualscoreOK) {
							if (bedRegions == null) {
								overlap = true;
							} else if (bedRegions != null && impqualscoreOK) {
								Feature variantFeature = variant.asFeature();

								for (Feature r : bedRegions) {
									if (r.overlaps(variantFeature)) {
										overlap = true;
										break;
									}
								}
							}

						}

						if (!impqualscoreOK || !overlap) {
							if (iter == 0) {
								logout.writeln(variant.getId()
										+ "\t" + variant.getChr()
										+ "\t" + variant.getPos()
										+ "\t" + Strings.concat(variant.getAlleles(), Strings.comma)
										+ "\t" + variant.getMinorAllele()
										+ "\t" + imputationqualityscore
										+ "\t" + null
										+ "\t" + overlap
										+ "\t" + null
										+ "\t" + false);
							}
						} else {
							// throw into a thread
							// TODO: conditional on dosages is separate from conditional on genotypes.. ?
							ArrayList<Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>>> conditional = new ArrayList<>();
							ArrayList<DoubleMatrix2D> conditionalDosages = new ArrayList<>();
							LRTestTask task = new LRTestTask(vcfLn,
									variant,
									iter,
									genotypeSamplesWithCovariatesAndDiseaseStatus,
									finalDiseaseStatus,
									finalCovariates,
									conditional,
									conditionalDosages,
									0,
									0,
									options);
							jobHandler.submit(task);
							submitted++;

							if (submitted % 1000 == 0) {
								clearQueue(logout, pvalout, iter, variants, jobHandler, assocResults);
							}
						}
					}
					variantCtr++;
					if (variantCtr % 1000 == 0) {
						String currentLowestDosagePValSNPId = null;
						if (currentLowestVariant != null) {
							currentLowestDosagePValSNPId = currentLowestVariant.getId();
						}
						System.out.println("Iteration: " + iter + "\tNr variants: " + variantCtr + " loaded.\tSubmitted to queue: " + submitted + "\tNr Tested: " + nrTested + "\tHighest P-val: " + highestLog10p + "\tSNP: " + currentLowestDosagePValSNPId);
					}
					vcfLn = vcfIn.readLine();
				}


				// clean up the queue
				clearQueue(logout, pvalout, iter, variants, jobHandler, assocResults);
				pvalout.close();
				logout.close();
				vcfIn.close();
			} else {
				System.out.println(variants.size() + " variants in total.");

				// run conditional per region..

				// get iter0 results
				// get best associated variant in region
				// testNormal conditional on best variant

				HashSet<Feature> visitedRegions = new HashSet<Feature>();
				for (VCFVariant v : variants) {
					for (int f = 0; f < bedRegions.size(); f++) {
						if (v.asFeature().overlaps(bedRegions.get(f))) {
							visitedRegions.add(bedRegions.get(f));
						}
					}
				}
				int iter = 1;
				ArrayList<Feature> remainingRegions = new ArrayList<>();
				remainingRegions.addAll(visitedRegions);


				ArrayList<ArrayList<VCFVariant>> conditionalVariants = new ArrayList<>();
				LRTestTask tasktmp = new LRTestTask();
				for (int i = iter; i < options.getMaxIter(); i++) {

					System.out.println("Output will be written here: " + options.getOutputdir() + "gwas-" + iter + ".txt");
					TextFile pvalout = new TextFile(options.getOutputdir() + "gwas-" + iter + ".txt", TextFile.W);
					TextFile modelsout = new TextFile(options.getOutputdir() + "gwas-" + iter + "-models.txt", TextFile.W);
					pvalout.writeln(header);

					// start new result array
					ArrayList<AssociationResult> assocResultsIter = new ArrayList<>();
					submitted = 0;
					for (int regionId = 0; regionId < remainingRegions.size(); regionId++) {
						// get variants in region
						ArrayList<VCFVariant> variantsInRegion = filterVariantsByRegion(variants, remainingRegions.get(regionId));
						VCFVariant bestVariantLastIter = getBestAssocForRegion(assocResults, remainingRegions.get(regionId), variantsInRegion);

						// get conditional variants for this region
						ArrayList<VCFVariant> conditionalVariantsForRegion = conditionalVariants.get(regionId);
						if (conditionalVariantsForRegion == null) {
							conditionalVariantsForRegion = new ArrayList<>();
						}
						conditionalVariantsForRegion.add(bestVariantLastIter);

						// now put the data in the correct data structure.. recode and whatever
						ArrayList<Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>>> conditional = new ArrayList<>();
						ArrayList<DoubleMatrix2D> conditionalDosages = new ArrayList<>();

						String model = remainingRegions.get(regionId).toBedString() + "\ty ~ intercept ";

						for (int c = 0; c < conditionalVariantsForRegion.size(); c++) {
							VCFVariant variant = conditionalVariantsForRegion.get(c);
							model += " + " + variant.getId();
							Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>> unfilteredGenotypeData = tasktmp.filterAndRecodeGenotypes(
									genotypeSamplesWithCovariatesAndDiseaseStatus,
									variant.getGenotypeAllelesAsMatrix2D(),
									finalDiseaseStatus,
									variant.getAlleles().length,
									finalCovariates.rows());
							conditional.add(unfilteredGenotypeData);
							DoubleMatrix2D dosages = variant.getImputedDosagesAsMatrix2D();
							conditionalDosages.add(dosages);
						}
						modelsout.writeln(model);

						int alleleOffsetGenotypes = 0;
						int alleleOffsetDosages = 0;

						for (Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>> c : conditional) {
							alleleOffsetGenotypes += c.getLeft().rows();
						}
						for (DoubleMatrix2D c : conditionalDosages) {
							alleleOffsetDosages += c.columns();
						}


						for (int v = 0; v < variantsInRegion.size(); v++) {
							// throw into a thread
							// TODO: conditional on dosages is separate from conditional on genotypes.. ?
							VCFVariant variant = variantsInRegion.get(v);
							LRTestTask task = new LRTestTask(null,
									variant,
									iter,
									genotypeSamplesWithCovariatesAndDiseaseStatus,
									finalDiseaseStatus,
									finalCovariates,
									conditional,
									conditionalDosages,
									alleleOffsetGenotypes,
									alleleOffsetDosages,
									options);
							jobHandler.submit(task);
							submitted++;

							if (submitted % 1000 == 0) {
								clearQueue(null, pvalout, iter, variants, jobHandler, assocResultsIter);
							}
						}
					}
					clearQueue(null, pvalout, iter, variants, jobHandler, assocResultsIter);
					assocResults = assocResultsIter;
					pvalout.close();
					modelsout.close();
				}
			}
			threadPool.shutdown();
		}
	}

	public void testExhaustivePairwise() throws IOException {
		boolean initialized = initialize();
		if (initialized) {

			// get a list of variants for the region
			TextFile vcfIn = new TextFile(options.getVcf(), TextFile.R);
			String vcfLn = vcfIn.readLine();
			while (vcfLn != null && vcfLn.startsWith("#")) {
				vcfLn = vcfIn.readLine();
			}

			HashMap<Feature, ArrayList<VCFVariant>> variantsInRegions = new HashMap<Feature, ArrayList<VCFVariant>>();
			HashSet<Feature> availableRegions = new HashSet<Feature>();

			LRTestTask tasktmp = new LRTestTask();
			int nrLoaded = 0;
			int nrParsed = 0;
			System.out.println("Parsing: " + options.getVcf());
			TextFile logout = new TextFile(options.getOutputdir() + "variantlog.txt", TextFile.W);
			System.out.println("Variant log is here: " + options.getOutputdir() + "-variantlog.txt.gz");
			while (vcfLn != null) {
				VCFVariant variant = new VCFVariant(vcfLn, VCFVariant.PARSE.HEADER);
				Double impqual = variant.getImputationQualityScore();
				if (impqual == null || impqual > options.getImputationqualitythreshold()) {
					boolean variantPassesQC = true;
					boolean variantHasBeenQCed = false;
					boolean variantIncluded = false;
					boolean overlap = false;
					double maf = 0;
					double hwep = 0;
					for (int f = 0; f < bedRegions.size() && variantPassesQC; f++) {
						Feature region = bedRegions.get(f);
						if (variant.asFeature().overlaps(region)) {

							overlap = true;
							if (!variantHasBeenQCed) {
								variant = new VCFVariant(vcfLn, VCFVariant.PARSE.ALL);
								// parse the genotype, do some QC checks
								Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>> unfilteredGenotypeData = tasktmp.filterAndRecodeGenotypes(
										genotypeSamplesWithCovariatesAndDiseaseStatus,
										variant.getGenotypeAllelesAsMatrix2D(),
										finalDiseaseStatus,
										variant.getAlleles().length,
										finalCovariates.rows());

								Triple<Integer, Double, Double> qc = unfilteredGenotypeData.getRight();

								maf = qc.getMiddle();
								hwep = qc.getRight();
								if (maf < options.getMafthresholdD() || hwep < options.getHWEPThreshold()) {
									variantPassesQC = false;
								}
								variantHasBeenQCed = true;
							}

							if (variantHasBeenQCed && variantPassesQC) {
								ArrayList<VCFVariant> variants = variantsInRegions.get(region);
								if (variants == null) {
									variants = new ArrayList<>();
								}
								variants.add(variant);
								variantsInRegions.put(region, variants);
								availableRegions.add(region);
								variantIncluded = true;
							}
						}
					}
					if (variantIncluded) {
						nrLoaded++;
					}
					logout.writeln(variant.getChr()
							+ "\t" + variant.getPos()
							+ "\t" + variant.getId()
							+ "\t" + maf
							+ "\t" + hwep
							+ "\t" + overlap
							+ "\t" + variantPassesQC
					);
				}
				vcfLn = vcfIn.readLine();
				nrParsed++;
				if (nrParsed % 1000 == 0) {
					System.out.print(nrParsed + " lines parsed.. " + nrLoaded + " loaded\r");
				}
			}
			vcfIn.close();
			logout.close();
			System.out.println("Done parsing. " + nrLoaded + " variants loaded in total");

			System.out.println(availableRegions.size() + " regions in file " + vcfIn.getFileName());

			// combinatorial madness!
			if (availableRegions.isEmpty() || variantsInRegions.isEmpty()) {
				// print some fancy error message;
				System.out.println("Sorry. No work.");
			} else {
				//			// boot up threadpool
				System.out.println("Setting up threadpool with: " + options.getNrThreads() + " threads..");
				ExecutorService threadPool = Executors.newFixedThreadPool(options.getNrThreads());
				CompletionService<AssociationResult> jobHandler = new ExecutorCompletionService<AssociationResult>(threadPool);

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
						pb.set(returned);


					} catch (InterruptedException e) {
						e.printStackTrace();
					} catch (ExecutionException e) {
						e.printStackTrace();
					}
				}

				pvalOut.close();
				pb.close();
				threadPool.shutdown();

			}


		}
	}

	private VCFVariant getBestAssocForRegion(ArrayList<AssociationResult> assocResults,
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
				if (v.getId().equals(topResult.getSnp().getName())) {
					return v;
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
								currentLowestVariant = variant;
							}
							pvalout.writeln(assoc.toString());
							if (associationResults != null) {
								associationResults.add(assoc);
							}
							if (variant != null && iter == 0 && options.getMaxIter() > 1) {
								variants.add(variant);
							}
							nrTested++;
						}

					}
				}
				returned++;
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
