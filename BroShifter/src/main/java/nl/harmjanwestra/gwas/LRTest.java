package nl.harmjanwestra.gwas;

import nl.harmjanwestra.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;

/**
 * Created by hwestra on 11/7/15.
 */
public class LRTest {

//	boolean useR = false;

	LRTestOptions options;
	private int submitted;
	private int returned;
	private double highestLog10P;
	private VCFVariant currentLowestVariant;
	private double highestLog10PProbs;
	private int nrTested;
	private int maxNrOfResultsInBuffer = 10000;

	public LRTest(LRTestOptions options) throws IOException {
		this.options = options;
		test();
	}

	public void test() throws IOException {

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
		HashMap<String, Integer> sampleToIntGenotypes = new HashMap<String, Integer>();
		int ctr = 0;
		boolean[] genotypesWithCovariatesAndDiseaseStatus = new boolean[vcfSamples.size()];
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
				genotypesWithCovariatesAndDiseaseStatus[vcfSampleCtr] = true;
				samplesIntersect.add(sample);
				ctr++;
			}
			vcfSampleCtr++;
		}

		// load bed regions, if any
		ArrayList<Feature> bedRegions = null;
		if (options.getBedfile() != null) {
			BedFileReader reader = new BedFileReader();
			bedRegions = reader.readAsList(options.getBedfile());
			System.out.println(bedRegions.size() + " regions loaded from: " + options.getBedfile());
			for (Feature region : bedRegions) {
				System.out.println(region.toString());
			}
			Collections.sort(bedRegions, new FeatureComparator(false));
		}

		// load list of variants to limit upon if any
		HashSet<String> snpLimit = null;
		if (options.getSnpLimitFile() != null) {
			TextFile excl = new TextFile(options.getSnpLimitFile(), TextFile.R);
			snpLimit = (HashSet<String>) excl.readAsSet(0, TextFile.tab);
			excl.close();
		}


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

		System.out.println(sampleToIntGenotypes.size() + " samples with disease status, covariates and genotypes");
		if (sampleToIntGenotypes.size() == 0) {
			System.out.println("Problem with matching samples...");
		} else {
			double[] finalDiseaseStatus = new double[sampleToIntGenotypes.size()];
			double[][] finalCovariates = new double[sampleToIntGenotypes.size()][covariates.columns()];

			TextFile sampleListOut = new TextFile(options.getOutputdir() + "samplelist.txt", TextFile.W);
			System.out.println(options.getOutputdir() + "samplelist.txt");
			for (int i = 0; i < finalCovariates.length; i++) {
				sampleListOut.writeln(samplesIntersect.get(i));
			}
			sampleListOut.close();

			System.out.println("Final covariate array size: " + finalCovariates.length);
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
						finalCovariates[id][col] = covariates.getElement(sid, col);
					}
					nrTotal++;
				}
			}

			System.out.println(nrCases + " cases ");
			System.out.println(nrControls + " controls ");
			System.out.println(nrUnknown + " unknown ");
			System.out.println(nrTotal + " with covariates");

			System.out.println("max/min status: " + maxstatus + "/" + minstatus);

			// TODO: other VCF info scores...

			// keep a list of genotypes to condition on
			int iter = 0;
			ArrayList<Triple<double[][], boolean[], Integer>> conditional = new ArrayList<>();
			ArrayList<double[][]> conditionalDosages = new ArrayList<>();
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
					String logoutheader = "SNP\tChr\tPos\tImputationQual\tMAF\tOverlapOK\tMAFOk\tImpQualOK";
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

				data.close();
				data.open();

				currentLowestVariant = null;

				int variantCtr = 0;
				nrTested = 0;
				submitted = 0;
				returned = 0;

				int alleleOffsetGenotypes = 0;
				int alleleOffsetDosages = 0;
				highestLog10P = 0;
				highestLog10PProbs = 0;

				for (Triple<double[][], boolean[], Integer> c : conditional) {
					alleleOffsetGenotypes += c.getLeft().length;
				}
				for (double[][] c : conditionalDosages) {
					alleleOffsetDosages += c[0].length;
				}


				while ((iter == 0 && data.hasNext()) || (iter > 0 && variantCtr < variants.size())) {

					VCFVariant variant = null;
					if (iter == 0) {
						variant = data.nextLoadHeader();
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
										+ "\t" + imputationqualityscore
										+ "\t" + null
										+ "\t" + overlap
										+ "\t" + null
										+ "\t" + false);
							}
						} else {

							// throw into a thread
							// TODO: conditional on dosages is separate from conditional on genotypes.. ?
							LRTestTask task = new LRTestTask(variant,
									iter,
									genotypesWithCovariatesAndDiseaseStatus,
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
								clearQueue(logout, pvalout, iter, variants, jobHandler);
							}
						}
					} // end snplimit.contains(snp)
					variantCtr++;

					if (variantCtr % 1000 == 0) {
						String currentLowestDosagePValSNPId = null;
						if (currentLowestVariant != null) {
							currentLowestDosagePValSNPId = currentLowestVariant.getId();
						}
						System.out.println("Iteration: " + iter + "\tNr variants: " + variantCtr + " loaded.\tSubmitted to queue: " + submitted + "\tNr Tested: " + nrTested + "\tHighest P-val: " + highestLog10P + "\tSNP: " + currentLowestDosagePValSNPId);
					}
				} // end data.hasnext


				if (submitted % maxNrOfResultsInBuffer == 0) {
					clearQueue(logout, pvalout, iter, variants, jobHandler);
				}

				String currentLowestDosagePValSNPId = null;
				if (currentLowestVariant != null) {
					currentLowestDosagePValSNPId = currentLowestVariant.getId();
				}
				System.out.println("Iteration: " + iter + " is done. \tNr variants: " + variantCtr + "\tNr Tested: " + nrTested + "\tHighest P-val: " + highestLog10P + "\tBy variant: " + currentLowestDosagePValSNPId);

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
								Triple<double[][], boolean[], Integer> unfilteredGenotypeData = tasktmp.filterAndRecodeGenotypes(
										genotypesWithCovariatesAndDiseaseStatus,
										variant.getGenotypeAlleles(),
										variant.getAlleles().length,
										finalCovariates.length);
								conditional.add(unfilteredGenotypeData);

								double[][] dosages = variant.getImputedDosages();
								conditionalDosages.add(dosages);
								conditionalVariantIds.add(variant.getId());
							}
						}
					}
				} else if (currentLowestVariant != null) {
					VCFVariant variant = currentLowestVariant;
					Triple<double[][], boolean[], Integer> unfilteredGenotypeData = tasktmp.filterAndRecodeGenotypes(
							genotypesWithCovariatesAndDiseaseStatus,
							variant.getGenotypeAlleles(),
							variant.getAlleles().length,
							finalCovariates.length);
					conditional.add(unfilteredGenotypeData);


					double[][] dosages = variant.getImputedDosages();
					conditional.add(unfilteredGenotypeData);
					conditionalDosages.add(dosages);
					conditionalVariantIds.add(variant.getId());
					System.out.println("Iteration: " + iter + "\tNr variants: " + variantCtr + "\tHighest P-val: " + highestLog10PProbs + "\tBy variant: " + variant.getId());

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


	private void clearQueue(TextFile logout, TextFile pvalout,
	                        int iter, ArrayList<VCFVariant> variants,
	                        CompletionService<Triple<String, AssociationResult, VCFVariant>> jobHandler) throws IOException {
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
						if (logStr != null) {
							logout.writeln(logStr);
						}
						if (assoc != null) {
							double p = assoc.getLog10Pval();
							if (p > highestLog10P) {
								highestLog10P = p;
								currentLowestVariant = variant;
							}
							pvalout.writeln(assoc.toString());
							nrTested++;
						}
						if (variant != null && iter == 0 && options.getMaxIter() > 1) {
							variants.add(variant);
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
