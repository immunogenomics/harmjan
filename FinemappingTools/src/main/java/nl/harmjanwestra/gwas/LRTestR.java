package nl.harmjanwestra.gwas;

import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import org.rosuda.REngine.REXP;
import org.rosuda.REngine.REXPLogical;
import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;
import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;

/**
 * Created by hwestra on 6/1/15.
 */
public class LRTestR implements Callable<Boolean> {

	private final static boolean lingpipe = false;
	private final String vcf;
	private final String outputdir;
	private final String diseaseStatusFile;
	private final String covariateFile;
	private final HashSet<String> covariatesToInclude;
	private final String samplesToExclude;
	private final boolean imputationqualityfilter;
	private final double imputationqualitythreshold;
	private final int minNObservedAllele;
	private final int threadnum;
	private final String famfile;
	private final HashSet<String> snpLimit;
	boolean skipMakingPseudoControls = false;
	boolean runIterative = false;
	private ArrayList<Feature> bedRegions;

	public LRTestR(String vcf,
				   String outputdir,
				   String diseaseStatusFile,
				   String covariateFile,
				   HashSet<String> covariatesToInclude,
				   HashSet<String> snpLimit,
				   String samplesToExclude,
				   String famfile,
				   boolean imputationqualityfilter,
				   double imputationqualitythreshold,
				   int minNObservedAllele,
				   int threadnum) {
		this.vcf = vcf;
		this.outputdir = outputdir;
		this.diseaseStatusFile = diseaseStatusFile;
		this.covariateFile = covariateFile;
		this.covariatesToInclude = covariatesToInclude;
		this.samplesToExclude = samplesToExclude;
		this.imputationqualityfilter = imputationqualityfilter;
		this.imputationqualitythreshold = imputationqualitythreshold;
		this.minNObservedAllele = minNObservedAllele;
		this.threadnum = threadnum;
		this.famfile = famfile;
		this.snpLimit = snpLimit;
	}

	public void setRunIterative(boolean b, ArrayList<Feature> bedRegions) {
		this.bedRegions = bedRegions;

		this.runIterative = b;
	}

	@Override
	public Boolean call() throws Exception {
		if (runIterative) {
			return runIterative();
		} else {
			return run();
		}
	}

	public boolean runIterative() throws Exception {


		int maxiterations = 5;

		// load covariates
		Gpio.createDir(outputdir);

		VCFVariant lastVariant = null;

		System.out.println("Assoc: " + vcf);
		System.out.println("Covar: " + covariateFile);
		System.out.println("Disease: " + diseaseStatusFile);
		System.out.println("Out: " + outputdir);
		// index covariate samples to vcf samples
		HashSet<String> excludeTheseSamples = new HashSet<String>();
		if (samplesToExclude != null) {
			TextFile excl = new TextFile(samplesToExclude, TextFile.R);
			excludeTheseSamples = (HashSet<String>) excl.readAsSet(0, TextFile.tab);
			excl.close();
		}
		HashMap<String, Integer> diseaseStatus = new HashMap<String, Integer>();

		TextFile tf = new TextFile(diseaseStatusFile, TextFile.R);
		String[] elems = tf.readLineElems(Strings.tab);
		while (elems != null) {

			String sample = elems[0];
			String statusStr = elems[1];

			Integer status = Integer.parseInt(statusStr);
			diseaseStatus.put(sample, status);

			elems = tf.readLineElems(Strings.tab);
		}
		tf.close();

		System.out.println(diseaseStatus.size() + " disease status samples loaded");

		//
		VCFGenotypeData data = new VCFGenotypeData(vcf);
//		VCFFunctions v = new VCFFunctions();
		ArrayList<String> vcfSamples = data.getSamples();

		HashSet<String> samplesWithDiseaseStatus = new HashSet<String>();
		for (String sample : vcfSamples) {
			if (diseaseStatus.containsKey(sample)) {
				if (excludeTheseSamples == null || !excludeTheseSamples.contains(sample)) {
					samplesWithDiseaseStatus.add(sample);
				}
			}
		}


		// assume rows are samples and columns are covariates
		// this loads only samples that also have a disease status loaded..
		DoubleMatrixDataset<String, String> covariates = DoubleMatrixDataset.loadSubsetOfTextDoubleData(covariateFile, "\t", samplesWithDiseaseStatus, covariatesToInclude);
		System.out.println("Covariate matrix: " + covariates.rows() + " samples " + covariates.columns() + " covariates");

		LinkedHashMap<String, Integer> covariateSampleHash = covariates.getHashRows();

		HashMap<String, Integer> sampleToIntGenotypes = new HashMap<String, Integer>();
		int ctr = 0;
		boolean[] genotypesWithCovariatesAndDiseaseStatus = new boolean[vcfSamples.size()];
		int vcfSampleCtr = 0;

		HashSet<String> alternatecovariateSamples = new HashSet<String>();
		Set<String> keys = covariateSampleHash.keySet();
		HashMap<String, String> altToSample = new HashMap<String, String>();
		for (String sample : keys) {
			alternatecovariateSamples.add(sample + "_" + sample);
			altToSample.put(sample, sample + "_" + sample);
		}


		ArrayList<String> samplesIntersect = new ArrayList<String>();
		for (String sample : vcfSamples) {
			if (covariateSampleHash.containsKey(sample) || alternatecovariateSamples.contains(sample)) {
				sampleToIntGenotypes.put(sample, ctr);
				genotypesWithCovariatesAndDiseaseStatus[vcfSampleCtr] = true;
				samplesIntersect.add(sample);
				ctr++;
			}
			vcfSampleCtr++;
		}


		String[] covariateColNames = covariates.getColObjects().toArray(new String[0]);

		System.out.println(sampleToIntGenotypes.size() + " samples with disease status, covariates and genotypes");


		if (sampleToIntGenotypes.size() == 0) {
			System.out.println("Problem with matching samples...");
		} else {
			try {
				double[] finalDiseaseStatus = null; // new double[sampleToIntGenotypes.size()];
				double[][] finalCovariates = null;  // new double[sampleToIntGenotypes.size()][covariates.columns()];

				Integer[][] sampleParents = new Integer[2][sampleToIntGenotypes.size()];
				ArrayList<Integer> kidsInTrios = new ArrayList<Integer>(); // lookup hash to check for each sample whether they are in a trio. Integers are relative to genotype
				HashMap<Integer, Integer> kidToPseudoCC = new HashMap<Integer, Integer>();

				HashSet<Integer> samplesToRemove = new HashSet<Integer>(); // relative to new ordering
				boolean arraysInit = false;
				if (famfile != null) {

					HashSet<Integer> samplesAlreadyAssignedToTrio = new HashSet<Integer>();
					HashSet<Integer> samplesAlreadyAssignedAsKid = new HashSet<Integer>();
					// get the trios (only trios where kid is a case)
					ArrayList<Pair<String, Triple<String, String, String>>> famData = getTrios(famfile); // format: familyname<kid, mom, dad>

					HashSet<String> familiesAlreadyUsed = new HashSet<String>();

					if (famData.size() > 0) {
						int nrCompleteTrios = 0;
						int nrIncompleteTrios = 0;
						int nrRelatedKids = 0;

						for (Pair<String, Triple<String, String, String>> family : famData) {
							String familyName = family.getLeft();
							Triple<String, String, String> familyMembers = family.getRight();

							String kidName = familyMembers.getLeft();
							Integer kidId = sampleToIntGenotypes.get(kidName);
							if (kidId == null) {
								kidId = sampleToIntGenotypes.get(altToSample.get(kidName));
							}

							String momName = familyMembers.getMiddle();
							Integer momId = sampleToIntGenotypes.get(momName);
							if (momId == null) {
								momId = sampleToIntGenotypes.get(altToSample.get(momName));
							}

							String dadName = familyMembers.getRight();
							Integer dadId = sampleToIntGenotypes.get(dadName);
							if (dadId == null) {
								dadId = sampleToIntGenotypes.get(altToSample.get(dadName));
							}

							if (kidId != null) {
								sampleParents[0][kidId] = dadId;
								sampleParents[1][kidId] = momId;
								if (dadId != null) {
									samplesToRemove.add(dadId);
								}
								if (momId != null) {
									samplesToRemove.add(momId);
								}
								// pick one case per family
								if (familiesAlreadyUsed.contains(familyName)) {
									samplesToRemove.add(kidId);
									nrRelatedKids++;
								} else if (dadId != null && momId != null && !familiesAlreadyUsed.contains(familyName)) {
									kidsInTrios.add(kidId);

									if (samplesAlreadyAssignedToTrio.contains(kidId)) {
										System.out.println("kid sample " + kidId + " for " + kidName + " already assigned?");
									}
									samplesAlreadyAssignedToTrio.add(kidId);
									samplesAlreadyAssignedAsKid.add(kidId);

									if (samplesAlreadyAssignedAsKid.contains(dadId)) {
										System.out.println("Father: " + dadId + " for " + dadName + " already assigned as kid.");
									}
									if (samplesAlreadyAssignedAsKid.contains(momId)) {
										System.out.println("Mother: " + momId + " for " + momName + " already assigned as kid.");
									}
									kidToPseudoCC.put(kidId, sampleToIntGenotypes.size() + nrCompleteTrios);
									familiesAlreadyUsed.add(familyName);
									nrCompleteTrios++;
								} else {
									// incomplete trio. mom and dad already removed, also remove kid, because, why not.
									samplesToRemove.add(kidId);
									nrIncompleteTrios++;
								}
							}
						}

						System.out.println(nrCompleteTrios + " full trios found in genotypes.. " + nrIncompleteTrios + " incomplete trios.." + nrRelatedKids + " kids related.." + samplesToRemove.size() + " samples to remove");
						// add vector to finalCovariates and finalDiseaseStatus
						if (!skipMakingPseudoControls) {
							if (nrCompleteTrios > 0) {

								finalDiseaseStatus = new double[sampleToIntGenotypes.size() + nrCompleteTrios];
								finalCovariates = new double[sampleToIntGenotypes.size() + nrCompleteTrios][covariates.columns()];

								// set pseudocontrol disease status
								for (int i = sampleToIntGenotypes.size(); i < finalCovariates.length; i++) {
									finalDiseaseStatus[i] = 0; // 1 == unaffected, 2 == affected --> will be remapped to 0/1
								}

								// copy covariate values from kids
								for (int i = 0; i < kidsInTrios.size(); i++) {
									Integer kidId = kidsInTrios.get(i);
									Integer pseudoCC = kidToPseudoCC.get(kidId);
									if (pseudoCC != null) {
										String kidName = samplesIntersect.get(kidId);
										samplesIntersect.add(kidName + "-PseudoControl");
										Integer covariateSampleId = covariateSampleHash.get(kidName);
										for (int col = 0; col < covariates.columns(); col++) {
											finalCovariates[pseudoCC][col] = covariates.getElement(covariateSampleId, col);
										}
									}
								}
								arraysInit = true;
							}
						}
					}
				}

				if (!arraysInit) {
					finalDiseaseStatus = new double[sampleToIntGenotypes.size()];
					finalCovariates = new double[sampleToIntGenotypes.size()][covariates.columns()];
				}


				TextFile sampleListOut = new TextFile(outputdir + "samplelist.txt", TextFile.W);
				System.out.println(outputdir + "samplelist.txt");

				for (int i = 0; i < finalCovariates.length; i++) {
					sampleListOut.writeln(samplesIntersect.get(i) + "\t" + samplesToRemove.contains(i));
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
				// make this load genotype probabilities as well...
				// TODO: genotype probabilities
				// TODO: other VCF info scores...
				// TODO: exclude samples
				// 1       1955373 rs4648791       C       T       .       PASS    AR2=0;DR2=0;AF=0.264    GT:DS:GP        0|0:0.526:0.543,0.388,0.069
			/*
			 ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated Allele Frequencies">
             ##INFO=<ID=AR2,Number=1,Type=Float,Description="Allelic R-Squared: estimated correlation between most probable ALT dose and true ALT dose">
             ##INFO=<ID=DR2,Number=1,Type=Float,Description="Dosage R-Squared: estimated correlation between estimated ALT dose [P(RA) + 2*P(AA)] and true ALT dose">
             ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
             ##FORMAT=<ID=DS,Number=1,Type=Float,Description="estimated ALT dose [P(RA) + P(AA)]">
             ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Estimated Genotype Probability">
             */

				RConnection rConnection = new RConnection();
				System.out.println("R server found: " + rConnection.getServerVersion());
				rConnection.assign("status", finalDiseaseStatus);
				System.out.println("Loaded disease status");


				assignAsRMatrix(rConnection, finalCovariates, "covariates", false);
				rConnection.assign("covariatenames", covariateColNames);
				rConnection.voidEval("colnames(covariates) <- covariatenames");
				System.out.println("Loaded covariates");


				//				rConnection.voidEval("print(head(status))");

//				rConnection.voidEval("print('" + threadnum + " thread environment set')");
				System.out.println("Set thread environment.. " + threadnum);

				int varctr = 0;
				TextFile logout = new TextFile(outputdir + "log.txt", TextFile.W);

				boolean previousVariantHasMissingGenotyes = false;
				boolean firstvariant = true;

				boolean continueTesting = true;

//				while (data.hasNext() && continueTesting) {
				ArrayList<double[][]> allGenotypes = new ArrayList<double[][]>();
				ArrayList<String> variantNames = new ArrayList<String>();
				HashMap<String, Integer> genotypeIndex = new HashMap<String, Integer>();


				ArrayList<Integer> variantPositions = new ArrayList<Integer>();

				HashMap<String, Double> variantMafs = new HashMap<String, Double>();
				while (data.hasNext()) {

					VCFVariant variant = data.next();

					String name = variant.getId();
					if (snpLimit == null || snpLimit.contains(name)) {

						lastVariant = variant; // debugging purposes :(

						Double imputationqualityscore = variant.getImputationQualityScore();
						boolean testvariant = false;
						if (imputationqualityfilter) {
							if (imputationqualityscore != null && imputationqualityscore >= imputationqualitythreshold) {
								testvariant = true;
							} else if (imputationqualityscore == null) {
								System.err.println("No imputaton quality score for variant: " + variant.getChr() + "-" + variant.getPos() + "-" + variant.getId());
								System.err.println("In path: " + vcf);
								logout.writeln("Imputation quality score below threshold:\t" + imputationqualityscore + "\t" + variant.getChr() + "-" + variant.getPos() + "-" + variant.getId());
							}
						} else {
							testvariant = true;
						}


						if (testvariant) {

							// recode the genotypes to the same ordering as the covariate table
							// generate pseudocontrol genotypes
							Triple<double[][], Double, Integer> recodedGenotypes = recodeGenotypes(
									variant.getChr() + "-" + variant.getPos() + "-" + variant.getId(),
									genotypesWithCovariatesAndDiseaseStatus,
									variant.getGenotypeAlleles(),
									variant.getAlleles().length,
									kidsInTrios,
									sampleParents,
									kidToPseudoCC,
									samplesToRemove,
									finalCovariates.length);


							Integer sampleswithgenotypes = recodedGenotypes.getRight();
							double maf = recodedGenotypes.getMiddle();

							double mafthresholdD = minNObservedAllele / (sampleswithgenotypes * 2);
							if (maf >= mafthresholdD) {
								double[][] genotypes = recodedGenotypes.getLeft();
								allGenotypes.add(genotypes);
								String variantname = variant.getChr().toString() + ":" + variant.getPos() + "-" + variant.getId();
								genotypeIndex.put(variantname, genotypeIndex.size());
								variantNames.add(variantname);
								variantPositions.add(variant.getPos());
								variantMafs.put(variantname, maf);
							} else {
								logout.writeln("variant skipped maf: " + maf + " below threshold " +
										variant.getId()
										+ "\t" + variant.getChr()
										+ "\t" + variant.getPos()
										+ "\t" + variant.getMAF()
								);
							}
						} else {
							logout.writeln("variant skipped rsq: " + imputationqualityscore + " below threshold " +
									variant.getId()
									+ "\t" + variant.getChr()
									+ "\t" + variant.getPos()
									+ "\t" + variant.getMAF()
							);
						}
					}
					varctr++;
					if (varctr % 1000 == 0) {

						System.out.print("thread:\t" + threadnum + "\tnrVars: " + varctr + "\n");

					}
				}
				logout.close();

				System.out.println("thread:\t" + threadnum + "\tAll data loaded: " + allGenotypes.size() + " variants loaded.");

				// now testNormal the variants
				System.out.println("thread:\t" + threadnum + "\tAbout to testNormal: " + bedRegions.size() + " regions");
				for (int r = 0; r < bedRegions.size(); r++) {
					// get a list of indexes of variants overlapping each bedRegion
					Pair<ArrayList<double[][]>, ArrayList<String>> overlappingVariants = overlapVariantsWithBedRegion(bedRegions.get(r),
							variantPositions,
							variantNames,
							allGenotypes);
					ArrayList<double[][]> regionGenotypes = overlappingVariants.getLeft();
					ArrayList<String> regionVariantNames = overlappingVariants.getRight();
					System.out.println("thread:\t" + threadnum + "\t" + bedRegions.get(r).toString() + "\t" + regionGenotypes.size() + " variants");
					int iteration = 0;
					double lowestP = 1;
					ArrayList<Integer> lowestPValGenotypes = new ArrayList<Integer>();
					boolean continuetesting = true;
					while (continuetesting) {
						System.out.println("Entering main loop");
						lowestP = 1;
						Integer lowestPIndex = null;

						TextFile pvalout = new TextFile(outputdir + bedRegions.get(r).toString() + "-gwas-" + iteration + ".txt", TextFile.W);
						System.out.println("thread:\t" + threadnum + "\twriting to: " + outputdir + bedRegions.get(r).toString() + "-gwas-" + iteration + ".txt");
						String modelStr = "#model: status ~ snp + covar ";
						if (!lowestPValGenotypes.isEmpty()) {
//							String modelStr = "#model: status ~ snp + covar ";
							for (int l = 0; l < lowestPValGenotypes.size(); l++) {
								// assign genotypes
								Integer index = lowestPValGenotypes.get(l);
								double[][] genotypes2 = regionGenotypes.get(index);
								String covarname = "genotypes" + l;
								assignAsRMatrix(rConnection, genotypes2, covarname, true);
								rConnection.voidEval("notnagenotypes" + l + "<-apply(!is.na(" + covarname + "),1,any)");
//								if (l == 0) {
								modelStr += " + " + regionVariantNames.get(index);
//								}
							}
							pvalout.writeln(modelStr);
						}

						if (finalCovariates[0].length == 0) {
							pvalout.writeln("VariantID" +
									"\tN" +
									"\tMAF" +
									"\tDevianceNull" +
									"\tDfNull" +
									"\tDevianceGeno" +
									"\tDfAlt" +
									"\tPvalIntersect" +
									"\t-Log10(pvalIntersect)" +
									"\tPval" +
									"\t-Log10(pval)");
						} else {
							pvalout.writeln("VariantID" +
									"\tN" +
									"\tMAF" +
									"\tDevianceNull" +
									"\tDfNull" +
									"\tDevianceGeno" +
									"\tDfAlt" +
									"\tBeta(Genotype)" +
									"\tSE(Genotype)" +
									"\tOR" +
									"\tOR-Hi" +
									"\tOR-Lo" +
									"\tPval" +
									"\t-Log10(pval)");
						}


						for (int variantId = 0; variantId < regionGenotypes.size(); variantId++) {
							double[][] genotypes = regionGenotypes.get(variantId);
							if (finalCovariates[0].length == 0) {
								// no need to run the null model as there are no covariates
								assignAsRMatrix(rConnection, genotypes, "genotypes", true);

								String model = "";
								if (!lowestPValGenotypes.isEmpty()) {
									model = "glm_geno <- glm(status[notnagenotypes] ~ genotypes[notnagenotypes,]";
									for (int l = 0; l < lowestPValGenotypes.size(); l++) {
										// assign genotypes
										String covarname = "genotypes" + l;
										model += " + " + covarname + "[notnagenotypes,]";
									}
									model += ", family=binomial(logit))";
								} else {
									model = "glm_geno <- glm(status[notnagenotypes] ~ genotypes[notnagenotypes,], family=binomial(logit))";
								}

								rConnection.voidEval(model);
								rConnection.voidEval("glm_summary <- summary(glm_geno)");

								REXP rexp = rConnection.eval("glm_summary$aliased[\"genotypes[notnagenotypes, ]\"]");
								boolean aliased = false;
								if (rexp.isLogical()) {
									boolean[] bools = ((REXPLogical) rexp).isTRUE();
									if (bools[0]) {
										aliased = true;
									}
								}

								if (!aliased) {
									Double pgeno = rConnection.eval("glm_summarycoefficients[\"genotypes\",4]").asDouble();
									Double pintersect = rConnection.eval("pintersect <- glm_summarycoefficients[1,4]").asDouble();


									String variantName = regionVariantNames.get(variantId);
									Double maf = variantMafs.get(variantName);
									int nrT = rConnection.eval("length(which(notnagenotypes == TRUE))").asInteger();
									String outstr = regionVariantNames.get(variantId)
											+ "\t" + nrT
											+ "\t" + maf
											+ "\t" + pintersect
											+ "\t" + (-Math.log10(pintersect))
											+ "\t" + pgeno
											+ "\t" + (-Math.log10(pgeno));
									pvalout.writeln(outstr);

									if (pgeno != null) {
										if (Double.isInfinite(pgeno) || Double.isNaN(pgeno) || pgeno < 0 || pgeno > 1) {

										} else {
											if (pgeno < lowestP) {
												lowestP = pgeno;
												lowestPIndex = variantId;
											}
										}
									}
								}
							} else {

								// avoid running the null model very often...
								// assign genotypes to a matrix
								assignAsRMatrix(rConnection, genotypes, "genotypes", true);

								if (!lowestPValGenotypes.isEmpty()) {

									rConnection.voidEval("notnagenotypes<-apply(!is.na(genotypes),1,any)");
									String model = "glm_geno <- glm(status[notnagenotypes] ~ genotypes[notnagenotypes,] + covariates[notnagenotypes,]";
									String nullmodel = "glm_null <- glm(status[notnagenotypes] ~ covariates[notnagenotypes,]";
									String natest = "notnagenotypes <- (notnagenotypes ";

									for (int l = 0; l < lowestPValGenotypes.size(); l++) {
										// assign genotypes
										String covarname = "genotypes" + l;
										String additionalCovar = " + " + covarname + "[notnagenotypes,]";
										model += additionalCovar;
										nullmodel += additionalCovar;
										natest += " & notnagenotypes" + l;
									}

									natest += ")";
									rConnection.voidEval(natest);
									model += ", family=binomial(logit))";
									nullmodel += ", family=binomial(logit))";

									rConnection.voidEval(model);
									rConnection.voidEval(nullmodel);
//
//									System.out.println("model: " + model);
//									System.out.println("nullmodel: " + nullmodel);
//									System.out.println("natest: " + natest);
//
//
//									System.out.println("samples left: " + nrT);


								} else {
									rConnection.voidEval("notnagenotypes<-apply(!is.na(genotypes),1,any)");
									rConnection.voidEval("glm_null <- glm(status[notnagenotypes] ~ covariates[notnagenotypes,], family=binomial(logit))");
									rConnection.voidEval("glm_geno <- glm(status[notnagenotypes] ~ genotypes[notnagenotypes,] + covariates[notnagenotypes,], family=binomial(logit))");
								}

								rConnection.voidEval("nullSummary<-summary(glm_null)");
								rConnection.voidEval("deviance_null <- summary(glm_null)$deviance");
								rConnection.voidEval("nulldf <- summary(glm_null)$df[1]");

								double devianceNull = rConnection.eval("deviance_null").asDouble();
								rConnection.voidEval("glm_summary <- summary(glm_geno)");

//								if (iteration > 0) {
//									String[] output = rConnection.eval("glm_summary").asStrings();
//									System.out.println("GLM summary: ");
//									for (String s : output) {
//										System.out.println(s);
//									}
//									System.out.println();
//
//
//									output = rConnection.eval("nullSummary").asStrings();
//									System.out.println("Null summary: ");
//									for (String s : output) {
//										System.out.println(s);
//									}
//									System.out.println();
//
//								}

								boolean aliased = false;
								if (genotypes.length == 1) {
									REXP rexp = rConnection.eval("glm_summary$aliased[\"genotypes[notnagenotypes, ]\"]");
									if (rexp.isLogical()) {
										boolean[] bools = ((REXPLogical) rexp).isTRUE();
										if (bools[0]) {
											aliased = true;
										}
									}
								} else {
									rConnection.voidEval("print(glm_summary$aliased)");
									boolean[] aliasPerAllele = new boolean[genotypes.length];
									for (int i = 1; i < genotypes.length + 1; i++) {
										REXP rexp = rConnection.eval("glm_summary$aliased[\"genotypes[notnagenotypes, ]" + i + "\"]");
										if (rexp.isLogical()) {
											boolean[] bools = ((REXPLogical) rexp).isTRUE();
											if (bools[0]) {
												aliasPerAllele[i - 1] = true;
											}
										} else {
											System.out.println("not logical");
										}
									}
									aliased = true;
									for (int i = 0; i < aliasPerAllele.length; i++) {
										if (!aliasPerAllele[i]) {
											aliased = false;
										}
										System.out.println(i + "\t" + aliasPerAllele[i]);
									}
								}

								if (!aliased) { // at least one genotype vector should remain...
									rConnection.voidEval("deviance_geno <- glm_summary$deviance");
									Double deviance_geno = rConnection.eval("deviance_geno").asDouble();

									rConnection.voidEval("deltaDeviance <- deviance_null - deviance_geno");
									rConnection.voidEval("altdf <- glm_summary$df[1]");
									double altdf = rConnection.eval("altdf").asDouble();
									double nulldf = rConnection.eval("nulldf").asDouble();
									rConnection.voidEval("dfdiff <- altdf - nulldf");

									Double pval = rConnection.eval("pchisq(deltaDeviance, df = dfdiff, lower.tail = FALSE, log.p = FALSE)").asDouble();

									if (devianceNull - deviance_geno < 0.0001) {
										pval = 1d;
									}
									double log10p = Math.abs((-Math.log10(pval)));

									int nrT = rConnection.eval("length(which(notnagenotypes == TRUE))").asInteger();

									if (genotypes.length == 1) {
										rConnection.voidEval("genocoeff <- as.vector(glm_summary$coefficients[2,])");
										Double beta = rConnection.eval("genocoeff[1]").asDouble();
										Double OR = Math.exp(beta); // rConnection.eval("exp(beta)").asDouble();
										Double betase = rConnection.eval("genocoeff[2]").asDouble();
										Double orLow = Math.exp(beta - 1.96 * betase);// rConnection.eval("exp(beta - 1.96 * betase)").asDouble();
										Double orHigh = Math.exp(beta + 1.96 * betase);//rConnection.eval("exp(beta + 1.96 * betase)").asDouble();

										String variantName = regionVariantNames.get(variantId);
										Double maf = variantMafs.get(variantName);

										String outstr = variantName
												+ "\t" + nrT
												+ "\t" + maf
												+ "\t" + devianceNull
												+ "\t" + nulldf
												+ "\t" + deviance_geno
												+ "\t" + altdf
												+ "\t" + beta
												+ "\t" + betase
												+ "\t" + OR
												+ "\t" + orLow
												+ "\t" + orHigh
												+ "\t" + pval
												+ "\t" + log10p;
//									System.out.println(outstr);
										pvalout.writeln(outstr);
//										rConnection.voidEval("save.image(path=\"" + outputdir + regionVariantNames.get(variantId) + ".rData\")");
//										System.out.println("Saved to: " + outputdir + regionVariantNames.get(variantId) + ".rData");
//										System.exit(-1);
									} else {
										String outstr = variantNames.get(variantId)
												+ "\t" + nrT
												+ "\t" + devianceNull
												+ "\t" + nulldf
												+ "\t" + deviance_geno
												+ "\t" + altdf
												+ "\tnull"
												+ "\tnull"
												+ "\tnull"
												+ "\tnull"
												+ "\tnull"
												+ "\t" + pval
												+ "\t" + log10p;
//									System.out.println(outstr);
										pvalout.writeln(outstr);


									}


									if (pval != null) {
										if (Double.isInfinite(pval) || Double.isNaN(pval) || pval < 0 || pval > 1) {
										} else {
											if (pval < lowestP) {
												lowestP = pval;
												lowestPIndex = variantId;
											}
										}
									}
								}
							}
							if (variantId % 50 == 0) {

								System.out.println("thread:\t" + threadnum
										+ "\titer: " + iteration
										+ "\tnrVars processed: " + variantId + "/" + regionGenotypes.size()
										+ "\t lowest P: " + lowestP
										+ "\tmodel: " + modelStr
								);

							}
						}

						if (lowestPIndex != null) {

							// calculate pairwise LD
							determinePairWiseLd(lowestPIndex,
									variantNames,
									regionGenotypes,
									outputdir + bedRegions.get(r).toString() + "-ld-" + iteration + ".txt");

							HashSet<Integer> q = new HashSet<Integer>();
							for (Integer z : lowestPValGenotypes) {
								q.add(z);
							}
							if (!q.contains(lowestPIndex)) {
								lowestPValGenotypes.add(lowestPIndex);
							}
							if (lowestP > 0.05) {
								continuetesting = false;
							}
						} else {
							continuetesting = false;
						}

						if (lowestPIndex != null) {
							System.out.println("thread: " + threadnum + "\titer: " + iteration + "\tlowestP: " + lowestP + "\t" + regionVariantNames.get(lowestPIndex));
						} else {
							System.out.println("thread: " + threadnum + "\titer: " + iteration + "\tlowestP: " + lowestP + "\t" + null);
						}

						iteration++;

						if (iteration >= maxiterations) {
							continuetesting = false;
						}

						pvalout.close();
					}

				}


				rConnection.close();


			} catch (RserveException ex) {
				System.err.println(ex.getMessage());
				System.err.println("Could not connect to RServe");
				System.err.println("ERRORRRRRR: " + vcf);
				System.err.println(lastVariant.getChr() + ":" + lastVariant.getPos() + "-" + lastVariant.getId());
				System.exit(0);
				return false;

			} catch (REngineException e) {
				e.printStackTrace();
				System.err.println("ERRORRRRRR: " + vcf);
				System.err.println(lastVariant.getChr() + ":" + lastVariant.getPos() + "-" + lastVariant.getId());
				System.exit(0);
				return false;
			} catch (REXPMismatchException e) {
				e.printStackTrace();
				System.err.println("ERRORRRRRR: " + vcf);
				System.err.println(lastVariant.getChr() + ":" + lastVariant.getPos() + "-" + lastVariant.getId());
				System.exit(0);
				return false;
			} catch (Exception e) {
				System.err.println("ERRORRRRRR: " + vcf);

				e.printStackTrace();
				System.err.println(lastVariant.getChr() + ":" + lastVariant.getPos() + "-" + lastVariant.getId());

				System.exit(0);
				return false;
			}

		}

		// load genotypes


		return true;
	}


	// thanks @ Lude Franke :) who wrote the original code below
	private void determinePairWiseLd(Integer lowestPIndex, ArrayList<String> variantNames, ArrayList<double[][]> regionGenotypes, String s) throws IOException {

		TextFile tf = new TextFile(s, TextFile.W);
		tf.writeln("Var1\tVar2\tR-squared\tD-prime");
		boolean print = false;
		double[][] genotypes1 = regionGenotypes.get(lowestPIndex); // [alleles][samples] (counts copies of genotypes...
		if (genotypes1.length < 2) {
			String name1 = variantNames.get(lowestPIndex);
			for (int q = 0; q < regionGenotypes.size(); q++) {
				String name2 = variantNames.get(q);
				double[][] genotypes2 = regionGenotypes.get(q);

				if (genotypes2.length < 2) {
					//Get genotypes:
					int[][] genotypes = new int[3][3];
					int nrCalledGenotypes = 0;

					for (int ind = 0; ind < genotypes1[0].length; ind++) {
//						if (individualsToInclude == INCLUDE_CASES_AND_CONTROLS
//								|| (genotypeData.getIsCase()[ind] != null && genotypeData.getIsCase()[ind] && individualsToInclude == INCLUDE_CASES)
//								|| (genotypeData.getIsCase()[ind] != null && !genotypeData.getIsCase()[ind] && individualsToInclude == INCLUDE_CONTROLS)) {
						double genotypeX = genotypes1[0][ind];
						double genotypeY = genotypes2[0][ind];
						if (!Double.isNaN(genotypeX) && !Double.isNaN(genotypeY)) {

							genotypes[(int) genotypeX][(int) genotypeY]++;
							nrCalledGenotypes++;
						}
//						}
					}
					//System.out.println("NrCalledGenotypes:\t" + nrCalledGenotypes);

					//Determine genotype frequencies:
					double[][] genotypesFreq = new double[3][3];
					for (int x = 0; x < 3; x++) {
						for (int y = 0; y < 3; y++) {
							genotypesFreq[x][y] = (double) genotypes[x][y] / (double) nrCalledGenotypes;
							if (print) {
								System.out.print(genotypes[x][y] + "\t");
							}
						}
						if (print) {
							System.out.println("");
						}
					}

					//Determine alle frequencies:
					double[][] alleleFreq = new double[2][2];
					//SNP X:
					alleleFreq[0][0] = (genotypesFreq[0][0] + genotypesFreq[0][1] + genotypesFreq[0][2]) + (genotypesFreq[1][0] + genotypesFreq[1][1] + genotypesFreq[1][2]) / 2d;
					alleleFreq[0][1] = (genotypesFreq[2][0] + genotypesFreq[2][1] + genotypesFreq[2][2]) + (genotypesFreq[1][0] + genotypesFreq[1][1] + genotypesFreq[1][2]) / 2d;
					//SNP Y:
					alleleFreq[1][0] = (genotypesFreq[0][0] + genotypesFreq[1][0] + genotypesFreq[2][0]) + (genotypesFreq[0][1] + genotypesFreq[1][1] + genotypesFreq[2][1]) / 2d;
					alleleFreq[1][1] = (genotypesFreq[0][2] + genotypesFreq[1][2] + genotypesFreq[2][2]) + (genotypesFreq[0][1] + genotypesFreq[1][1] + genotypesFreq[2][1]) / 2d;

					if (print) {
						System.out.println("Allele freq:");
						System.out.println(alleleFreq[0][0]);
						System.out.println(alleleFreq[0][1]);
						System.out.println(alleleFreq[1][0]);
						System.out.println(alleleFreq[1][1]);
					}

					//Precalculate triangles of non-double heterozygote:
					double[][] genotypesTriangleFreq = new double[3][3];
					genotypesTriangleFreq[0][0] = 2d * genotypesFreq[0][0] + genotypesFreq[1][0] + genotypesFreq[0][1];
					genotypesTriangleFreq[2][0] = 2d * genotypesFreq[2][0] + genotypesFreq[1][0] + genotypesFreq[2][1];
					genotypesTriangleFreq[0][2] = 2d * genotypesFreq[0][2] + genotypesFreq[0][1] + genotypesFreq[1][2];
					genotypesTriangleFreq[2][2] = 2d * genotypesFreq[2][2] + genotypesFreq[1][2] + genotypesFreq[2][1];
					if (print) {
						System.out.println("Triangle freq:");
						System.out.println(genotypesTriangleFreq[0][0]);
						System.out.println(genotypesTriangleFreq[0][2]);
						System.out.println(genotypesTriangleFreq[2][0]);
						System.out.println(genotypesTriangleFreq[2][2]);
					}

					//Calculate expected genotypes, assuming equilibrium, take this as start:
					double h11 = alleleFreq[0][0] * alleleFreq[1][0];
					double h12 = alleleFreq[0][0] * alleleFreq[1][1];
					double h21 = alleleFreq[0][1] * alleleFreq[1][0];
					double h22 = alleleFreq[0][1] * alleleFreq[1][1];

					//Calculate the frequency of the two double heterozygotes:
					double x12y12 = h11 * h22 / (h11 * h22 + h12 * h21) * genotypesFreq[1][1];
					double x12y21 = h12 * h21 / (h11 * h22 + h12 * h21) * genotypesFreq[1][1];

					if (print) {
						System.out.println(h11 + "\t" + h12 + "\t" + h21 + "\t" + h22 + "\t\t" + x12y12 + "\t" + x12y21);
					}

					//Perform iterations using EM algorithm:
					for (int itr = 0; itr < 25; itr++) {

						h11 = (x12y12 + genotypesTriangleFreq[0][0]) / 2;
						h12 = (x12y21 + genotypesTriangleFreq[0][2]) / 2;
						h21 = (x12y21 + genotypesTriangleFreq[2][0]) / 2;
						h22 = (x12y12 + genotypesTriangleFreq[2][2]) / 2;

						x12y12 = h11 * h22 / (h11 * h22 + h12 * h21) * genotypesFreq[1][1];
						x12y21 = h12 * h21 / (h11 * h22 + h12 * h21) * genotypesFreq[1][1];

						if (print) {
							System.out.println(h11 + "\t" + h12 + "\t" + h21 + "\t" + h22 + "\t\t" + x12y12 + "\t" + x12y21);
						}

					}

					double d = h11 - (alleleFreq[0][0] * alleleFreq[1][0]);

//        if (returnType == RETURN_R_SQUARED) {
					double rSquared = d * d / (alleleFreq[0][0] * alleleFreq[0][1] * alleleFreq[1][0] * alleleFreq[1][1]);
					//return rSquared;
//        } else {
					double dMax = 0;
					if (d < 0) {
						double a = alleleFreq[0][1] * alleleFreq[1][1];
						if (alleleFreq[0][0] > alleleFreq[1][0]) {
							a = alleleFreq[0][0] * alleleFreq[1][0];
						}
						double b = alleleFreq[0][0] * alleleFreq[1][0];
						if (alleleFreq[0][0] > alleleFreq[1][0]) {
							b = alleleFreq[0][1] * alleleFreq[1][1];
						}
						dMax = Math.min(a, b);
					} else {
						double a = alleleFreq[0][1] * alleleFreq[1][0];
						if (alleleFreq[0][0] > alleleFreq[1][0]) {
							a = alleleFreq[0][0] * alleleFreq[1][1];
						}
						double b = alleleFreq[0][0] * alleleFreq[1][1];
						if (alleleFreq[0][0] > alleleFreq[1][0]) {
							b = alleleFreq[0][1] * alleleFreq[1][0];
						}
						dMax = Math.min(a, b);
					}
					double dPrime = Math.abs(d / dMax);
					dPrime = Math.min(1, dPrime);

					tf.writeln(name1 + "\t" + name2 + "\t" + rSquared + "\t" + dPrime);
				}
			}



		/*
		 if (dPrime>1.01) {
         System.out.println("");
         System.out.println(genotypes[0][0] + "\t" + genotypes[0][1] + "\t" + genotypes[0][2]);
         System.out.println(genotypes[1][0] + "\t" + genotypes[1][1] + "\t" + genotypes[1][2]);
         System.out.println(genotypes[2][0] + "\t" + genotypes[2][1] + "\t" + genotypes[2][2]);
         System.out.println(alleleFreq[0][0] + "\t" + alleleFreq[0][1] + "\t" + alleleFreq[1][0] + "\t" + alleleFreq[1][1]);
         System.out.println(h11 + "\t" + h12 + "\t" + h21 + "\t" + h22);
         System.out.println(d + "\t" + dMax + "\t" + dPrime);
         }
         */
//			return new Pair<Double, Double>(Math.min(1, dPrime), rSquared);

		}

		tf.close();

	}

	private Pair<ArrayList<double[][]>, ArrayList<String>> overlapVariantsWithBedRegion(Feature region, ArrayList<Integer> variantPositions, ArrayList<String> variantNames, ArrayList<double[][]> genotypes) {

		ArrayList<double[][]> output = new ArrayList<>();
		ArrayList<String> variantNamesOut = new ArrayList<>();

		for (int j = 0; j < variantPositions.size(); j++) {
			int varPos = variantPositions.get(j);
			Feature f = new Feature(region.getChromosome(), varPos, (varPos + 1));
			if (region.overlaps(f)) {
				output.add(genotypes.get(j));
				variantNamesOut.add(variantNames.get(j));
			}
		}

		return new Pair<ArrayList<double[][]>, ArrayList<String>>(output, variantNamesOut);
	}

	public boolean run() throws Exception {
		Gpio.createDir(outputdir);

		VCFVariant lastVariant = null;

		System.out.println("Assoc: " + vcf);
		System.out.println("Covar: " + covariateFile);
		System.out.println("Disease: " + diseaseStatusFile);
		System.out.println("Out: " + outputdir);
		// index covariate samples to vcf samples
		HashSet<String> excludeTheseSamples = new HashSet<String>();
		if (samplesToExclude != null) {
			TextFile excl = new TextFile(samplesToExclude, TextFile.R);
			excludeTheseSamples = (HashSet<String>) excl.readAsSet(0, TextFile.tab);
			excl.close();
		}
		HashMap<String, Integer> diseaseStatus = new HashMap<String, Integer>();

		TextFile tf = new TextFile(diseaseStatusFile, TextFile.R);
		String[] elems = tf.readLineElems(Strings.tab);
		while (elems != null) {

			String sample = elems[0];
			String statusStr = elems[1];

			Integer status = Integer.parseInt(statusStr);
			diseaseStatus.put(sample, status);

			elems = tf.readLineElems(Strings.tab);
		}
		tf.close();

		System.out.println(diseaseStatus.size() + " disease status samples loaded");

		//
		VCFGenotypeData data = new VCFGenotypeData(vcf);
//		VCFFunctions v = new VCFFunctions();
		ArrayList<String> vcfSamples = data.getSamples();

		HashSet<String> samplesWithDiseaseStatus = new HashSet<String>();
		for (String sample : vcfSamples) {
			if (diseaseStatus.containsKey(sample)) {
				if (excludeTheseSamples == null || !excludeTheseSamples.contains(sample)) {
					samplesWithDiseaseStatus.add(sample);
				}
			}
		}


		// assume rows are samples and columns are covariates
		// this loads only samples that also have a disease status loaded..
		DoubleMatrixDataset<String, String> covariates = DoubleMatrixDataset.loadSubsetOfTextDoubleData(covariateFile, "\t", samplesWithDiseaseStatus, covariatesToInclude);
		System.out.println("Covariate matrix: " + covariates.rows() + " samples " + covariates.columns() + " covariates");

		LinkedHashMap<String, Integer> covariateSampleHash = covariates.getHashRows();

		HashMap<String, Integer> sampleToIntGenotypes = new HashMap<String, Integer>();
		int ctr = 0;
		boolean[] genotypesWithCovariatesAndDiseaseStatus = new boolean[vcfSamples.size()];
		int vcfSampleCtr = 0;

		HashSet<String> alternatecovariateSamples = new HashSet<String>();
		Set<String> keys = covariateSampleHash.keySet();
		HashMap<String, String> altToSample = new HashMap<String, String>();
		for (String sample : keys) {
			alternatecovariateSamples.add(sample + "_" + sample);
			altToSample.put(sample, sample + "_" + sample);
		}


		ArrayList<String> samplesIntersect = new ArrayList<String>();
		for (String sample : vcfSamples) {
			if (covariateSampleHash.containsKey(sample) || alternatecovariateSamples.contains(sample)) {
				sampleToIntGenotypes.put(sample, ctr);
				genotypesWithCovariatesAndDiseaseStatus[vcfSampleCtr] = true;
				samplesIntersect.add(sample);
				ctr++;
			}
			vcfSampleCtr++;
		}


		String[] covariateColNames = covariates.getColObjects().toArray(new String[0]);

		System.out.println(sampleToIntGenotypes.size() + " samples with disease status, covariates and genotypes");

		if (sampleToIntGenotypes.size() == 0) {
			System.out.println("Problem with matching samples...");
		} else {
			try {
				double[] finalDiseaseStatus = null; // new double[sampleToIntGenotypes.size()];
				double[][] finalCovariates = null;  // new double[sampleToIntGenotypes.size()][covariates.columns()];

				Integer[][] sampleParents = new Integer[2][sampleToIntGenotypes.size()];
				ArrayList<Integer> kidsInTrios = new ArrayList<Integer>(); // lookup hash to check for each sample whether they are in a trio. Integers are relative to genotype
				HashMap<Integer, Integer> kidToPseudoCC = new HashMap<Integer, Integer>();

				HashSet<Integer> samplesToRemove = new HashSet<Integer>(); // relative to new ordering
				boolean arraysInit = false;
				if (famfile != null) {

					HashSet<Integer> samplesAlreadyAssignedToTrio = new HashSet<Integer>();
					HashSet<Integer> samplesAlreadyAssignedAsKid = new HashSet<Integer>();
					// get the trios (only trios where kid is a case)
					ArrayList<Pair<String, Triple<String, String, String>>> famData = getTrios(famfile); // format: familyname<kid, mom, dad>

					HashSet<String> familiesAlreadyUsed = new HashSet<String>();

					if (famData.size() > 0) {
						int nrCompleteTrios = 0;
						int nrIncompleteTrios = 0;
						int nrRelatedKids = 0;

						for (Pair<String, Triple<String, String, String>> family : famData) {
							String familyName = family.getLeft();
							Triple<String, String, String> familyMembers = family.getRight();

							String kidName = familyMembers.getLeft();
							Integer kidId = sampleToIntGenotypes.get(kidName);
							if (kidId == null) {
								kidId = sampleToIntGenotypes.get(altToSample.get(kidName));
							}

							String momName = familyMembers.getMiddle();
							Integer momId = sampleToIntGenotypes.get(momName);
							if (momId == null) {
								momId = sampleToIntGenotypes.get(altToSample.get(momName));
							}

							String dadName = familyMembers.getRight();
							Integer dadId = sampleToIntGenotypes.get(dadName);
							if (dadId == null) {
								dadId = sampleToIntGenotypes.get(altToSample.get(dadName));
							}

							if (kidId != null) {
								sampleParents[0][kidId] = dadId;
								sampleParents[1][kidId] = momId;
								if (dadId != null) {
									samplesToRemove.add(dadId);
								}
								if (momId != null) {
									samplesToRemove.add(momId);
								}
								// pick one case per family
								if (familiesAlreadyUsed.contains(familyName)) {
									samplesToRemove.add(kidId);
									nrRelatedKids++;
								} else if (dadId != null && momId != null && !familiesAlreadyUsed.contains(familyName)) {
									kidsInTrios.add(kidId);

									if (samplesAlreadyAssignedToTrio.contains(kidId)) {
										System.out.println("kid sample " + kidId + " for " + kidName + " already assigned?");
									}
									samplesAlreadyAssignedToTrio.add(kidId);
									samplesAlreadyAssignedAsKid.add(kidId);

									if (samplesAlreadyAssignedAsKid.contains(dadId)) {
										System.out.println("Father: " + dadId + " for " + dadName + " already assigned as kid.");
									}
									if (samplesAlreadyAssignedAsKid.contains(momId)) {
										System.out.println("Mother: " + momId + " for " + momName + " already assigned as kid.");
									}
									kidToPseudoCC.put(kidId, sampleToIntGenotypes.size() + nrCompleteTrios);
									familiesAlreadyUsed.add(familyName);
									nrCompleteTrios++;
								} else {
									// incomplete trio. mom and dad already removed, also remove kid, because, why not.
									samplesToRemove.add(kidId);
									nrIncompleteTrios++;
								}
							}
						}

						System.out.println(nrCompleteTrios + " full trios found in genotypes.. " + nrIncompleteTrios + " incomplete trios.." + nrRelatedKids + " kids related.." + samplesToRemove.size() + " samples to remove");
						// add vector to finalCovariates and finalDiseaseStatus
						if (!skipMakingPseudoControls) {
							if (nrCompleteTrios > 0) {

								finalDiseaseStatus = new double[sampleToIntGenotypes.size() + nrCompleteTrios];
								finalCovariates = new double[sampleToIntGenotypes.size() + nrCompleteTrios][covariates.columns()];

								// set pseudocontrol disease status
								for (int i = sampleToIntGenotypes.size(); i < finalCovariates.length; i++) {
									finalDiseaseStatus[i] = 0; // 1 == unaffected, 2 == affected --> will be remapped to 0/1
								}

								// copy covariate values from kids
								for (int i = 0; i < kidsInTrios.size(); i++) {
									Integer kidId = kidsInTrios.get(i);
									Integer pseudoCC = kidToPseudoCC.get(kidId);
									if (pseudoCC != null) {
										String kidName = samplesIntersect.get(kidId);
										samplesIntersect.add(kidName + "-PseudoControl");
										Integer covariateSampleId = covariateSampleHash.get(kidName);
										for (int col = 0; col < covariates.columns(); col++) {
											finalCovariates[pseudoCC][col] = covariates.getElement(covariateSampleId, col);
										}
									}
								}
								arraysInit = true;
							}
						}
					}
				}

				if (!arraysInit) {
					finalDiseaseStatus = new double[sampleToIntGenotypes.size()];
					finalCovariates = new double[sampleToIntGenotypes.size()][covariates.columns()];
				}


				TextFile sampleListOut = new TextFile(outputdir + "samplelist.txt", TextFile.W);
				System.out.println(outputdir + "samplelist.txt");

				for (int i = 0; i < finalCovariates.length; i++) {
					sampleListOut.writeln(samplesIntersect.get(i) + "\t" + samplesToRemove.contains(i));
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
				// make this load genotype probabilities as well...
				// TODO: genotype probabilities
				// TODO: other VCF info scores...
				// TODO: exclude samples
				// 1       1955373 rs4648791       C       T       .       PASS    AR2=0;DR2=0;AF=0.264    GT:DS:GP        0|0:0.526:0.543,0.388,0.069
			/*
			 ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated Allele Frequencies">
             ##INFO=<ID=AR2,Number=1,Type=Float,Description="Allelic R-Squared: estimated correlation between most probable ALT dose and true ALT dose">
             ##INFO=<ID=DR2,Number=1,Type=Float,Description="Dosage R-Squared: estimated correlation between estimated ALT dose [P(RA) + 2*P(AA)] and true ALT dose">
             ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
             ##FORMAT=<ID=DS,Number=1,Type=Float,Description="estimated ALT dose [P(RA) + P(AA)]">
             ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Estimated Genotype Probability">
             */

				RConnection rConnection = new RConnection();
				System.out.println("R server found: " + rConnection.getServerVersion());
				rConnection.assign("status", finalDiseaseStatus);
				System.out.println("Loaded disease status");


				assignAsRMatrix(rConnection, finalCovariates, "covariates", false);
				rConnection.assign("covariatenames", covariateColNames);
				rConnection.voidEval("colnames(covariates) <- covariatenames");
				System.out.println("Loaded covariates");


//				rConnection.voidEval("print(head(status))");

//				rConnection.voidEval("print('" + threadnum + " thread environment set')");
				System.out.println("Set thread environment.. " + threadnum);

//				System.exit(-1);

				// testNormal without genotypes
				// glm_null <-  glm(status ~ condmatrix.min + collect + gender, family=binomial(logit))
				// glm.null <- glm(status ~ collect + gender, family=binomial(logit))
				// deviance_null <- summary(glm_null)$deviance

				int varctr = 0;
				TextFile logout = new TextFile(outputdir + "log.txt", TextFile.W);
				TextFile pvalout = new TextFile(outputdir + "gwas.txt", TextFile.W);
				if (finalCovariates[0].length == 0) {

					/*

					 */
					pvalout.writeln("VariantID" +
							"\tVariantChr" +
							"\tVariantPos" +
							"\tmaf" +
							"\tn" +
							"\tDevianceNull" +
							"\tDfNull" +
							"\tDevianceGeno" +
							"\tDfAlt" +
							"\tPvalIntersect" +
							"\t-Log10(pvalIntersect)" +
							"\tPval" +
							"\t-Log10(pval)");
				} else {
					pvalout.writeln("VariantID" +
							"\tVariantChr" +
							"\tVariantPos" +
							"\tmaf" +
							"\tn" +
							"\tDevianceNull" +
							"\tDfNull" +
							"\tDevianceGeno" +
							"\tDfAlt" +
							"\tPval" +
							"\t-Log10(pval)");
				}
				boolean previousVariantHasMissingGenotyes = false;
				boolean firstvariant = true;

				boolean continueTesting = true;

//				while (data.hasNext() && continueTesting) {
				while (data.hasNext()) {

					// ref <- refAllele[vcol]
					// vdose <-unlist(dosage[,vcol])
					//
					//
					VCFVariant variant = data.next();

					String name = variant.getId();
					if (snpLimit == null || snpLimit.contains(name)) {

						lastVariant = variant; // debugging purposes :(

						Double imputationqualityscore = variant.getImputationQualityScore();
						boolean testvariant = false;
						if (imputationqualityfilter) {
							if (imputationqualityscore != null && imputationqualityscore >= imputationqualitythreshold) {
								testvariant = true;
							} else if (imputationqualityscore == null) {
								System.err.println("No imputaton quality score for variant: " + variant.getChr() + "-" + variant.getPos() + "-" + variant.getId());
								System.err.println("In path: " + vcf);
								logout.writeln("Imputation quality score below threshold:\t" + imputationqualityscore + "\t" + variant.getChr() + "-" + variant.getPos() + "-" + variant.getId());
							}
						} else {
							testvariant = true;
						}

						if (testvariant) {

							// recode the genotypes to the same ordering as the covariate table
							// generate pseudocontrol genotypes
							Triple<double[][], Double, Integer> recodedGenotypes = recodeGenotypes(
									variant.getChr() + "-" + variant.getPos() + "-" + variant.getId(),
									genotypesWithCovariatesAndDiseaseStatus,
									variant.getGenotypeAlleles(),
									variant.getAlleles().length,
									kidsInTrios,
									sampleParents,
									kidToPseudoCC,
									samplesToRemove,
									finalCovariates.length);


							Integer sampleswithgenotypes = recodedGenotypes.getRight();
							double maf = recodedGenotypes.getMiddle();

							double mafthresholdD = minNObservedAllele / (sampleswithgenotypes * 2);
							if (maf >= mafthresholdD) {

								double[][] genotypes = recodedGenotypes.getLeft();

								if (finalCovariates[0].length == 0) {
									// no need to run the null model as there are no covariates
									assignAsRMatrix(rConnection, genotypes, "genotypes", true);
									rConnection.voidEval("glm_geno <- glm(status[notnagenotypes] ~ genotypes[notnagenotypes,], family=binomial(logit))");
									rConnection.voidEval("glm_summary <- summary(glm_geno)");

									REXP rexp = rConnection.eval("glm_summary$aliased[\"genotypes[notnagenotypes, ]\"]");
									boolean aliased = false;
									if (rexp.isLogical()) {
										boolean[] bools = ((REXPLogical) rexp).isTRUE();
										if (bools[0]) {
											aliased = true;
										}
									}

									if (!aliased) {
										Double pgeno = rConnection.eval("glm_summarycoefficients[\"genotypes\",4]").asDouble();
										Double pintersect = rConnection.eval("pintersect <- glm_summarycoefficients[1,4]").asDouble();
										String outstr = variant.getId()
												+ "\t" + variant.getChr()
												+ "\t" + variant.getPos()
												+ "\t" + maf
												+ "\t" + sampleswithgenotypes
												+ "\t" + null
												+ "\t" + null
												+ "\t" + null
												+ "\t" + null
												+ "\t" + pintersect
												+ "\t" + (-Math.log10(pintersect))
												+ "\t" + pgeno
												+ "\t" + (-Math.log10(pgeno));
										pvalout.writeln(outstr);
									} else {
										logout.writeln("variant genotypes aliased: " + imputationqualityscore + " below threshold " + mafthresholdD
												+ "\t" + variant.getId()
												+ "\t" + variant.getChr()
												+ "\t" + variant.getPos()
												+ "\t" + variant.getMAF()
										);
									}
								} else {
									// avoid running the null model very often...
									assignAsRMatrix(rConnection, genotypes, "genotypes", true);
									if (!sampleswithgenotypes.equals(genotypes[0].length)) {
										rConnection.voidEval("notnagenotypes<-apply(!is.na(genotypes),1,any)");
										rConnection.voidEval("glm_null <- glm(status[notnagenotypes] ~ covariates[notnagenotypes,], family=binomial(logit))");
										rConnection.voidEval("nullSummary<-summary(glm_null)");
//										rConnection.voidEval("print(nullSummary)");
										rConnection.voidEval("deviance_null <- summary(glm_null)$deviance");
										rConnection.voidEval("nulldf <- summary(glm_null)$df[1]");
										previousVariantHasMissingGenotyes = true;
									} else {
										if (previousVariantHasMissingGenotyes || firstvariant) {
											rConnection.voidEval("notnagenotypes<-apply(!is.na(genotypes),1,any)");
											rConnection.voidEval("glm_null <- glm(status[notnagenotypes] ~ covariates[notnagenotypes,], family=binomial(logit))");
											rConnection.voidEval("nullSummary<-summary(glm_null)");
//											rConnection.voidEval("print(nullSummary)");
											rConnection.voidEval("deviance_null <- summary(glm_null)$deviance");
											rConnection.voidEval("nulldf <- summary(glm_null)$df[1]");
											firstvariant = false;
										}
										previousVariantHasMissingGenotyes = false;
									}

									// write out the null model summary
									double devianceNull = rConnection.eval("deviance_null").asDouble();


									rConnection.voidEval("glm_geno <- glm(status[notnagenotypes] ~ genotypes[notnagenotypes,] + covariates[notnagenotypes,], family=binomial(logit))");

									rConnection.voidEval("glm_summary <- summary(glm_geno)");
									// rConnection.voidEval("print(glm_summary)");

									// this depends on the number of alleles (!!)

									boolean aliased = false;
									if (genotypes.length == 1) {
										REXP rexp = rConnection.eval("glm_summary$aliased[\"genotypes[notnagenotypes, ]\"]");
										if (rexp.isLogical()) {
											boolean[] bools = ((REXPLogical) rexp).isTRUE();
											if (bools[0]) {
												aliased = true;
											}
										}
									} else {
										rConnection.voidEval("print(glm_summary$aliased)");
										boolean[] aliasPerAllele = new boolean[genotypes.length];
										for (int i = 1; i < genotypes.length + 1; i++) {
											REXP rexp = rConnection.eval("glm_summary$aliased[\"genotypes[notnagenotypes, ]" + i + "\"]");
											if (rexp.isLogical()) {
												boolean[] bools = ((REXPLogical) rexp).isTRUE();
												if (bools[0]) {
													aliasPerAllele[i - 1] = true;
												}
											} else {
												System.out.println("not logical");
											}
										}
										aliased = true;
										for (int i = 0; i < aliasPerAllele.length; i++) {
											if (!aliasPerAllele[i]) {
												aliased = false;
											}
											System.out.println(i + "\t" + aliasPerAllele[i]);
										}
									}


//							REXP rexp = rConnection.eval("glm_summary$aliased[\"genotypes1\"] || glm_summary$aliased[\"genotypes2\"]");

// rConnection.voidEval("print(glm_summary$aliased[\"genotypes1\"] || glm_summary$aliased[\"genotypes2\"])");
//							if (rexp.isLogical()) {
//								boolean[] bools = ((REXPLogical) rexp).isTRUE();
//								if (!bools[0]) {

//							Double beta = rConnection.eval("glm_summary$coefficients[\"genotypes1\", 1]").asDouble();
//							Double OR = rConnection.eval("exp(beta)").asDouble();
//							rConnection.voidEval("betase<-glm_summary$coefficients[\"genotypes\", 2]");
//							Double betase = rConnection.eval("betase").asDouble();
//							Double orLow = rConnection.eval("exp(beta - 1.96 * betase)").asDouble();
//							Double orHigh = rConnection.eval("exp(beta + 1.96 * betase)").asDouble();


									if (!aliased) { // at least one genotype vector should remain...
										rConnection.voidEval("deviance_geno <- glm_summary$deviance");
										Double deviance_geno = rConnection.eval("deviance_geno").asDouble();

										rConnection.voidEval("deltaDeviance <- deviance_null - deviance_geno");
										rConnection.voidEval("altdf <- glm_summary$df[1]");
										double altdf = rConnection.eval("altdf").asDouble();
										double nulldf = rConnection.eval("nulldf").asDouble();
										rConnection.voidEval("dfdiff <- altdf - nulldf");

										//rConnection.voidEval("print(dfdiff)");

										Double pval = rConnection.eval("pchisq(deltaDeviance, df = dfdiff, lower.tail = FALSE, log.p = FALSE)").asDouble();

										if (pval == null || Double.isNaN(pval) || Double.isInfinite(pval)) {
											pval = 1d;
										}
										if (Math.abs(devianceNull - deviance_geno) < 0.0001) {
											pval = 1d;
										}


										double pvallog = Math.abs(-Math.log10(pval));


//							Double pval = rConnection.eval("pchisq(deltaDeviance, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)").asDouble();

//									Double caseFreq = rConnection.eval("mean(genotypes[status == 1]) / 2").asDouble();
//									Double controlFreq = rConnection.eval("mean(vdose[status == 0]) / 2").asDouble();

										String outstr = variant.getId()
												+ "\t" + variant.getChr()
												+ "\t" + variant.getPos()
												+ "\t" + maf
												+ "\t" + sampleswithgenotypes
												+ "\t" + devianceNull
												+ "\t" + nulldf
												+ "\t" + deviance_geno
												+ "\t" + altdf
												+ "\t" + pval
												+ "\t" + pvallog;
//								System.out.println(outstr);
										pvalout.writeln(outstr);
									} else {
										logout.writeln("variant genotypes aliased: " + imputationqualityscore + " below threshold " + mafthresholdD
												+ "\t" + variant.getId()
												+ "\t" + variant.getChr()
												+ "\t" + variant.getPos()
												+ "\t" + variant.getMAF()

										);
									}
								}


//								}
//							if (variant.getId().equals("rs7919913")) {
//								System.out.println(outputdir + "/rs7919913.txt");
//								rConnection.voidEval("sink('" + outputdir + "/rs7919913.txt')");
//
//								rConnection.voidEval("print('Alleles (ref,alt): " + Strings.concat(variant.getAlleles(), Strings.comma) + "')");
//								rConnection.voidEval("print(glm_summary)");
//								rConnection.voidEval("print(genotypes)");
//								rConnection.voidEval("sink()");
//								System.exit(0);
//							}
//							}

							} else {
								logout.writeln("variant skipped maf: " + maf + " below threshold " +
										variant.getId()
										+ "\t" + variant.getChr()
										+ "\t" + variant.getPos()
										+ "\t" + variant.getMAF()

								);
							}
						} else {
							logout.writeln("variant skipped rsq: " + imputationqualityscore + " below threshold " +
									variant.getId()
									+ "\t" + variant.getChr()
									+ "\t" + variant.getPos()
									+ "\t" + variant.getMAF()

							);

						}

					}


					// encode genotypes
					// testNormal with genotype (how to handle R2 < 0.8)?
					// glm_test <- glm(status ~ condmatrix.min + collect + gender + vdose, family=binomial(logit))

                    /*
					 if (!summary(glm.out)$aliased["vdose"]) {

                     }
                     */
					varctr++;
					if (varctr == 250) {
						continueTesting = false;
					}
					if (varctr % 10 == 0) {
						if (varctr % 100 == 0) {
							System.out.print(varctr + "\n");
						} else {
							System.out.print(".");
						}


					}
				}

				logout.close();
				pvalout.close();
				rConnection.close();


			} catch (RserveException ex) {
				System.err.println(ex.getMessage());
				System.err.println("Could not connect to RServe");
				System.err.println("ERRORRRRRR: " + vcf);
				System.err.println(lastVariant.getChr() + ":" + lastVariant.getPos() + "-" + lastVariant.getId());
				System.exit(0);
				return false;

			} catch (REngineException e) {
				e.printStackTrace();
				System.err.println("ERRORRRRRR: " + vcf);
				System.err.println(lastVariant.getChr() + ":" + lastVariant.getPos() + "-" + lastVariant.getId());
				System.exit(0);
				return false;
			} catch (REXPMismatchException e) {
				e.printStackTrace();
				System.err.println("ERRORRRRRR: " + vcf);
				System.err.println(lastVariant.getChr() + ":" + lastVariant.getPos() + "-" + lastVariant.getId());
				System.exit(0);
				return false;
			} catch (Exception e) {
				System.err.println("ERRORRRRRR: " + vcf);

				e.printStackTrace();
				System.err.println(lastVariant.getChr() + ":" + lastVariant.getPos() + "-" + lastVariant.getId());

				System.exit(0);
				return false;
			}

		}
		return true;
	}

	private double[][] transpose(double[][] genotypes) {
		double[][] output = new double[genotypes[0].length][genotypes.length];
		for (int i = 0; i < genotypes.length; i++) {
			for (int j = 0; j < genotypes[i].length; j++) {
				output[j][i] = genotypes[i][j];
			}
		}
		return output;
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
//				System.out.println(kid + "\t" + control);
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

	public void assignAsRMatrix(RConnection rEngine, double[][] sourceArray, String nameToAssignOn, boolean transposeafter) throws REngineException, REXPMismatchException {

//		System.out.println("assigning matrix: " + sourceArray.length + "x" + sourceArray[0].length + " to var: " + nameToAssignOn);

//		rEngine.voidEval("if(!exists('" + nameToAssignOn + "') || !identical(as.numeric(dim(" + nameToAssignOn + ")+0), c(" + sourceArray.length + "," + sourceArray[0].length + "))){" +
//				"print('assigning new matrix.')\n" +
//				nameToAssignOn + " <- matrix(,nrow=" + sourceArray.length + ",ncol=" + sourceArray[0].length + ")" +
//				"}");


		if (transposeafter) {
			rEngine.voidEval("if(exists('" + nameToAssignOn + "')){" + nameToAssignOn + " <- t(" + nameToAssignOn + ")}");
		}
		// rEngine.voidEval("if(!exists('" + nameToAssignOn + "')){" + nameToAssignOn + " <- matrix(,nrow=" + sourceArray.length + ",ncol=" + sourceArray[0].length + ")}");

		rEngine.voidEval(nameToAssignOn + " <- matrix(,nrow=" + sourceArray.length + ",ncol=" + sourceArray[0].length + ")");

		for (int i = 0; i < sourceArray.length; i++) {
			rEngine.assign("temp", sourceArray[i]);
			rEngine.voidEval(nameToAssignOn + "[" + (i + 1) + ",] <- temp");
//			if (i % 10000 == 0) {
//				System.out.println(i + " rows processed");
//			}
		}

		if (transposeafter) {
			rEngine.voidEval(nameToAssignOn + " <- t(" + nameToAssignOn + ")");
		}

		//rEngine.voidEval("print(head(" + nameToAssignOn + "))");


	}

	private Triple<double[][], Double, Integer> recodeGenotypes(
			String variant,
			boolean[] includeGenotype,
			double[][] genotypeAlleles,
			int nrAlleles,
			ArrayList<Integer> kidsInTrios,
			Integer[][] sampleParents,
			HashMap<Integer, Integer> kidToPseudoCC,
			HashSet<Integer> samplesToRemove,
			int nrsamples) {


		double[][] outputGenotypes = new double[nrAlleles - 1][nrsamples];
		int individualCounter = 0;

		// first iterate the genotyped samples to load the genotypes
		for (int i = 0; i < genotypeAlleles[0].length; i++) {
			if (includeGenotype[i]) { // this is set to false if the genotype doesn't have a disease status or covariate data.
				byte b1 = (byte) genotypeAlleles[0][i];
				byte b2 = (byte) genotypeAlleles[1][i];
				if (b1 == -1) {
					for (int q = 0; q < outputGenotypes.length; q++) {
						outputGenotypes[q][individualCounter] = Double.NaN;
					}
				} else {
					if (b1 == b2) {
						// homozygote
						if (b1 == 0) {
							// do nothing
						} else {
							int allele = b1 - 1;
							if (allele >= 0) {
								outputGenotypes[allele][individualCounter] = 2;
							}
						}
					} else {
						int allele1 = b1 - 1;
						int allele2 = b2 - 1;
						if (allele1 >= 0) {
							outputGenotypes[allele1][individualCounter] = 1;
						}
						if (allele2 >= 0) {
							outputGenotypes[allele2][individualCounter] = 1;
						}
					}
				}
				individualCounter++;
			}
		}

		HashSet<Integer> idsToNa = new HashSet<Integer>();
		idsToNa.addAll(samplesToRemove);

		if (!skipMakingPseudoControls) {
			System.out.println(individualCounter + " individuals processed before PseudoCC");
			// now make the pseudocontrols
			int nrMendelianErrors = 0;
			int pseudocontrolsSet = 0;
			for (int i = 0; i < kidsInTrios.size(); i++) {
//			System.out.println("generating pseudoCC genotypes: " + nrAlleles + " alleles");
				Integer kidId = kidsInTrios.get(i);
				Integer dadId = sampleParents[0][kidId];
				Integer momId = sampleParents[1][kidId];

				Integer pseudoControlNr = kidToPseudoCC.get(kidId);

				if (pseudoControlNr != null) {

					double[] kidGt = new double[nrAlleles - 1];
					double[] momGt = new double[nrAlleles - 1];
					double[] dadGt = new double[nrAlleles - 1];
					for (int j = 0; j < nrAlleles - 1; j++) {
						kidGt[j] = outputGenotypes[j][kidId];
						momGt[j] = outputGenotypes[j][momId];
						dadGt[j] = outputGenotypes[j][dadId];
					}
					if (Double.isNaN(kidGt[0])) {
						// parents can stay
						// set the genotypes of the pseudocontrol missing
						for (int j = 0; j < nrAlleles - 1; j++) {
							idsToNa.add(pseudoControlNr);
						}
					} else {
						// remove mom and dad from the equation (later on)
						idsToNa.add(momId);
						idsToNa.add(dadId);

						if (Double.isNaN(momGt[0]) || Double.isNaN(dadGt[0])) {
							// broken trio due to missing genotypes.. remove the pseudocontrol
							// and the parents
							idsToNa.add(pseudoControlNr);
							idsToNa.add(momId);
							idsToNa.add(dadId);
						} else {
							// mendelian error check
							int[] allelesInDad = recode(dadId, genotypeAlleles);
							int[] allelesInMom = recode(momId, genotypeAlleles);
							int[] allelesInKid = recode(kidId, genotypeAlleles);

							HashSet<Pair<Integer, Integer>> allowedGenotypes = new HashSet<Pair<Integer, Integer>>();
							ArrayList<Integer> dadAlleles = new ArrayList<Integer>();
							ArrayList<Integer> momAlleles = new ArrayList<Integer>();
							HashMap<Integer, Integer> uniqueAlleles = new HashMap<Integer, Integer>();
							// count the occurrence of each allele
							for (int a = 0; a < allelesInDad.length; a++) {
								dadAlleles.add(allelesInDad[a]);
								Integer nr = uniqueAlleles.get(allelesInDad[a]);
								if (nr == null) {
									nr = 1;
								} else {
									nr++;
								}
								uniqueAlleles.put(allelesInDad[a], nr);
								momAlleles.add(allelesInMom[a]);
								nr = uniqueAlleles.get(allelesInMom[a]);
								if (nr == null) {
									nr = 1;
								} else {
									nr++;
								}
								uniqueAlleles.put(allelesInMom[a], nr);

								// cross with mom: generate all possible alleles
								for (int b = 0; b < allelesInMom.length; b++) {
									allowedGenotypes.add(new Pair<Integer, Integer>(allelesInDad[a], allelesInMom[b]));
									allowedGenotypes.add(new Pair<Integer, Integer>(allelesInMom[b], allelesInDad[a]));
								}
							}

							Pair<Integer, Integer> allelesKid = new Pair<Integer, Integer>(allelesInKid[0], allelesInKid[1]);
							if (!allowedGenotypes.contains(allelesKid)) {
								// mendelian error.
								// remove the ids for this variant.. all of them
								idsToNa.add(pseudoControlNr);
								idsToNa.add(kidId);
								idsToNa.add(momId);
								idsToNa.add(dadId);

								nrMendelianErrors++;
							} else {
								pseudocontrolsSet++;
								// remove the kids alleles from the list of available alleles
								Set<Integer> keys = uniqueAlleles.keySet();
								for (int a = 0; a < allelesInKid.length; a++) {
									int nr = uniqueAlleles.get(allelesInKid[a]); // get count for each allele
//								System.out.println("kid is taking one allele " + allelesInKid[a] + " - " + nr);
									nr--; // claim one
									uniqueAlleles.put(allelesInKid[a], nr); // put it back
								}

								// make pseudo from remaining alleles
								for (Integer k : keys) { // key to an allele
									int nr = uniqueAlleles.get(k); // the number of this type of allele that is left
//								System.out.println(k + " - " + nr);
									if (nr == 2) { // 2 alleles of this type left: pseudoCC must be homozygote for this allele
										if (k.equals(0)) { // this is the reference allele, which we don't set
//										pseudo = k + "/" + k;
//										System.out.println("pseudo hom: " + pseudo);
										} else {
											outputGenotypes[k - 1][pseudoControlNr] = 2; //
//										pseudo = k + "/" + k;
//										System.out.println("pseudo hom: " + pseudo);
											break;
										}
									} else if (nr != 0) {
										if (k.equals(0)) {
											// do nothing
//										System.out.println();
//										System.out.println("pseudo allele: " + 0);
										} else {
											outputGenotypes[k - 1][pseudoControlNr] = 1; // only one allele of this type left. PseudoCC must be het
//										System.out.println("pseudo allele: " + k);
										}
									}
								}

								// remove mom and dad..
								idsToNa.add(momId);
								idsToNa.add(dadId);
							}
						}
					}
					individualCounter++;

				}
			}

			System.out.println(individualCounter + " individuals processed after PseudoCC");
			System.out.println(pseudocontrolsSet + " pseudocontrols set");
			System.out.println(idsToNa.size() + " individuals to remove");

			if (nrMendelianErrors > 0) {
				System.out.println(variant + " has mendelian errors: " + nrMendelianErrors);
			}


		}

		// remove the samples that were selected for removal
		int nrNaN = 0;
		Integer[] meh = idsToNa.toArray(new Integer[0]);
		for (int z = 0; z < meh.length; z++) {
			int id = meh[z];

			nrNaN++;
			for (int j = 0; j < nrAlleles - 1; j++) {
				outputGenotypes[j][id] = Double.NaN;
			}
		}
//		System.out.println(nrNaN + " individuals removed after (relatedness and such)");


		int[] nrAllelesPresent = new int[nrAlleles];
		int called = 0;

		// reassess allele frequency
		int nrSamplesWGenotypes = 0;
		for (int sample = 0; sample < nrsamples; sample++) {
			if (!Double.isNaN(outputGenotypes[0][sample])) {
				called += 2;
				nrSamplesWGenotypes++;
				int nrAllelesLeft = 2;
				for (int j = 0; j < nrAlleles - 1; j++) {
					if (outputGenotypes[j][sample] == 2d) {
						nrAllelesPresent[j + 1] += 2;
						nrAllelesLeft -= 2;
					} else if (outputGenotypes[j][sample] == 1d) {
						nrAllelesPresent[j + 1] += 1;
						nrAllelesLeft -= 1;
					}
				}
				nrAllelesPresent[0] += nrAllelesLeft;
			}
		}

//		System.out.println(nrSamplesWGenotypes + " samples with genotypes remain.");

		// recalculate maf
		double maf = 1;
		for (int i = 0; i < nrAlleles; i++) {
			double d = (double) nrAllelesPresent[i] / called;
			if (d < maf) {
				maf = d;
			}
		}


		return new Triple<double[][], Double, Integer>(outputGenotypes, maf, called / 2);
	}

	int[] recode(int id, double[][] gt) {
		int[] alleles = new int[2];


		alleles[0] = (int) gt[0][id];
		alleles[1] = (int) gt[1][id];
		return alleles;
	}

	boolean isHet(int[] alleles) {
		return (alleles[0] == alleles[1]);
	}

	public void skipMakingPseudoControls(boolean makePseudoControls) {
		this.skipMakingPseudoControls = makePseudoControls;
	}
}
