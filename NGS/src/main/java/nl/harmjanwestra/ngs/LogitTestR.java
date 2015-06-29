package nl.harmjanwestra.ngs;

import nl.harmjanwestra.ngs.GenotypeFormats.VCF.VCFGenotypeData;
import nl.harmjanwestra.ngs.GenotypeFormats.VCFVariant;
import org.rosuda.REngine.REXP;
import org.rosuda.REngine.REXPLogical;
import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;
import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;

/**
 * Created by hwestra on 6/1/15.
 */
public class LogitTestR implements Callable<Boolean> {

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

	public LogitTestR(String vcf,
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

	@Override
	public Boolean call() throws Exception {

		VCFVariant lastVariant = null;
		try {
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
			DoubleMatrixDataset<String, String> covariates = DoubleMatrixDataset.loadSubsetOfTextDoubleData(covariateFile, "\t", samplesWithDiseaseStatus, covariatesToInclude);
			System.out.println("Covariate matrix: " + covariates.rows() + " samples " + covariates.columns() + " covariates");

			LinkedHashMap<String, Integer> covariateSamples = covariates.getHashRows();

			HashMap<String, Integer> sampleToIntGenotypes = new HashMap<String, Integer>();
			int ctr = 0;
			boolean[] includeGenotype = new boolean[vcfSamples.size()];
			int vcfSampleCtr = 0;

			HashSet<String> alternatecovariateSamples = new HashSet<String>();
			Set<String> keys = covariateSamples.keySet();
			HashMap<String, String> altToSample = new HashMap<String, String>();
			for (String sample : keys) {
				alternatecovariateSamples.add(sample + "_" + sample);
				altToSample.put(sample, sample + "_" + sample);
			}

			for (String sample : vcfSamples) {
				if (covariateSamples.containsKey(sample) || alternatecovariateSamples.contains(sample)) {
					sampleToIntGenotypes.put(sample, ctr);
					ctr++;
					includeGenotype[vcfSampleCtr] = true;
				}
				vcfSampleCtr++;
			}


			String[] covariateColNames = covariates.getColObjects().toArray(new String[0]);

			System.out.println(sampleToIntGenotypes.size() + " samples with disease status, covariates and genotypes");

			if (sampleToIntGenotypes.size() == 0) {
				System.out.println("Problem with matching samples...");
			} else {


				double[] finalDiseaseStatus = null; // new double[sampleToIntGenotypes.size()];
				double[][] finalCovariates = null;  // new double[sampleToIntGenotypes.size()][covariates.columns()];

				Integer[][] sampleParents = new Integer[2][sampleToIntGenotypes.size()];
				boolean[] sampleIsParent = new boolean[sampleToIntGenotypes.size()]; // we need this to set the genotype to N/A after making the pseudocontrol
				ArrayList<Integer> kidsInTrios = new ArrayList<Integer>(); // lookup hash to check for each sample whether they are in a trio.
				int nrTrios = 0;
				HashSet<Integer> samplesAlreadyAssignedToTrio = new HashSet<Integer>();
				if (famfile != null) {
					ArrayList<Triple<String, String, String>> trios = getTrios(famfile); // format: kid, dad, mom
					int nrSamplesToRemove = 0;
					if (trios.size() > 0) {
						// make pseudocontrols

						for (Triple<String, String, String> triple : trios) {
							String kidName = triple.getLeft();
							Integer kidId = sampleToIntGenotypes.get(kidName);


							String dadName = triple.getMiddle();
							Integer dadId = sampleToIntGenotypes.get(dadName);


							String momName = triple.getRight();
							Integer momId = sampleToIntGenotypes.get(momName);

							if (kidId != null) {

								if (dadId != null && momId != null) {
									// no problem, just include the case

//									if (samplesAlreadyAssignedToTrio.contains(kidId)) {
//										System.out.println("kid sample " + kidId + " for " + kidName + " already assigned?");
//									}
//									samplesAlreadyAssignedToTrio.add(kidId);
//									if (samplesAlreadyAssignedToTrio.contains(dadId)) {
//										System.out.println("dad sample " + dadId + " for " + dadName + " already assigned?");
//									}
//									samplesAlreadyAssignedToTrio.add(dadId);
//									if (samplesAlreadyAssignedToTrio.contains(momId)) {
//										System.out.println("dad sample " + momId + " for " + momName + " already assigned?");
//									}
//									samplesAlreadyAssignedToTrio.add(momId);

									sampleIsParent[dadId] = true;
									sampleIsParent[momId] = true;
									nrSamplesToRemove += 2;

									sampleParents[0][kidId] = dadId;
									sampleParents[1][kidId] = momId;
									kidsInTrios.add(kidId);
									nrTrios++;
								} else {
									if (dadId != null) {
										sampleIsParent[dadId] = true;
										nrSamplesToRemove++;
									}
									if (momId != null) {
										sampleIsParent[momId] = true;
										nrSamplesToRemove++;
									}
								}
							}
						}


						System.out.println(nrSamplesToRemove + " samples will be removed due to relatedness.. " + nrTrios + " full trios found in genotypes");

						// return an array of booleans: set these people to N/A because they are parents to one of the cases

						// return an int[samples][parents] with the positions of the parents genotypes

						// add vector to finalCovariates and finalDiseaseStatus

						if (nrTrios > 0) {

							finalDiseaseStatus = new double[sampleToIntGenotypes.size() + nrTrios];
							finalCovariates = new double[sampleToIntGenotypes.size() + nrTrios][covariates.columns() + 1];
							String[] covariateColNamesTmp = new String[covariateColNames.length + 1];
							System.arraycopy(covariateColNames, 0, covariateColNamesTmp, 0, covariateColNames.length);
							covariateColNamesTmp[covariateColNamesTmp.length - 1] = "PseudoCC";
							covariateColNames = covariateColNamesTmp;
							System.out.println("also adding extra covariate: " + (covariates.columns() + 1) + " covariates will be set");
							for (int i = sampleToIntGenotypes.size(); i < finalCovariates.length; i++) {
								finalCovariates[i][finalCovariates[i].length - 1] = 1;
								finalDiseaseStatus[i] = 0;
							}

							// copy covariate values
							for (int i = 0; i < kidsInTrios.size(); i++) {
								Integer kidId = kidsInTrios.get(i);
								int pseudoCC = sampleToIntGenotypes.size() + i;
								for (int j = 0; j < finalCovariates[finalCovariates.length - 1].length - 1; j++) {
									finalCovariates[pseudoCC][j] = finalCovariates[kidId][j];
								}
							}
						} else {
							finalDiseaseStatus = new double[sampleToIntGenotypes.size()];
							finalCovariates = new double[sampleToIntGenotypes.size()][covariates.columns()];
						}
					} else {
						finalDiseaseStatus = new double[sampleToIntGenotypes.size()];
						finalCovariates = new double[sampleToIntGenotypes.size()][covariates.columns()];
					}

				} else {
					finalDiseaseStatus = new double[sampleToIntGenotypes.size()];
					finalCovariates = new double[sampleToIntGenotypes.size()][covariates.columns()];
				}


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

				assignAsRMatrix(rConnection, finalCovariates, "covariates", false);

				rConnection.assign("covariatenames", covariateColNames);
//				rConnection.voidEval("print(covariatenames)");
				rConnection.voidEval("colnames(covariates) <- covariatenames");

//				rConnection.voidEval("print('loaded covariates')");
				System.out.println("Loaded covariates");

//				rConnection.voidEval("print(head(covariates))");

				rConnection.assign("status", finalDiseaseStatus);
				System.out.println("Loaded disease status");
//				rConnection.voidEval("print(head(status))");

//				rConnection.voidEval("print('" + threadnum + " thread environment set')");
				System.out.println("Set thread environment.. " + threadnum);

//				System.exit(-1);

				// test without genotypes
				// glm_null <-  glm(status ~ condmatrix.min + collect + gender, family=binomial(logit))
				// glm.null <- glm(status ~ collect + gender, family=binomial(logit))
				// deviance_null <- summary(glm_null)$deviance

				int varctr = 0;
				TextFile logout = new TextFile(outputdir + "log.txt", TextFile.W);
				TextFile pvalout = new TextFile(outputdir + "assoc.txt", TextFile.W);
				pvalout.writeln("VariantID" +
						"\tVariantChr" +
						"\tVariantPos" +
						"\tmaf" +
						"\tn" +
						"\tDevianceNull\tDfNull\tDevianceGeno\tDfAlt\tPval\t-Log10(pval)");

				boolean previousVariantHasMissingGenotyes = false;
				boolean firstvariant = true;
				while (data.hasNext()) {

					// ref <- refAllele[vcol]
					// vdose <-unlist(dosage[,vcol])
					//
					//
					VCFVariant variant = data.next();

					String name = variant.getId();
					if (snpLimit == null || snpLimit.contains(name)) {

						lastVariant = variant; // debugging purposes :(

//					double[][] probs = variant.getGenotypeProbs();
//					double[] dosage = variant.getGenotypeDosages();
//					HashMap<String, Double> info = variant.getInfo();

//					byte[][] genotypeAlleles = variant.getGenotypeAlleles();
//
//					String[] availableAlleles = variant.getAlleles();

						Double imputationqualityscore = variant.getInfo().get("AR2");
						boolean testvariant = false;
						if (imputationqualityfilter) {
							if (imputationqualityscore != null && imputationqualityscore >= imputationqualitythreshold) {
								testvariant = true;
							} else if (imputationqualityscore == null) {
								System.err.println("No imputaton quality score for variant: " + variant.getChr() + "-" + variant.getPos() + "-" + variant.getId());
								System.err.println("In file: " + vcf);
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
									includeGenotype,
									variant.getGenotypeAllelesNew(),
									variant.getAlleles().length,
									kidsInTrios,
									sampleParents,
									finalCovariates.length);


							Integer sampleswithgenotypes = recodedGenotypes.getRight();
							double maf = recodedGenotypes.getMiddle();

							double mafthresholdD = minNObservedAllele / (sampleswithgenotypes * 2);
							if (maf >= mafthresholdD) {

								double[][] genotypes = recodedGenotypes.getLeft();

								assignAsRMatrix(rConnection, genotypes, "genotypes", true);


								// avoid running the null model very often...
								if (!sampleswithgenotypes.equals(genotypes[0].length)) {
									rConnection.voidEval("notnagenotypes<-apply(!is.na(genotypes),1,any)");
									rConnection.voidEval("glm_null <- glm(status[notnagenotypes] ~ covariates[notnagenotypes,], family=binomial(logit))");
									rConnection.voidEval("nullSummary<-summary(glm_null)");
									rConnection.voidEval("deviance_null <- summary(glm_null)$deviance");
									rConnection.voidEval("nulldf <- summary(glm_null)$df[1]");
									previousVariantHasMissingGenotyes = true;
								} else {
									if (previousVariantHasMissingGenotyes || firstvariant) {
										rConnection.voidEval("notnagenotypes<-apply(!is.na(genotypes),1,any)");
										rConnection.voidEval("glm_null <- glm(status[notnagenotypes] ~ covariates[notnagenotypes,], family=binomial(logit))");
										rConnection.voidEval("nullSummary<-summary(glm_null)");
										rConnection.voidEval("deviance_null <- summary(glm_null)$deviance");
										rConnection.voidEval("nulldf <- summary(glm_null)$df[1]");
										firstvariant = false;
									}
									previousVariantHasMissingGenotyes = false;
								}


//							rConnection.voidEval("print(notnagenotypes)");
//
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
											+ "\t" + (-Math.log10(pval));
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
					// test with genotype (how to handle R2 < 0.8)?
					// glm_test <- glm(status ~ condmatrix.min + collect + gender + vdose, family=binomial(logit))

                    /*
					 if (!summary(glm.out)$aliased["vdose"]) {

                     }
                     */
					varctr++;
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

			}
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


		return true;
	}

	private ArrayList<Triple<String, String, String>> getTrios(String famfile) throws IOException {
		System.out.println("Loading trios from FAM file: " + famfile);
		ArrayList<Triple<String, String, String>> output = new ArrayList<Triple<String, String, String>>();
		TextFile tf = new TextFile(famfile, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {

			String kid = elems[1];
			String dad = elems[2];
			String mom = elems[3];

			if (!dad.equals("0") && !mom.equals("0")) {
				Integer control = Integer.parseInt(elems[5]);
//				System.out.println(kid + "\t" + control);
				if (control.equals(2)) {
					output.add(new Triple<String, String, String>(kid, dad, mom));
				}
			}
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();
		System.out.println(output.size() + " trios found in FAM file");
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
			byte[][] genotypeAlleles,
			int nrAlleles,
			ArrayList<Integer> kidsInTrios,

			Integer[][] sampleParents,
			int nrsamples
	) {


		double[][] outputGenotypes = new double[nrAlleles - 1][nrsamples];


		int pseudoControlCounter = 0;


		for (int i = 0; i < genotypeAlleles[0].length; i++) {
			if (includeGenotype[i]) {
				byte b1 = genotypeAlleles[0][i];
				byte b2 = genotypeAlleles[1][i];
				if (b1 != -1) {
					if (b1 == b2) {
						// homozygote
						if (b1 == 0) {
							// do nothing
						} else {
							int allele = b1 - 1;
							if (allele >= 0) {
								outputGenotypes[allele][pseudoControlCounter] = 2;
							}
						}


					} else {
						int allele1 = b1 - 1;
						int allele2 = b2 - 1;
						if (allele1 >= 0) {
							outputGenotypes[allele1][pseudoControlCounter] = 1;
						}

						if (allele2 >= 0) {
							outputGenotypes[allele2][pseudoControlCounter] = 1;
						}

					}

				} else {
					for (int q = 0; q < outputGenotypes.length; q++) {
						outputGenotypes[q][pseudoControlCounter] = Double.NaN;
					}

				}
				pseudoControlCounter++;


			}
		}

		// make the pseudocontrols
		ArrayList<Integer> idsToNa = new ArrayList<Integer>();
		int nrMendelianErrors = 0;
		for (int i = 0; i < kidsInTrios.size(); i++) {
//			System.out.println("generating pseudoCC genotypes: " + nrAlleles + " alleles");
			Integer kidId = kidsInTrios.get(i);
			Integer dadId = sampleParents[0][kidId];
			Integer momId = sampleParents[1][kidId];
			double[] kidGt = new double[nrAlleles - 1];

			double[] momGt = new double[nrAlleles - 1];
			double[] dadGt = new double[nrAlleles - 1];
			for (int j = 0; j < nrAlleles - 1; j++) {
				kidGt[j] = outputGenotypes[j][kidId];
				momGt[j] = outputGenotypes[j][momId];
				dadGt[j] = outputGenotypes[j][dadId];
			}


			if (includeGenotype[kidId]) {
//				System.out.println(kidId + "\t" + momId + "\t" + dadId);
//				System.out.println("kid: " + Strings.concat(kidGt, Strings.forwardslash) + "\t" + includeGenotype[kidId]);
//				System.out.println("mom: " + Strings.concat(momGt, Strings.forwardslash) + "\t" + includeGenotype[momId]);
//				System.out.println("dad: " + Strings.concat(dadGt, Strings.forwardslash) + "\t" + includeGenotype[dadId]);
				if (Double.isNaN(kidGt[0])) {
					// parents can stay
					// set the genotypes of the pseudocontrol missing
					for (int j = 0; j < nrAlleles - 1; j++) {
						outputGenotypes[j][pseudoControlCounter] = Double.NaN;
					}


				} else {
					// remove mom and dad from the equation
//					for (int j = 0; j < nrAlleles - 1; j++) {
//						outputGenotypes[j][momId] = Double.NaN;
//						outputGenotypes[j][dadId] = Double.NaN;
//					}
					idsToNa.add(momId);
					idsToNa.add(dadId);

					if (!Double.isNaN(momGt[0]) && !Double.isNaN(dadGt[0])) {
						// mendelian error check
						int[] allelesInDad = recode(dadId, genotypeAlleles);
						int[] allelesInMom = recode(momId, genotypeAlleles);
						int[] allelesInKid = recode(kidId, genotypeAlleles);

//						System.out.println("recoded kid: " + Strings.concat(allelesInKid, Strings.forwardslash));
//						System.out.println("recoded mom: " + Strings.concat(allelesInMom, Strings.forwardslash));
//						System.out.println("recoded dad: " + Strings.concat(allelesInDad, Strings.forwardslash));


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
//							for (int j = 0; j < nrAlleles - 1; j++) {
//								outputGenotypes[j][ctr] = Double.NaN;
//								outputGenotypes[j][kidId] = Double.NaN;
//							}
							idsToNa.add(pseudoControlCounter);
							idsToNa.add(kidId);
//							System.out.println("Mendelian error for kid: ");
//							System.out.println("recoded kid: " + Strings.concat(allelesInKid, Strings.forwardslash) + "\t" + outputGenotypes[0][kidId]);
//							System.out.println("recoded mom: " + Strings.concat(allelesInMom, Strings.forwardslash) + "\t" + outputGenotypes[0][momId]);
//							System.out.println("recoded dad: " + Strings.concat(allelesInDad, Strings.forwardslash) + "\t" + outputGenotypes[0][dadId]);

							// remove the ids for this variant
							idsToNa.add(kidId);
							idsToNa.add(momId);
							idsToNa.add(dadId);

							nrMendelianErrors++;
							//System.exit(-1);
						} else {
							// remove the kids alleles from the list of available alleles
							Set<Integer> keys = uniqueAlleles.keySet();
							for (int a = 0; a < allelesInKid.length; a++) {
								int nr = uniqueAlleles.get(allelesInKid[a]); // get count for each allele
//								System.out.println("kid is taking one allele " + allelesInKid[a] + " - " + nr);
								nr--; // claim one
								uniqueAlleles.put(allelesInKid[a], nr); // put it back
							}

							// make pseudo from remaining alleles
							String pseudo = "";
							for (Integer k : keys) { // key to an allele
								int nr = uniqueAlleles.get(k); // the number of this type of allele that is left
//								System.out.println(k + " - " + nr);
								if (nr == 2) { // 2 alleles of this type left: pseudoCC must be homozygote for this allele
									if (k.equals(0)) { // this is the reference allele, which we don't set
//										pseudo = k + "/" + k;
//										System.out.println("pseudo hom: " + pseudo);
									} else {
										outputGenotypes[k - 1][pseudoControlCounter] = 2; //
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
										outputGenotypes[k - 1][pseudoControlCounter] = 1; // only one allele of this type left. PseudoCC must be het
//										System.out.println("pseudo allele: " + k);
									}
								}
							}


						}

					} else {
						// broken trio due to missing genotypes.. remove the pseudocontrol
						// and the parents
//						for (int j = 0; j < nrAlleles - 1; j++) {
//							outputGenotypes[j][ctr] = Double.NaN;
						idsToNa.add(pseudoControlCounter);
						idsToNa.add(momId);
						idsToNa.add(dadId);

//						}
					}
				}


			}


			pseudoControlCounter++;


		}

		if (nrMendelianErrors > 0) {
			System.out.println(variant + " has mendelian errors: " + nrMendelianErrors);
		}

		// remove the samples that were selected for removal
		for (Integer id : idsToNa) {
			for (int j = 0; j < nrAlleles - 1; j++) {
				outputGenotypes[j][id] = Double.NaN;
			}
		}

		int[] nrAllelesPresent = new int[nrAlleles];
		int called = 0;

		// reassess allele frequency
		for (int sample = 0; sample < nrsamples; sample++) {
			if (!Double.isNaN(outputGenotypes[0][sample])) {
				called += 2;
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

	int[] recode(int id, byte[][] gt) {
		int[] alleles = new int[2];


		alleles[0] = gt[0][id];
		alleles[1] = gt[1][id];
		return alleles;
	}

	boolean isHet(int[] alleles) {
		return (alleles[0] == alleles[1]);
	}

}
