package nl.harmjanwestra.assoc;

import JSci.maths.ArrayMath;
import nl.harmjanwestra.assoc.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.math.LogisticRegression;
import nl.harmjanwestra.utilities.math.LogisticRegressionResult;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;
import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.ChiSquare;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.*;

/**
 * Created by hwestra on 11/7/15.
 */
public class LRTest {

	boolean useR = false;

	LRTestOptions options;

	public LRTest(LRTestOptions options) throws IOException {
		this.options = options;
		test();
	}

	public void test() throws IOException {

		// load covariates
		VCFVariant lastVariant = null;

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

		//
		HashSet<String> covariatesToInclude = null;
		if (options.getCovariatesToInclude() != null) {
			// load hashet
			TextFile excl = new TextFile(options.getCovariatesToInclude(), TextFile.R);
			covariatesToInclude = (HashSet<String>) excl.readAsSet(0, TextFile.tab);
			excl.close();
			System.out.println(covariatesToInclude + " covariates selected from: " + options.getCovariatesToInclude());
		}

		DoubleMatrixDataset<String, String> covariates = DoubleMatrixDataset.loadSubsetOfTextDoubleData(options.getCovariateFile(), "\t", vcfSamplesWithDiseaseStatus, covariatesToInclude);
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

		// index the final list of samples
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

		System.out.println(sampleToIntGenotypes.size() + " samples with disease status, covariates and genotypes");

		if (sampleToIntGenotypes.size() == 0) {
			System.out.println("Problem with matching samples...");
		} else {
			try {
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
				RConnection rConnection = null;
				if (useR) {
					rConnection = new RConnection();
					System.out.println("R server found: " + rConnection.getServerVersion());
				}


				int iter = 0;
				boolean continueTesting = true;

				// keep a list of genotypes to condition on
				ArrayList<Triple<double[][], boolean[], Integer>> conditional = new ArrayList<>();
				ArrayList<double[][]> conditionalDosages = new ArrayList<>();
				ArrayList<String> conditionalVariantIds = new ArrayList<String>();

				AssociationFile associationFile = new AssociationFile();
				String header = associationFile.getHeader();

				HashSet<String> snpLimit = null;
				if (options.getSnpLimitFile() != null) {
					TextFile excl = new TextFile(options.getSnpLimitFile(), TextFile.R);
					snpLimit = (HashSet<String>) excl.readAsSet(0, TextFile.tab);
					excl.close();
				}
//options.getOutputdir()
				while (iter < options.getMaxIter() && continueTesting) {

					TextFile logout = new TextFile(options.getOutputdir() + "log-iter" + iter + ".txt", TextFile.W);
					TextFile pvalout = new TextFile(options.getOutputdir() + "assoc-" + iter + ".txt", TextFile.W);
					pvalout.writeln(header);
//					TextFile pvalout2 = new TextFile(outputdir + "assoc-prob-" + iter + ".txt", TextFile.W);
//					pvalout2.writeln(header);

					data.close();
					boolean debugbreak = false;
					data.open();

					// variables containing SNP with lowest pval in region
					Triple<double[][], boolean[], Integer> currentLowestPValSNP = null;
					double[][] currentLowestDosagePValSNP = null;
					String currentLowestDosagePValSNPId = null;

					int alleleOffsetGenotypes = 0;
					int alleleOffsetDosages = 0;
					double highestLog10P = 0;
					double highestLog10PProbs = 0;

					for (Triple<double[][], boolean[], Integer> c : conditional) {
						alleleOffsetGenotypes += c.getLeft().length;
					}
					for (double[][] c : conditionalDosages) {
						alleleOffsetDosages += c[0].length;
					}

					if (iter == 0) {
						System.out.println("Iteration " + iter + " starting. Model: y ~ SNP + covar. Output: " + pvalout.getFileName());

					} else {
						System.out.println("Iteration " + iter + " starting. Model: y ~ " + Strings.concat(conditionalVariantIds, Strings.semicolon) + " + SNP + covar. Output: " + pvalout.getFileName());
					}


					int variantCtr = 0;
					while (data.hasNext() && !debugbreak) {

						VCFVariant variant = data.next();

						String name = variant.getId();
						if (snpLimit == null || snpLimit.contains(name)) {

							Double imputationqualityscore = variant.getInfo().get("AR2");
							if (imputationqualityscore == null) {
								imputationqualityscore = Double.NaN;
							}

							boolean testvariant = false;

							if ((Double.isNaN(imputationqualityscore) && options.isTestVariantsWithoutImputationQuality())
									|| (!Double.isNaN(imputationqualityscore) && imputationqualityscore >= options.getImputationqualitythreshold())) {
								testvariant = true;
							} else if (Double.isNaN(imputationqualityscore)) {
								System.err.println("No imputaton quality score for variant: " + variant.getChr() + "-" + variant.getPos() + "-" + variant.getId());
								System.err.println("In file: " + options.getVcf());
								logout.writeln("Imputation quality score below threshold:\t" + imputationqualityscore + "\t" + variant.getChr() + "-" + variant.getPos() + "-" + variant.getId());
							}


							if (!testvariant) {
								logout.writeln("variant skipped rsq: " + imputationqualityscore + " below threshold " +
										variant.getId()
										+ "\t" + variant.getChr()
										+ "\t" + variant.getPos()
										+ "\t" + variant.getMAF()
								);
							} else {
								// TODO: switch this around: make the ordering of the covariate table the same as the genotype file...
								// recode the genotypes to the same ordering as the covariate table
								Triple<double[][], boolean[], Integer> unfilteredGenotypeData = filterAndRecodeGenotypes(
										genotypesWithCovariatesAndDiseaseStatus,
										variant.getGenotypeAllelesNew(),
										variant.getAlleles().length,
										finalCovariates.length);

								// generate pseudocontrol genotypes
								Pair<Pair<double[][], double[]>,
										Pair<Double, Double>> recodeGenotypes = prepareGenotypeMatrix(
										unfilteredGenotypeData,
										conditional,
										finalDiseaseStatus,
										finalCovariates,
										finalCovariates.length);

								Pair<double[][], double[]> genotypedata = recodeGenotypes.getLeft();
								Pair<Double, Double> summary = recodeGenotypes.getRight();

								double maf = summary.getLeft();
								if (maf < options.getMafthresholdD()) {
									logout.writeln("variant skipped maf: " + maf + " below threshold " +
											variant.getId()
											+ "\t" + variant.getChr()
											+ "\t" + variant.getPos()
											+ "\t" + variant.getMAF()
											+ "\t" + imputationqualityscore
									);
								} else {

									double[][] probs = variant.getGenotypeProbsNew(); // [probs][inds]
									int nrAlleles = variant.getAlleles().length;
									double[] y = genotypedata.getRight();
									double[][] x = genotypedata.getLeft();

									if (probs == null) {
										double p = pruneTestAndWrite(x, y, nrAlleles, alleleOffsetGenotypes, rConnection, pvalout, variant, maf);
										if (p > highestLog10P) {
											highestLog10P = p;
											currentLowestPValSNP = unfilteredGenotypeData;
										}
									} else {
										// recode to dosages
										double[][] dosages = variant.getImputedDosages();
										double[][] dosageMat = prepareDosageMatrix(genotypesWithCovariatesAndDiseaseStatus, dosages, finalCovariates, conditionalDosages);
										double p2 = pruneTestAndWrite(dosageMat, y, nrAlleles, alleleOffsetDosages, rConnection, pvalout, variant, maf);
										if (p2 > highestLog10PProbs) {
											highestLog10PProbs = p2;
											currentLowestDosagePValSNP = dosages;
											currentLowestDosagePValSNPId = variant.getChr().toString() + ":" + variant.getPos() + "-" + variant.getId();
										}

									}
									variantCtr++;

									if (variantCtr % 10 == 0) {
										System.out.println("Iteration: " + iter + "\tnr variants: " + variantCtr + "\tHighest P-val: " + highestLog10PProbs);
									}
								} // end maf > threshold
							}
						} // end snplimit.contains(snp)
					} // end data.hasnext

					if (currentLowestPValSNP != null) {
						conditional.add(currentLowestPValSNP);
					}
					if (currentLowestDosagePValSNP != null) {
						System.out.println("Iteration: " + iter + "\tnr variants: " + variantCtr + "\tHighest P-val: " + highestLog10PProbs);
						conditionalDosages.add(currentLowestDosagePValSNP);
						conditionalVariantIds.add(currentLowestDosagePValSNPId);
					}

					logout.close();
					pvalout.close();
//					pvalout2.close();
					iter++;
				}
				if (useR) {
					rConnection.close();
				}
			} catch (RserveException ex) {
				System.err.println(ex.getMessage());
				System.err.println("Could not connect to RServe");
				System.err.println("ERRORRRRRR: " + options.getVcf());
				System.err.println(lastVariant.getChr() + ":" + lastVariant.getPos() + "-" + lastVariant.getId());
				System.exit(0);

			} catch (REngineException e) {
				e.printStackTrace();
				System.err.println("ERRORRRRRR: " + options.getVcf());
				System.err.println(lastVariant.getChr() + ":" + lastVariant.getPos() + "-" + lastVariant.getId());
				System.exit(0);

			} catch (REXPMismatchException e) {
				e.printStackTrace();
				System.err.println("ERRORRRRRR: " + options.getVcf());
				System.err.println(lastVariant.getChr() + ":" + lastVariant.getPos() + "-" + lastVariant.getId());
				System.exit(0);

			} catch (Exception e) {
				System.err.println("ERRORRRRRR: " + options.getVcf());

				e.printStackTrace();
				System.err.println(lastVariant.getChr() + ":" + lastVariant.getPos() + "-" + lastVariant.getId());

				System.exit(0);
			}
		}
	}

	private double pruneTestAndWrite(double[][] x,
									 double[] y,
									 int nrAlleles,
									 int alleleOffset,
									 RConnection rConnection,
									 TextFile out,
									 VCFVariant variant,
									 double maf) throws REXPMismatchException, REngineException, IOException {
		Pair<double[][], boolean[]> pruned = removeCollinearVariables(x);
		x = pruned.getLeft(); // x is now probably shorter than original X
		boolean[] notaliased = pruned.getRight(); // length of original X

		// check if the alleles we put in are aliased
		int firstAllele = 1 + alleleOffset; // intercept + allleles we conditioned on (original X indexing)
		int lastAllele = 1 + alleleOffset + (nrAlleles - 1); // intercept + allleles we conditioned on + alleles for this variant (original X indexing)


		// index new position in X
		int nrRemaining = 0;
		int[] alleleIndex = new int[nrAlleles - 1];
		boolean[] colsToRemove = new boolean[x[0].length];
		for (int i = 0; i < alleleIndex.length; i++) {
			alleleIndex[i] = -1;
		}

		int newXIndex = 0;
		for (int i = 0; i < notaliased.length; i++) {
			if (notaliased[i]) {
				if (i >= firstAllele && i < lastAllele) {
					alleleIndex[i - firstAllele] = newXIndex;
					colsToRemove[newXIndex] = true;
					nrRemaining++;
				}
				if (i == 0) {
					colsToRemove[newXIndex] = false;
				}
				newXIndex++;
			}
		}

		double log10p = 0;
		AssociationResult result = new AssociationResult();
		Feature snp = new Feature(Chromosome.parseChr(variant.getChr()), variant.getPos(), variant.getPos());
		snp.setName(variant.getId());
		result.setSnp(snp);
		result.setN(x.length);
		result.setMaf(maf);

		Double imputationqualityscore = variant.getInfo().get("AR2");
		result.setImputationQualScore(imputationqualityscore);

		if (nrRemaining > 0) {
			// perform test on full model
			// remove genotypes and run test on reduced model
			if (useR) {
				// debug: test with R
				rConnection.assign("y", y);
				assignAsRMatrix(rConnection, x, "x", false);
				double[][] covarsOnly = removeGenotypes(x, colsToRemove);
				assignAsRMatrix(rConnection, covarsOnly, "z", false);

				rConnection.voidEval("glm1 <- glm(y ~ x, family=binomial(logit))");
				rConnection.voidEval("glm2 <- glm(y ~ z, family=binomial(logit))");
				rConnection.voidEval("sum1 <- summary(glm1)");
				rConnection.voidEval("sum2 <- summary(glm2)");
				rConnection.voidEval("deva <- sum1$deviance");
				rConnection.voidEval("devn <- sum2$deviance");
				double devr = rConnection.eval("sum1$deviance").asDouble();
				rConnection.voidEval("deltaDeviance <- devn - deva");

				rConnection.voidEval("altdf <- sum1$df[1]");
				rConnection.voidEval("nulldf <- sum2$df[1]");
				rConnection.voidEval("dfdiff <- altdf - nulldf");
				Double pvalR = rConnection.eval("pchisq(deltaDeviance, df = dfdiff, lower.tail = FALSE, log.p = FALSE)").asDouble();
				Double dfR = rConnection.eval("dfdiff").asDouble();
				Double deltaR = rConnection.eval("deltaDeviance").asDouble();

				double[] betas = new double[nrRemaining];

				int ctr = 0;
				for (int i = 0; i < alleleIndex.length; i++) {
					int idx = alleleIndex[i];
					if (idx != -1) {
						int idxplusone = idx + 1;
						betas[ctr] = rConnection.eval("as.vector(sum1$coefficients[" + idxplusone + ",])[1]").asDouble();
						ctr++;
					}
				}

				String betaR = Strings.concat(betas, Strings.semicolon);


			} else {

				LogisticRegression reg = new LogisticRegression();
				LogisticRegressionResult resultX = reg.univariate(y, x);
				double devx = resultX.getDeviance();
				double[][] covarsOnly = removeGenotypes(x, colsToRemove);
				LogisticRegressionResult resultCovars = reg.univariate(y, covarsOnly);
				double devnull = resultCovars.getDeviance();
				double[] betasmlelr = new double[nrRemaining];
				double[] stderrsmlelr = new double[nrRemaining];
				double[] or = new double[nrRemaining];
				double[] orhi = new double[nrRemaining];
				double[] orlo = new double[nrRemaining];

				int ctr = 0;

				for (int i = 0; i < alleleIndex.length; i++) {
					int idx = alleleIndex[i];
					if (idx != -1) {
						double beta = -resultX.getBeta()[idx];
						double se = resultX.getStderrs()[idx];
						betasmlelr[ctr] = beta;
						stderrsmlelr[ctr] = se;

						double OR = Math.exp(beta);
						double orLow = Math.exp(beta - 1.96 * se);
						double orHigh = Math.exp(beta + 1.96 * se);
						or[ctr] = OR;
						orhi[ctr] = orHigh;
						orlo[ctr] = orLow;
						ctr++;
					}
				}

				double deltaDeviance = devnull - devx;
				int df = x[0].length - covarsOnly[0].length;
				double p = ChiSquare.getP(df, deltaDeviance);
				log10p = Math.abs((-Math.log10(p)));


				result.setDevianceNull(devnull);
				result.setDevianceGeno(devx);
				result.setDf(df);

				result.setBeta(betasmlelr);
				result.setSe(stderrsmlelr);
				result.setPval(p);

//StringBuilder builder = new StringBuilder();
//				builder.append(variant.getChr());
//				builder.append("\t").append(variant.getPos());
//				builder.append("\t").append(variant.getId());
//				builder.append("\t").append(variant.getChr().toString()).append(":").append(variant.getPos()).append("-").append(variant.getId());
//				builder.append("\t").append(x.length);
//				builder.append("\t").append(maf);
//				builder.append("\t").append(devnull);
//				builder.append("\t").append(devx);
//				builder.append("\t").append(df);
//				builder.append("\t").append(Strings.concat(betasmlelr, Strings.semicolon));
//				builder.append("\t").append(Strings.concat(stderrsmlelr, Strings.semicolon));
//				builder.append("\t").append(Strings.concat(or, Strings.semicolon));
//				builder.append("\t").append(Strings.concat(orhi, Strings.semicolon));
//				builder.append("\t").append(Strings.concat(orlo, Strings.semicolon));
//				builder.append("\t").append(p);
//				builder.append("\t").append(log10p);
				out.writeln(result.toString());
			}
//			System.out.println(betaMLE + "\t" + betaR + "\t" + deltaDeviance + "\t" + deltaR + "\t" + df + "\t" + dfR + "\t" + p + "\t" + pvalR);
		} else {
			// result is null..
//			StringBuilder builder = new StringBuilder();
//			builder.append(variant.getChr());
//			builder.append("\t").append(variant.getPos());
//			builder.append("\t").append(variant.getId());
//			builder.append("\t").append(variant.getChr().toString()).append(":").append(variant.getPos()).append("-").append(variant.getId());
//			builder.append("\t").append(x.length);
//			builder.append("\t").append(maf);
//			builder.append("\t").append(0);
//			builder.append("\t").append(0);
//			builder.append("\t").append(0);
//			builder.append("\t").append(0);
//			builder.append("\t").append(0);
//			builder.append("\t").append(1);
//			builder.append("\t").append(1);
//			builder.append("\t").append(1);
//			builder.append("\t").append(1);
//			builder.append("\t").append(0);
			out.writeln(result.toString());

		}
		return log10p;
	}

	public double[][] removeGenotypes(double[][] x, boolean[] colsToRemove) {

		int nrToRemove = 0;
		for (int i = 0; i < colsToRemove.length; i++) {
			if (colsToRemove[i]) {
				nrToRemove++;
			}
		}

		double[][] output = new double[x.length][colsToRemove.length - nrToRemove];
		int ctr = 0;
		for (int j = 0; j < colsToRemove.length; j++) {
			if (!colsToRemove[j]) {
				for (int i = 0; i < x.length; i++) {
					output[i][ctr] = x[i][j];
				}
				ctr++;
			}
		}
		return output;
	}

	// remove variables with zero variance and perfect correlation
	private Pair<double[][], boolean[]> removeCollinearVariables(double[][] mat) {

		boolean[] includeCol = new boolean[mat[0].length];
		for (int c = 0; c < includeCol.length; c++) {
			includeCol[c] = true;
		}

		includeCol[0] = true; // intercept
		for (int j = 1; j < mat[0].length; j++) {
			if (includeCol[j]) {
				double[] vals = new double[mat.length];
				for (int i = 0; i < mat.length; i++) {
					vals[i] = mat[i][j];
				}

				double variance = ArrayMath.variance(vals);
				if (variance == 0d) {
					includeCol[j] = false;
				} else {
					for (int j2 = j + 1; j2 < mat[0].length; j2++) {
						if (includeCol[j2]) {
							double[] vals2 = new double[mat.length];
							for (int i = 0; i < mat.length; i++) {
								vals2[i] = mat[i][j2];
							}
							// correlate
							double variance2 = ArrayMath.variance(vals2);
							if (variance2 == 0d) {
								includeCol[j2] = false;
							} else {
								double corr = ArrayMath.correlation(vals, vals2);
								if (Math.abs(corr) == 1d) {
									includeCol[j2] = false;
								}
							}

						}
					}
				}
			}
		}

		int nrRemaining = 0;
		for (int c = 0; c < includeCol.length; c++) {
			if (includeCol[c])
				nrRemaining++;
		}

		// remove correlated variables
		double[][] matOut = new double[mat.length][nrRemaining];
		int ctr = 0;
		for (int j = 0; j < mat[0].length; j++) {
			if (includeCol[j]) {
				for (int i = 0; i < mat.length; i++) {
					matOut[i][ctr] = mat[i][j];
				}
				ctr++;
			}
		}


		return new Pair<>(matOut, includeCol);
	}


	private double[][] prepareDosageMatrix(boolean[] includeGenotype, double[][] dosages, double[][] covariates, ArrayList<double[][]> conditional) {

		int nrSamples = 0;
		for (int i = 0; i < dosages.length; i++) {
			if (includeGenotype[i]) {
				nrSamples++;
			}
		}
		int offset = 1;

		int extraAlleles = 0;
		for (int i = 0; i < conditional.size(); i++) {
			double[][] data = conditional.get(i);
			extraAlleles += data[0].length;
		}

		double[][] matrix = new double[nrSamples][dosages[0].length + covariates[0].length + extraAlleles + offset];

		int nrGenotypeSamples = dosages.length;

		// copy conditional alleles
		int addedAlleles = 0;
		for (int z = 0; z < conditional.size(); z++) {
			double[][] data = conditional.get(z);
			int ctr = 0;
			for (int i = 0; i < nrGenotypeSamples; i++) {
				if (includeGenotype[i]) {
					for (int j = 0; j < data[0].length; j++) {
						matrix[ctr][offset + addedAlleles + j] = data[i][j];
					}
					ctr++;
				}
			}
			addedAlleles += data[0].length;
		}

		int ctr = 0;
		int nrAlleles = dosages[0].length;
		for (int i = 0; i < nrGenotypeSamples; i++) {
			if (includeGenotype[i]) {
				matrix[ctr][0] = 1; // intercept

				for (int j = 0; j < nrAlleles; j++) {
					matrix[ctr][j + extraAlleles + offset] = dosages[i][j];
				}

				for (int j = 0; j < covariates[ctr].length; j++) {
					matrix[ctr][j + extraAlleles + offset + nrAlleles] = covariates[ctr][j];
				}

				ctr++;
			}
		}
		return matrix;
	}

	public void save(String variantname, double[][] yout, double[][] cov, String outdir) throws IOException {

		TextFile out = new TextFile(outdir + variantname + ".txt", TextFile.W);

		String header = "";
		for (int y = 0; y < yout[0].length; y++) {
			if (y == 0) {
				header += "y" + y;
			} else {
				header += "\ty" + y;
			}

		}
		for (int y = 0; y < cov[0].length; y++) {
			header += "\tcov" + y;
		}
		out.writeln(header);
		for (int i = 0; i < yout.length; i++) {
			String str = Strings.concat(yout[i], Strings.tab) + "\t" + Strings.concat(cov[i], Strings.tab);
			out.writeln(str);
		}
		out.close();
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
		System.out.println(output.size() + " trios found in FAM file");
		return output;
	}

	public void assignAsRMatrix(RConnection rEngine, double[][] sourceArray, String nameToAssignOn,
								boolean transposeafter) throws REngineException, REXPMismatchException {

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

	private Triple<double[][], boolean[], Integer> filterAndRecodeGenotypes(
			boolean[] includeGenotype,
			byte[][] genotypeAlleles,
			int nrAlleles,
			int nrsamples) {

		double[][] tmpgenotypes = new double[nrAlleles - 1][nrsamples];
		int individualCounter = 0;

		// first iterate the genotyped samples to load the genotypes
		int nrWithMissingGenotypes = 0;
		boolean[] genotypeMissing = new boolean[nrsamples];
		for (int i = 0; i < genotypeAlleles[0].length; i++) {
			if (includeGenotype[i]) { // this is set to false if the genotype doesn't have a disease status or covariate data.
				byte b1 = genotypeAlleles[0][i];
				byte b2 = genotypeAlleles[1][i];
				if (b1 == -1) {
					for (int q = 0; q < tmpgenotypes.length; q++) {
						tmpgenotypes[q][individualCounter] = Double.NaN;
					}
					nrWithMissingGenotypes++;
					genotypeMissing[individualCounter] = true;
				} else {
					if (b1 == b2) {
						// homozygote
						if (b1 == 0) {
							// do nothing
						} else {
							int allele = b1 - 1;
							if (allele >= 0) {
								tmpgenotypes[allele][individualCounter] = 2;
							}
						}
					} else {
						int allele1 = b1 - 1;
						int allele2 = b2 - 1;
						if (allele1 >= 0) {
							tmpgenotypes[allele1][individualCounter] = 1;
						}
						if (allele2 >= 0) {
							tmpgenotypes[allele2][individualCounter] = 1;
						}
					}
				}
				individualCounter++;
			}
		}

		return new Triple<>(tmpgenotypes, genotypeMissing, nrWithMissingGenotypes);
	}

	// output format genotypes+covariates, matched disease status, maf, cr
	private Pair<Pair<double[][], double[]>, Pair<Double, Double>> prepareGenotypeMatrix(
			Triple<double[][], boolean[], Integer> genotypeData,
			ArrayList<Triple<double[][], boolean[], Integer>> conditional,
			double[] diseaseStatus,
			double[][] covariates,
			int nrsamples) {

		int nrWithMissingGenotypes = genotypeData.getRight();
		boolean[] genotypeMissing = genotypeData.getMiddle();
		double[][] genotypes = genotypeData.getLeft();


		double[][] tmpgenotypes = genotypes;

		if (!conditional.isEmpty()) {
			// count extra columns from conditional genotypes
			int extraAlleles = 0;
			for (Triple<double[][], boolean[], Integer> t : conditional) {
				double[][] data = t.getLeft();
				int nrCols = data.length;
				extraAlleles += nrCols;
			}

			tmpgenotypes = new double[genotypes.length + extraAlleles][genotypes[0].length];

			// copy the data to the new matrix
			int ctr = 0;
			for (Triple<double[][], boolean[], Integer> t : conditional) {
				double[][] data = t.getLeft();
				boolean[] missingGT = t.getMiddle();
				int nrAllelesForSNP = data.length;
				for (int i = 0; i < nrAllelesForSNP; i++) {
					tmpgenotypes[i + ctr] = data[i];
				}

				for (int i = 0; i < missingGT.length; i++) {
					genotypeMissing[i] = genotypeMissing[i] && missingGT[i];
				}
				ctr += nrAllelesForSNP;
			}


// copy the original genotypes to the matrix
			for (int i = 0; i < genotypes.length; i++) {
				tmpgenotypes[i + ctr] = genotypes[i];
			}

			nrWithMissingGenotypes = 0;
			for (int i = 0; i < genotypeMissing.length; i++) {
				if (genotypeMissing[i]) {
					nrWithMissingGenotypes++;
				}
			}
		}

		// determine MAF
		int[] nrAllelesPresent = new int[genotypes.length + 1];
		int called = 0;
		for (int i = 0; i < genotypes[0].length; i++) { // individuals
			int nrAllelesLeft = 2;
			for (int j = 0; j < genotypes.length; j++) { // alleles
				if (!genotypeMissing[i]) {
					called += 2;
					if (genotypes[j][i] == 2d) {
						nrAllelesPresent[j + 1] += 2;
						nrAllelesLeft -= 2;
					} else if (genotypes[j][i] == 1d) {
						nrAllelesPresent[j + 1] += 1;
						nrAllelesLeft -= 1;
					}
				}
			}
			nrAllelesPresent[0] += nrAllelesLeft;
		}

		double maf = 1;
		for (int i = 0; i < nrAllelesPresent.length; i++) {
			double d = (double) nrAllelesPresent[i] / called;
			if (d < maf) {
				maf = d;
			}
		}

		// merge with covariate matrix, and remove missing genotypes
		int nrAlleles = tmpgenotypes.length;
		int nrCovariates = covariates[0].length;
		double[] outputdiseaseStatus = new double[nrsamples - nrWithMissingGenotypes];
		double[][] outputmatrixwgenotypes = new double[nrsamples - nrWithMissingGenotypes][nrAlleles + nrCovariates + 1];

		int ctr = 0;

		for (int i = 0; i < covariates.length; i++) {
			if (!genotypeMissing[i]) {

				// add the intercept
				outputmatrixwgenotypes[ctr][0] = 1;

				// copy the genotypes
				for (int a = 0; a < tmpgenotypes.length; a++) {
					double tmp = tmpgenotypes[a][i];
					outputmatrixwgenotypes[ctr][a + 1] = tmp;
				}

				// copy the covariates
				for (int a = 0; a < covariates[0].length; a++) {
					outputmatrixwgenotypes[ctr][nrAlleles + a + 1] = covariates[i][a];
				}
				outputdiseaseStatus[ctr] = diseaseStatus[i];
				ctr++;

			}
		}

		return new Pair<>(
				new Pair<>(outputmatrixwgenotypes, outputdiseaseStatus),
				new Pair<>(maf, (double) called / 2));
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


	private Pair<ArrayList<double[][]>, ArrayList<String>> overlapVariantsWithBedRegion(Feature
																								region, ArrayList<Integer> variantPositions, ArrayList<String> variantNames, ArrayList<double[][]>
																								genotypes) {

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

	public void determineSamplesToRemoveFromFam
			(HashSet<String> samplesToRemove, HashSet<String> samplesWithGenotypeAndDiseaseStatus) throws IOException {

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
