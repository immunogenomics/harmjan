package nl.harmjanwestra.gwas;

import JSci.maths.ArrayMath;
import nl.harmjanwestra.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.math.LogisticRegression;
import nl.harmjanwestra.utilities.math.LogisticRegressionResult;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.math.stats.ChiSquare;

import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.Callable;

/**
 * Created by Harm-Jan on 04/20/16.
 */

// callable for easy multithreading..
public class LRTestTask implements Callable<Triple<String, AssociationResult, VCFVariant>> {

	private VCFVariant variant;
	private int iter;
	private boolean[] genotypesWithCovariatesAndDiseaseStatus;
	private double[] finalDiseaseStatus;
	private double[][] finalCovariates;
	private ArrayList<Triple<double[][], boolean[], Integer>> conditional;
	private ArrayList<double[][]> conditionalDosages;
	private int alleleOffsetGenotypes;
	private int alleleOffsetDosages;
	private LRTestOptions options;
	private String vcfLn;


	public LRTestTask() {
	}

	public LRTestTask(String vcfLn,
					  VCFVariant variant,
					  int iter,
					  boolean[] genotypesWithCovariatesAndDiseaseStatus,
					  double[] finalDiseaseStatus,
					  double[][] finalCovariates,
					  ArrayList<Triple<double[][], boolean[], Integer>> conditional, ArrayList<double[][]> conditionalDosages,
					  int alleleOffsetGenotypes,
					  int alleleOffsetDosages,
					  LRTestOptions options) {
		this.variant = variant;
		this.vcfLn = vcfLn;
		this.iter = iter;
		this.genotypesWithCovariatesAndDiseaseStatus = genotypesWithCovariatesAndDiseaseStatus;
		this.finalDiseaseStatus = finalDiseaseStatus;
		this.finalCovariates = finalCovariates;
		this.conditional = conditional;
		this.conditionalDosages = conditionalDosages;
		this.alleleOffsetGenotypes = alleleOffsetGenotypes;
		this.alleleOffsetDosages = alleleOffsetDosages;
		this.options = options;
	}

	@Override
	public Triple<String, AssociationResult, VCFVariant> call() throws Exception {
		if (iter == 0) {
			// do some more parsing if this is the first time we're seeing this variant...
			variant = new VCFVariant(vcfLn, VCFVariant.PARSE.ALL);
			vcfLn = null;
		}

		// TODO: switch this around: make the ordering of the covariate table the same as the genotype file...
		// recode the genotypes to the same ordering as the covariate table
		Triple<double[][], boolean[], Integer> unfilteredGenotypeData = filterAndRecodeGenotypes(
				genotypesWithCovariatesAndDiseaseStatus,
				variant.getGenotypeAlleles(),
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
			if (iter == 0) {
//								SNP	Chr	Pos	ImputationQual	MAF	OverlapOK	MAFOk	ImpQualOK
				String output = variant.getId()
						+ "\t" + variant.getChr()
						+ "\t" + variant.getPos()
						+ "\t" + variant.getImputationQualityScore()
						+ "\t" + true
						+ "\t" + true
						+ "\t" + false
						+ "\t" + true;

				return new Triple<>(output, null, null);
//					logout.writeln(
				//);
			} else {
				return null;
			}
		} else {

//				if (iter == 0 && maxNrIter > 1) {
//					variants.add(variant);
//				}

			double[][] probs = variant.getGenotypeProbabilies(); // [probs][inds]
			int nrAlleles = variant.getAlleles().length;
			double[] y = genotypedata.getRight();
			double[][] x = genotypedata.getLeft();

			AssociationResult result = null;
			if (probs == null) {
				result = pruneAndTest(x, y, nrAlleles, alleleOffsetGenotypes, variant, maf);
//
			} else {
				// recode to dosages
				double[][] dosages = variant.getImputedDosages();
				double[][] dosageMat = prepareDosageMatrix(genotypesWithCovariatesAndDiseaseStatus, dosages, finalCovariates, conditionalDosages);
				result = pruneAndTest(dosageMat, y, nrAlleles, alleleOffsetDosages, variant, maf);

			}
//				nrTested++;
			//								SNP	Chr	Pos	ImputationQual	MAF	OverlapOK	MAFOk	ImpQualOK
			String output = variant.getId()
					+ "\t" + variant.getChr()
					+ "\t" + variant.getPos()
					+ "\t" + variant.getImputationQualityScore()
					+ "\t" + true
					+ "\t" + true
					+ "\t" + true
					+ "\t" + true;
			return new Triple<>(output, result, variant);
		} // end maf > threshold
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

	private double[][] removeGenotypes(double[][] x, boolean[] colsToRemove) {

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

	Triple<double[][], boolean[], Integer> filterAndRecodeGenotypes(
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

	private AssociationResult pruneAndTest(double[][] x,
										   double[] y,
										   int nrAlleles,
										   int alleleOffset,
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

		Double imputationqualityscore = variant.getImputationQualityScore();
		result.setImputationQualScore(imputationqualityscore);

		if (nrRemaining > 0) {
			// perform test on full model
			// remove genotypes and run test on reduced model
//				if (useR) {
//					// debug: test with R
//					rConnection.assign("y", y);
//					assignAsRMatrix(rConnection, x, "x", false);
//					double[][] covarsOnly = removeGenotypes(x, colsToRemove);
//					assignAsRMatrix(rConnection, covarsOnly, "z", false);
//
//					rConnection.voidEval("glm1 <- glm(y ~ x, family=binomial(logit))");
//					rConnection.voidEval("glm2 <- glm(y ~ z, family=binomial(logit))");
//					rConnection.voidEval("sum1 <- summary(glm1)");
//					rConnection.voidEval("sum2 <- summary(glm2)");
//					rConnection.voidEval("deva <- sum1$deviance");
//					rConnection.voidEval("devn <- sum2$deviance");
//					double devr = rConnection.eval("sum1$deviance").asDouble();
//					rConnection.voidEval("deltaDeviance <- devn - deva");
//
//					rConnection.voidEval("altdf <- sum1$df[1]");
//					rConnection.voidEval("nulldf <- sum2$df[1]");
//					rConnection.voidEval("dfdiff <- altdf - nulldf");
//					Double pvalR = rConnection.eval("pchisq(deltaDeviance, df = dfdiff, lower.tail = FALSE, log.p = FALSE)").asDouble();
//					Double dfR = rConnection.eval("dfdiff").asDouble();
//					Double deltaR = rConnection.eval("deltaDeviance").asDouble();
//
//					double[] betas = new double[nrRemaining];
//
//					int ctr = 0;
//					for (int i = 0; i < alleleIndex.length; i++) {
//						int idx = alleleIndex[i];
//						if (idx != -1) {
//							int idxplusone = idx + 1;
//							betas[ctr] = rConnection.eval("as.vector(sum1$coefficients[" + idxplusone + ",])[1]").asDouble();
//							ctr++;
//						}
//					}
//
//					String betaR = Strings.concat(betas, Strings.semicolon);
//
//
//				} else {

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
//				log10p = Math.abs((-Math.log10(p)));


			result.setDevianceNull(devnull);
			result.setDevianceGeno(devx);
			result.setDf(df);

			result.setBeta(betasmlelr);
			result.setSe(stderrsmlelr);
			result.setPval(p);
			return result;
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
//					out.writeln(result.toString());
//				}
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
//				out.writeln(result.toString());
			return result;

		}
	}
}

