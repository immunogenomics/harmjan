package nl.harmjanwestra.gwas;

import JSci.maths.ArrayMath;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import nl.harmjanwestra.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.math.LogisticRegressionOptimized;
import nl.harmjanwestra.utilities.math.LogisticRegressionResult;
import nl.harmjanwestra.utilities.matrix.ByteMatrix2D;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.math.stats.ChiSquare;
import umcg.genetica.text.Strings;

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
	private DiseaseStatus[] finalDiseaseStatus;
	private DoubleMatrix2D finalCovariates;
	private ArrayList<Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>>> conditional;
	private ArrayList<DoubleMatrix2D> conditionalDosages;
	private int alleleOffsetGenotypes;
	private int alleleOffsetDosages;
	private LRTestOptions options;


	public LRTestTask() {
	}

	public LRTestTask(String vcfLn,
					  VCFVariant variant,
					  int iter,
					  boolean[] genotypesWithCovariatesAndDiseaseStatus,
					  DiseaseStatus[] finalDiseaseStatus,
					  DoubleMatrix2D finalCovariates,
					  ArrayList<Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>>> conditional,
					  ArrayList<DoubleMatrix2D> conditionalDosages,
					  int alleleOffsetGenotypes,
					  int alleleOffsetDosages,
					  LRTestOptions options) {
		this.variant = variant;

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

		// TODO: switch this around: make the ordering of the covariate table the same as the genotype file...
		// recode the genotypes to the same ordering as the covariate table
		Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>> unfilteredGenotypeData = filterAndRecodeGenotypes(
				genotypesWithCovariatesAndDiseaseStatus,
				variant.getGenotypeAllelesAsMatrix2D(),
				finalDiseaseStatus,
				variant.getAlleles().length,
				finalCovariates.rows());

		Triple<Integer, Double, Double> stats = unfilteredGenotypeData.getRight();
		double maf = stats.getMiddle();
		double hwep = stats.getRight();

		if (maf < options.getMafthresholdD() && hwep < options.getHWEPThreshold()) {

			if (iter == 0) {
				String output = variant.getId()
						+ "\t" + variant.getChr()
						+ "\t" + variant.getPos()
						+ "\t" + Strings.concat(variant.getAlleles(), Strings.comma)
						+ "\t" + variant.getMinorAllele()
						+ "\t" + variant.getImputationQualityScore()
						+ "\t" + maf
						+ "\t" + hwep
						+ "\t" + false;

				return new Triple<>(output, null, null);
			} else {
				return new Triple<>(null, null, null);
			}
		} else {
			// generate pseudocontrol genotypes
			Pair<DoubleMatrix2D, double[]> genotypedata = prepareGenotypeMatrix(
					unfilteredGenotypeData,
					conditional,
					finalDiseaseStatus,
					finalCovariates,
					finalCovariates.rows());

			int nrAlleles = variant.getAlleles().length;
			double[] y = genotypedata.getRight();
			DoubleMatrix2D x = genotypedata.getLeft();

			AssociationResult result = null;
			if (!variant.hasImputationDosages()) {
				result = pruneAndTest(x, y, nrAlleles, alleleOffsetGenotypes, variant, maf);
			} else {
				// recode to dosages
				DoubleMatrix2D dosageMat = prepareDosageMatrix(genotypesWithCovariatesAndDiseaseStatus,
						variant.getImputedDosagesAsMatrix2D(),
						finalCovariates,
						conditionalDosages);
				result = pruneAndTest(dosageMat, y, nrAlleles, alleleOffsetDosages, variant, maf);
			}

			if (result == null) {
				return new Triple<>(null, null, null);
			}
			result.setHWEP(hwep);
			//								SNP	Chr	Pos	ImputationQual	MAF	OverlapOK	MAFOk	ImpQualOK
			String output = variant.getId()
					+ "\t" + variant.getChr()
					+ "\t" + variant.getPos()
					+ "\t" + Strings.concat(variant.getAlleles(), Strings.comma)
					+ "\t" + variant.getMinorAllele()
					+ "\t" + variant.getImputationQualityScore()
					+ "\t" + maf
					+ "\t" + hwep
					+ "\t" + true;

			return new Triple<>(output, result, variant);
		} // end maf > threshold
	}

	// remove variables with zero variance and perfect correlation
	public Pair<DoubleMatrix2D, boolean[]> removeCollinearVariables(DoubleMatrix2D mat) {

		boolean[] includeCol = new boolean[mat.columns()];
		for (int c = 0; c < includeCol.length; c++) {
			includeCol[c] = true;
		}

		includeCol[0] = true; // intercept
		for (int j = 1; j < mat.columns(); j++) {
			if (includeCol[j]) {
				double[] vals = new double[mat.rows()];
				for (int i = 0; i < mat.rows(); i++) {
					vals[i] = mat.getQuick(i, j);
				}

				double variance = ArrayMath.variance(vals);
				if (variance == 0d) {
					includeCol[j] = false;
				} else {
					for (int j2 = j + 1; j2 < mat.columns(); j2++) {
						if (includeCol[j2]) {
							double[] vals2 = new double[mat.rows()];
							for (int i = 0; i < mat.rows(); i++) {
								vals2[i] = mat.getQuick(i, j2);
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
		DoubleMatrix2D matOut = new DenseDoubleMatrix2D(mat.rows(), nrRemaining);
		int ctr = 0;
		for (int j = 0; j < mat.columns(); j++) {
			if (includeCol[j]) {
				for (int i = 0; i < mat.rows(); i++) {
					matOut.setQuick(i, ctr, mat.getQuick(i, j));
				}
				ctr++;
			}
		}


		return new Pair<>(matOut, includeCol);
	}

	public DoubleMatrix2D removeGenotypes(DoubleMatrix2D x, boolean[] colsToRemove) {

		int nrToRemove = 0;
		for (int i = 0; i < colsToRemove.length; i++) {
			if (colsToRemove[i]) {
				nrToRemove++;
			}
		}

		DoubleMatrix2D output = new DenseDoubleMatrix2D(x.rows(), colsToRemove.length - nrToRemove);
		int ctr = 0;
		for (int j = 0; j < colsToRemove.length; j++) {
			if (!colsToRemove[j]) {
				for (int i = 0; i < x.rows(); i++) {
					output.setQuick(i, ctr, x.getQuick(i, j));
				}
				ctr++;
			}
		}
		return output;
	}

	Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>> filterAndRecodeGenotypes(
			boolean[] includeGenotype,
			ByteMatrix2D genotypeAlleles,
			DiseaseStatus[] diseaseStatus,
			int nrAlleles,
			int nrsamples) {

		DoubleMatrix2D tmpgenotypes = new DenseDoubleMatrix2D(nrAlleles - 1, nrsamples);
		int individualCounter = 0;

		// first iterate the genotyped samples to load the genotypes
		// calculate maf and hwep

		// initialize hwep stuff
		int nrCombinations = (nrAlleles * (nrAlleles + 1)) / 2;
		int[] obs = new int[nrCombinations];
		int[] nrHomozygous = new int[nrAlleles];
		int nrCalled = 0;
		int[][] index = new int[nrAlleles][nrAlleles];

		int ctr = 0;
		for (int i = 0; i < nrAlleles; i++) {
			for (int j = i; j < nrAlleles; j++) {
				index[i][j] = ctr;
				index[j][i] = ctr;
				ctr++;
			}
		}

		int nrWithMissingGenotypes = 0;
		boolean[] genotypeMissing = new boolean[nrsamples];
		for (int i = 0; i < genotypeAlleles.columns(); i++) {
			if (includeGenotype[i]) { // this is set to false if the genotype doesn't have a disease status or covariate data.
				byte b1 = genotypeAlleles.getQuick(0, i);
				byte b2 = genotypeAlleles.getQuick(1, i);

				if (b1 == -1) {
					for (int q = 0; q < tmpgenotypes.rows(); q++) {
						tmpgenotypes.setQuick(q, individualCounter, Double.NaN);
					}
					nrWithMissingGenotypes++;
					genotypeMissing[individualCounter] = true;
				} else {
					if (diseaseStatus[individualCounter].equals(DiseaseStatus.CONTROL)) { // controls
						if (b1 == b2) {
							nrHomozygous[b1]++;
						}

						int id = index[b1][b2];
						obs[id]++;
						nrCalled++;
					}
					if (b1 == b2) {
						// homozygote
						if (b1 == 0) {
							// do nothing
						} else {
							int allele = b1 - 1;
							if (allele >= 0) {
								tmpgenotypes.set(allele, individualCounter, 2);
							}
						}
					} else {
						int allele1 = b1 - 1;
						int allele2 = b2 - 1;
						if (allele1 >= 0) {
							tmpgenotypes.set(allele1, individualCounter, 1);
						}
						if (allele2 >= 0) {
							tmpgenotypes.set(allele2, individualCounter, 1);
						}
					}
				}
				individualCounter++;
			}
		}

		double[] freqs = new double[nrAlleles];
		for (int i = 0; i < nrAlleles; i++) {
			freqs[i] = Math.sqrt((double) nrHomozygous[i] / nrCalled);
		}

		ctr = 0;
		double chisq = 0;
		for (int i = 0; i < nrAlleles; i++) {
			for (int j = i; j < nrAlleles; j++) {
				double expectedFreq;
				if (i == j) {
					expectedFreq = (freqs[i] * freqs[i]) * nrCalled; // homozygote freq
				} else {
					expectedFreq = (2 * (freqs[i] * freqs[j])) * nrCalled; // heterozygote freq
				}
				double observedFreq = obs[ctr];
				double oe = (observedFreq - expectedFreq);
				if (oe > 0 && expectedFreq > 0) {
					chisq += ((oe * oe) / expectedFreq);
				}
				ctr++;
			}
		}

		int df = (nrCombinations - nrAlleles);
		double hwep = umcg.genetica.math.stats.ChiSquare.getP(df, chisq);

		int[] nrAllelesPresent = new int[tmpgenotypes.rows() + 1];
		int called = 0;
		for (int i = 0; i < tmpgenotypes.columns(); i++) { // individuals
			int nrAllelesLeft = 2;
			for (int j = 0; j < tmpgenotypes.rows(); j++) { // alleles
				if (!genotypeMissing[i]) {
					double gji = tmpgenotypes.getQuick(j, i);
					called += 2;

					if (gji == 2d) {
						nrAllelesPresent[j + 1] += 2;
						nrAllelesLeft -= 2;
					} else if (gji == 1d) {
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


		return new Triple<>(tmpgenotypes, genotypeMissing, new Triple<>(nrWithMissingGenotypes, maf, hwep));
	}

	// output format genotypes+covariates, matched disease status, maf, cr
	public Pair<DoubleMatrix2D, double[]> prepareGenotypeMatrix(
			Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>> genotypeData,
			ArrayList<Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>>> conditional,
			DiseaseStatus[] diseaseStatus,
			DoubleMatrix2D covariates,
			int nrsamples) {

		int nrWithMissingGenotypes = genotypeData.getRight().getLeft();
		boolean[] genotypeMissing = genotypeData.getMiddle();
		DoubleMatrix2D genotypes = genotypeData.getLeft();
		DoubleMatrix2D tmpgenotypes = genotypes;

		if (!conditional.isEmpty()) {
			// count extra columns from conditional genotypes
			int extraAlleles = 0;
			for (Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>> t : conditional) {
				DoubleMatrix2D data = t.getLeft();
				int nrCols = data.rows();
				extraAlleles += nrCols;
			}

			tmpgenotypes = new DenseDoubleMatrix2D(genotypes.rows() + extraAlleles, genotypes.columns());

			// copy the data to the new matrix
			int ctr = 0;
			for (Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>> t : conditional) {
				DoubleMatrix2D data = t.getLeft();
				boolean[] missingGT = t.getMiddle();
				int nrAllelesForSNP = data.rows();
				for (int i = 0; i < nrAllelesForSNP; i++) {
					for (int j = 0; j < data.columns(); j++) {
						tmpgenotypes.setQuick(i + ctr, j, data.getQuick(i, j));
					}
				}

				for (int i = 0; i < missingGT.length; i++) {
					genotypeMissing[i] = genotypeMissing[i] && missingGT[i];
				}
				ctr += nrAllelesForSNP;
			}


// copy the original genotypes to the matrix
			for (int i = 0; i < genotypes.rows(); i++) {
				for (int j = 0; j < genotypes.columns(); j++) {
					tmpgenotypes.setQuick(i + ctr, j, genotypes.getQuick(i, j));
				}
			}

			nrWithMissingGenotypes = 0;
			for (int i = 0; i < genotypeMissing.length; i++) {
				if (genotypeMissing[i]) {
					nrWithMissingGenotypes++;
				}
			}
		}

		// merge with covariate matrix, and remove missing genotypes
		int nrAlleles = tmpgenotypes.rows();
		int nrCovariates = covariates.columns();
		double[] outputdiseaseStatus = new double[nrsamples - nrWithMissingGenotypes];
		DoubleMatrix2D outputmatrixwgenotypes = new DenseDoubleMatrix2D(nrsamples - nrWithMissingGenotypes, nrAlleles + nrCovariates + 1);

		int ctr = 0;

		for (int i = 0; i < covariates.rows(); i++) {
			if (!genotypeMissing[i]) {
				// add the intercept
				outputmatrixwgenotypes.setQuick(ctr, 0, 1);

				// copy the genotypes
				for (int a = 0; a < tmpgenotypes.rows(); a++) {
					double tmp = tmpgenotypes.getQuick(a, i);
					outputmatrixwgenotypes.setQuick(ctr, a + 1, tmp);
				}

				// copy the covariates
				for (int a = 0; a < covariates.columns(); a++) {
					outputmatrixwgenotypes.setQuick(ctr, nrAlleles + a + 1, covariates.getQuick(i, a));
				}
				outputdiseaseStatus[ctr] = diseaseStatus[i].getNumber();
				ctr++;
			}
		}


//
		return new Pair<>(outputmatrixwgenotypes, outputdiseaseStatus);

	}

//	private double calculateHWEP(DoubleMatrix2D genotypes, double[] diseaseStatus, int nrAlleles) {
//
//
//		return umcg.genetica.math.stats.ChiSquare.getP(df, chisq);
//	}

	public DoubleMatrix2D prepareDosageMatrix(boolean[] includeGenotype,
											  DoubleMatrix2D dosages,
											  DoubleMatrix2D covariates,
											  ArrayList<DoubleMatrix2D> conditional) {

		int nrSamples = 0;
		for (int i = 0; i < dosages.rows(); i++) {
			if (includeGenotype[i]) {
				nrSamples++;
			}
		}
		int offset = 1;

		int extraAlleles = 0;
		for (int i = 0; i < conditional.size(); i++) {
			DoubleMatrix2D data = conditional.get(i);
			extraAlleles += data.columns();
		}

		DoubleMatrix2D matrix = new DenseDoubleMatrix2D(nrSamples, dosages.columns() + covariates.columns() + extraAlleles + offset);

		int nrGenotypeSamples = dosages.rows();

		// copy conditional alleles
		int addedAlleles = 0;
		for (int z = 0; z < conditional.size(); z++) {
			DoubleMatrix2D data = conditional.get(z);
			int ctr = 0;
			for (int i = 0; i < nrGenotypeSamples; i++) {
				if (includeGenotype[i]) {
					for (int j = 0; j < data.columns(); j++) {
						matrix.setQuick(ctr, offset + addedAlleles + j, data.getQuick(i, j));
					}
					ctr++;
				}
			}
			addedAlleles += data.columns();
		}

		int ctr = 0;
		int nrAlleles = dosages.columns();
		for (int i = 0; i < nrGenotypeSamples; i++) {
			if (includeGenotype[i]) {
				matrix.setQuick(ctr, 0, 1); // intercept

				for (int j = 0; j < nrAlleles; j++) {
					matrix.setQuick(ctr, j + extraAlleles + offset, dosages.getQuick(i, j));
				}

				for (int j = 0; j < covariates.columns(); j++) {
					matrix.setQuick(ctr, j + extraAlleles + offset + nrAlleles, covariates.getQuick(ctr, j));
				}
				ctr++;
			}
		}
		return matrix;
	}

	private AssociationResult pruneAndTest(DoubleMatrix2D x,
										   double[] y,
										   int nrAlleles,
										   int alleleOffset,
										   VCFVariant variant,
										   double maf) throws REXPMismatchException, REngineException, IOException {
		Pair<DoubleMatrix2D, boolean[]> pruned = removeCollinearVariables(x);
		x = pruned.getLeft(); // x is now probably shorter than original X
		boolean[] notaliased = pruned.getRight(); // length of original X

		// check if the alleles we put in are aliased
		int firstAllele = 1 + alleleOffset; // intercept + allleles we conditioned on (original X indexing)
		int lastAllele = 1 + alleleOffset + (nrAlleles - 1); // intercept + allleles we conditioned on + alleles for this variant (original X indexing)


		// index new position in X
		int nrRemaining = 0;
		int[] alleleIndex = new int[nrAlleles - 1];
		boolean[] colsToRemove = new boolean[x.columns()];
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
		SNPFeature snp = new SNPFeature(Chromosome.parseChr(variant.getChr()), variant.getPos(), variant.getPos());
		snp.setName(variant.getId());
		result.setSnp(snp);
		result.setN(x.rows());
		result.setMaf(maf);

		Double imputationqualityscore = variant.getImputationQualityScore();
		snp.setImputationQualityScore(imputationqualityscore);
		snp.setAlleles(variant.getAlleles());
		snp.setMinorAllele(variant.getMinorAllele());


		if (nrRemaining > 0) {

			LogisticRegressionOptimized reg = new LogisticRegressionOptimized();
			LogisticRegressionResult resultX = reg.univariate(y, x);
			if (resultX == null) {

				// try once more with some more iterations
				LogisticRegressionOptimized reg2 = new LogisticRegressionOptimized(1000);
				resultX = reg2.univariate(y, x);
				if (resultX == null) {
					System.err.println("ERROR: did not converge.");
					System.err.println("Variant: " + snp.getChromosome().toString()
							+ "\t" + snp.getStart()
							+ "\t" + snp.getName()
							+ "\t" + Strings.concat(snp.getAlleles(), Strings.comma)
							+ "\t" + snp.getMinorAllele()
							+ "\t" + maf);
					System.err.println("-----");
					return null;
				}
			}

			double devx = resultX.getDeviance();
			DoubleMatrix2D covarsOnly = removeGenotypes(x, colsToRemove);
			LogisticRegressionResult resultCovars = reg.univariate(y, covarsOnly);
			if (resultCovars == null) {
				System.err.println("ERROR: covariate regression did not converge. ");
				return null;
			}
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
			int df = x.columns() - covarsOnly.columns();
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

