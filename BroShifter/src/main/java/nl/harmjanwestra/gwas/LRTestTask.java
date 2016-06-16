package nl.harmjanwestra.gwas;

import JSci.maths.ArrayMath;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import nl.harmjanwestra.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.math.LogisticRegressionOptimized;
import nl.harmjanwestra.utilities.math.LogisticRegressionResult;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.math.stats.ChiSquare;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.Callable;

/**
 * Created by Harm-Jan on 04/20/16.
 */

// callable for easy multithreading..
public class LRTestTask implements Callable<Triple<String, AssociationResult, VCFVariant>> {

	private VCFVariant variant;
	private int iter;
	private DiseaseStatus[] finalDiseaseStatus;
	private DoubleMatrix2D finalCovariates;
	private ArrayList<Pair<VCFVariant, Triple<int[], boolean[], Triple<Integer, Double, Double>>>> conditional;
	private int alleleOffsetGenotypes;
	private LRTestOptions options;
	private DenseDoubleAlgebra dda = new DenseDoubleAlgebra();
	private LogisticRegressionResult resultCovars;
	private int nrCovars;


	public LRTestTask() {
	}

	public LRTestTask(VCFVariant variant,
	                  int iter,
	                  DiseaseStatus[] finalDiseaseStatus,
	                  DoubleMatrix2D finalCovariates,
	                  ArrayList<Pair<VCFVariant, Triple<int[], boolean[], Triple<Integer, Double, Double>>>> conditional,
	                  int alleleOffsetGenotypes,
	                  LRTestOptions options) {
		this.variant = variant;
		this.iter = iter;
		this.finalDiseaseStatus = finalDiseaseStatus;
		this.finalCovariates = finalCovariates;
		this.conditional = conditional;
		this.alleleOffsetGenotypes = alleleOffsetGenotypes;
		this.options = options;
	}

	@Override
	public Triple<String, AssociationResult, VCFVariant> call() {
		// recode the genotypes to the same ordering as the covariate table
		LRTestVariantQCTask lrq = new LRTestVariantQCTask();
		Triple<int[], boolean[], Triple<Integer, Double, Double>> qcdata = lrq.filterAndRecodeGenotypes(
				variant.getGenotypeAllelesAsMatrix2D(),
				finalDiseaseStatus,
				variant.getAlleles().length,
				finalCovariates.rows());

		Triple<Integer, Double, Double> stats = qcdata.getRight();
		double maf = stats.getMiddle();
		double hwep = stats.getRight();

		// generate pseudocontrol genotypes
		Pair<DoubleMatrix2D, double[]> xandy = prepareMatrices(
				variant,
				qcdata.getLeft(),
				conditional,
				finalDiseaseStatus,
				finalCovariates
		);

		int nrAlleles = variant.getAlleles().length;
		double[] y = xandy.getRight(); // get the phenotypes for all non-missing genotypes
		DoubleMatrix2D x = xandy.getLeft();


		AssociationResult result = pruneAndTest(x, y, 1, 1 + (nrAlleles - 1), variant, maf);

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
//		} // end maf > threshold
	}

	public Pair<DoubleMatrix2D, double[]> prepareMatrices(VCFVariant variant, int[] left, ArrayList<Pair<VCFVariant, Triple<int[], boolean[], Triple<Integer, Double, Double>>>> conditional, DiseaseStatus[] finalDiseaseStatus, DoubleMatrix2D finalCovariates) {

		// prepare genotype matrix
		DoubleMatrix2D x = variant.getDosagesAsMatrix2D();
		HashSet<Integer> missingGenotypeIds = new HashSet<Integer>();
		for (int i : left) {
			missingGenotypeIds.add(i);
		}

		for (int c = 0; c < conditional.size(); c++) {
			Pair<VCFVariant, Triple<int[], boolean[], Triple<Integer, Double, Double>>> conditionalData = conditional.get(c);
			VCFVariant var2 = conditionalData.getLeft();
			int[] left2 = conditionalData.getRight().getLeft();

			DoubleMatrix2D dosages = var2.getDosagesAsMatrix2D();

			for (int i : left2) {
				missingGenotypeIds.add(i);
			}
			x = DoubleFactory2D.dense.appendColumns(x, dosages);
		}

		// now append covariates
		x = DoubleFactory2D.dense.appendColumns(x, finalCovariates);

		// add intercept column
		DoubleMatrix2D vector = DoubleFactory2D.dense.make(x.rows(), 1);
		vector.assign(1);
		x = DoubleFactory2D.dense.appendColumns(vector, x);

		if (missingGenotypeIds.isEmpty()) {
			double[] y = new double[finalDiseaseStatus.length];
			for (int i = 0; i < finalDiseaseStatus.length; i++) {
				y[i] = finalDiseaseStatus[i].getNumber();
			}
			return new Pair<>(x, y);
		} else {
			// filter missing samples
			Integer[] rowIndexesArr = missingGenotypeIds.toArray(new Integer[0]);
			int[] rowIndexes = Primitives.toPrimitiveArr(rowIndexesArr);

			x = dda.subMatrix(x, rowIndexes, 0, x.columns() - 1);

			double[] y = new double[finalDiseaseStatus.length - missingGenotypeIds.size()];
			int ctr = 0;
			for (int i = 0; i < finalDiseaseStatus.length; i++) {
				if (!missingGenotypeIds.contains(i)) {
					y[ctr] = finalDiseaseStatus[i].getNumber();
					ctr++;
				}
			}
			return new Pair<>(x, y);
		}
	}


	// remove variables with zero variance and perfect correlation
	public Pair<DoubleMatrix2D, boolean[]> removeCollinearVariables(DoubleMatrix2D mat) {

		DoubleMatrix2D corrmat = new DenseDoubleMatrix2D(mat.columns(), mat.columns());
		for (int i = 0; i < mat.columns(); i++) {
			for (int j = i; j < mat.columns(); j++) {
				double c = Correlation.correlate(mat.viewColumn(i).toArray(), mat.viewColumn(j).toArray());
				if (Double.isNaN(c)) {
					corrmat.setQuick(i, j, 0);
				} else {
					corrmat.setQuick(i, j, c);
				}
			}
		}


		boolean[] includeCol = new boolean[mat.columns()];
		for (int c = 0; c < includeCol.length; c++) {
			includeCol[c] = true;
		}

		for (int j = 1; j < mat.columns(); j++) {
			double[] col1 = mat.viewColumn(j).toArray();
			if (ArrayMath.variance(col1) == 0) {
				includeCol[j] = false;
			} else {
				for (int j2 = j + 1; j2 < mat.columns(); j2++) {
					if (includeCol[j2]) {
						double c = corrmat.getQuick(j, j2);
						double[] col2 = mat.viewColumn(2).toArray();
						if (ArrayMath.variance(col2) == 0) {
							includeCol[j2] = false;
						} else if (Math.abs(c) == 1d) {
							includeCol[j] = false;  // prefer to remove the earlier column
						}
					}
				}
			}
		}

		ArrayList<Integer> remaining = new ArrayList<>();
		for (int c = 0; c < includeCol.length; c++) {
			if (includeCol[c]) {
				remaining.add(c);
			}
		}

		int[] colindexes = Primitives.toPrimitiveArr(remaining.toArray(new Integer[0]));
		// remove correlated variables
		DoubleMatrix2D matOut = dda.subMatrix(mat, 0, mat.rows() - 1, colindexes);
		return new Pair<>(matOut, includeCol);
	}

	private AssociationResult pruneAndTest(DoubleMatrix2D x,
	                                       double[] y,
	                                       int firstColumnToRemove,
	                                       int lastColumnToRemove,
	                                       VCFVariant variant,
	                                       double maf) {



		Pair<DoubleMatrix2D, boolean[]> pruned = removeCollinearVariables(x);
		x = pruned.getLeft(); // x is now probably shorter than original X

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

		if (nrRemaining == 0) {
			result.setDevianceNull(0);
			result.setDevianceGeno(0);
			result.setDf(0);
			result.setBeta(new double[]{0d});
			result.setSe(new double[]{0d});
			result.setPval(1);
			return result;
		} else {
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

			if (resultCovars == null) {
				DoubleMatrix2D xprime = dda.subMatrix(x, 0, x.rows() - 1, Primitives.toPrimitiveArr(colIndexArr.toArray(new Integer[0])));
				resultCovars = reg.univariate(y, xprime);
				if (resultCovars == null) {
					System.err.println("ERROR: null-model regression did not converge. ");
					System.err.println("Variant: " + snp.getChromosome().toString()
							+ "\t" + snp.getStart()
							+ "\t" + snp.getName()
							+ "\t" + Strings.concat(snp.getAlleles(), Strings.comma)
							+ "\t" + snp.getMinorAllele()
							+ "\t" + maf);
					System.err.println(x.rows() + "\t" + x.columns());
					System.err.println(xprime.rows() + "\t" + xprime.columns());
					System.err.println("-----");
					return null;
				}
				nrCovars = xprime.columns();
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
			int df = x.columns() - nrCovars;
			double p = ChiSquare.getP(df, deltaDeviance);
//				log10p = Math.abs((-Math.log10(p)));


			result.setDevianceNull(devnull);
			result.setDevianceGeno(devx);
			result.setDf(df);
			result.setDfalt(x.columns());
			result.setDfnull(nrCovars);

			result.setBeta(betasmlelr);
			result.setSe(stderrsmlelr);
			result.setPval(p);
			return result;
		}
	}

	public void setResultNullmodel(LogisticRegressionResult r, int nrvars) {
		this.resultCovars = r;
		this.nrCovars = nrvars;

	}
}

