package nl.harmjanwestra.finemappingtools.gwas.tasks;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import nl.harmjanwestra.finemappingtools.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.association.AssociationResultPairwise;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.math.LogisticRegressionOptimized;
import nl.harmjanwestra.utilities.math.LogisticRegressionResult;
import nl.harmjanwestra.utilities.vcf.SampleAnnotation;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.math.stats.ChiSquare;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.Callable;

/**
 * Created by hwestra on 5/24/16.
 */
public class LRTestExhaustiveTask implements Callable<AssociationResultPairwise> {

	private final boolean[] genotypesWithCovariatesAndDiseaseStatus;

	private final LRTestOptions options;
	private final SampleAnnotation sampleAnnotation;
	private int snpid1;
	private int snpid2;
	private ArrayList<VCFVariant> variants;
	private DenseDoubleAlgebra dda = new DenseDoubleAlgebra();
	private LogisticRegressionResult resultCovars;
	private int nrCovars;

	public LRTestExhaustiveTask(ArrayList<VCFVariant> variants, int i, int j,
								boolean[] genotypesWithCovariatesAndDiseaseStatus,
								SampleAnnotation sampleAnnotation,
								LRTestOptions options) {
		this.variants = variants;
		this.snpid1 = i;
		this.snpid2 = j;
		this.genotypesWithCovariatesAndDiseaseStatus = genotypesWithCovariatesAndDiseaseStatus;
		this.sampleAnnotation = sampleAnnotation;
		this.options = options;
	}

	@Override
	public AssociationResultPairwise call() throws Exception {

		VCFVariant variant1 = variants.get(snpid1);
		VCFVariant variant2 = variants.get(snpid2);
		LRTestTask taskObj = new LRTestTask(sampleAnnotation);

		LRTestVariantQCTask lrq = new LRTestVariantQCTask();


		Triple<int[], boolean[], Integer> qcdata2 = lrq.determineMissingGenotypes(
				variant2.getGenotypeAllelesAsMatrix2D(),
				sampleAnnotation.getCovariates().rows());


		double maf1 = variant1.getMAFControls();
		double hwep1 = variant1.getHwepControls();
		double maf2 = variant2.getMAFControls();
		double hwep2 = variant2.getHwepControls();


		ArrayList<Pair<VCFVariant, Triple<int[], boolean[], Integer>>> conditional = new ArrayList<>();
		conditional.add(new Pair<>(variant2, qcdata2));

		Pair<DoubleMatrix2D, double[]> xandy = taskObj.prepareMatrices(
				variant1,
				conditional
		);

		int numberOfColumns = variant1.getAlleles().length - 1 + variant2.getAlleles().length - 1;
		double[] y = xandy.getRight();
		DoubleMatrix2D x = xandy.getLeft();

		AssociationResultPairwise output = pruneAndTest(x, y, 1, 1 + numberOfColumns, taskObj);

		if (output == null) {
			return null;
		}

		SNPFeature snp = new SNPFeature(Chromosome.parseChr(variant1.getChr()), variant1.getPos(), variant1.getPos());

		Double imputationqualityscore = variant1.getImputationQualityScore();
		if (imputationqualityscore == null) {
			imputationqualityscore = 0d;
		}
		snp.setName(variant1.getId());
		output.setSnp(snp);
		output.setN(x.rows());
		snp.setMaf(maf1);
		snp.setHwep(hwep1);
		snp.setImputationQualityScore(imputationqualityscore);
		snp.setAlleles(variant1.getAlleles());
		snp.setMinorAllele(variant1.getMinorAllele());

		SNPFeature snp2 = new SNPFeature(Chromosome.parseChr(variant2.getChr()), variant2.getPos(), variant2.getPos());
		Double imputationqualityscore2 = variant2.getImputationQualityScore();
		if (imputationqualityscore2 == null) {
			imputationqualityscore2 = 0d;
		}
		snp2.setName(variant2.getId());
		output.setSnp2(snp2);
		snp2.setMaf(maf2);
		snp2.setHwep(hwep2);
		snp2.setImputationQualityScore(imputationqualityscore2);
		snp2.setAlleles(variant2.getAlleles());
		snp2.setMinorAllele(variant2.getMinorAllele());

		// calculate the ld between the variants :)
		DetermineLD ldcalc = new DetermineLD();
		if (variant1.getNrAlleles() == 2 && variant2.getNrAlleles() == 2) {
			Pair<Double, Double> ld = ldcalc.getLD(variant1, variant2);
			output.setLDRSquared(ld.getRight());
			output.setLdDprime(ld.getLeft());
		} else {
			output.setLDRSquared(Double.NaN);
			output.setLdDprime(Double.NaN);
		}

		return output;
	}

	private AssociationResultPairwise pruneAndTest(DoubleMatrix2D x,
												   double[] y,
												   int firstGenotypeColumn,
												   int lastGenotypeColumn,
												   LRTestTask testObj) throws IOException {
		LRTestTask lrt = new LRTestTask(sampleAnnotation, options);
		if (options.debug) {
			System.out.println(x.columns() + " before pruning");
		}

		Pair<DoubleMatrix2D, boolean[]> pruned = lrt.removeCollinearVariables(x);
		if (pruned == null) {
			System.err.println("Error pruning " + variants.get(0).getId() + " and " + variants.get(1).getId());
			return null;
		}
		x = pruned.getLeft(); // x is now probably shorter than original X

		if (options.debug) {
			System.out.println(x.columns() + " after pruning...");
		}

		// figure out which columns are genotypes, and which ones are covariates
		boolean[] notaliased = pruned.getRight(); // length of original X

		ArrayList<Integer> remainingGenotypeColumns = new ArrayList<>();
		ArrayList<Integer> remainingCovariateColumns = new ArrayList<>();
		int ctr = 0;
		int[] alleleIndex = new int[lastGenotypeColumn - firstGenotypeColumn];

		boolean conditionalOk = true;

		int allelectr = 0;
		for (int i = 0; i < notaliased.length; i++) {
			if (notaliased[i]) {
				if (i >= firstGenotypeColumn && i < lastGenotypeColumn) {
					remainingGenotypeColumns.add(ctr);
					alleleIndex[allelectr] = ctr;
					allelectr++;
				} else {
					remainingCovariateColumns.add(ctr);
				}

				ctr++;
			} else {
				if (i >= firstGenotypeColumn && i < lastGenotypeColumn) {
					alleleIndex[allelectr] = -1;
					allelectr++;
				}
			}
		}

		AssociationResultPairwise result = new AssociationResultPairwise();

		if (remainingGenotypeColumns.isEmpty()) {
			return result;
		} else {
			LogisticRegressionOptimized reg = new LogisticRegressionOptimized();
			if (resultCovars == null) {
				DoubleMatrix2D xprime = dda.subMatrix(x, 0, x.rows() - 1, Primitives.toPrimitiveArr(remainingCovariateColumns.toArray(new Integer[0])));
				resultCovars = reg.univariate(y, xprime);
				if (resultCovars == null) {
					System.err.println("ERROR: null-model regression did not converge. ");
					System.err.println(x.rows() + "\t" + x.columns());
					System.err.println(xprime.rows() + "\t" + xprime.columns());
					System.err.println("-----");
					return null;
				}
				nrCovars = xprime.columns();
			}

			// perform testNormal on full model
			// remove genotypes and run testNormal on reduced model

			LogisticRegressionResult resultX = reg.univariate(y, x);
			if (resultX == null) {
				System.err.println("ERROR: did not converge. ");
				VCFVariant variant1 = variants.get(snpid1);
				VCFVariant variant2 = variants.get(snpid2);
				System.err.println("Variant1: " + variant1.getChr()
						+ "\t" + variant1.getPos()
						+ "\t" + variant1.getId()
						+ "\t" + Strings.concat(variant1.getAlleles(), Strings.comma)
						+ "\t" + variant1.getMinorAllele()
						+ "\t" + variant1.getMAF());
				System.err.println("Variant2: " + variant2.getChr()
						+ "\t" + variant2.getPos()
						+ "\t" + variant2.getId()
						+ "\t" + Strings.concat(variant2.getAlleles(), Strings.comma)
						+ "\t" + variant2.getMinorAllele()
						+ "\t" + variant2.getMAF());
				DetermineLD d = new DetermineLD();
				Pair<Double, Double> ld = d.getLD(variant1, variant2);
				System.err.println("Distance:" + Math.abs(variant1.getPos() - variant2.getPos()));
				System.err.println("LD dp:" + ld.getLeft() + " / rsq " + ld.getRight());
				System.err.println("-----");
				return null;
			}

			double devx = resultX.getDeviance();
//			DoubleMatrix2D xprime = dda.subMatrix(x, 0, x.rows() - 1, Primitives.toPrimitiveArr(colIndexArr.toArray(new Integer[0])));
//			LogisticRegressionResult resultCovars = reg.univariate(y, xprime);
//			if (resultCovars == null) {
//				System.err.println("ERROR: covariate regression did not converge. ");
//				return null;
//			}
			int nrRemaining = remainingGenotypeColumns.size();
			double devnull = resultCovars.getDeviance();
			double[] betasmlelr = new double[nrRemaining];
			double[] stderrsmlelr = new double[nrRemaining];
			double[] or = new double[nrRemaining];
			double[] orhi = new double[nrRemaining];
			double[] orlo = new double[nrRemaining];

			ctr = 0;

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
