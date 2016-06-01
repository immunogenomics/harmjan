package nl.harmjanwestra.gwas;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.math.LogisticRegressionOptimized;
import nl.harmjanwestra.utilities.math.LogisticRegressionResult;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.math.stats.ChiSquare;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.Callable;

/**
 * Created by hwestra on 5/24/16.
 */
public class LRTestExhaustiveTask implements Callable<AssociationResult> {

	private final boolean[] genotypesWithCovariatesAndDiseaseStatus;
	private final DiseaseStatus[] finalDiseaseStatus;
	private final DoubleMatrix2D finalCovariates;
	private final LRTestOptions options;
	private int snpid1;
	private int snpid2;
	private ArrayList<VCFVariant> variants;

	public LRTestExhaustiveTask(ArrayList<VCFVariant> variants, int i, int j,
								boolean[] genotypesWithCovariatesAndDiseaseStatus,
								DiseaseStatus[] finalDiseaseStatus,
								DoubleMatrix2D finalCovariates,
								LRTestOptions options) {
		this.variants = variants;
		this.snpid1 = i;
		this.snpid2 = j;
		this.genotypesWithCovariatesAndDiseaseStatus = genotypesWithCovariatesAndDiseaseStatus;
		this.finalDiseaseStatus = finalDiseaseStatus;
		this.finalCovariates = finalCovariates;
		this.options = options;
	}

	@Override
	public AssociationResult call() throws Exception {

		VCFVariant variant1 = variants.get(snpid1);
		VCFVariant variant2 = variants.get(snpid2);
		LRTestTask taskObj = new LRTestTask();

		Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>> unfilteredGenotypeData1 = taskObj.filterAndRecodeGenotypes(
				genotypesWithCovariatesAndDiseaseStatus,
				variant1.getGenotypeAllelesAsMatrix2D(),
				finalDiseaseStatus,
				variant1.getAlleles().length,
				finalCovariates.rows());

		Triple<Integer, Double, Double> stats1 = unfilteredGenotypeData1.getRight();
		double maf1 = stats1.getMiddle();
		double hwep1 = stats1.getRight();

		Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>> unfilteredGenotypeData2 = taskObj.filterAndRecodeGenotypes(
				genotypesWithCovariatesAndDiseaseStatus,
				variant2.getGenotypeAllelesAsMatrix2D(),
				finalDiseaseStatus,
				variant2.getAlleles().length,
				finalCovariates.rows());

		Triple<Integer, Double, Double> stats2 = unfilteredGenotypeData2.getRight();
		double maf2 = stats2.getMiddle();
		double hwep2 = stats2.getRight();

		ArrayList<Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>>> conditional = new ArrayList<>();
		ArrayList<DoubleMatrix2D> conditionalDosages = new ArrayList<>();
		conditional.add(unfilteredGenotypeData2);
		DoubleMatrix2D dosages2 = variant2.getImputedDosagesAsMatrix2D();
		conditionalDosages.add(dosages2);

		Pair<DoubleMatrix2D, double[]> genotypedata = taskObj.prepareGenotypeMatrix(
				unfilteredGenotypeData1,
				conditional,
				finalDiseaseStatus,
				finalCovariates,
				finalCovariates.rows());

		int nrAlleles = variant1.getAlleles().length + variant2.getAlleles().length;
		double[] y = genotypedata.getRight();
		DoubleMatrix2D x = genotypedata.getLeft();

		AssociationResult output = null;
		if (!variant1.hasImputationDosages() || !variant2.hasImputationDosages()) {
			output = pruneAndTest(x, y, nrAlleles, taskObj);
		} else {
			DoubleMatrix2D dosageMat = taskObj.prepareDosageMatrix(genotypesWithCovariatesAndDiseaseStatus,
					variant1.getImputedDosagesAsMatrix2D(),
					finalCovariates,
					conditionalDosages);
			output = pruneAndTest(dosageMat, y, nrAlleles, taskObj);
		}

		if (output == null) {
			return null;
		}

		SNPFeature snp = new SNPFeature(Chromosome.parseChr(variant1.getChr()), variant1.getPos(), variant1.getPos());

		Double imputationqualityscore = variant1.getImputationQualityScore();
		snp.setName(variant1.getId());
		output.setSnp(snp);
		output.setN(x.rows());
		output.setMaf(maf1);
		output.setHWEP(hwep1);
		snp.setImputationQualityScore(imputationqualityscore);
		snp.setAlleles(variant1.getAlleles());
		snp.setMinorAllele(variant1.getMinorAllele());

		SNPFeature snp2 = new SNPFeature(Chromosome.parseChr(variant2.getChr()), variant2.getPos(), variant2.getPos());
		Double imputationqualityscore2 = variant2.getImputationQualityScore();
		snp2.setName(variant2.getId());
		output.setSnp2(snp2);
		output.setMaf2(maf2);
		output.setHWEP2(hwep2);
		snp2.setImputationQualityScore(imputationqualityscore2);
		snp2.setAlleles(variant2.getAlleles());
		snp2.setMinorAllele(variant2.getMinorAllele());

		output.setPairWise(true);


		// calculate the ld between the variants :)
//		DetermineLD ldcalc = new DetermineLD();
//		Pair<Double, Double> ld = ldcalc.getLD(variant1, variant2);
		output.setLDRSquared(0d);
		output.setLdDprime(0d);

		return output;
	}

	private AssociationResult pruneAndTest(DoubleMatrix2D x,
										   double[] y,
										   int nrAlleles,
										   LRTestTask testObj) throws IOException {
		Pair<DoubleMatrix2D, boolean[]> pruned = testObj.removeCollinearVariables(x);
		x = pruned.getLeft(); // x is now probably shorter than original X
		boolean[] notaliased = pruned.getRight(); // length of original X

		// check if the alleles we put in are aliased
		int firstAllele = 1; // intercept (first column) + allleles we conditioned on (original X indexing)
		int lastAllele = 1 + (nrAlleles - 2); // intercept + allleles we conditioned on + alleles for this variant (original X indexing)

		// index new position in X
		int nrRemaining = 0;
		int[] alleleIndex = new int[nrAlleles - 2];
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


		if (nrRemaining > 0) {
			// perform testNormal on full model
			// remove genotypes and run testNormal on reduced model

			LogisticRegressionOptimized reg = new LogisticRegressionOptimized();
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
			DoubleMatrix2D covarsOnly = testObj.removeGenotypes(x, colsToRemove);
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

			result.setDevianceNull(devnull);
			result.setDevianceGeno(devx);
			result.setDf(df);
			result.setBeta(betasmlelr);
			result.setSe(stderrsmlelr);
			result.setPval(p);
			return result;
		} else {
			return result;
		}
	}
}
