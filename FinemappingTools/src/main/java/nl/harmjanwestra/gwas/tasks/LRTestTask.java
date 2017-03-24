package nl.harmjanwestra.gwas.tasks;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import nl.harmjanwestra.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.math.LogisticRegressionOptimized;
import nl.harmjanwestra.utilities.math.LogisticRegressionResult;
import nl.harmjanwestra.utilities.vcf.SampleAnnotation;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.math.stats.ChiSquare;
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

	SampleAnnotation sampleAnnotation;
	private VCFVariant variant;
	private int iter;
	private ArrayList<Pair<VCFVariant, Triple<int[], boolean[], Integer>>> conditional;
	private int alleleOffsetGenotypes;
	private LRTestOptions options;
	private DenseDoubleAlgebra dda = new DenseDoubleAlgebra();
	private LogisticRegressionResult resultCovars;
	private int nrCovars;

	public LRTestTask(SampleAnnotation sampleAnnotation) {
		this.sampleAnnotation = sampleAnnotation;
	}

	public LRTestTask(VCFVariant variant,
					  int iter,
					  ArrayList<Pair<VCFVariant, Triple<int[], boolean[], Integer>>> conditional,
					  int alleleOffsetGenotypes,
					  SampleAnnotation sampleAnnotation,
					  LRTestOptions options) {
		this.variant = variant;
		this.iter = iter;
		this.conditional = conditional;
		this.alleleOffsetGenotypes = alleleOffsetGenotypes;
		this.options = options;
		this.sampleAnnotation = sampleAnnotation;
	}

	@Override
	public Triple<String, AssociationResult, VCFVariant> call() {
		// recode the genotypes to the same ordering as the covariate table

		// generate pseudocontrol genotypes
		Pair<DoubleMatrix2D, double[]> xandy = prepareMatrices(
				variant,
				conditional
		);

		int nrAlleles = variant.getAlleles().length;
		double[] y = xandy.getRight(); // get the phenotypes for all non-missing genotypes
		DoubleMatrix2D x = xandy.getLeft();

		if (x.rows() != y.length) {
			System.err.println("Unequal length for X and Y:  x: " + x.rows() + "x" + x.columns() + "\ty: " + y.length + "\tcallrate " + variant.getCallrate());
			System.exit(-1);
		}

		double mafControls = variant.getMAFControls();
		double hwep = variant.getHwepControls();

		AssociationResult result = pruneAndTest(x, y, 1, 1 + (nrAlleles - 1), variant, mafControls);

		if (result == null) {
			return new Triple<>(null, null, null);
		}

		result.getSnp().setHwep(hwep);
		if (variant.getAlleleFrequenciesCases() != null) {
			result.getSnp().setAFCases(variant.getAlleleFrequenciesCases()[0]);
			result.getSnp().setAFControls(variant.getAlleleFrequenciesControls()[0]);
		}


		//								SNP	Chr	Pos	ImputationQual	MAF	OverlapOK	MAFOk	ImpQualOK
		String logoutput = variant.getId()
				+ "\t" + variant.getChr()
				+ "\t" + variant.getPos()
				+ "\t" + Strings.concat(variant.getAlleles(), Strings.comma)
				+ "\t" + variant.getMinorAllele()
				+ "\t" + variant.getImputationQualityScore()
				+ "\t" + mafControls
				+ "\t" + hwep
				+ "\t" + true;

		return new Triple<>(logoutput, result, variant);
//		} // end maf > threshold
	}

	public Pair<DoubleMatrix2D, double[]> prepareMatrices(VCFVariant variant,
														  ArrayList<Pair<VCFVariant, Triple<int[], boolean[], Integer>>> conditional) {

		LRTestVariantQCTask lrq = new LRTestVariantQCTask();
		Triple<int[], boolean[], Integer> qcdata = lrq.determineMissingGenotypes(
				variant.getGenotypeAllelesAsMatrix2D(),
				sampleAnnotation.getCovariates().rows());

		int[] left = qcdata.getLeft();

		// prepare genotype matrix
		DoubleMatrix2D x = variant.getDosagesAsMatrix2D();
		HashSet<Integer> missingGenotypeIds = new HashSet<Integer>();
		for (int i : left) {
			missingGenotypeIds.add(i);
		}

		if (conditional != null && !conditional.isEmpty()) {
			for (int c = 0; c < conditional.size(); c++) {
				Pair<VCFVariant, Triple<int[], boolean[], Integer>> conditionalData = conditional.get(c);
				VCFVariant var2 = conditionalData.getLeft();
				int[] left2 = conditionalData.getRight().getLeft();

				DoubleMatrix2D dosages = var2.getDosagesAsMatrix2D();

				for (int i : left2) {
					missingGenotypeIds.add(i);
				}
				x = DoubleFactory2D.dense.appendColumns(x, dosages);
			}
		}

		// now append covariates
		DoubleMatrix2D finalCovariates = sampleAnnotation.getCovariates();
		DiseaseStatus[][] finalDiseaseStatus = sampleAnnotation.getSampleDiseaseStatus();

		x = DoubleFactory2D.dense.appendColumns(x, finalCovariates);

		// add intercept column
		DoubleMatrix2D vector = DoubleFactory2D.dense.make(x.rows(), 1);
		vector.assign(1);
		x = DoubleFactory2D.dense.appendColumns(vector, x);

		if (missingGenotypeIds.isEmpty()) {
			double[] y = new double[finalDiseaseStatus.length];
			for (int i = 0; i < finalDiseaseStatus.length; i++) {
				y[i] = finalDiseaseStatus[i][0].getNumber();
			}
			return new Pair<>(x, y);
		} else {
			// filter missing samples
			int[] rowIndexes = new int[x.rows() - missingGenotypeIds.size()];
			double[] y = new double[finalDiseaseStatus.length - missingGenotypeIds.size()];

			int q = 0;
			for (int i = 0; i < x.rows(); i++) {
				if (!missingGenotypeIds.contains(i)) {
					rowIndexes[q] = i;
					q++;
				}
			}

			x = dda.subMatrix(x, rowIndexes, 0, x.columns() - 1);
			if (x.rows() != y.length) {
				System.out.println("Something weird has happened");
			}

			int ctr = 0;
			for (int i = 0; i < finalDiseaseStatus.length; i++) {
				if (!missingGenotypeIds.contains(i)) {
					y[ctr] = finalDiseaseStatus[i][0].getNumber();
					ctr++;
				}
			}
			return new Pair<>(x, y);
		}
	}


	// remove variables with zero variance and perfect correlation
	public Pair<DoubleMatrix2D, boolean[]> removeCollinearVariables(DoubleMatrix2D mat) {

//		DoubleMatrix2D corrmat = new DenseDoubleMatrix2D(mat.columns(), mat.columns());
//		for (int i = 0; i < mat.columns(); i++) {
//			for (int j = i; j < mat.columns(); j++) {
//				double c = Math.abs(Correlation.correlate(mat.viewColumn(i).toArray(), mat.viewColumn(j).toArray()));
//				if (Double.isNaN(c)) {
//					corrmat.setQuick(i, j, 0);
//				} else {
//					corrmat.setQuick(i, j, c);
//				}
//			}
//		}

//		System.out.println(mat.columns() + " columns to start with...");
		// TODO: calculate variance inflation factor using OLS
		OLSMultipleLinearRegression olsMultipleLinearRegression = new OLSMultipleLinearRegression();
		// skip first column, because it is the intercept
		ArrayList<Integer> columns = new ArrayList<>();

		for (int i = 1; i < mat.columns(); i++) {
			columns.add(i);
		}
//		try {
		int iter = 0;
		int nrColinear = mat.columns();

		while (nrColinear > 0) {
//				TextFile vifs = new TextFile("/Data/tmp/sh2b3fix/vif-" + variant.getId() + "-" + iter + ".txt", TextFile.W);
			int[] colidx = new int[mat.columns() - 2];
			nrColinear = 0;

			ArrayList<Integer> noncolinear = new ArrayList<Integer>();
			ArrayList<Integer> colinear = new ArrayList<Integer>();
			for (int i : columns) {
				double[] y = mat.viewColumn(i).toArray();
				int ctr = 0;
				for (int j : columns) {
					if (j != i) {
						colidx[ctr] = j;
						ctr++;
					}
				}
				DoubleMatrix2D othercols = dda.subMatrix(mat, 0, mat.rows() - 1, colidx);
				olsMultipleLinearRegression.newSampleData(y, othercols.toArray());
				double rsq = olsMultipleLinearRegression.calculateAdjustedRSquared();
				double vif = 1 / (1 - rsq);

				if (Double.isInfinite(vif) || vif > 10) {
					nrColinear++;
					colinear.add(i);
				} else {
					noncolinear.add(i);
				}
//					vifs.writeln(i + "\t" + rsq + "\t" + vif);
//					vifs.flush();
			}

			// if there are colinear columns, remove one
			if (nrColinear > 0) {
				// add all other columns, except for last colinear one
				ArrayList<Integer> currcolumns = noncolinear;
				for (int q = 0; q < colinear.size() - 1; q++) {
					currcolumns.add(colinear.get(q));
				}
				columns = currcolumns;
			}

//				vifs.close();
			iter++;
		}


//			TextFile out = new TextFile("/Data/tmp/sh2b3fix/cor-" + variant.getId() + ".txt", TextFile.W);
//			String header = "-";
//			for (int j = 0; j < corrmat.columns(); j++) {
//				header += "\tvar" + j;
//			}
//			out.writeln(header);
//			for (int i = 0; i < corrmat.rows(); i++) {
//				String ln = "var" + i;
//				for (int j = 0; j < corrmat.columns(); j++) {
//					ln += "\t" + corrmat.getQuick(i, j);
//				}
//				out.writeln(ln);
//			}
//			out.close();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

//
//		boolean[] includeCol = new boolean[mat.columns()];
//		for (int c = 0; c < includeCol.length; c++) {
//			includeCol[c] = true;
//		}
//
//		for (int j = 1; j < mat.columns(); j++) {
//			double[] col1 = mat.viewColumn(j).toArray();
//			if (ArrayMath.variance(col1) == 0) {
//				includeCol[j] = false;
//			} else {
//				for (int j2 = j + 1; j2 < mat.columns(); j2++) {
//					if (includeCol[j2]) {
//						double c = corrmat.getQuick(j, j2);
//						double[] col2 = mat.viewColumn(j2).toArray();
//						if (ArrayMath.variance(col2) == 0) {
//							includeCol[j2] = false;
//						} else if (c == 1d) {
//							includeCol[j] = false;  // prefer to remove the earlier column
//						}
//					}
//				}
//			}
//		}
////
//		ArrayList<Integer> remaining = new ArrayList<>();
//		for (int c = 0; c < includeCol.length; c++) {
//			if (includeCol[c]) {
//				remaining.add(c);
//			}
//		}

		ArrayList<Integer> cols = new ArrayList<>();
		cols.add(0);
		cols.addAll(columns);
		columns = cols;
		boolean[] includeCol = new boolean[mat.columns()];
		includeCol[0] = true;
		for (int i : columns) {
			includeCol[i] = true;
		}

		int[] colindexes = Primitives.toPrimitiveArr(columns.toArray(new Integer[0]));
		// remove correlated variables

		DoubleMatrix2D matOut = dda.subMatrix(mat, 0, mat.rows() - 1, colindexes);
		return new Pair<>(matOut, includeCol);
//		int[] colindexes = Primitives.toPrimitiveArr(remaining.toArray(new Integer[0]));
//		// remove correlated variables
//		DoubleMatrix2D matOut = dda.subMatrix(mat, 0, mat.rows() - 1, colindexes);
//		return new Pair<>(matOut, includeCol);
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


		ArrayList<Integer> remainingAlleles = new ArrayList<>();
		int[] alleleIndex = new int[lastColumnToRemove - firstColumnToRemove];
		for (int i = 0; i < alleleIndex.length; i++) {
			alleleIndex[i] = -1;
		}

		ArrayList<Integer> remainingCovariates = new ArrayList<>();
		remainingCovariates.add(0);
		int nrRemaining = 0;
//		int allelectr = 1; // assume the first allele is now column 1
//		for (int i = firstAllele; i < lastAllele; i++) {
//			if (notaliased[i]) {
//				remainingAlleles.add(i);
//				alleleIndex[] =
//				allelectr++;
//				nrRemaining++;
//			}
//		}
		for (int i = lastAllele; i < x.columns(); i++) {
			if (notaliased[i]) {
				remainingCovariates.add(i);
			}
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
		snp.setMaf(maf);

		Double imputationqualityscore = variant.getImputationQualityScore();
		if (imputationqualityscore != null) {
			snp.setImputationQualityScore(imputationqualityscore);
		}
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

			if (resultCovars == null) {
				DoubleMatrix2D xprime = null;
				try {
					xprime = dda.subMatrix(x, 0, x.rows() - 1, Primitives.toPrimitiveArr(remainingCovariates.toArray(new Integer[0])));
//					System.out.println(variant.getId() + "\t" + x.columns() + "\t" + xprime.columns());
					resultCovars = reg.univariate(y, xprime);
				} catch (IndexOutOfBoundsException q) {
					System.out.println(variant.getId() + "\tcols:" + x.columns());
					for (int i = 0; i < remainingCovariates.size(); i++) {

						System.out.println(i + "\t" + remainingCovariates.get(i));

					}
					System.out.println();
					for (int z = 0; z < notaliased.length; z++) {
						System.out.println(z + "\t" + notaliased[z]);
					}
					System.exit(-1);
				}


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


			if (nrRemaining > 1 && options.testMultiAllelicVariantsIndependently) {
				// multi allelic variant...

				AssociationResult[] subresults = new AssociationResult[nrRemaining];
				System.out.println(nrRemaining + " results remaining");
				int nrAlleles = variant.getNrAlleles() - 1;
				int allelectr = 0;
				for (int a = 0; a < nrAlleles; a++) {
					// get index of allele
					int idx = alleleIndex[a];
					if (idx != -1) {
						System.out.println("Allele " + a);

						// make a new x-matrix using the original x-matrix
						ArrayList<Integer> colIndexArrAllele = new ArrayList<>(x.columns());
						newXIndex = 0;
						int suballeleindex = -1;
						for (int i = 0; i < notaliased.length; i++) {
							if (notaliased[i]) {
								if (i >= firstAllele && i < lastAllele) {
									if (newXIndex == idx) {
										colIndexArrAllele.add(newXIndex);
										suballeleindex = newXIndex;
									}
									nrRemaining++;
								} else {
									colIndexArrAllele.add(newXIndex);
								}
								newXIndex++;
							}
						}

						DoubleMatrix2D xallele = dda.subMatrix(x, 0, x.rows() - 1, Primitives.toPrimitiveArr(colIndexArrAllele.toArray(new Integer[0])));
						LogisticRegressionResult resultX = reg.univariate(y, xallele);
						AssociationResult subresult = new AssociationResult();

						SNPFeature subsnp = new SNPFeature();
						subsnp.setChromosome(snp.getChromosome());
						subsnp.setStart(snp.getStart());
						subsnp.setStop(snp.getStop());
						subsnp.setMinorAllele("NA");
						subsnp.setAlleles(new String[]{snp.getAlleles()[0], snp.getAlleles()[a + 1]});
						subsnp.setName(snp.getName() + "_" + snp.getAlleles()[a + 1]);
						subresult.setSnp(subsnp);
						subresult.setN(x.rows());
						subsnp.setMaf(maf);

						double devx = resultX.getDeviance();
						double devnull = resultCovars.getDeviance();
						double[] betasmlelr = new double[1];
						double[] stderrsmlelr = new double[1];
						double[] or = new double[1];
						double[] orhi = new double[1];
						double[] orlo = new double[1];

						int ctr = 0;

						double beta = -resultX.getBeta()[suballeleindex];
						double se = resultX.getStderrs()[suballeleindex];
						betasmlelr[ctr] = beta;
						stderrsmlelr[ctr] = se;

						double OR = Math.exp(beta);
						double orLow = Math.exp(beta - 1.96 * se);
						double orHigh = Math.exp(beta + 1.96 * se);
						or[ctr] = OR;
						orhi[ctr] = orHigh;
						orlo[ctr] = orLow;

						double deltaDeviance = devnull - devx;
						int df = xallele.columns() - nrCovars;
						double p = ChiSquare.getP(df, deltaDeviance);

						subresult.setDevianceNull(devnull);
						subresult.setDevianceGeno(devx);
						subresult.setDf(df);
						subresult.setDfalt(xallele.columns());
						subresult.setDfnull(nrCovars);

						subresult.setBeta(betasmlelr);
						subresult.setSe(stderrsmlelr);
						subresult.setPval(p);

						subresults[allelectr] = subresult;
						allelectr++;
					}
				}

				result.setSubresults(subresults);
			}


			LogisticRegressionResult resultX = reg.univariate(y, x);
			if (resultX == null) {
				// try once more with some more iterations
				LogisticRegressionOptimized reg2 = new LogisticRegressionOptimized(1000);
				resultX = reg2.univariate(y, x);
				System.out.println("Rerunning model...");
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
			double devnull = resultCovars.getDeviance();
			double[] betasmlelr = new double[nrRemaining];
			double[] stderrsmlelr = new double[nrRemaining];
			double[] or = new double[nrRemaining];
			double[] orhi = new double[nrRemaining];
			double[] orlo = new double[nrRemaining];

			int ctr = 0;

//			DoubleMatrix2D covariates = sampleAnnotation.getCovariates();
//			DiseaseStatus[][] disease = sampleAnnotation.getSampleDiseaseStatus();

//			if (covariates.rows() == x.rows()) {
//				try {
//					TextFile tfout = new TextFile("/Data/tmp/sh2b3fix/" + variant.getId() + ".txt", TextFile.W);
//
//					String header = "sample\tpheno\tgt";
//					for (int c = 0; c < covariates.columns(); c++) {
//						header += "\tcov" + c;
//					}
//
//					tfout.writeln(header);
//
//					for (int r = 0; r < covariates.rows(); r++) {
//						String ln = "sample" + r + "\t" + disease[r][0].getNumber() + "\t" + x.getQuick(r, 1);
//						for (int c = 0; c < covariates.columns(); c++) {
//							ln += "\t" + covariates.getQuick(r, c);
//						}
//						tfout.writeln(ln);
//					}
//					tfout.close();
//				} catch (IOException e) {
//					e.printStackTrace();
//				}
//			} else {
//			try {
//				TextFile tfout = new TextFile("/Data/tmp/sh2b3fix/" + variant.getId() + "-x.txt", TextFile.W);
//
//				String header = "sample\tpheno\tpheno2";
//				for (int c = 0; c < x.columns(); c++) {
//					header += "\tvar" + c;
//				}
//
//				tfout.writeln(header);
//
//				for (int r = 0; r < x.rows(); r++) {
//					String ln = "sample" + r + "\t" + disease[r][0].getNumber() + "\t" + y[r];
//					for (int c = 0; c < x.columns(); c++) {
//						ln += "\t" + x.getQuick(r, c);
//					}
//					tfout.writeln(ln);
//				}
//				tfout.close();
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//			}


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

