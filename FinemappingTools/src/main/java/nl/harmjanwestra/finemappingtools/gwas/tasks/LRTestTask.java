package nl.harmjanwestra.finemappingtools.gwas.tasks;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import nl.harmjanwestra.finemappingtools.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Triple;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.ChiSquare;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.Descriptives;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import nl.harmjanwestra.utilities.legacy.genetica.util.Primitives;
import nl.harmjanwestra.utilities.math.LogisticRegressionOptimized;
import nl.harmjanwestra.utilities.math.LogisticRegressionResult;
import nl.harmjanwestra.utilities.vcf.SampleAnnotation;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
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
	
	public LRTestTask(SampleAnnotation sampleAnnotation, LRTestOptions options) {
		this.sampleAnnotation = sampleAnnotation;
		this.options = options;
	}
	
	@Override
	public Triple<String, AssociationResult, VCFVariant> call() {
		// recode the genotypes to the same ordering as the covariate table
		
		// generate pseudocontrol genotypes
		Pair<DoubleMatrix2D, double[][]> xandy = prepareMatrices(
				variant,
				conditional
		);
		
		int nrAlleles = variant.getAlleles().length;
		double[][] y = xandy.getRight(); // get the phenotypes for all non-missing genotypes
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
	
	public Pair<DoubleMatrix2D, double[][]> prepareMatrices(VCFVariant variant,
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
			double[][] y = new double[finalDiseaseStatus.length][finalDiseaseStatus[0].length];
			for (int i = 0; i < finalDiseaseStatus.length; i++) {
				for (int j = 0; j < finalDiseaseStatus[i].length; j++) {
					y[i][j] = finalDiseaseStatus[i][j].getNumber();
				}
			}
			return new Pair<>(x, y);
		} else {
			// filter missing samples
			int[] rowIndexes = new int[x.rows() - missingGenotypeIds.size()];
			double[][] y = new double[finalDiseaseStatus.length - missingGenotypeIds.size()][finalDiseaseStatus[0].length];
			
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
					for (int j = 0; j < finalDiseaseStatus[i].length; j++) {
						y[ctr][j] = finalDiseaseStatus[i][j].getNumber();
					}
					
					ctr++;
				}
			}
			return new Pair<>(x, y);
		}
	}
	
	
	// remove variables with zero variance and perfect correlation
	public Pair<DoubleMatrix2D, boolean[]> removeCollinearVariables(DoubleMatrix2D mat) {
		
		OLSMultipleLinearRegression olsMultipleLinearRegression = new OLSMultipleLinearRegression();
		// skip first column, because it is the intercept
		ArrayList<Integer> columns = new ArrayList<>();
		
		// check mean and variance
		for (int i = 1; i < mat.columns(); i++) {
			double[] col = mat.viewColumn(i).toArray();
			if (Descriptives.variance(col) > 0) {
				columns.add(i);
			} else if (options.debug) {
				System.out.println("Column has zero variance: " + i + "\t" + Descriptives.variance(col));
			}
		}
//		try {
		int iter = 0;
		int nrColinear = mat.columns();
		
		
		try {
			if (options.debug) {
				System.out.println("Finding colinear covariates from: " + columns.size() + " cols");
			}
			
			while (nrColinear > 0) {
				TextFile vifs = null;
				if (options.debug && variant != null) {
					vifs = new TextFile(options.getOutputdir() + "vif-" + variant.getId() + "-" + iter + ".txt", TextFile.W);
				}
				nrColinear = 0;
				
				ArrayList<Integer> noncolinear = new ArrayList<Integer>();
				ArrayList<Integer> colinear = new ArrayList<Integer>();
				for (int i : columns) {
					double[] zscores = null;
					if (options.debug) {
						zscores = new double[columns.size()];
					}
					
					int[] colidx = new int[columns.size() - 1];
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
					double rsq = 1;
					try {
						rsq = olsMultipleLinearRegression.calculateAdjustedRSquared();
						if (options.debug) {
							double[] betas = olsMultipleLinearRegression.estimateRegressionParameters();
							double[] ses = olsMultipleLinearRegression.estimateRegressionParametersStandardErrors();
							for (int q = 0; q < betas.length; q++) {
								zscores[q] = betas[q] / ses[q];
							}
						}
						
					} catch (SingularMatrixException e) {
						
						if (options.debug) {
							System.out.println(e.getMessage());
							System.out.println();
							System.out.println("Error testing variant: " + variant.getId());
							System.out.println("Singular matrix detected when testing for multi-collinearity.");
							System.out.println("Iter: " + iter + "\tcol:" + i);
							System.out.println(columns.size() + " columns left out of " + mat.columns());
							
							for (int c = 0; c < mat.columns(); c++) {
								double[] col = mat.viewColumn(c).toArray();
								System.out.println(c + "\t" + Descriptives.mean(col) + "\t" + Descriptives.variance(col));
							}
							
							System.exit(-1);
						} else {
							return null;
						}
					}
					
					if (rsq >= options.collinearitythreshold) {
						nrColinear++;
						colinear.add(i);
					} else {
						noncolinear.add(i);
					}
//
					if (options.debug && variant != null) {
						String ln = i + "\t" + rsq + "\t" + Descriptives.variance(y);

//						for (int q = 0; q < zscores.length; q++) {
//							ln += "\t" + zscores[q];
//						}
						
						vifs.writeln(ln);
						vifs.flush();
					}
				}
				
				// if there are colinear columns, remove one
				if (nrColinear > 0) {
					// add all other columns, except for last colinear one
					ArrayList<Integer> currcolumns = noncolinear;
					for (int q = 0; q < colinear.size() - 1; q++) {
						currcolumns.add(colinear.get(q));
					}
					Collections.sort(currcolumns);
					if (options.debug) {
						System.out.println("Removing column: " + colinear.get(colinear.size() - 1));
					}
					columns = currcolumns;
				}
				if (options.debug && variant != null) {
					vifs.close();
				}
				iter++;
			}
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
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
	}
	
	private AssociationResult pruneAndTest(DoubleMatrix2D x,
										   double[][] y,
										   int firstGenotypeColumn,
										   int lastGenotypeColumn,
										   VCFVariant variant,
										   double maf) {
		
		
		if (options.debug) {
			System.out.println(x.columns() + " before pruning");
		}
		
		Pair<DoubleMatrix2D, boolean[]> pruned = removeCollinearVariables(x);
		if (pruned == null) {
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
		
		if (conditional != null && !conditional.isEmpty()) {
			// count the number of conditional columns that is not aliased
			
			int nrConditionalAlleleCols = 0;
			for (int q = 0; q < conditional.size(); q++) {
				VCFVariant v = conditional.get(q).getLeft();
				nrConditionalAlleleCols += v.getAlleles().length - 1;
			}
			
			int lastConditionalCol = lastGenotypeColumn + nrConditionalAlleleCols; // 2 normal allles:1,2  3 conditional alleles: 2,3,4
			int nrConditionalAlleleColsNotAliased = 0;
			for (int q = lastGenotypeColumn; q < lastConditionalCol; q++) {
				if (notaliased[q]) {
					nrConditionalAlleleColsNotAliased++;
				}
			}
			conditionalOk = (nrConditionalAlleleCols == nrConditionalAlleleColsNotAliased);
		}
		
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
		
		if (options.debug) {
			System.out.println(remainingGenotypeColumns.size() + " remaining genotype cols: ");
			for (int col : remainingGenotypeColumns) {
				System.out.println(col);
			}
			for (int i = 0; i < alleleIndex.length; i++) {
				System.out.println("Allele: " + i + "\tnew idx: " + alleleIndex[i]);
			}
			System.out.println(remainingCovariateColumns.size() + " remaining covariate cols");
			for (int col : remainingCovariateColumns) {
				System.out.println(col);
			}
		}
		
		
		AssociationResult result = new AssociationResult();
		SNPFeature snp = new SNPFeature(Chromosome.parseChr(variant.getChr()), variant.getPos(), variant.getPos());
		snp.setName(variant.getId());
		result.setSnp(snp);
		result.setN(x.rows());
		snp.setMaf(maf);
		snp.setCrCases(variant.getCallrateCases());
		snp.setCrControls(variant.getCallrateControls());
		snp.setMissingnessP(variant.getDiffMissingnessP());
		
		Double imputationqualityscore = variant.getImputationQualityScore();
		if (imputationqualityscore != null) {
			snp.setImputationQualityScore(imputationqualityscore);
		}
		snp.setAlleles(variant.getAlleles());
		snp.setMinorAllele(variant.getMinorAllele());
		
		
		int nrDiseases = y[0].length;
		int nrAlleles = snp.getAlleles().length - 1;
		
		if (remainingGenotypeColumns.isEmpty() || !conditionalOk) {
			result.setDevianceNull(0);
			result.setDevianceGeno(0);
			result.setDf(0);
			
			double[][] nullMat = new double[nrDiseases][nrAlleles];
			result.setBeta(nullMat);
			result.setSe(nullMat);
			result.setPval(1);
			return result;
		} else {
			LogisticRegressionOptimized reg = new LogisticRegressionOptimized();
			
			if (resultCovars == null) {
				DoubleMatrix2D xprime = null;
				try {
					xprime = dda.subMatrix(x, 0, x.rows() - 1, Primitives.toPrimitiveArr(remainingCovariateColumns.toArray(new Integer[0])));
//					System.out.println(variant.getId() + "\t" + x.columns() + "\t" + xprime.columns());
					if (options.debug) {
						System.out.println("xprime size: " + xprime.rows() + "x" + xprime.columns());
					}
					
					resultCovars = reg.multinomial(y, xprime);
				} catch (IndexOutOfBoundsException q) {
					System.out.println(variant.getId() + "\tcols:" + x.columns());
					for (int i = 0; i < remainingCovariateColumns.size(); i++) {
						
						System.out.println(i + "\t" + remainingCovariateColumns.get(i));
						
					}
					System.out.println();
					for (int z = 0; z < notaliased.length; z++) {
						System.out.println(z + "\t" + notaliased[z]);
					}
					System.exit(-1);
				}
				
				
				if (options.debug) {
					// debug shizzle
					
					reg.debug = options.debug;
					reg.printIters = true;
					reg.flipCoding = false;
					reg.output = options.getOutputdir() + "debugdata-" + snp.getName() + ".txt";
					resultCovars = reg.multinomial(y, xprime);
					
					System.out.println("---");
					System.out.println("Debugging full model..");
					LogisticRegressionResult r2 = reg.multinomial(y, x);
					
//					double deltadeviance = r2.getDeviance() - resultCovars.getDeviance();
//					System.out.println("Delta deviance: " + deltadeviance);

//					System.out.println("Flipping signs");
//					reg.flipCoding = true;
//					resultCovars = reg.multinomial(y, xprime);
					
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
			
			
			if (remainingGenotypeColumns.size() > 1 && options.testMultiAllelicVariantsIndependently) {
				// multi allelic variant...
				if (options.debug) {
					System.out.println("Testing alleles independently");
				}
				
				int nrRemaining = remainingGenotypeColumns.size();
				nrAlleles = lastGenotypeColumn - firstGenotypeColumn;
				AssociationResult[] subresults = new AssociationResult[nrRemaining];
				System.out.println(nrRemaining + " results remaining");
				
				allelectr = 0;
				for (int a = 0; a < nrAlleles; a++) {
					// get index of allele
					int idx = alleleIndex[a];
					if (idx != -1) {
						System.out.println("Allele " + a);
						
						// make a new x-matrix using the original x-matrix
						ArrayList<Integer> colIndexArrAllele = new ArrayList<>(x.columns());
						int newXIndex = 0;
						int suballeleindex = -1;
						for (int i = 0; i < notaliased.length; i++) {
							if (notaliased[i]) {
								if (i >= firstGenotypeColumn && i < lastGenotypeColumn) {
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
						LogisticRegressionResult resultX = reg.multinomial(y, xallele);
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
						double[][] betasmlelr = new double[nrDiseases][1];
						double[][] stderrsmlelr = new double[nrDiseases][1];
//						double[] or = new double[1];
//						double[] orhi = new double[1];
//						double[] orlo = new double[1];
						
						ctr = 0;
						for (int d = 0; d < nrDiseases; d++) {
							double beta = resultX.getBeta()[d][suballeleindex];
							double se = resultX.getStderrs()[d][suballeleindex];
							betasmlelr[d][ctr] = beta;
							stderrsmlelr[d][ctr] = se;
						}


//						double OR = Math.exp(beta);
//						double orLow = Math.exp(beta - 1.96 * se);
//						double orHigh = Math.exp(beta + 1.96 * se);
//						or[ctr] = OR;
//						orhi[ctr] = orHigh;
//						orlo[ctr] = orLow;
						
						double deltaDeviance = devnull - devx;

//						int df = xallele.columns() - nrCovars;
						
						int totalDf = xallele.columns() * y[0].length;
						int deltaDf = totalDf - (nrCovars * y[0].length);
						
						double p = ChiSquare.getP(deltaDf, deltaDeviance);
						subresult.setDevianceNull(devnull);
						subresult.setDevianceGeno(devx);
						subresult.setDf(deltaDf);
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
			
			
			if (options.debug) {
				System.out.println("x size: " + x.rows() + "x" + x.columns());
			}
			
			LogisticRegressionResult resultX = reg.multinomial(y, x);
			if (resultX == null) {
				// try once more with some more iterations
				LogisticRegressionOptimized reg2 = new LogisticRegressionOptimized(1000);
				resultX = reg2.multinomial(y, x);
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
			
			
			nrAlleles = lastGenotypeColumn - firstGenotypeColumn;
			double devx = resultX.getDeviance();
			double devnull = resultCovars.getDeviance();
			double[][] betasmlelr = new double[nrDiseases][nrAlleles];
			double[][] stderrsmlelr = new double[nrDiseases][nrAlleles];
//			double[][] or = new double[nrAlleles];
//			double[] orhi = new double[nrAlleles];
//			double[] orlo = new double[nrAlleles];
//
			ctr = 0;
			
			for (int i = 0; i < alleleIndex.length; i++) {
				int idx = alleleIndex[i];
				if (idx != -1) {
					for (int d = 0; d < nrDiseases; d++) {
						double beta = resultX.getBeta()[d][idx];
						double se = resultX.getStderrs()[d][idx];
						betasmlelr[d][i] = beta;
						stderrsmlelr[d][i] = se;
					}

//					double OR = Math.exp(beta);
//					double orLow = Math.exp(beta - 1.96 * se);
//					double orHigh = Math.exp(beta + 1.96 * se);
//					or[i] = OR;
//					orhi[i] = orHigh;
//					orlo[i] = orLow;
					ctr++;
				}
			}
			
			double deltaDeviance = devnull - devx;
			
			int totalDf = x.columns() * y[0].length;
			int deltaDf = totalDf - (nrCovars * y[0].length);
			if (options.debug) {
				System.out.println("df: " + deltaDf);
			}
			double p = ChiSquare.getP(deltaDf, deltaDeviance);
			
			result.setDevianceNull(devnull);
			result.setDevianceGeno(devx);
			result.setDf(deltaDf);
			result.setDfalt(totalDf);
			result.setDfnull(nrCovars * y[0].length);
			
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

