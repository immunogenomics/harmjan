package nl.harmjanwestra.utilities.vcf;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;

import java.util.concurrent.*;

/**
 * Created by hwestra on 7/1/16.
 * Blatantly copied from Lude Franke's WGAViewerModel.java
 */
public class GeneticSimilarity {


	// provide a list of variants (within dataset comparison)
	// assume variants are pruned.
	// assumes variants have same strandedness...
	// format: variants[nrdatasets][nrvariants]

	public Pair<DoubleMatrix2D, DoubleMatrix2D> calculate2(VCFVariant[][] variants, ExecutorService service) {

		int nrInds1 = variants[0][0].getNrSamples();
		int nrInds2 = variants[1][0].getNrSamples();
		int nrVariants = variants[0].length;

		DoubleMatrix2D geneticSimilarity = new DenseDoubleMatrix2D(nrInds1, nrInds2);
		DoubleMatrix2D geneticSimilaritySameGenotypes = new DenseDoubleMatrix2D(nrInds1, nrInds2);

		DoubleMatrix2D genotypes1 = new DenseDoubleMatrix2D(nrInds1, nrVariants);
		DoubleMatrix2D genotypes2 = new DenseDoubleMatrix2D(nrInds1, nrVariants);
		double[] alleleFreqs = new double[nrVariants];
		double[] callrates = new double[nrVariants];

		ProgressBar pb1 = new ProgressBar(nrVariants, "Calculating allelefrequencies");

		for (int snp = 0; snp < nrVariants; snp++) {
			VCFVariant variant1 = variants[0][snp];
			VCFVariant variant2 = variants[1][snp];
			DoubleMatrix2D genotypeAlleles1 = variant1.getGenotypeAllelesAsMatrix2D();
			double alleleFreq0 = 0;
			int call = 0;
			for (int i = 0; i < nrInds1; i++) {
				if (genotypeAlleles1.getQuick(i, 0) != -1) {
					double genotype0I = 0;
					if (0 == genotypeAlleles1.getQuick(i, 0)) genotype0I += .5;
					if (0 == genotypeAlleles1.getQuick(i, 1)) genotype0I += .5;
					alleleFreq0 += genotype0I;
					genotypes1.setQuick(i, snp, genotype0I);
//					callrates[i]++;
					call++;
				} else {
					genotypes1.setQuick(i, snp, -1);
				}
			}

			DoubleMatrix2D genotypeAlleles2 = variant2.getGenotypeAllelesAsMatrix2D();

			for (int i = 0; i < nrInds2; i++) {
				if (genotypeAlleles2.getQuick(i, 0) != -1) {
					double genotype0I = 0;
					if (0 == genotypeAlleles2.getQuick(i, 0)) genotype0I += .5;
					if (0 == genotypeAlleles2.getQuick(i, 1)) genotype0I += .5;
					alleleFreq0 += genotype0I;
					genotypes2.setQuick(i, snp, genotype0I);
//					callRate[i]++;
					call++;
				} else {
					genotypes2.setQuick(i, snp, -1);
				}
			}

			double snpCallRate = call / (nrInds1 + nrInds2);
			callrates[snp] = snpCallRate;
			alleleFreqs[snp] = alleleFreq0 / call;
			pb1.set(snp);
		}
		pb1.close();


		long nrSubmitted = 0;
		CompletionService<Triple<Integer, Integer, Pair<Double, Double>>> jobHandler = new ExecutorCompletionService<>(service);
		ProgressBar pb2 = new ProgressBar(nrInds1 * nrInds2, "Calculating distances");
		long returned = 0;
		for (int i = 0; i < nrInds1; i++) {
			for (int j = 0; j < nrInds2; j++) {
				DetermineGeneticSimilarityTask t = new DetermineGeneticSimilarityTask(
						i,
						j,
						genotypes1,
						genotypes2,
						nrVariants,
						alleleFreqs,
						callrates
				);
				jobHandler.submit(t);

				nrSubmitted++;
				if (nrSubmitted % 1000000 == 0) {
					System.out.println("Clearing buffer..... " + nrSubmitted + " / " + (nrInds1 * nrInds2));
					while (returned < nrSubmitted) {

						try {
							Triple<Integer, Integer, Pair<Double, Double>> future = jobHandler.take().get();

							int i1 = future.getLeft();
							int j1 = future.getMiddle();
							double similarity = future.getRight().getLeft();
							double simgenotypes = future.getRight().getRight();

							geneticSimilarity.setQuick(i1, j1, similarity);
							geneticSimilaritySameGenotypes.setQuick(i1, j1, simgenotypes);

							returned++;

						} catch (InterruptedException e) {
							e.printStackTrace();
						} catch (ExecutionException e) {
							e.printStackTrace();
						}

					}

					System.out.println("Done clearing buffer.. " + returned);
					pb2.set(returned);
					pb2.print();
				}
			}
		}
		System.out.println("Done submitting..");

		while (returned < nrSubmitted) {
			try {
				Triple<Integer, Integer, Pair<Double, Double>> future = jobHandler.take().get();

				int i1 = future.getLeft();
				int j1 = future.getMiddle();
				double similarity = future.getRight().getLeft();
				double simgenotypes = future.getRight().getRight();

				geneticSimilarity.setQuick(i1, j1, similarity);
				geneticSimilaritySameGenotypes.setQuick(i1, j1, simgenotypes);
				pb2.iterate();
				returned++;

			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}
		pb2.close();

		System.out.println("Done calculating genetic distance.");

		return new Pair<DoubleMatrix2D, DoubleMatrix2D>(geneticSimilarity, geneticSimilaritySameGenotypes);
	}

	public class DetermineGeneticSimilarityTask implements Callable<Triple<Integer, Integer, Pair<Double, Double>>> {

		double[] alleleFreqs;
		double[] callrates;

		int i;
		int j;
		int nrVariants;
		DoubleMatrix2D genotypes1;
		DoubleMatrix2D genotypes2;

		public DetermineGeneticSimilarityTask(int i, int j, DoubleMatrix2D genotypes1, DoubleMatrix2D genotypes2, int nrVariants, double[] alleleFreqs, double[] callrates) {
			this.i = i;
			this.j = j;
			this.genotypes1 = genotypes1;
			this.genotypes2 = genotypes2;
			this.nrVariants = nrVariants;
			this.alleleFreqs = alleleFreqs;
			this.callrates = callrates;
		}


		@Override
		public Triple<Integer, Integer, Pair<Double, Double>> call() throws Exception {

			int nrIdentical = 0;
			int nrGenotypes = 0;
			double sum = 0;
			for (int snp = 0; snp < nrVariants; snp++) {
				double af = alleleFreqs[snp];
				double cr = callrates[snp];

				if (!Double.isNaN(af) && cr >= 0.10) {
					double denominator = af * (1.0d - af);
					if (!Double.isNaN(denominator) && denominator > 0) {
						double genotype0I = genotypes1.getQuick(i, snp);
						double genotype0J = genotypes2.getQuick(j, snp);
						if (genotype0I != -1 && genotype0J != -1) {
							sum += ((genotype0I - af) * (genotype0J - af) / denominator);
							if (genotype0I == genotype0J) {
								nrIdentical++;
							}
							nrGenotypes++;
						}
					}
				}
			}
			return new Triple<Integer, Integer, Pair<Double, Double>>(i, j, new Pair<Double, Double>(sum / nrGenotypes, (double) nrIdentical / nrGenotypes));
		}
	}


	// provide a list of variants (within dataset comparison)
	// assume variants are pruned.
	public Pair<double[][], double[][]> calculate(VCFVariant[] variants) {

		int nrInds = variants[0].getNrSamples();
		int nrVariants = variants.length;

		double[] callRate = new double[nrInds];
		double[][] geneticSimilarity = new double[nrInds][nrInds];
		double[][] geneticSimilaritySameGenotypes = new double[nrInds][nrInds];
		int[][] geneticSimilarityCalled = new int[nrInds][nrInds];


		for (int snpID = 0; snpID < nrVariants; snpID++) {

			VCFVariant variant = variants[snpID];
			double[][] genotypeAlleles = variant.getGenotypeAlleles();
			double alleleFreq0 = 0;
			double call = 0;
			for (int i = 0; i < nrInds; i++) {
				if (genotypeAlleles[i][0] != -1) {
					double genotype0I = 0;
					if (0 == genotypeAlleles[i][0]) genotype0I += .5;
					if (0 == genotypeAlleles[i][1]) genotype0I += .5;
					alleleFreq0 += genotype0I;
					callRate[i]++;
					call++;
				}
			}

			alleleFreq0 /= call;
			double snpCallRate = call / nrInds;

			if (!Double.isNaN(alleleFreq0) && snpCallRate >= 0.10) {
				double denominator = alleleFreq0 * (1.0d - alleleFreq0);
				if (!Double.isNaN(denominator) && denominator > 0) {
					for (int i = 0; i < nrInds; i++) {
						if (genotypeAlleles[i][0] != -1) {
							double genotype0I = 0;
							if (0 == genotypeAlleles[i][0]) genotype0I += .5;
							if (0 == genotypeAlleles[i][1]) genotype0I += .5;
							for (int j = i + 1; j < nrInds; j++) {
								if (genotypeAlleles[j][0] != -1) {
									double genotype0J = 0;
									if (0 == genotypeAlleles[j][0]) genotype0J += .5;
									if (0 == genotypeAlleles[j][1]) genotype0J += .5;
									geneticSimilarity[i][j] += (genotype0I - alleleFreq0) * (genotype0J - alleleFreq0) / denominator;
									if (genotype0I == genotype0J) geneticSimilaritySameGenotypes[i][j]++;
									geneticSimilarityCalled[i][j]++;
								}
							}
						}
					}
				}
			}
		}

		for (int i = 0; i < nrInds; i++) {
			for (int j = i + 1; j < nrInds; j++) {
				geneticSimilarity[i][j] /= geneticSimilarityCalled[i][j];
				geneticSimilaritySameGenotypes[i][j] /= geneticSimilarityCalled[i][j];
			}
		}

		return new Pair<double[][], double[][]>(geneticSimilarity, geneticSimilaritySameGenotypes);

	}


}
