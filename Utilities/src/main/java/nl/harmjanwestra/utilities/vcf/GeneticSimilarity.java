package nl.harmjanwestra.utilities.vcf;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
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

	public Triple<DoubleMatrix2D, DoubleMatrix2D, DoubleMatrix2D> calculateBetweenDatasets(VCFVariant[][] variants, ExecutorService service) {

		int nrInds1 = variants[0][0].getNrSamples();
		int nrInds2 = variants[1][0].getNrSamples();
		int nrVariants = variants[0].length;

		DoubleMatrix2D geneticSimilarity = new DenseDoubleMatrix2D(nrInds1, nrInds2);
		DoubleMatrix2D geneticSimilaritySameGenotypes = new DenseDoubleMatrix2D(nrInds1, nrInds2);
		DoubleMatrix2D geneticSimilarityCorrelation = new DenseDoubleMatrix2D(nrInds1, nrInds2);

		Pair<Pair<double[], double[]>, Pair<DoubleMatrix2D, DoubleMatrix2D>> data = calculateAlleleFreqAndCallRate(null, variants);
		DoubleMatrix2D genotypes1 = data.getRight().getLeft();
		DoubleMatrix2D genotypes2 = data.getRight().getRight();
		double[] callrates = data.getLeft().getLeft();
		double[] alleleFreqs = data.getLeft().getRight();


		long nrSubmitted = 0;
		CompletionService<Triple<Integer, Integer, Triple<Double, Double, Double>>> jobHandler = new ExecutorCompletionService<>(service);
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
							Triple<Integer, Integer, Triple<Double, Double, Double>> future = jobHandler.take().get();

							int i1 = future.getLeft();
							int j1 = future.getMiddle();
							double similarity = future.getRight().getLeft();
							double simgenotypes = future.getRight().getMiddle();
							double correlation = future.getRight().getRight();
							geneticSimilarity.setQuick(i1, j1, similarity);
							geneticSimilaritySameGenotypes.setQuick(i1, j1, simgenotypes);
							geneticSimilarityCorrelation.setQuick(i1, j1, correlation);


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
				Triple<Integer, Integer, Triple<Double, Double, Double>> future = jobHandler.take().get();

				int i1 = future.getLeft();
				int j1 = future.getMiddle();
				double similarity = future.getRight().getLeft();
				double simgenotypes = future.getRight().getMiddle();
				double correlation = future.getRight().getRight();
				geneticSimilarity.setQuick(i1, j1, similarity);
				geneticSimilaritySameGenotypes.setQuick(i1, j1, simgenotypes);
				geneticSimilarityCorrelation.setQuick(i1, j1, correlation);
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

		return new Triple<DoubleMatrix2D, DoubleMatrix2D, DoubleMatrix2D>(geneticSimilarity, geneticSimilaritySameGenotypes, geneticSimilarityCorrelation);
	}

	// provide a list of variants (within dataset comparison)
	// assume variants are pruned.
	public Triple<DoubleMatrix2D, DoubleMatrix2D, DoubleMatrix2D> calculateWithinDataset(VCFVariant[] variants, ExecutorService service) {

		int nrInds1 = variants[0].getNrSamples();

		int nrVariants = variants.length;

		DoubleMatrix2D geneticSimilarity = new DenseDoubleMatrix2D(nrInds1, nrInds1);
		DoubleMatrix2D geneticSimilaritySameGenotypes = new DenseDoubleMatrix2D(nrInds1, nrInds1);
		DoubleMatrix2D geneticSimilarityCorrelation = new DenseDoubleMatrix2D(nrInds1, nrInds1);


		Pair<Pair<double[], double[]>, Pair<DoubleMatrix2D, DoubleMatrix2D>> data = calculateAlleleFreqAndCallRate(variants, null);
		DoubleMatrix2D genotypes1 = data.getRight().getLeft();
		double[] callrates = data.getLeft().getLeft();
		double[] alleleFreqs = data.getLeft().getRight();


		long nrSubmitted = 0;
		CompletionService<Triple<Integer, Integer, Triple<Double, Double, Double>>> jobHandler = new ExecutorCompletionService<>(service);
		ProgressBar pb2 = new ProgressBar(((nrInds1 * nrInds1) - nrInds1) / 2, "Calculating distances");
		long returned = 0;
		for (int i = 0; i < nrInds1; i++) {
			for (int j = i; j < nrInds1; j++) {
				DetermineGeneticSimilarityTask t = new DetermineGeneticSimilarityTask(
						i,
						j,
						genotypes1,
						null,
						nrVariants,
						alleleFreqs,
						callrates
				);
				jobHandler.submit(t);

				nrSubmitted++;
				if (nrSubmitted % 1000000 == 0) {

					while (returned < nrSubmitted) {

						try {
							Triple<Integer, Integer, Triple<Double, Double, Double>> future = jobHandler.take().get();

							int i1 = future.getLeft();
							int j1 = future.getMiddle();
							double similarity = future.getRight().getLeft();
							double simgenotypes = future.getRight().getMiddle();
							double correlation = future.getRight().getRight();

							geneticSimilarity.setQuick(i1, j1, similarity);
							geneticSimilarity.setQuick(j1, i1, similarity);
							geneticSimilaritySameGenotypes.setQuick(i1, j1, simgenotypes);
							geneticSimilaritySameGenotypes.setQuick(j1, i1, simgenotypes);
							geneticSimilarityCorrelation.setQuick(i1, j1, correlation);
							geneticSimilarityCorrelation.setQuick(j1, i1, correlation);

							returned++;
							pb2.iterate();
						} catch (InterruptedException e) {
							e.printStackTrace();
						} catch (ExecutionException e) {
							e.printStackTrace();
						}

					}
				}
			}
		}
		System.out.println("Done submitting..");

		while (returned < nrSubmitted) {
			try {
				Triple<Integer, Integer, Triple<Double, Double, Double>> future = jobHandler.take().get();

				int i1 = future.getLeft();
				int j1 = future.getMiddle();
				double similarity = future.getRight().getLeft();
				double simgenotypes = future.getRight().getMiddle();
				double correlation = future.getRight().getRight();

				geneticSimilarity.setQuick(i1, j1, similarity);
				geneticSimilarity.setQuick(j1, i1, similarity);
				geneticSimilaritySameGenotypes.setQuick(i1, j1, simgenotypes);
				geneticSimilaritySameGenotypes.setQuick(j1, i1, simgenotypes);
				geneticSimilarityCorrelation.setQuick(i1, j1, correlation);
				geneticSimilarityCorrelation.setQuick(j1, i1, correlation);
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

		return new Triple<DoubleMatrix2D, DoubleMatrix2D, DoubleMatrix2D>(geneticSimilarity, geneticSimilaritySameGenotypes, geneticSimilarityCorrelation);
	}

	public Pair<Pair<double[], double[]>, Pair<DoubleMatrix2D, DoubleMatrix2D>> calculateAlleleFreqAndCallRate(VCFVariant[] variants, VCFVariant[][] variantMatrix) {

		DoubleMatrix2D genotypes1;
		DoubleMatrix2D genotypes2;

		int nrVariants = 0;
		int nrInds1 = 0;
		int nrInds2 = 0;
		if (variants == null) {
			nrVariants = variantMatrix[0].length;
			nrInds1 = variantMatrix[0][0].getNrSamples();
			nrInds2 = variantMatrix[1][0].getNrSamples();
			genotypes1 = new DenseDoubleMatrix2D(nrInds1, nrVariants);
			genotypes2 = new DenseDoubleMatrix2D(nrInds2, nrVariants);
		} else {
			nrVariants = variants.length;
			nrInds1 = variants[0].getNrSamples();
			genotypes1 = new DenseDoubleMatrix2D(nrInds1, nrVariants);
			genotypes2 = null;
		}

		ProgressBar pb1 = new ProgressBar(nrVariants, "Calculating allelefrequencies");


		double[] alleleFreqs = new double[nrVariants];
		double[] callrates = new double[nrVariants];
		for (int snp = 0; snp < nrVariants; snp++) {

			VCFVariant variant1;
			VCFVariant variant2;

			if (variants == null) {
				variant1 = variantMatrix[0][snp];
				variant2 = variantMatrix[1][snp];
			} else {
				variant1 = variants[snp];
				variant2 = null;
			}

			DoubleMatrix2D genotypeAlleles1 = variant1.getGenotypeAllelesAsMatrix2D();
			double alleleFreq0 = 0;
			int call = 0;
			for (int i = 0; i < nrInds1; i++) {
				if (genotypeAlleles1.getQuick(i, 0) != -1) {
					double gt = genotypeAlleles1.getQuick(i, 0) + genotypeAlleles1.getQuick(i, 1);
					alleleFreq0 += gt;
					genotypes1.setQuick(i, snp, gt);
					call++;
				} else {
					genotypes1.setQuick(i, snp, -1);
				}
			}

			if (variant2 != null) {
				DoubleMatrix2D genotypeAlleles2 = variant2.getGenotypeAllelesAsMatrix2D();

				for (int i = 0; i < nrInds2; i++) {
					if (genotypeAlleles2.getQuick(i, 0) != -1) {
						double gt = genotypeAlleles2.getQuick(i, 0) + genotypeAlleles2.getQuick(i, 1);
						alleleFreq0 += gt;
						genotypes2.setQuick(i, snp, gt);
						call++;
					} else {
						genotypes2.setQuick(i, snp, -1);
					}
				}
			}

			double snpCallRate = call / (nrInds1 + nrInds2);
			callrates[snp] = snpCallRate;
			alleleFreqs[snp] = alleleFreq0 / (call * 2);
			pb1.set(snp);
		}
		pb1.close();

		return new Pair<Pair<double[], double[]>, Pair<DoubleMatrix2D, DoubleMatrix2D>>(
				new Pair<double[], double[]>(callrates, alleleFreqs),
				new Pair<DoubleMatrix2D, DoubleMatrix2D>(genotypes1, genotypes2)
		);
	}


	public class DetermineGeneticSimilarityTask implements Callable<Triple<Integer, Integer, Triple<Double, Double, Double>>> {

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
		public Triple<Integer, Integer, Triple<Double, Double, Double>> call() throws Exception {

			int nrIdentical = 0;
			int nrGenotypes = 0;
			double sum = 0;

			double[] x = new double[nrVariants];
			double[] y = new double[nrVariants];
			int nrMissing = 0;
			for (int snp = 0; snp < nrVariants; snp++) {
				double af = alleleFreqs[snp];
				double cr = callrates[snp];

				double genotype0I = genotypes1.getQuick(i, snp);
				double genotype0J = -1;
				if (genotypes2 == null) {
					genotype0J = genotypes1.getQuick(j, snp);
				} else {
					genotype0J = genotypes2.getQuick(j, snp);
				}

				if (genotype0I != -1 && genotype0J != -1) {
					x[snp] = genotype0I;
					y[snp] = genotype0J;
				} else {
					x[snp] = -1;
					y[snp] = -1;
					nrMissing++;
				}

				if (!Double.isNaN(af) && cr >= 0.10) {
					double denominator = 2 * (af * (1.0d - af));
					if (!Double.isNaN(denominator) && denominator > 0) {


						if (genotype0I != -1 && genotype0J != -1) {
							double g01 = genotype0I - (2 * af);
							double g02 = genotype0J - (2 * af);

							double fhat = (g01 * g02) / denominator;
							sum += fhat;
							if (genotype0I == genotype0J) {
								nrIdentical++;
							}
							nrGenotypes++;
						}
					}
				}
			}

			// include correlation here as well?
			double[] xn = new double[nrVariants - nrMissing];
			double[] yn = new double[nrVariants - nrMissing];

			int ctr = 0;
			for (int s = 0; s < nrVariants; s++) {
				if (x[s] != -1) {
					xn[ctr] = x[s];
					yn[ctr] = y[s];
					ctr++;
				}
			}

			PearsonsCorrelation c = new PearsonsCorrelation();
			double corr = c.correlation(x, y);
			return new Triple<Integer, Integer, Triple<Double, Double, Double>>(i, j, new Triple<>(sum / nrGenotypes, (double) nrIdentical / nrGenotypes, corr));
		}
	}


}
