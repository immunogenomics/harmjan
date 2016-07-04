package nl.harmjanwestra.utilities.vcf;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;

/**
 * Created by hwestra on 7/1/16.
 * Blatantly copied from Lude Franke's WGAViewerModel.java
 */
public class GeneticSimilarity {


	// provide a list of variants (within dataset comparison)
	// assume variants are pruned.
	// assumes variants have same strandedness...
	// format: variants[nrdatasets][nrvariants]
	public Pair<double[][], double[][]> calculate(VCFVariant[][] variants) {

		int nrInds1 = variants[0][0].getNrSamples();
		int nrInds2 = variants[1][0].getNrSamples();
		int nrVariants = variants[0].length;

		double[] callRate = new double[nrInds1];
		double[][] geneticSimilarity = new double[nrInds1][nrInds2];
		double[][] geneticSimilaritySameGenotypes = new double[nrInds1][nrInds2];
		int[][] geneticSimilarityCalled = new int[nrInds1][nrInds2];

		ProgressBar pb = new ProgressBar(nrVariants, "Calculating genetic similarity");
		for (int snpID = 0; snpID < nrVariants; snpID++) {

			VCFVariant variant1 = variants[0][snpID];
			VCFVariant variant2 = variants[1][snpID];


			double alleleFreq0 = 0;
			double call = 0;
			double[][] genotypeAlleles1 = variant1.getGenotypeAlleles();
			for (int i = 0; i < nrInds1; i++) {
				if (genotypeAlleles1[i][0] != -1) {
					double genotype0I = 0;
					if (0 == genotypeAlleles1[i][0]) genotype0I += .5;
					if (0 == genotypeAlleles1[i][1]) genotype0I += .5;
					alleleFreq0 += genotype0I;
					callRate[i]++;
					call++;
				}
			}

			double[][] genotypeAlleles2 = variant2.getGenotypeAlleles();
			for (int i = 0; i < nrInds2; i++) {
				if (genotypeAlleles2[i][0] != -1) {
					double genotype0I = 0;
					if (0 == genotypeAlleles2[i][0]) genotype0I += .5;
					if (0 == genotypeAlleles2[i][1]) genotype0I += .5;
					alleleFreq0 += genotype0I;
					callRate[i]++;
					call++;
				}
			}

			alleleFreq0 /= call;
			double snpCallRate = call / (nrInds1 + nrInds2);

			if (!Double.isNaN(alleleFreq0) && snpCallRate >= 0.10) {
				double denominator = alleleFreq0 * (1.0d - alleleFreq0);
				if (!Double.isNaN(denominator) && denominator > 0) {
					for (int i = 0; i < nrInds1; i++) {
						if (genotypeAlleles1[i][0] != -1) {
							double genotype0I = 0;
							if (0 == genotypeAlleles1[i][0]) genotype0I += .5;
							if (0 == genotypeAlleles1[i][1]) genotype0I += .5;
							for (int j = 0; j < nrInds2; j++) {
								if (genotypeAlleles2[j][0] != -1) {
									double genotype0J = 0;
									if (0 == genotypeAlleles2[j][0]) genotype0J += .5;
									if (0 == genotypeAlleles2[j][1]) genotype0J += .5;
									geneticSimilarity[i][j] += (genotype0I - alleleFreq0) * (genotype0J - alleleFreq0) / denominator;
									if (genotype0I == genotype0J) geneticSimilaritySameGenotypes[i][j]++;
									geneticSimilarityCalled[i][j]++;
								}
							}
						}
					}
				}
			}
			pb.set(snpID);
		}
		pb.close();

		for (int i = 0; i < nrInds1; i++) {
			for (int j = 0; j < nrInds2; j++) {
				geneticSimilarity[i][j] /= geneticSimilarityCalled[i][j];
				geneticSimilaritySameGenotypes[i][j] /= geneticSimilarityCalled[i][j];
			}
		}

		return new Pair<double[][], double[][]>(geneticSimilarity, geneticSimilaritySameGenotypes);
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
