package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.genotypes.GenotypeTools;
import nl.harmjanwestra.utilities.vcf.VCFFunctions;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by Harm-Jan on 02/15/16.
 */
public class MixupTest {

	public void test(String vcf1, String vcf2, String out) throws IOException {

		VCFGenotypeData data1 = new VCFGenotypeData(vcf1);
		ArrayList<String> samples1 = data1.getSamples();

		VCFGenotypeData data2 = new VCFGenotypeData(vcf2);
		ArrayList<String> samples2 = data2.getSamples();


		HashSet<String> samples1Hash = new HashSet<String>();
		samples1Hash.addAll(samples1);


		ArrayList<String> overlapSamples = new ArrayList<>();
		int ctr = 0;
		for (String sample : samples2) {
			if (samples1Hash.contains(sample)) {
				overlapSamples.add(sample);
			}
		}


		int[] index1 = index(samples1, overlapSamples);
		int[] index2 = index(samples2, overlapSamples);

		HashMap<String, VCFVariant> variantMap = new HashMap<>();

		while (data1.hasNext()) {
			VCFVariant variant1 = data1.next();
			variantMap.put(variant1.toString(), variant1);
		}

		VCFMerger merger = new VCFMerger();

		ArrayList<VCFVariant> variantsOverlap1 = new ArrayList<>();
		ArrayList<VCFVariant> variantsOverlap2 = new ArrayList<>();
		while (data2.hasNext()) {
			VCFVariant variant2 = data2.next();
			VCFVariant variant1 = variantMap.get(variant2.toString());

			if (variant1 != null) {

				// TODO: this should really be put in a separate class.
				// check if alleles are comparable
				String[] refAlleles = variant1.getAlleles();
				String refMinorAllele = variant1.getMinorAllele();
				String[] testVariantAlleles = variant2.getAlleles();
				String testVariantMinorAllele = variant2.getMinorAllele();


				int nridenticalalleles = merger.countIdenticalAlleles(refAlleles, testVariantAlleles);
				GenotypeTools gtools = new GenotypeTools();

				boolean complement = false;
				if (nridenticalalleles == 0) {
					// try complement
					complement = true;
					String[] complementAlleles2 = gtools.convertToComplement(testVariantAlleles);
					testVariantMinorAllele = gtools.getComplement(testVariantMinorAllele);
					nridenticalalleles = merger.countIdenticalAlleles(refAlleles, complementAlleles2);
				}

				VCFFunctions t = new VCFFunctions();

				String logoutputln = "";
				boolean flipped = false;
				if (variant1.getAlleles().length == 2 && variant2.getAlleles().length == 2) {
					// simple case: both are biallelic..
					// check if the minor alleles are equal. else, skip the variant.
					if (nridenticalalleles == 2) {
						if (testVariantMinorAllele.equals(refMinorAllele) || (variant2.getMAF() > 0.45 && variant1.getMAF() > 0.45)) {
							// check whether the reference allele is equal
							String[] tmpAlleles = testVariantAlleles;
							if (complement) {
								variant2.convertAllelesToComplement();
								tmpAlleles = variant2.getAlleles();
							}

							if (!refAlleles[0].equals(tmpAlleles[0])) {
								variant2.flipReferenceAlelele();
							}
							variantsOverlap1.add(variant1);
							variantsOverlap2.add(variant2);
						} else {
							// different minor allele
						}
					} else {
						// write to log?
						// incompatible alleles
					}
				}

			}
		}

		// clear some memory
		variantMap = null;

		// now put everything in a matrix

		byte[][] genotypematrix1 = new byte[overlapSamples.size()][variantsOverlap1.size()];
		byte[][] genotypematrix2 = new byte[overlapSamples.size()][variantsOverlap1.size()];
		double[][] corrmat = new double[overlapSamples.size()][overlapSamples.size()];
		for (int i = 0; i < variantsOverlap1.size(); i++) {
			VCFVariant variant1 = variantsOverlap1.get(i);
			VCFVariant variant2 = variantsOverlap2.get(i);
			putInMatrix(genotypematrix1, variant1, index1, i);
			putInMatrix(genotypematrix2, variant2, index2, i);
		}

		// now correlate the samples
		for (int i = 0; i < corrmat.length; i++) {
			for (int j = 0; j < corrmat.length; j++) {
				corrmat[i][j] = correlate(genotypematrix1[i], genotypematrix2[j]);
			}
		}

		// save the correlation matrix
		DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<>();
		ds.setMatrix(corrmat);
		try {
			ds.setColObjects(overlapSamples);
			ds.setRowObjects(overlapSamples);

		} catch (Exception e) {

		}

		ds.save(out);

	}

	private double correlate(byte[] b1, byte[] b2) {
		return Correlation.correlate(toDouble(b1), toDouble(b2));

	}

	public double[] toDouble(byte[] b) {
		double[] o = new double[b.length];
		for (int i = 0; i < b.length; i++) {
			o[i] = b[i];
		}
		return o;
	}

	private void putInMatrix(byte[][] genotypematrix1, VCFVariant variant1, int[] index1, int variantId) {
		byte[][] alleles = variant1.getGenotypeAlleles();
		for (int i = 0; i < index1.length; i++) {
			int origSample = index1[i];
			byte genotype = (byte) (alleles[0][origSample] + alleles[1][origSample]);
			genotypematrix1[i][variantId] = genotype;
		}
	}

	private int[] index(ArrayList<String> samples, ArrayList<String> overlap) {
		HashMap<String, Integer> tmp = new HashMap<String, Integer>();
		for (int i = 0; i < samples.size(); i++) {
			tmp.put(samples.get(i), i);
		}

		int[] index = new int[overlap.size()];

		for (int i = 0; i < overlap.size(); i++) {
			index[i] = tmp.get(overlap.get(i));
		}

		return index;
	}

}
