package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.genotypes.GenotypeTools;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.text.Strings;

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
		System.out.println(samples1.size() + " samples in " + vcf1);
		VCFGenotypeData data2 = new VCFGenotypeData(vcf2);
		ArrayList<String> samples2 = data2.getSamples();
		System.out.println(samples2.size() + " samples in " + vcf2);

		// get sample overlap
		HashMap<String, Integer> overlapSampleMap = new HashMap<>();
		ArrayList<String> overlappingSamples = new ArrayList<>();
		HashSet<String> samples1Hash = new HashSet<String>();
		samples1Hash.addAll(samples1);
		int ctr = 0;
		for (String sample : samples2) {
			if (samples1Hash.contains(sample) && !overlapSampleMap.containsKey(sample)) {
				overlapSampleMap.put(sample, ctr);
				overlappingSamples.add(sample);
				ctr++;
			}
		}

		Pair<int[], ArrayList<String>> p1 = reorderSamples(samples1, overlapSampleMap);
		ArrayList<String> samples1Reordered = new ArrayList<>();
		samples1Reordered.addAll(overlappingSamples);
		samples1Reordered.addAll(p1.getRight());
		int[] sample1Order = p1.getLeft();

		Pair<int[], ArrayList<String>> p2 = reorderSamples(samples2, overlapSampleMap);
		ArrayList<String> samples2Reordered = new ArrayList<>();
		samples2Reordered.addAll(overlappingSamples);
		samples2Reordered.addAll(p2.getRight());
		int[] sample2Order = p2.getLeft();

		System.out.println(overlapSampleMap.size() + " samples shared between VCFs");
		HashMap<String, VCFVariant> variantMap = new HashMap<>();
		while (data1.hasNext()) {
			VCFVariant variant1 = data1.next();
			Chromosome chr = Chromosome.parseChr(variant1.getChr());
			if (chr.isAutosome() && variant1.getMAF() > 0.01) {
				variantMap.put(variant1.toString(), variant1);
			}
		}

		System.out.println(variantMap.size() + " variants loaded from: " + vcf1);

		VCFMerger merger = new VCFMerger();

		ArrayList<VCFVariant> variantsOverlap1 = new ArrayList<>();
		ArrayList<VCFVariant> variantsOverlap2 = new ArrayList<>();
		int wrongalleles = 0;
		while (data2.hasNext()) {
			VCFVariant variant2 = data2.next();

			VCFVariant variant1 = variantMap.get(variant2.toString());

			if (variant1 != null) {

				// TODO: this should really be put in a separate class.
				// check if alleles are comparable
				String[] alleles1 = variant1.getAlleles();
				String minorAllele1 = variant1.getMinorAllele();
				String[] alleles2 = variant2.getAlleles();
				String minorAllele2 = variant2.getMinorAllele();

				if (minorAllele1 != null && minorAllele2 != null) {

					int nridenticalalleles = merger.countIdenticalAlleles(alleles1, alleles2);
					GenotypeTools gtools = new GenotypeTools();

					boolean complement = false;
					if (variant1.getAlleles().length == 2 && variant2.getAlleles().length == 2) {
						if (nridenticalalleles == 0) {
							// try complement
							System.out.println("Before complement:\t\t" + Strings.concat(alleles1, Strings.forwardslash) + "\t" + minorAllele1
									+ "\t" + Strings.concat(alleles2, Strings.forwardslash) + "\t" + minorAllele2 + "\t" + nridenticalalleles);

							complement = true;
							String[] complementAlleles2 = gtools.convertToComplement(alleles2);
							minorAllele2 = gtools.getComplement(minorAllele2);
							nridenticalalleles = merger.countIdenticalAlleles(alleles1, complementAlleles2);

							System.out.println("After complement:\t\t" + Strings.concat(alleles1, Strings.forwardslash) + "\t" + minorAllele1
									+ "\t" + Strings.concat(complementAlleles2, Strings.forwardslash) + "\t" + minorAllele2 + "\t" + nridenticalalleles);


						}

						// simple case: both are biallelic..
						// check if the minor alleles are equal. else, skip the variant.
						if (nridenticalalleles == 2) {
							if (minorAllele2.equals(minorAllele1) || (variant2.getMAF() > 0.45 && variant1.getMAF() > 0.45)) {
								// check whether the reference allele is equal
								String[] tmpAlleles = alleles2;
								if (complement) {
									variant2.convertAllelesToComplement();
									tmpAlleles = variant2.getAlleles();
								}

								if (!alleles1[0].equals(tmpAlleles[0])) {
									variant2.flipReferenceAlelele();

									alleles2 = variant2.getAlleles();
									minorAllele2 = variant2.getMinorAllele();
									System.out.println("After flip:\t\t" + Strings.concat(alleles1, Strings.forwardslash) + "\t" + minorAllele1
											+ "\t" + Strings.concat(alleles2, Strings.forwardslash) + "\t" + minorAllele2);

								}
								variantsOverlap1.add(variant1);
								variantsOverlap2.add(variant2);
							}
						} else {
							System.out.println("Incompatible alleles:\t\t" + Strings.concat(alleles1, Strings.forwardslash) + "\t" + minorAllele1
									+ "\t" + Strings.concat(alleles2, Strings.forwardslash) + "\t" + minorAllele2);
						}
					}
				}


			}
		}

		// clear some memory
		variantMap = null;

		// now put everything in a matrix
		System.out.println(variantsOverlap1.size() + " variants overlap..." + wrongalleles + " have wrong alleles?");
		byte[][] genotypematrix1 = new byte[samples1.size()][variantsOverlap1.size()];
		byte[][] genotypematrix2 = new byte[samples2.size()][variantsOverlap1.size()];
		double[][] corrmat = new double[samples1.size()][samples2.size()];
		for (int i = 0; i < variantsOverlap1.size(); i++) {
			VCFVariant variant1 = variantsOverlap1.get(i);
			VCFVariant variant2 = variantsOverlap2.get(i);
			putInMatrix(genotypematrix1, variant1, sample1Order, i);
			putInMatrix(genotypematrix2, variant2, sample2Order, i);
		}

		// now correlate the samples
		for (int i = 0; i < samples1.size(); i++) {
			for (int j = 0; j < samples2.size(); j++) {
				corrmat[i][j] = correlate(genotypematrix1[i], genotypematrix2[j]);
			}
		}

		// save the correlation matrix
		DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<>();
		ds.setMatrix(corrmat);
		try {
			ds.setRowObjects(samples1Reordered);
			ds.setColObjects(samples2Reordered);
		} catch (Exception e) {

		}

		ds.save(out + "-matrix.txt.gz");

		TextFile outf = new TextFile(out + "samples1vs2.txt", TextFile.W);
		outf.writeln("Sample\tCorrSelf\tBestMatch\tBestMatchCorr\tIdenticalSample");
		// get highest absolute correlation per individual
		for (int i = 0; i < samples1Reordered.size(); i++) {
			double maxCorr = 0;
			int jmax = -1;
			for (int j = 0; j < samples2Reordered.size(); j++) {
				if (corrmat[i][j] > maxCorr) {
					maxCorr = corrmat[i][j];
					jmax = j;
				}
			}

			String sample = samples1Reordered.get(i);
			if (overlapSampleMap.containsKey(sample)) {
				outf.writeln(samples1Reordered.get(i)
						+ "\t" + corrmat[i][i]
						+ "\t" + samples2Reordered.get(jmax)
						+ "\t" + maxCorr
						+ "\t" + (samples1Reordered.get(i).equals(samples2Reordered.get(jmax)))
				);
			} else {
				outf.writeln(samples1Reordered.get(i)
						+ "\t" + null
						+ "\t" + samples2Reordered.get(jmax)
						+ "\t" + maxCorr
						+ "\t" + (samples1Reordered.get(i).equals(samples2Reordered.get(jmax)))
				);
			}
		}

		outf.close();


		outf = new TextFile(out + "samples2vs1.txt", TextFile.W);
		// get highest absolute correlation per individual
		for (int i = 0; i < samples2Reordered.size(); i++) {
			double maxCorr = 0;
			int jmax = -1;
			for (int j = 0; j < samples1Reordered.size(); j++) {
				if (corrmat[j][i] > maxCorr) {
					maxCorr = corrmat[j][i];
					jmax = j;
				}
			}

			String sample = samples2Reordered.get(i);
			if (overlapSampleMap.containsKey(sample)) {
				outf.writeln(samples2Reordered.get(i)
						+ "\t" + corrmat[i][i]
						+ "\t" + samples1Reordered.get(jmax)
						+ "\t" + maxCorr
						+ "\t" + (samples1Reordered.get(i).equals(samples2Reordered.get(jmax)))
				);
			} else {
				outf.writeln(samples2Reordered.get(i)
						+ "\t" + null
						+ "\t" + samples1Reordered.get(jmax)
						+ "\t" + maxCorr
						+ "\t" + (samples1Reordered.get(i).equals(samples2Reordered.get(jmax)))
				);
			}
			outf.writeln(samples2Reordered.get(i) + "\t" + samples1Reordered.get(jmax) + "\t" + maxCorr);

		}

		outf.close();


	}

	private Pair<int[], ArrayList<String>> reorderSamples(ArrayList<String> samples, HashMap<String, Integer> overlapSampleMap) {

		int[] sampleIndex = new int[samples.size()];
		for (int i = 0; i < sampleIndex.length; i++) {
			sampleIndex[i] = -1;
		}
		int ctr = 0;

		ArrayList<String> nonoverlapping = new ArrayList<>();
		for (int i = 0; i < samples.size(); i++) {
			String sample = samples.get(i);
			Integer index = overlapSampleMap.get(sample);
			if (index == null) {
				index = overlapSampleMap.size() + ctr;
				nonoverlapping.add(sample);
				ctr++;
			}
			sampleIndex[i] = index;
		}

		return new Pair<int[], ArrayList<String>>(sampleIndex, nonoverlapping);

	}

	private double correlate(byte[] b1, byte[] b2) {

		// remove missing genotypes

		int nrmissing = 0;
		for (int i = 0; i < b1.length; i++) {
			if (b1[i] == -1 || b2[i] == -1) {
				nrmissing++;
			}
		}
		if (nrmissing > 0) {
			byte[] tmpb1 = new byte[b1.length - nrmissing];
			byte[] tmpb2 = new byte[b1.length - nrmissing];
			int ctr = 0;
			for (int i = 0; i < b1.length; i++) {
				if (b1[i] != -1 && b2[i] != -1) {
					tmpb1[ctr] = b1[i];
					tmpb2[ctr] = b2[i];
					ctr++;
				}
			}
			b1 = tmpb1;
			b2 = tmpb2;
		}

		return Correlation.correlate(toDouble(b1), toDouble(b2));

	}

	public double[] toDouble(byte[] b) {
		double[] o = new double[b.length];
		for (int i = 0; i < b.length; i++) {
			o[i] = b[i];
		}
		return o;
	}

	private void putInMatrix(byte[][] genotypematrix1, VCFVariant variant1, int[] sampleOrder, int variantId) {
		byte[][] alleles = variant1.getGenotypeAlleles();
		for (int i = 0; i < alleles[0].length; i++) {
			if (alleles[0][i] == -1) {
				genotypematrix1[i][variantId] = -1;
			} else {
				byte genotype = (byte) (alleles[0][i] + alleles[1][i]);
				int sampleIndex = sampleOrder[i];

				if (sampleIndex >= genotypematrix1.length) {
					System.err.println("ERROR: sample index higher than number of shared samples?");
					for (int q = 0; q < sampleOrder.length; q++) {
						if ((sampleOrder[q] >= genotypematrix1.length)) {
							System.out.println(q + "\t" + sampleOrder[q] + "\t" + (sampleOrder[q] >= genotypematrix1.length));
						}
					}
					System.exit(-1);
				}
				if (sampleIndex != -1) {
					genotypematrix1[sampleIndex][variantId] = genotype;
				}
			}
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
