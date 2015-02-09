package nl.harmjanwestra.ngs;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.VCFVariantType;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.HWE;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 2/4/15.
 */
public class VCFTools {


	public static void main(String[] args) {
		try {
			VCFTools t = new VCFTools();
			String vcf = args[0];
			String plink = args[1];
			String out = args[2];
			t.compareVCFGenotypesToPedAndMap(vcf, plink, out);
		} catch (IOException e) {
			e.printStackTrace();

		}

	}

	public void compareVCFGenotypesToPedAndMap(String vcf, String plinkfile, String out) throws IOException {

		// get samples from plink data

		ArrayList<String> famSamples = parseFam(plinkfile + ".fam");
		System.out.println(famSamples.size() + " samples loaded from PED/MAP");
		// get samples from VCF
		ArrayList<String> vcfSamples = getVCFSamples(vcf);
		System.out.println(vcfSamples.size() + " samples loaded from VCF");

		// intersect samples
		HashSet<String> intersectedSamples = intersectSamples(famSamples, vcfSamples);

		System.out.println(intersectedSamples.size() + " samples shared");
		writeHash(intersectedSamples, out + "sharedSamples.txt");

		// get variants from VCF
		ArrayList<Feature> variantsOnVCF = getVariantsFromVCF(vcf);
		System.out.println(variantsOnVCF.size() + " variants on VCF");

		// get variants from PLINK file that have same position
		ArrayList<Feature> variantsOnMap = getVariantsFromMapFile(plinkfile + ".map");
		System.out.println(variantsOnMap.size() + " variants on MAP");

		// intersect variants
		HashSet<Feature> intersectedVariants = intersectVariants(variantsOnMap, variantsOnVCF);
		System.out.println(intersectedVariants.size() + " shared variants");
		writeFeatureHash(intersectedVariants, out + "sharedVariants.txt");

		ArrayList<String> samples = new ArrayList<String>();
		samples.addAll(intersectedSamples);
		Collections.sort(samples);

		HashMap<String, Integer> sampleMap = new HashMap<String, Integer>();
		int s = 0;

		for (String sample : samples) {
			sampleMap.put(sample, s);
			s++;
		}

		HashMap<Feature, Integer> variantMap = new HashMap<Feature, Integer>();

		ArrayList<String> variants = new ArrayList<String>();
		HashMap<String, Feature> strToFeat = new HashMap<String, Feature>();
		for (Feature f : intersectedVariants) {
			variants.add(f.getChromosome().getName() + ":" + f.getStart());
			strToFeat.put(f.getChromosome().getName() + ":" + f.getStart(), f);
		}
		Collections.sort(variants);

		int v = 0;
		for (String var : variants) {
			Feature f = strToFeat.get(var);
			variantMap.put(f, v);
			v++;
		}


		// load genotypes for VCF
		// format: byte[sample][genotype] String[nrVariants];
		Pair<byte[][], String[]> vcfGenotypesPair = loadVCFGenotypes(vcf, sampleMap, variantMap);
		writeGenotypes(vcfGenotypesPair, samples, variants, out + "vcfGenotypes.txt");


		// load genotypes for PLINK
		Pair<byte[][], String[]> pedGenotypesPair = loadPedGenotypes(plinkfile + ".ped", plinkfile + ".map", sampleMap, variantMap);
		writeGenotypes(pedGenotypesPair, samples, variants, out + "pedGenotypes.txt");


		// compare alleles

		String[] allelesVCF = vcfGenotypesPair.getRight();
		String[] allelesPED = pedGenotypesPair.getRight();

		boolean[] excludeVariant = new boolean[allelesPED.length];
		boolean[] flipGenotypes = new boolean[allelesPED.length];

		System.out.println("vcf vcf correlations");
		correlateRows(vcfGenotypesPair, vcfGenotypesPair, samples, out + "vcfSampleCorrelations.txt", excludeVariant, flipGenotypes);
		System.out.println("ped ped correlations");
		correlateRows(pedGenotypesPair, pedGenotypesPair, samples, out + "pedSampleCorrelations.txt", excludeVariant, flipGenotypes);

		Pair<boolean[], boolean[]> flipvariants = determineAlleleFlips(allelesPED, allelesVCF, variants, out + "vcfpedAlleleFlips.txt");
		excludeVariant = flipvariants.getLeft();
		flipGenotypes = flipvariants.getRight();

		// correlate samples
		System.out.println("vcf ped sample correlations");
		double[][] correlations = correlateRows(vcfGenotypesPair, pedGenotypesPair, samples, out + "vcfpedSampleCorrelations.txt", excludeVariant, flipGenotypes);


		writeBestCorrelations(correlations, samples, out + "vcfpedBestSampleCorrelations.txt");

		// correlate variants
		System.out.println("vcf ped genotype correlations");
		double[][] variantCorrelations = correlateColumns(vcfGenotypesPair, pedGenotypesPair, variants, out + "vcfpedGenotypeCorrelations.txt", excludeVariant, flipGenotypes);
		writeBestCorrelations(variantCorrelations, variants, out + "vcfpedBestGenotypeCorrelations.txt");

		// determine genetic similarity in VCF

		// determine call rates per sample in VCF


	}

	private void writeBestCorrelations(double[][] correlations, ArrayList<String> labels, String out) throws IOException {
		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln("Label\tCorrelation\tBestMatch\tBestCorrelation\tMatch");
		for (int i = 0; i < labels.size(); i++) {

			String ln = labels.get(i);
			ln += "\t" + correlations[i][i];
			double max = -1;
			double origMax = -1;
			int maxS = i;
			for (int j = 0; j < labels.size(); j++) {
				double corr = correlations[i][j];
				double absCorr = Math.abs(corr);
				if (absCorr > max) {
					max = absCorr;
					origMax = corr;
					maxS = j;
				}
			}
			ln += "\t" + labels.get(maxS) + "\t" + origMax + "\t" + (maxS == i);
			outf.writeln(ln);
		}
		outf.close();
	}

	private Pair<boolean[], boolean[]> determineAlleleFlips(String[] allelesPED, String[] allelesVCF, ArrayList<String> variants, String out) throws IOException {
		boolean[] excludeVariant = new boolean[allelesPED.length];
		boolean[] flipGenotypes = new boolean[allelesPED.length];

		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln("variant\tallelesPED\tallelesVCF\tflip\tIncompatibility\tExcludeVariant");


		for (int a = 0; a < allelesPED.length; a++) {
			String allelePED = allelesPED[a];
			String alleleVCF = allelesVCF[a];
			String allele1ped = allelePED.split("/")[0];
			String allele1vcf = alleleVCF.split("/")[0];
			String allele2ped = allelePED.split("/")[1];
			String allele2vcf = alleleVCF.split("/")[1];
			Boolean flipAlleles = BaseAnnot.flipalleles(allelePED, allele1ped, alleleVCF, allele1vcf);
			if (flipAlleles == null) {


				if (allele2ped == null || allele1ped == null || allele1vcf == null || allele2vcf == null) {
					excludeVariant[a] = true;
//					System.err.println("ERROR in alleles : ped: " + allelePED + "\t" + alleleVCF + " for variant: " + variants.get(a));
					outf.writeln(variants.get(a) + "\t" + allelesPED[a] + "\t" + allelesVCF[a] + "\t" + false + "\tNull\ttrue");
				} else {
					if (isComplement(allele1ped, allele2ped) || isComplement(allele1vcf, allele2vcf)) {
//						System.err.println("Complimentary alleles : ped: " + allelePED + "\t" + alleleVCF + " for variant: " + variants.get(a));
						outf.writeln(variants.get(a) + "\t" + allelesPED[a] + "\t" + allelesVCF[a] + "\t" + false + "\tCompliment\tfalse");
					} else {
						// System.err.println("Incompatible alleles : ped: " + allelePED + "\t" + alleleVCF + " for variant: " + variants.get(a));
						outf.writeln(variants.get(a) + "\t" + allelesPED[a] + "\t" + allelesVCF[a] + "\t" + false + "\tIncompatible\ttrue");
						excludeVariant[a] = true;
					}
				}


			} else {
				if (flipAlleles) {
//					System.out.println("Alleles should be flipped: ped: " + allelePED + "\t" + alleleVCF + " for variant: " + variants.get(a));
					outf.writeln(variants.get(a) + "\t" + allelesPED[a] + "\t" + allelesVCF[a] + "\t" + flipAlleles + "\tOK\tfalse");
				} else {
//					System.out.println("Alleles equal: ped: " + allelePED + "\t" + alleleVCF + " for variant: " + variants.get(a));
					outf.writeln(variants.get(a) + "\t" + allelesPED[a] + "\t" + allelesVCF[a] + "\t" + flipAlleles + "\tOK\tfalse");
				}
				flipGenotypes[a] = flipAlleles;
			}
		}
		outf.close();
		return new Pair<boolean[], boolean[]>(excludeVariant, flipGenotypes);
	}

	private boolean isComplement(String allele1ped, String allele2ped) {
		if ((allele1ped.equals("A") && allele2ped.equals("T")) || (allele1ped.equals("T") && allele2ped.equals("A"))) {
			return true;
		} else if ((allele1ped.equals("G") && allele2ped.equals("C")) || (allele1ped.equals("C") && allele2ped.equals("G"))) {
			return true;
		} else {
			return false;
		}
	}

	private double[][] correlateRows(Pair<byte[][], String[]> genotypesPair1,
									 Pair<byte[][], String[]> genotypesPair2, ArrayList<String> samples, String s,
									 boolean[] excludeCols, boolean[] flipAlleles) throws IOException {

		double[][] correlations = new double[samples.size()][samples.size()];

		byte[][] genotypes1 = genotypesPair1.getLeft();
		byte[][] genotypes2 = genotypesPair2.getLeft();

		for (int i = 0; i < samples.size(); i++) {

			for (int j = 0; j < samples.size(); j++) {
				ArrayList<Double> d1 = new ArrayList<Double>();
				ArrayList<Double> d2 = new ArrayList<Double>();
				for (int g = 0; g < genotypes1[i].length; g++) {
					boolean flip = flipAlleles[g];
					if (genotypes1[i][g] != -1 && genotypes2[j][g] != -1 && !excludeCols[g]) {
						d1.add((double) genotypes1[i][g]);
						double gt = (double) genotypes2[j][g];
						if (flip) {
							gt = Math.abs(gt - 2);
						}
						d2.add(gt);

					}
				}

				double[] g1 = Primitives.toPrimitiveArr(d1.toArray(new Double[0]));
				double[] g2 = Primitives.toPrimitiveArr(d2.toArray(new Double[0]));
				correlations[i][j] = Correlation.correlate(g1, g2);
				correlations[j][i] = correlations[i][j];
			}
		}

		writeCorrelationMatrix(correlations, samples, s);

		return correlations;
	}

	private double[][] correlateColumns(Pair<byte[][], String[]> genotypesPair1,
										Pair<byte[][], String[]> genotypesPair2, ArrayList<String> variants, String s,
										boolean[] excludeCols, boolean[] flipAlleles) throws IOException {

		double[][] correlations = new double[variants.size()][variants.size()];

		byte[][] genotypes1 = genotypesPair1.getLeft();
		byte[][] genotypes2 = genotypesPair2.getLeft();

		for (int i = 0; i < variants.size(); i++) {

			for (int j = 0; j < variants.size(); j++) {
				ArrayList<Double> d1 = new ArrayList<Double>();
				ArrayList<Double> d2 = new ArrayList<Double>();

				boolean flip2 = flipAlleles[j];

				for (int k = 0; k < genotypes1.length; k++) {

					double gt1 = (double) genotypes2[k][i];
					double gt2 = (double) genotypes2[k][j];

					if (gt1 != -1 && gt2 != -1 && !excludeCols[i] && !excludeCols[j]) {
						d1.add(gt1);

						if (flip2) {
							gt2 = Math.abs(gt2 - 2);
						}
						d2.add(gt2);
					}
				}

				double[] g1 = Primitives.toPrimitiveArr(d1.toArray(new Double[0]));
				double[] g2 = Primitives.toPrimitiveArr(d2.toArray(new Double[0]));

				double c = Correlation.correlate(g1, g2);
				if (d1.isEmpty()) {
					if (!excludeCols[i] && !excludeCols[j]) {
						System.err.println("ERROR: no genotypes for variants: " + variants.get(i) + " or " + variants.get(j));
					} else {
						c = 0;
					}
				}

				correlations[i][j] = c;
				correlations[j][i] = correlations[i][j];

			}
		}


		writeCorrelationMatrix(correlations, variants, s);

		return correlations;
	}

	private void writeCorrelationMatrix(double[][] correlations, ArrayList<String> labels, String s) throws IOException {
		TextFile out = new TextFile(s, TextFile.W);
		String header = "-";
		for (int i = 0; i < labels.size(); i++) {
			header += "\t" + labels.get(i);
		}
		out.writeln(header);
		for (int i = 0; i < labels.size(); i++) {
			String ln = labels.get(i);
			for (int j = 0; j < correlations.length; j++) {
				ln += "\t" + correlations[i][j];
			}
			out.writeln(ln);
		}
		out.close();
	}


	private void writeGenotypes(Pair<byte[][], String[]> vcfGenotypesPair, ArrayList<String> samples, ArrayList<String> variants, String s) throws IOException {
		TextFile out = new TextFile(s, TextFile.W);

		byte[][] genotypes = vcfGenotypesPair.getLeft();
		String[] alleles = vcfGenotypesPair.getRight();
		String header = "variant\talleles";
		for (int i = 0; i < genotypes.length; i++) {
			header += "\t" + samples.get(i);
		}
		out.writeln(header);
		for (int i = 0; i < genotypes[0].length; i++) {
			String ln = variants.get(i) + "\t" + alleles[i];
			for (int j = 0; j < genotypes.length; j++) {
				ln += "\t" + genotypes[j][i];
			}
			out.writeln(ln);
		}
		out.close();

	}

	private void writeHash(HashSet<String> set, String s) throws IOException {
		TextFile out = new TextFile(s, TextFile.W);
		for (String v : set) {
			out.writeln(v);
		}
		out.close();
	}

	private void writeFeatureHash(HashSet<Feature> set, String s) throws IOException {
		TextFile out = new TextFile(s, TextFile.W);
		for (Feature v : set) {
			out.writeln(v.getChromosome().getName() + "\t" + v.getStart());
		}
		out.close();
	}



	private ArrayList<String> parseFam(String famfile) throws IOException {
		ArrayList<String> famSamples = new ArrayList<String>();
		TextFile tf = new TextFile(famfile, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {

			String sample = elems[1];
			famSamples.add(sample);
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();

		return famSamples;
	}

	private ArrayList<String> getVCFSamples(String vcf) throws IOException {

		ArrayList<String> samples = new ArrayList<String>();
		TextFile tf = new TextFile(vcf, TextFile.R);
		String[] lnElems = tf.readLineElems(TextFile.tab);
		boolean headerFound = false;
		while (lnElems != null && !headerFound) {

			if (lnElems[0].startsWith("##")) {
				// metadata
			} else if (lnElems[0].startsWith("#CHROM")) {

				// header
				for (int i = 9; i < lnElems.length; i++) {
					// T1D/RA specific stuff
					String[] sampleNameElems = lnElems[i].split("/");
					if (sampleNameElems.length == 1) {
						samples.add(sampleNameElems[0]);
					} else if (sampleNameElems.length == 2) {

						// split further
						String tmpSampleName = sampleNameElems[1];
						String[] tmpSampleNameElems = tmpSampleName.split("_");
						samples.add(tmpSampleNameElems[0]);

					} else {

						System.err.println("Error: unexpected sample name");
					}
				}
				headerFound = true;


			}

			lnElems = tf.readLineElems(TextFile.tab);
		}

		tf.close();
		return samples;
	}

	private HashSet<String> intersectSamples(ArrayList<String> s1, ArrayList<String> s2) {
		HashSet<String> intersect = new HashSet<String>();
		HashSet<String> set1 = new HashSet<String>();
		set1.addAll(s1);
		for (String s : s2) {
			if (set1.contains(s)) {
				intersect.add(s);
			}
		}
		return intersect;
	}


	private ArrayList<Feature> getVariantsFromVCF(String vcf) throws IOException {
		TextFile tf = new TextFile(vcf, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		ArrayList<Feature> output = new ArrayList<Feature>();
		while (elems != null) {

			Chromosome chr = Chromosome.NA;
			if (elems[0].startsWith("##")) {
				// superheader

			} else if (elems[0].startsWith("#CHROM")) {
				// header
				// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1   SAMPLE
			} else {
				// line
				String ref = elems[3];
				String alt = elems[4];
				chr = Chromosome.parseChr(elems[0]);
				String[] alternates = alt.split(",");
				String[] alleles = new String[1 + alternates.length];
				alleles[0] = ref.intern();
				for (int i = 0; i < alternates.length; i++) {
					alleles[i + 1] = alternates[i].intern();
				}


				VCFVariantType type = VCFVariantType.parseType(ref, alt);
				boolean biAllelic = type.isBiallelic();
				if (biAllelic) {
					// determine MAF
					Integer pos = Integer.parseInt(elems[1]);
					Feature f = new Feature();
					f.setChromosome(chr);
					f.setStart(pos);
					f.setStop(pos);
					output.add(f);
				}

			}


			elems = tf.readLineElems(TextFile.tab);
		}


		tf.close();
		return output;
	}

	private ArrayList<Feature> getVariantsFromMapFile(String mapfile) throws IOException {
		ArrayList<Feature> output = new ArrayList<Feature>();

		TextFile tf = new TextFile(mapfile, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {

			Chromosome chr = Chromosome.parseChr(elems[0]);

			Feature f = new Feature();
			f.setChromosome(chr);
			Integer position = -1;
			try {
				position = Integer.parseInt(elems[3]);
			} catch (NumberFormatException e) {

			}


			f.setStart(position);
			f.setStop(position);
			output.add(f);
			elems = tf.readLineElems(Strings.whitespace);
		}

		tf.close();
		return output;
	}

	private HashSet<Feature> intersectVariants(ArrayList<Feature> f1, ArrayList<Feature> f2) {
		HashSet<Feature> shared = new HashSet<Feature>();
		for (Feature feat1 : f1) {
			for (Feature feat2 : f2) {
				if (feat1.overlaps(feat2)) {
					shared.add(feat1);
				}
			}
		}
		return shared;

	}


	private Pair<byte[][], String[]> loadPedGenotypes(String ped,
													  String map,
													  HashMap<String, Integer> sampleMap,
													  HashMap<Feature, Integer> variantMap) throws IOException {

		// load list of variants to include
		ArrayList<Feature> mapVariants = getVariantsFromMapFile(map);


		int[] variantToNewVariant = new int[mapVariants.size()];
		for (int v = 0; v < mapVariants.size(); v++) {
			Integer newVariantId = variantMap.get(mapVariants.get(v));
			if (newVariantId != null) {
				variantToNewVariant[v] = newVariantId;
			} else {
				variantToNewVariant[v] = -1;
			}
		}

		byte[][] genotypes = new byte[sampleMap.size()][variantMap.size()];
		String[] alleles1 = new String[variantMap.size()];
		String[] alleles2 = new String[variantMap.size()];

		// iterate ped file
		/*
		Family ID
     Individual ID
     Paternal ID
     Maternal ID
     Sex (1=male; 2=female; other=unknown)
     Phenotype
		 */

		TextFile tf = new TextFile(ped, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);
		int nrSamplesLoaded = 0;
		while (elems != null) {
			String sample = elems[1];
			Integer newSampleId = sampleMap.get(sample);

			// check length of the line
			if ((elems.length - 6) / 2 != mapVariants.size()) {
				System.err.println("ERROR: different number of variants in ped than in map! Found " + elems.length + "-6 columns in PED, but found " + mapVariants.size() + " in map file.");
				System.exit(-1);
			} else {
				if (newSampleId != null) {

					for (int i = 6; i < elems.length; i += 2) {

						int newVariantId = variantToNewVariant[(i - 6) / 2];

						if (newVariantId != -1) {
							String allele1 = elems[i];
							String allele2 = elems[i + 1];
							if (allele1.equals("A") || allele1.equals("C") || allele1.equals("T") || allele1.equals("G")) {
								if (allele2.equals("A") || allele2.equals("C") || allele2.equals("T") || allele2.equals("G")) {
									if (alleles1[newVariantId] == null) {
										alleles1[newVariantId] = allele1;
										if (!allele1.equals(allele2)) {
											alleles2[newVariantId] = allele2;
										}
									} else {

										if (alleles2[newVariantId] == null) {
											String prevAllele1 = alleles1[newVariantId];
											if (!allele1.equals(prevAllele1)) {
												alleles2[newVariantId] = allele1;
											} else {
												if (!allele2.equals(prevAllele1)) {
													alleles2[newVariantId] = allele2;
												}
											}
										}

										if (alleles1[newVariantId] != null && alleles2[newVariantId] != null) {
											// do a quick check to see if we're doing OK
											String prevAllele1 = alleles1[newVariantId];
											String prevAllele2 = alleles2[newVariantId];
											boolean variantOK = true;
											if (!prevAllele1.equals(allele1) && !prevAllele2.equals(allele1)) {
												Feature f = mapVariants.get(i - 6);
												System.out.println("Found new allele for variant: " + f.getChromosome().getName() + ":" + f.getStart() + " " + allele1 + " found where " + alleles1[newVariantId] + "/" + alleles2[newVariantId] + " expected");
												variantOK = false;
											}
											if (!prevAllele1.equals(allele2) && !prevAllele2.equals(allele2)) {
												Feature f = mapVariants.get(i - 6);
												System.out.println("Found new allele for variant: " + f.getChromosome().getName() + ":" + f.getStart() + " " + allele2 + " found where " + alleles1[newVariantId] + "/" + alleles2[newVariantId] + " expected");
												variantOK = false;
											}
											if (!variantOK) {
												System.exit(-1);
											}
										}
									}
								}
							}
						}
					}
					nrSamplesLoaded++;
					if (nrSamplesLoaded % 50 == 0) {
						System.out.println(nrSamplesLoaded + " samples loaded");
					}
				}
			}
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();


		tf.open();
		elems = tf.readLineElems(Strings.whitespace);
		nrSamplesLoaded = 0;
		while (elems != null) {
			String sample = elems[1];
			Integer newSampleId = sampleMap.get(sample);

			// check length of the line
			if ((elems.length - 6) / 2 != mapVariants.size()) {
				System.err.println("ERROR: different number of variants in ped than in map! Found " + elems.length + "-6 columns in PED, but found " + mapVariants.size() + " in map file.");
				System.exit(-1);
			} else {
				if (newSampleId != null) {

					for (int i = 6; i < elems.length; i += 2) {

						int newVariantId = variantToNewVariant[(i - 6) / 2];

						if (newVariantId != -1) {
							String allele1 = elems[i];
							String allele2 = elems[i + 1];
							if (allele1.equals("A") || allele1.equals("C") || allele1.equals("T") || allele1.equals("G")) {
								if (allele2.equals("A") || allele2.equals("C") || allele2.equals("T") || allele2.equals("G")) {

									String prevAllele1 = alleles1[newVariantId];
									String prevAllele2 = alleles2[newVariantId];

									byte genotype = -1;
									if (allele1.equals(prevAllele1) && allele2.equals(prevAllele1)) {
										genotype = 0;
									} else if (allele1.equals(prevAllele2) && allele2.equals(prevAllele2)) {
										genotype = 2;
									} else {
										genotype = 1;
									}
									genotypes[newSampleId][newVariantId] = genotype;


								} else {
									// System.out.println("Unknown allele: " + allele1 + " - " + allele2);
									genotypes[newSampleId][newVariantId] = -1;
								}
							} else {
								// System.out.println("Unknown allele: " + allele1 + " - " + allele2);
								genotypes[newSampleId][newVariantId] = -1;
							}
						}
					}
					nrSamplesLoaded++;
					if (nrSamplesLoaded % 50 == 0) {
						System.out.println(nrSamplesLoaded + " samples loaded");
					}
				}
			}
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();


		String[] alleles = new String[alleles1.length];

		for (int v = 0; v < alleles1.length; v++) {
			alleles[v] = alleles1[v] + "/" + alleles2[v];
		}
		return new Pair<byte[][], String[]>(genotypes, alleles);
	}

	private Pair<byte[][], String[]> loadVCFGenotypes(String vcf, HashMap<String, Integer> sampleMap, HashMap<Feature, Integer> variantMap) throws IOException {

		TextFile tf = new TextFile(vcf, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);

		int[] colToNewSample = null;

		int nrVariantsFound = 0;

		byte[][] genotypes = new byte[sampleMap.size()][variantMap.size()];
		String[] alleles = new String[variantMap.size()];

		while (elems != null) {

			Chromosome chr = Chromosome.NA;
			if (elems[0].startsWith("##")) {
				// superheader

			} else if (elems[0].startsWith("#CHROM")) {
				// header
				// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1   SAMPLE
				// header
				colToNewSample = new int[elems.length];
				for (int i = 9; i < elems.length; i++) {
					// T1D/RA specific stuff
					String[] sampleNameElems = elems[i].split("/");
					String sample = "";
					if (sampleNameElems.length == 1) {

						sample = sampleNameElems[0];

					} else if (sampleNameElems.length == 2) {

						// split further
						String tmpSampleName = sampleNameElems[1];
						String[] tmpSampleNameElems = tmpSampleName.split("_");
						sample = tmpSampleNameElems[0];

					} else {

						System.err.println("Error: unexpected sample name");
					}


					Integer id = sampleMap.get(sample);
					if (id == null) {
						colToNewSample[i] = -1;
					} else {
						colToNewSample[i] = id;
					}
				}
			} else {
				// line
				String ref = elems[3];
				String alt = elems[4];
				VCFVariantType type = VCFVariantType.parseType(ref, alt);
				boolean biAllelic = type.isBiallelic();
				if (biAllelic) {
					chr = Chromosome.parseChr(elems[0]);

					// determine MAF
					Integer pos = Integer.parseInt(elems[1]);
					Feature f = new Feature();
					f.setChromosome(chr);
					f.setStart(pos);
					f.setStop(pos);
					Integer newVariantId = variantMap.get(f);
					if (newVariantId == null) {
						// don't do anything just yet
					} else {
						nrVariantsFound++;

						// load the genotypes
						String[] alternates = alt.split(",");
						String[] currentalleles = new String[1 + alternates.length];
						currentalleles[0] = ref.intern();
						for (int i = 0; i < alternates.length; i++) {
							currentalleles[i + 1] = alternates[i].intern();
						}

						String[] format = elems[8].split(":");
						int gtCol = -1; // genotype
						int adCol = -1; // Allelic depths for the ref and alt alleles in the order listed
						int dpCol = -1; // Approximate read depth (reads with MQ=255 or with bad mates are filtered)
						int gqCol = -1; // Genotype Quality
						int plCol = -1; // Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification
						for (int c = 0; c < format.length; c++) {
							if (format[c].equals("GT")) {
								gtCol = c;
							} else if (format[c].equals("AD")) {
								adCol = c;
							} else if (format[c].equals("DP")) {
								dpCol = c;
							} else if (format[c].equals("GQ")) {
								gqCol = c;
							} else if (format[c].equals("PL")) {
								plCol = c;
							}
						}

						alleles[newVariantId] = ref + "/" + alt;

						for (int e = 9; e < elems.length; e++) {
							int newSampleId = colToNewSample[e];
							if (newSampleId > -1) {
								String sampleColumn = elems[e];
								if (sampleColumn.equals("./.")) {
									genotypes[newSampleId][newVariantId] = -1;
								} else {
									String[] sampleElems = sampleColumn.split(":");

									if (gtCol != -1) {
										String gt = sampleElems[gtCol];
										String[] gtElems = gt.split("/");

										if (gtElems[0].equals(".") || gtElems[1].equals(".")) {
											genotypes[newSampleId][newVariantId] = -1;
										} else {
											Byte g1 = Byte.parseByte(gtElems[0]);
											Byte g2 = Byte.parseByte(gtElems[1]);

											byte dosage = g1;
											dosage += g2;
											genotypes[newSampleId][newVariantId] = dosage;
										}


									}
								}
							}
						}
					}
				}
			}


			elems = tf.readLineElems(TextFile.tab);
		}

		System.out.println("genotypes loaded for " + nrVariantsFound + " variants using VCF");
		tf.close();
		return new Pair<byte[][], String[]>(genotypes, alleles);
	}


	public void summarizeVCF(String vcffile, String outputFileLoc) throws IOException {

		TextFile tf = new TextFile(vcffile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);

		String[] samples = null;

		HashMap<String, Integer> infoToColumn = new HashMap<String, Integer>();

		int novelSNPs = 0;
		int novelIndels = 0;
		int knownIndels = 0;
		int knownSNPs = 0;

		int nrVariants = 0;
		ArrayList<String> columnNames = new ArrayList<String>();

		while (elems != null) {

			if (elems[0].startsWith("##")) {
				// superheader

			} else if (elems[0].startsWith("#CHROM")) {
				// header
				// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1   SAMPLE2
			} else {
				// line
				String[] infoColumns = elems[7].split(";");
				for (String infoColumn : infoColumns) {
					// column name separated by  =
					String[] infocolumnelems = infoColumn.split("=");
					String infocolumnName = infocolumnelems[0];
					if (!infoToColumn.containsKey(infocolumnName)) {
						infoToColumn.put(infocolumnName, infoToColumn.size());
						columnNames.add(infocolumnName);
					}
				}

				nrVariants++;
			}
			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();


		System.out.println(nrVariants + " variants and  " + infoToColumn.size() + " pieces of info per variant");

		double[][] values = new double[nrVariants][infoToColumn.size()];
		boolean[] biAllelic = new boolean[nrVariants];
		boolean[] isKnown = new boolean[nrVariants];

		int variant = 0;
		tf.open();
		elems = tf.readLineElems(TextFile.tab);
		TextFile outputFileWriter = new TextFile(outputFileLoc, TextFile.W);
		String header = "Chrom\tPos\tId\tRef\tAlt\tQual\tFilter\tType\tisBiallelic\tCallRate\tHWEP\tMAF\tminorAllele\tnrAlleles\t" + Strings.concat(columnNames, Strings.tab);

		outputFileWriter.writeln(header);
		while (elems != null) {

			Chromosome chr = Chromosome.NA;
			if (elems[0].startsWith("##")) {
				// superheader

			} else if (elems[0].startsWith("#CHROM")) {
				// header
				// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1   SAMPLE2
				samples = new String[elems.length - 9];

				System.arraycopy(elems, 9, samples, 0, elems.length - 9);

			} else {
				// line
				String ref = elems[3];
				String alt = elems[4];
				chr = Chromosome.parseChr(elems[0]);
				String[] alternates = alt.split(",");
				String[] alleles = new String[1 + alternates.length];
				alleles[0] = ref.intern();
				for (int i = 0; i < alternates.length; i++) {
					alleles[i + 1] = alternates[i].intern();
				}


				VCFVariantType type = VCFVariantType.parseType(ref, alt);
				biAllelic[variant] = type.isBiallelic();
				String[] infoColumns = elems[7].split(";");
				for (String infoColumn : infoColumns) {
					// column name separated by  =
					String[] infocolumnelems = infoColumn.split("=");
					String infocolumnName = infocolumnelems[0];
					if (infocolumnelems.length == 1) {
						Integer id = infoToColumn.get(infocolumnName);
						values[variant][id] = 1d;

					} else {
						String value = infocolumnelems[1];

						Integer id = infoToColumn.get(infocolumnName);

						if (infocolumnName.equals("culprit")) {
							Integer id2 = infoToColumn.get(value);
							values[variant][id] = id2;
						} else {

							try {
								values[variant][id] = Double.parseDouble(value);
							} catch (NumberFormatException e) {
								String[] valueElems = value.split(",");
								if (valueElems.length > 1 && infocolumnName.equals("MLEAC") ||
										infocolumnName.equals("MLEAF") ||
										infocolumnName.equals("AC") ||
										infocolumnName.equals("AF")) {
									values[variant][id] = -1; // now we know that these are multi-allelic
								} else {
									System.out.println("Error: could not parse " + value + " for infoElem: " + infocolumnName + " and variant: " + elems[0] + "\t" + elems[1] + "\t" + elems[2]);
								}

							}

						}
					}
				}

				String[] format = elems[8].split(":");
				int gtCol = -1; // genotype
				int adCol = -1; // Allelic depths for the ref and alt alleles in the order listed
				int dpCol = -1; // Approximate read depth (reads with MQ=255 or with bad mates are filtered)
				int gqCol = -1; // Genotype Quality
				int plCol = -1; // Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification
				for (int c = 0; c < format.length; c++) {
					if (format[c].equals("GT")) {
						gtCol = c;
					} else if (format[c].equals("AD")) {
						adCol = c;
					} else if (format[c].equals("DP")) {
						dpCol = c;
					} else if (format[c].equals("GQ")) {
						gqCol = c;
					} else if (format[c].equals("PL")) {
						plCol = c;
					}
				}

				byte[] sampleGenotypes = new byte[elems.length - 9];

				byte[][] sampleAlleles = new byte[elems.length - 9][2];

				int nrHets = 0;
				int nrHomAA = 0;
				int nrHomBB = 0;

				int nrTotal = 0;
				int nrCalled = 0;
				for (int e = 9; e < elems.length; e++) {
					int samplePos = e - 9;
					nrTotal += 2;
					String sampleColumn = elems[e];
					if (sampleColumn.equals("./.")) {
						// not called
						sampleGenotypes[samplePos] = -1;
					} else {

						String[] sampleElems = sampleColumn.split(":");

						if (gtCol != -1) {
							String gt = sampleElems[gtCol];
							String[] gtElems = gt.split("/");
							String allele1 = gtElems[0];
							String allele2 = gtElems[1];
							byte allele1b = -1;
							byte allele2b = -1;
							if (allele1.equals(".")) {
								// missing
								sampleAlleles[samplePos][0] = -1;
								System.out.println("allele 1 mising");
							} else {
								nrCalled++;
								allele1b = Byte.parseByte(allele1);

							}
							if (allele2.equals(".")) {
								// missing
								sampleAlleles[samplePos][1] = -1;
								System.out.println("allele 2 mising");
							} else {
								nrCalled++;
								allele2b = Byte.parseByte(allele2);
							}


							if (type.isBiallelic()) {
								// this should change for chr X: there we should only count alleles for females
								if (allele1b != -1 && allele2b != -1) {
									if (allele1b == allele2b) {
										sampleGenotypes[samplePos] = allele1b; // homozygote AA
										if (allele1b == 1) {
											sampleGenotypes[samplePos]++; // homozygote BB
											nrHomBB++;
										} else {
											nrHomAA++;
										}
									} else {
										sampleGenotypes[samplePos] = 1; // heterozygote AB
										nrHets++;
									}

								}

							} else {
								// leave it for now..
								sampleGenotypes[samplePos] = -1;
							}

						}

						if (adCol != -1) {
							// Allelic depths for the ref and alt alleles in the order listed
							String ad = sampleElems[adCol];
							String[] adElems = ad.split(",");

						}

						if (dpCol != -1) {
							// Approximate read depth (reads with MQ=255 or with bad mates are filtered)
							String dp = sampleElems[dpCol];
							String[] dpElems = dp.split(",");
						}

						if (gqCol != -1) {
							// Genotype Quality
							String gq = sampleElems[gqCol];
							String[] gqElems = gq.split(",");
						}

						if (plCol != -1) {
							// Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification
							String pl = sampleElems[plCol];
							String[] plElems = pl.split(",");
						}


					}
				}

				// determine MAF, HWEP, CR
				double maf = 0;
				double hwep = 1;
				double cr = (double) nrCalled / nrTotal;
				String minorAllele = alleles[0];
				if (type.isBiallelic()) {
					int freqAlleleA = 2 * nrHomAA + nrHets; // genotypeFreq[0] + genotypeFreq[1];
					int freqAlleleB = 2 * nrHomBB + nrHets; // genotypeFreq[2] + genotypeFreq[1];

					maf = (double) (freqAlleleA) / (nrCalled);
					minorAllele = alleles[0];
					if (freqAlleleA > freqAlleleB) {
						minorAllele = alleles[1];
						maf = 1 - maf;
					}
					hwep = HWE.calculateExactHWEPValue(nrHets, nrHomAA, nrHomBB);

				}

				// write the info elems to disk
				outputFileWriter.append(elems[0] + "\t" +
						elems[1] + "\t" +
						elems[2] + "\t" +
						elems[3] + "\t" +
						elems[4] + "\t" +
						elems[5] + "\t" +
						elems[6]);
				outputFileWriter.append("\t");
				outputFileWriter.append(type.toString());
				outputFileWriter.append("\t");
				outputFileWriter.append("" + type.isBiallelic());
				outputFileWriter.append("\t");
				outputFileWriter.append("" + cr);
				outputFileWriter.append("\t");
				outputFileWriter.append("" + hwep);
				outputFileWriter.append("\t");
				outputFileWriter.append("" + maf);
				outputFileWriter.append("\t");
				outputFileWriter.append("" + minorAllele);
				outputFileWriter.append("\t");
				outputFileWriter.append("" + alleles.length);

				for (int i = 0; i < values[variant].length; i++) {
					outputFileWriter.append("\t");
					outputFileWriter.append("" + values[variant][i]);

				}
				outputFileWriter.append("\n");

				variant++;
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		outputFileWriter.close();

		//

	}
}
