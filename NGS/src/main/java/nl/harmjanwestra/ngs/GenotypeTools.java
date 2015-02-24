package nl.harmjanwestra.ngs;

import com.xeiam.xchart.BitmapEncoder;
import com.xeiam.xchart.Chart;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.VCFVariantType;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.HWE;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.*;

/**
 * Created by hwestra on 2/4/15.
 */
public class GenotypeTools {


	public static void main(String[] args) {
		try {
			GenotypeTools t = new GenotypeTools();
//
//			t.compareFamFileWithSampleList("/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC.fam",
//					"/Data/Projects/2014-FR-Reseq/RASequencingSamples-ReWrite-seqIdToImmunoChip.txt");
//
//			String mapfile = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/filtered/ra.map";
//			String mapfileout = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/filtered/ra_liftover.map";
//			String liftoverbed = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/filtered/ra_liftover.bed";
////			t.convertPostLiftOverMAP(mapfile, mapfileout, liftoverbed);
//			t.rewriteRSNamesUsingDBSNP();
			String vcf = args[0];
			String plink = args[1];
			String out = args[2];
			t.summarizeVCF(vcf, out + "Summary.txt");
			t.determineVCFSummaryStatistics(vcf, out + "SampleCallrate.txt");
			t.compareVCFGenotypesToPedAndMap(vcf, plink, out);

//			String seqIdToImmunoChip = "/Data/Projects/2014-FR-Reseq/RASequencingSamples-ReWrite-seqIdToImmunoChip.txt";
//			String famfileIn = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/filtered/ra.ped";
//			String famfileOut = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/filtered/ra_rewrite.ped";
//			t.rewriteFamFileSampleNames(seqIdToImmunoChip, famfileIn, famfileOut);


		} catch (IOException e) {
			e.printStackTrace();

		}

	}

	private void rewriteFamFileSampleNames(String seqIdToImmunoChip, String famfileIn, String famfileOut) throws IOException {
		TextFile tf = new TextFile(seqIdToImmunoChip, TextFile.R);
		HashMap<String, String> idToId = new HashMap<String, String>();
		String[] elems = tf.readLineElems(Strings.tab);
		while (elems != null) {
			idToId.put(elems[1], elems[0]);
			elems = tf.readLineElems(Strings.tab);
		}
		tf.close();

		TextFile tf2 = new TextFile(famfileIn, TextFile.R);
		TextFile tf2out = new TextFile(famfileOut, TextFile.W);
		elems = tf2.readLineElems(Strings.tab);
		while (elems != null) {


			String sample = elems[1];
			String newId = idToId.get(sample);
			if (newId != null) {
				elems[1] = newId;
			}

			tf2out.writeln(Strings.concat(elems, Strings.tab));

			elems = tf2.readLineElems(Strings.tab);
		}
		tf2out.close();
		tf2.close();

	}


	private void rewriteRSNamesUsingDBSNP() throws IOException {
		String rsMergeFile = "/Data/dbSNP/b142/RsMergeArch.bcp";
		String rsRemoveFile = "";
		String mapfileIn = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/filtered/ra_liftover.map";
		String mapfileOut = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/filtered/ra_liftover_newRSIds.map";
		int dbsnpB = 138;
		HashMap<String, String> rsToNewRs = new HashMap<String, String>();

		TextFile tfmap = new TextFile(mapfileIn, TextFile.R);
		String[] elems = tfmap.readLineElems(Strings.whitespace);
		while (elems != null) {
			String rs = elems[1];
			// System.out.println("Found a snp: " + rs);
			rsToNewRs.put(rs, null);
			elems = tfmap.readLineElems(Strings.whitespace);
		}
		tfmap.close();
		System.out.println("Done parsing map file: " + mapfileIn);
		HashSet<String> removedRs = new HashSet<String>();

//			TextFile tfrem = new TextFile(rsRemoveFile, TextFile.R);
//			elems = tfrem.readLineElems(TextFile.tab);
//			while (elems != null) {
//
//
//
//				elems = tfrem.readLineElems(TextFile.tab);
//			}
//			tfrem.close();


		System.out.println("Parsing: " + rsMergeFile);
		TextFile tfRsMerge = new TextFile(rsMergeFile, TextFile.R);
		elems = tfRsMerge.readLineElems(TextFile.tab);

		HashMap<String, String> lowToHigh = new HashMap<String, String>();
		HashSet<String> visitedLow = new HashSet<String>();
		int ln = 0;
		while (elems != null) {
			String high = "rs" + elems[0]; // original ID (possibly in our MAP file)
			String low = "rs" + elems[1]; // new ID
			String build = elems[2];
			Integer bInt = Integer.parseInt(build);
			if (bInt <= dbsnpB) {
				if (rsToNewRs.containsKey(high)) { // check whether this SNP is relevant
					if (lowToHigh.containsKey(low)) { // check whether we've already assigned this low RS to another variant in our list
						System.err.println("WARNING: introducing duplicate: " + high + "\tand " + lowToHigh.get(low) + " share " + low);

						int i = 1;
						String newLow = low + "_dup" + i;
						while (lowToHigh.containsKey(newLow)) {
							newLow = low + "_dup" + i;
							i++;
						}

						lowToHigh.put(newLow, high);
						rsToNewRs.put(high, newLow);
					} else {

						if (rsToNewRs.containsKey(low)) {
							System.err.println("WARNING: replacing: " + high + " with " + low + " but " + low + " is already in the map file.");

							int i = 1;
							String newLow = low + "_dup" + i;
							while (lowToHigh.containsKey(newLow)) {
								newLow = low + "_dup" + i;
								i++;
							}

							rsToNewRs.put(high, newLow);
							lowToHigh.put(newLow, high);
							System.out.println("New Name: " + newLow);
						} else {
							lowToHigh.put(low, high);
							rsToNewRs.put(high, low);
						}
					}


				}
			}
			elems = tfRsMerge.readLineElems(TextFile.tab);
			ln++;
			if (ln % 1000000 == 0) {
				System.out.println(ln + " lns parsed.");
			}
		}
		tfRsMerge.close();

		System.out.println("Parsing stuff");
		TextFile outmap = new TextFile(mapfileOut, TextFile.W);
		tfmap.open();
		elems = tfmap.readLineElems(Strings.whitespace);
		while (elems != null) {
			String rs = elems[1];
			String newStr = rsToNewRs.get(rs);


			if (newStr == null) {
				// just output

			} else {
				System.out.println("Replacing: " + rs + " by " + newStr);
				elems[1] = newStr;

			}
			String outln = Strings.concat(elems, Strings.tab);
			outmap.writeln(outln);

			elems = tfmap.readLineElems(Strings.whitespace);
		}
		tfmap.close();
		outmap.close();

		System.out.println("Done");
	}

	private void convertPostLiftOverMAP(String mapfile, String mapfileout, String liftoverbed) throws IOException {


		HashMap<String, Integer> variantToPos = new HashMap<String, Integer>();
		TextFile tf1 = new TextFile(liftoverbed, TextFile.R);
		String[] bedelems = tf1.readLineElems(TextFile.tab);
		while (bedelems != null) {
			String variant = bedelems[3];
			Integer pos = Integer.parseInt(bedelems[1]);
			variantToPos.put(variant, pos);
			bedelems = tf1.readLineElems(TextFile.tab);
		}
		tf1.close();

		TextFile tf = new TextFile(mapfile, TextFile.R);
		TextFile out = new TextFile(mapfileout, TextFile.W);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String variant = elems[1];
			Integer newPosition = variantToPos.get(variant);
			if (newPosition == null) {
				for (int i = 0; i < elems.length; i++) {
					elems[i] = "hg18_" + elems[i];
				}
			} else {
				elems[3] = "" + newPosition;
			}

			out.writeln(Strings.concat(elems, Strings.tab));
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		out.close();
	}

	private HashMap<String, Integer> index(String[] s) {
		int c = 0;
		HashMap<String, Integer> index = new HashMap<String, Integer>();
		for (String s1 : s) {
			index.put(s1, c);
			c++;
		}
		return index;
	}

	private void determineVCFSummaryStatistics(String vcf, String out) throws IOException {
		ArrayList<String> vcfSamples = getVCFSamples(vcf);
		System.out.println(vcfSamples.size() + " samples loaded from VCF");

		HashMap<String, Integer> sampleMap = index(vcfSamples.toArray(new String[0]));
// get variants from VCF
		ArrayList<Feature> variantsOnVCF = getVariantsFromVCF(vcf);
		System.out.println(variantsOnVCF.size() + " variants on VCF");

		ArrayList<String> autosomal = new ArrayList<String>();
		HashMap<Feature, Integer> variantMap = new HashMap<Feature, Integer>();
		int ctr = 0;
		for (Feature f : variantsOnVCF) {
			if (!f.getChromosome().equals(Chromosome.X)) {
				variantMap.put(f, ctr);
				autosomal.add(f.getChromosome().getName() + ":" + f.getStart());
				ctr++;
			}
		}
		Pair<byte[][], String[]> genotypePair = loadVCFGenotypes(vcf, sampleMap, variantMap);
		// byte[sample][variant]

		writeGenotypes(genotypePair, vcfSamples, autosomal, out + "allVCFGenotypes.txt");
		byte[][] genotypes = genotypePair.getLeft();
		String[] alleles = genotypePair.getRight();
		Triple<Pair<double[], String[]>, double[], double[]> params = getGenotypeParams(genotypePair);
		double[] mafs = params.getLeft().getLeft();


		TextFile outf = new TextFile(out, TextFile.W);
		double[] xvals = new double[genotypes.length];
		double[] sampleCR = new double[genotypes.length];
		for (int sample = 0; sample < genotypes.length; sample++) {
			int nrTotal = 0;
			int nrCalled = 0;
			for (int gt = 0; gt < genotypes[sample].length; gt++) {
				//if (mafs[gt] > 0.01) {
				if (genotypes[sample][gt] == -1) {

				} else {
					nrCalled++;
				}
				nrTotal++;
				//}
			}
			double cr = (double) nrCalled / nrTotal;
			sampleCR[sample] = cr;
			outf.writeln(vcfSamples.get(sample) + "\t" + cr + "\t" + nrCalled + "\t" + nrTotal);
		}
		outf.close();

		Arrays.sort(sampleCR);

		Coverage c = new Coverage();
		Chart barchart = c.initBarChart("Call-rate per sample", "Sample", "Call-rate", 1200, 800, true, 50);
		barchart.addSeries("Series1", xvals, sampleCR);
		barchart.getStyleManager().setLegendVisible(false);
		BitmapEncoder.saveBitmapWithDPI(barchart, out + "SampleCallRate", BitmapEncoder.BitmapFormat.PNG, 300);


	}


	/*
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
	 */

	private Triple<Pair<double[], String[]>, double[], double[]> getGenotypeParams(Pair<byte[][], String[]> genotypePair) {
		byte[][] genotypes = genotypePair.getLeft();
		String[] alleles = genotypePair.getRight();

		String[] minorAlleles = new String[alleles.length];
		double[] mafs = new double[alleles.length];
		double[] crs = new double[alleles.length];
		double[] hwes = new double[alleles.length];
		for (int gt = 0; gt < genotypes[0].length; gt++) {
			String[] gtAlleles = alleles[gt].split("/");
			String allele1 = gtAlleles[0];
			String allele2 = gtAlleles[1];
			int nrCalled = 0;
			int nrTotal = 0;
			int nrHomAA = 0;
			int nrHomBB = 0;
			int nrHets = 0;
			for (int sample = 0; sample < genotypes.length; sample++) {
				if (genotypes[sample][gt] == -1) {
					// missing
				} else if (genotypes[sample][gt] == 0) {
					nrHomAA++;
					nrCalled += 2;
				} else if (genotypes[sample][gt] == 1) {
					nrHets++;
					nrCalled += 2;
				} else {
					nrHomBB++;
					nrCalled += 2;
				}
				nrTotal += 2;
			}
			int freqAlleleA = 2 * nrHomAA + nrHets; // genotypeFreq[0] + genotypeFreq[1];
			int freqAlleleB = 2 * nrHomBB + nrHets; // genotypeFreq[2] + genotypeFreq[1];
			double maf = (double) (freqAlleleA) / (nrCalled);
			String minorAllele = allele1;
			if (freqAlleleA > freqAlleleB) {
				maf = 1 - maf;
				minorAllele = allele2;
			}
			minorAlleles[gt] = minorAllele;
			double cr = (double) nrCalled / nrTotal;
			double hwep = HWE.calculateExactHWEPValue(nrHets, nrHomAA, nrHomBB);
			mafs[gt] = maf;
			crs[gt] = cr;
			hwes[gt] = hwep;
			// System.out.println(gt + "\t" + cr + "\t" + maf);
		}
		return new Triple<Pair<double[], String[]>, double[], double[]>(new Pair<double[], String[]>(mafs, minorAlleles), crs, hwes);
	}


	public void compareFamFileWithSampleList(String famfile, String sampleList) throws IOException {
		TextFile tf = new TextFile(sampleList, TextFile.R);
		String[] listofSamples = tf.readAsArray(1, TextFile.tab);
		tf.close();

		tf = new TextFile(sampleList, TextFile.R);
		String[] listofSamples2 = tf.readAsArray(0, TextFile.tab);
		tf.close();


		ArrayList<Pair<String, String>> famSamples = parseFam(famfile);
		HashMap<String, Pair<String, String>> samplesInFam = new HashMap<String, Pair<String, String>>();
		for (Pair<String, String> p : famSamples) {
			samplesInFam.put(p.getRight(), p);
		}


		for (int num = 0; num < listofSamples.length; num++) {
			String s = listofSamples[num];
			String s2 = listofSamples2[num];
			if (samplesInFam.containsKey(s)) {
				Pair<String, String> p = samplesInFam.get(s);
				System.out.println("Found:\t" + p.getLeft() + "\t" + p.getRight());
			} else {
				System.out.println("Not found:\t" + s2 + "\t" + s);
			}
		}

	}

	public void compareVCFGenotypesToPedAndMap(String vcf, String plinkfile, String out) throws IOException {

		// get samples from plink data

		ArrayList<Pair<String, String>> famSamples = parseFam(plinkfile + ".fam");
		System.out.println(famSamples.size() + " samples loaded from PED/MAP");
		// get samples from VCF
		ArrayList<String> vcfSamples = getVCFSamples(vcf);
		System.out.println(vcfSamples.size() + " samples loaded from VCF");

		// intersect samples
		ArrayList<String> samplesInFam = new ArrayList<String>();
		for (Pair<String, String> p : famSamples) {
			samplesInFam.add(p.getRight());
		}
		HashSet<String> intersectedSamples = intersectSamples(samplesInFam, vcfSamples);

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


//		System.out.println("vcf vcf correlations");
//		correlateSamples(vcfGenotypesPair, vcfGenotypesPair, samples, out + "vcfSampleCorrelations.txt", excludeVariant, flipGenotypes);
//		System.out.println("ped ped correlations");
//		correlateSamples(pedGenotypesPair, pedGenotypesPair, samples, out + "pedSampleCorrelations.txt", excludeVariant, flipGenotypes);


		Triple<Pair<double[], String[]>, double[], double[]> vcfParams = getGenotypeParams(vcfGenotypesPair);
		Triple<Pair<double[], String[]>, double[], double[]> pedParams = getGenotypeParams(pedGenotypesPair);
		double mafThreshold = 0.01;
		double crThreshold = 0.05;
		Pair<boolean[], boolean[]> flipvariants = determineAlleleFlipsAndExcludeVariants(allelesPED, pedParams.getLeft(), allelesVCF, vcfParams.getLeft(), variants, out + "vcfpedAlleleFlips.txt",
				mafThreshold, crThreshold);
		boolean[] excludeVariant = flipvariants.getLeft();
		boolean[] flipGenotypes = flipvariants.getRight();

		// correlate samples
		System.out.println("vcf ped sample correlations");


		Pair<double[][], double[][]> correlationOutputSamples = correlateSamples(vcfGenotypesPair, pedGenotypesPair, samples, out + "vcfpedSampleCorrelations.txt", excludeVariant, flipGenotypes);


		boolean[] sampleSwapped = writeBestCorrelations(correlationOutputSamples, samples, out + "vcfpedBestSampleCorrelations.txt");

		// correlate variants
		System.out.println("vcf ped genotype correlations");


		Pair<double[][], double[][]> correlationOutputVariants = correlateVariants(vcfGenotypesPair, pedGenotypesPair, variants, out + "vcfpedGenotypeCorrelations.txt", excludeVariant, flipGenotypes, sampleSwapped);
		writeBestCorrelations(correlationOutputVariants, variants, out + "vcfpedBestGenotypeCorrelations.txt");

		// write discordantGenotypes
		writeDiscordantGenotypes(vcfGenotypesPair, pedGenotypesPair, variants, out + "discordantGenotypes.txt", excludeVariant, flipGenotypes, sampleSwapped, samples);


	}

	private void writeDiscordantGenotypes(Pair<byte[][], String[]> vcfGenotypesPair, Pair<byte[][], String[]> pedGenotypesPair, ArrayList<String> variants, String outputFileName, boolean[] excludeVariant, boolean[] flipGenotypes, boolean[] sampleSwapped, ArrayList<String> samples) throws IOException {

		TextFile out = new TextFile(outputFileName, TextFile.W);

		String header = "Variant\tAllelesVCF\tAllelesPED (after flip)\tFlipped\tSample\tGenotypeVCF\tGenotypePED (after flip)\tGenotypeVCFAlleles\tGenotypePEDAlleles (after flip)";
		out.writeln(header);

		byte[][] genotypesVCF = vcfGenotypesPair.getLeft();
		String[] allelesVCF = vcfGenotypesPair.getRight();
		byte[][] genotypesPED = pedGenotypesPair.getLeft();
		String[] allelesPED = pedGenotypesPair.getRight();

		for (int v = 0; v < variants.size(); v++) {
			if (!excludeVariant[v]) {
				String variantOut = variants.get(v);

				String variantAllelesVCF = allelesVCF[v];
				String[] allelesVCFArr = variantAllelesVCF.split("/");
				String variantAllelesPED = allelesPED[v];
				String[] allelesPEDArr = variantAllelesPED.split("/");

				boolean flipped = flipGenotypes[v];
				if (flipped) {
					String all1 = allelesPEDArr[0];
					String all2 = allelesPEDArr[1];
					allelesPEDArr[0] = all2;
					allelesPEDArr[1] = all1;
				}


				variantOut += "\t" + variantAllelesVCF + "\t" + allelesPEDArr[0] + "/" + allelesPEDArr[1] + "\t" + flipped;

				for (int s = 0; s < genotypesPED.length; s++) {
					if (!sampleSwapped[s]) {

						int gt1 = genotypesVCF[s][v];
						int gt2 = genotypesPED[s][v];

						if (flipped) {
							gt2 = Math.abs(gt2 - 2);
							if (gt2 == 3) {
								gt2 = -1;
							}
						}


						if (gt1 != gt2 || variants.get(v).equals("Chr15:79152806")) {
							String sample = samples.get(s);

							String ln = variantOut + "\t" +
									sample + "\t" +
									gt1 + "\t" +
									gt2 + "\t" +
									gtToAlleles(allelesVCFArr, gt1) + "\t" +
									gtToAlleles(allelesPEDArr, gt2);
							out.writeln(ln);
						}

					}
				}
			}
		}

		out.close();

	}

	private String gtToAlleles(String[] alleles, int gt) {
		if (gt == 0) {
			return alleles[0] + "/" + alleles[0];
		} else if (gt == 1) {
			return alleles[0] + "/" + alleles[1];
		} else if (gt == 2) {
			return alleles[1] + "/" + alleles[1];
		} else {
			return "NA";
		}
	}

	private boolean[] writeBestCorrelations(Pair<double[][], double[][]> correlationInput, ArrayList<String> labels, String out) throws IOException {
		double[][] correlations = correlationInput.getLeft();
		double[][] n = correlationInput.getRight();
		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln("Label\t" +
				"Correlation\t" +
				"n\t" +
				"BestMatch\t" +
				"BestCorrelation\t" +
				"n\t" +
				"Match");
		boolean[] sampleSwapped = new boolean[labels.size()];
		for (int i = 0; i < labels.size(); i++) {

			String ln = labels.get(i) + "\t"
					+ correlations[i][i]
					+ "\t" + n[i][i];

			double max = -1;
			double origMax = -1;
			int maxS = i;
			for (int j = 0; j < labels.size(); j++) {
				double corr = correlations[i][j];
				if (!Double.isNaN(corr)) {
					double absCorr = Math.abs(corr);
					if (absCorr > max) {
						max = absCorr;
						origMax = corr;
						maxS = j;
					}
				}
			}
			if (max == -1) {
				ln += "\tUnknown" +
						"\tNaN" +
						"\t0" +
						"\tUnknown";
				sampleSwapped[i] = true;
			} else {
				ln += "\t" + labels.get(maxS) + "\t"
						+ origMax + "\t"
						+ n[i][maxS] + "\t"
						+ (maxS == i);
				sampleSwapped[i] = !(maxS == i);
			}
			outf.writeln(ln);
		}
		outf.close();
		return sampleSwapped;
	}

	private Pair<boolean[], boolean[]> determineAlleleFlipsAndExcludeVariants(String[] allelesPED, Pair<double[], String[]> pedMinorAlleles,
																			  String[] allelesVCF, Pair<double[], String[]> vcfMinorAlleles,
																			  ArrayList<String> variants, String out,
																			  double mafThreshold, double crThreshold) throws IOException {
		boolean[] excludeVariant = new boolean[allelesPED.length];
		boolean[] flipGenotypes = new boolean[allelesPED.length];

		double[] mafsPED = pedMinorAlleles.getLeft();
		String[] minorAllelesPED = pedMinorAlleles.getRight();

		double[] mafsVCF = vcfMinorAlleles.getLeft();
		String[] minorAllelesVCF = vcfMinorAlleles.getRight();

		/*
		for (int variant = 0; variant < excludeVariant.length; variant++) {
			if (mafs[variant] < mafThreshold || cr[variant] < crTheshold)) {
				excludeVariant[variant] = true;
				System.out.println("excluding variant: " + variants.get(variant) + "\t" + mafs[variant] + "\t" + cr[variant]);
			}
		}
		 */
		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln("variant\t" +
				"allelesPED\tminorAllelePED\tmafPED\t" +
				"allelesVCF\tminorAlleleVCF\tmafVCF\t" +
				"deltaMAF\tflip\tReason\tExcludeVariant");


		for (int a = 0; a < allelesPED.length; a++) {
			String allelePED = allelesPED[a];
			String alleleVCF = allelesVCF[a];
			String allele1ped = allelePED.split("/")[0];
			String allele1vcf = alleleVCF.split("/")[0];
			String allele2ped = allelePED.split("/")[1];
			String allele2vcf = alleleVCF.split("/")[1];


			String output = variants.get(a) + "\t" +
					allelesPED[a] + "\t" + minorAllelesPED[a] + "\t" + mafsPED[a] + "\t" +
					allelesVCF[a] + "\t" + minorAllelesVCF[a] + "\t" + mafsVCF[a] + "\t" + Math.abs(mafsVCF[a] - mafsPED[a]);


			if (allele2ped.equals("null") || allele1ped.equals("null") || allele1vcf.equals("null") || allele2vcf.equals("null")) {
				output += "\tfalse\tNullAllele\ttrue";
				excludeVariant[a] = true;
			} else {
				if (allele1ped.equals(BaseAnnot.getComplement(allele2ped)) && allele1vcf.equals(BaseAnnot.getComplement(allele2vcf))) {


					// A/T or G/C SNP
					// try to figure out allele flip on the basis of minor alleles

					String minorAllelePED = minorAllelesPED[a];
					String minorAlleleVCF = minorAllelesVCF[a];
					boolean flipAlleles = false;
					if ((allele1ped.equals(minorAllelePED) && allele1vcf.equals(minorAlleleVCF)) || (allele2ped.equals(minorAllelePED) && allele2vcf.equals(minorAlleleVCF))) {
						// same direction

					} else {
						flipAlleles = true;
					}

					output += "\t" + flipAlleles + "\tComplement\t" + false;
					flipGenotypes[a] = flipAlleles;

					excludeVariant[a] = false;
				} else {
					Boolean flipAlleles = BaseAnnot.flipalleles(allelePED, allele1ped, alleleVCF, allele1vcf);
					if (flipAlleles == null) {
						excludeVariant[a] = true;
						output += "\tfalse\tIncompatibleAlleles\ttrue";
					} else {
						output += "\t" + flipAlleles + "\tOK\tfalse";
						flipGenotypes[a] = flipAlleles;
					}
				}
			}
			outf.writeln(output);
		}
		outf.close();
		return new Pair<boolean[], boolean[]>(excludeVariant, flipGenotypes);
	}

	private Pair<double[][], double[][]> correlateSamples(Pair<byte[][], String[]> genotypesPair1,
														  Pair<byte[][], String[]> genotypesPair2, ArrayList<String> samples, String s,
														  boolean[] excludeCols, boolean[] flipAlleles) throws IOException {

		double[][] correlations = new double[samples.size()][samples.size()];

		byte[][] genotypes1 = genotypesPair1.getLeft();
		byte[][] genotypes2 = genotypesPair2.getLeft();
		double[][] n = new double[samples.size()][samples.size()];

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
				n[i][j] = d1.size();
				n[j][i] = d1.size();
			}
		}

		writeCorrelationMatrix(correlations, samples, s);
		writeCorrelationMatrix(n, samples, s + "-n.txt");

		return new Pair<double[][], double[][]>(correlations, n);
	}

	private Pair<double[][], double[][]> correlateVariants(Pair<byte[][], String[]> genotypesPair1,
														   Pair<byte[][], String[]> genotypesPair2, ArrayList<String> variants, String s,
														   boolean[] excludeVariants, boolean[] flipAlleles, boolean[] sampleSwapped) throws IOException {

		double[][] correlations = new double[variants.size()][variants.size()];
		double[][] n = new double[variants.size()][variants.size()];
		byte[][] genotypes1 = genotypesPair1.getLeft();
		byte[][] genotypes2 = genotypesPair2.getLeft();

		for (int variantId1 = 0; variantId1 < variants.size(); variantId1++) {

			if (!excludeVariants[variantId1]) {
				for (int variantId2 = 0; variantId2 < variants.size(); variantId2++) {
					if (!excludeVariants[variantId2]) {


						ArrayList<Double> d1 = new ArrayList<Double>();
						ArrayList<Double> d2 = new ArrayList<Double>();

						boolean flip2 = flipAlleles[variantId2];

						for (int sample = 0; sample < genotypes1.length; sample++) {

							double gt1 = (double) genotypes1[sample][variantId1];
							double gt2 = (double) genotypes2[sample][variantId2];

							if (!sampleSwapped[sample] && gt1 != -1 && gt2 != -1) {
								d1.add(gt1);

								if (flip2) {
									gt2 = Math.abs(gt2 - 2);
								}
								d2.add(gt2);
							}
						}

						double[] g1 = Primitives.toPrimitiveArr(d1.toArray(new Double[0]));
						double[] g2 = Primitives.toPrimitiveArr(d2.toArray(new Double[0]));

						double c = JSci.maths.ArrayMath.correlation(g1, g2);
						if (d1.isEmpty()) {
							if (!excludeVariants[variantId1] && !excludeVariants[variantId2]) {
								System.err.println("ERROR: no genotypes for variants: " + variants.get(variantId1) + " or " + variants.get(variantId2));
								c = Double.NaN;
							} else {
								c = Double.NaN;
							}
						}
						correlations[variantId1][variantId2] = c;
						n[variantId1][variantId2] = d1.size();
						n[variantId2][variantId1] = d1.size();
						correlations[variantId2][variantId1] = c;

						if (variants.get(variantId1).equals("Chr6:167433729") && variants.get(variantId2).equals("Chr6:167433729")) {
							for (int q = 0; q < d1.size(); q++) {
								System.out.println(q + "\t" + d1.get(q) + "\t" + d2.get(q));
							}

							System.out.println(c);
						}


					} else {
						correlations[variantId1][variantId2] = Double.NaN;
						correlations[variantId2][variantId1] = correlations[variantId1][variantId2];
					}
				}
			} else {
				for (int variantId2 = 0; variantId2 < variants.size(); variantId2++) {
					correlations[variantId1][variantId2] = Double.NaN;
					correlations[variantId2][variantId1] = correlations[variantId1][variantId2];
				}
			}

		}


		writeCorrelationMatrix(correlations, variants, s);
		writeCorrelationMatrix(n, variants, s + "-n.txt");

		return new Pair<double[][], double[][]>(correlations, n);
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


	private ArrayList<Pair<String, String>> parseFam(String famfile) throws IOException {
		ArrayList<Pair<String, String>> famSamples = new ArrayList<Pair<String, String>>();
		TextFile tf = new TextFile(famfile, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {

			String sample = elems[1];
			famSamples.add(new Pair<String, String>(elems[0], sample));
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
					if (!chr.equals(Chromosome.X)) {
						f.setChromosome(chr);
						f.setStart(pos);
						f.setStop(pos);
						output.add(f);
					}
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

	private Pair<byte[][], String[]> loadVCFGenotypes(String vcf,
													  HashMap<String, Integer> sampleMap,
													  HashMap<Feature, Integer> variantMap) throws IOException {

		TextFile tf = new TextFile(vcf, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);

		int[] colToNewSample = null;

		int nrVariantsFound = 0;

		byte[][] genotypes = new byte[sampleMap.size()][variantMap.size()];
		String[] alleles = new String[variantMap.size()];

		int nrSamplesFound = 0;
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
					// System.out.println(sample + "\t" + id);
					if (id == null) {
						colToNewSample[i] = -1;
					} else {
						colToNewSample[i] = id;
						nrSamplesFound++;
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

		System.out.println("genotypes loaded for " + nrVariantsFound + " variants using VCF for " + nrSamplesFound + " samples");
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
								// System.out.println("allele 1 mising");
							} else {
								nrCalled++;
								allele1b = Byte.parseByte(allele1);

							}
							if (allele2.equals(".")) {
								// missing
								sampleAlleles[samplePos][1] = -1;
								// System.out.println("allele 2 mising");
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

//						if (adCol != -1) {
//							// Allelic depths for the ref and alt alleles in the order listed
//							String ad = sampleElems[adCol];
//							String[] adElems = ad.split(",");
//
//						}
//
//						if (dpCol != -1&& dpCol < sampleElems.length) {
//							// Approximate read depth (reads with MQ=255 or with bad mates are filtered)
//							String dp = sampleElems[dpCol];
//							String[] dpElems = dp.split(",");
//						}
//
						if (gqCol != -1 && gqCol < sampleElems.length) {
							// Genotype Quality
//							String gq = sampleElems[gqCol];
//							String[] gqElems = gq.split(",");
							if (gqCol >= sampleElems.length) {
								System.err.println("ERROR: GQCol>length: " + gqCol + "\t" + Strings.concat(sampleElems, Strings.comma));
							}
						}
//
//						if (plCol != -1) {
//							// Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification
//							String pl = sampleElems[plCol];
//							String[] plElems = pl.split(",");
//						}


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

		// determine call rates per sample

	}
}
