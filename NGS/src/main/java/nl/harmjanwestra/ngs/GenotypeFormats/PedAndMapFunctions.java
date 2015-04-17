package nl.harmjanwestra.ngs.GenotypeFormats;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by hwestra on 4/8/15.
 */
public class PedAndMapFunctions {

	public void filterFAM(String individuals, String famfile, String out) throws IOException {

		TextFile tf = new TextFile(individuals, TextFile.R);

		Set<String> set = tf.readAsSet(0, TextFile.tab);
		tf.close();

		TextFile tf2 = new TextFile(famfile, TextFile.R);
		TextFile tfout = new TextFile(out, TextFile.W);

		String[] elems = tf2.readLineElems(Strings.whitespace);

		HashSet<String> detectedSamples = new HashSet<String>();
		while (elems != null) {

			if (set.contains(elems[1])) {
				tfout.writeln(elems[0] + "\t" + elems[1]);
				detectedSamples.add(elems[1]);
			}
			elems = tf2.readLineElems(Strings.whitespace);
		}

		tfout.close();
		tf2.close();

		for (String s : set) {
			if (!detectedSamples.contains(s)) {
				System.out.println("sample " + s + " not found");
			}
		}

		// iterate ped file
		/*
		Family ID
     Individual ID
     Paternal ID
     Maternal ID
     Sex (1=male; 2=female; other=unknown)
     Phenotype
		 */

	}


	public void deduplicateMAP(String map1, String map2) throws IOException {
		TextFile tf1 = new TextFile(map1, TextFile.R);

		TextFile out = new TextFile(map2, TextFile.W);
		String[] elems = tf1.readLineElems(TextFile.tab);
		HashSet<String> visitedSNPs = new HashSet<String>();
		while (elems != null) {

			String snp = elems[1];
			int ctr = 0;
			while (visitedSNPs.contains(snp)) {
				snp = elems[1] + "_" + ctr;
				ctr++;
			}

			elems[1] = snp;
			visitedSNPs.add(snp);

			out.writeln(Strings.concat(elems, Strings.tab));
			elems = tf1.readLineElems(TextFile.tab);
		}
		out.close();

		tf1.close();

	}

	public void comparePlinkAlleleFrequencies(String freqfile1, String freqfile2, String outf) throws IOException {

		HashMap<String, String> snpToMinorAllele = new HashMap<String, String>();
		HashMap<String, String> snpToMAF = new HashMap<String, String>();
		HashMap<String, String> snpToMissingnessssss = new HashMap<String, String>();
		HashMap<String, String> snpToAlternate = new HashMap<String, String>();


		TextFile tf = new TextFile(freqfile1, TextFile.R);
		tf.readLine();
		String ln = tf.readLine();
		while (ln != null) {

			while (ln.contains("  ")) {
				ln = ln.replaceAll("  ", " ");
			}
			String[] lnElems = ln.split(" ");
//			int i = 0;
//			for (String s : lnElems) {
//				System.out.println(i + "\t" + s);
//				i++;
//
//			}
//			System.exit(-1);
			String snp = lnElems[2];
			String minorallele = lnElems[3];
			String alternate = lnElems[4];
			String maf = lnElems[5];
			String missingnessssss = lnElems[6];
			snpToMinorAllele.put(snp, minorallele);
			snpToMAF.put(snp, maf);
			snpToMissingnessssss.put(snp, missingnessssss);
			snpToAlternate.put(snp, alternate);
			ln = tf.readLine();
		}
		tf.close();

		TextFile tf2 = new TextFile(freqfile2, TextFile.R);
		tf2.readLine();
		TextFile out = new TextFile(outf, TextFile.W);
		out.writeln("SNP\tMinorAlleleA\tMafA\tMinorAlleleB\tMafB" +
				"\tMissingnessA\tMissingnessB\tAbsDeltaMaf\tMinorAlleleAEqualsB\tSimpleStrandFlip?\tChange?");
		TextFile outlf = new TextFile(outf + "-lowfreq.txt", TextFile.W);
		ln = tf2.readLine();

		while (ln != null) {

			while (ln.contains("  ")) {
				ln = ln.replaceAll("  ", " ");
			}
			String[] lnElems = ln.split(" ");

			String snp = lnElems[2];
			String minorallele = lnElems[3];
			String maf = lnElems[5];
			String alternate = lnElems[4];
			String missingnessssss = lnElems[6];
			String otherMinorAllele = snpToMinorAllele.get(snp);
			String otherMAF = snpToMAF.get(snp);
			String othermissingnesssss = snpToMissingnessssss.get(snp);
			String otheralternate = snpToAlternate.get(snp);
			Double mafx = 0d;
			Double mafx2 = 0d;
			try {
				if (otherMAF != null) {
					mafx = Double.parseDouble(otherMAF);
				}
				mafx2 = Double.parseDouble(maf);

				boolean include = false;
				String mod = "";
				if (mafx == 0 || mafx2 == 0) {
					include = false;
					mod = "monomorph";
				} else if (minorallele.equals(otherMinorAllele)) {
					include = true;
					mod = "none";
				} else {
					if (mafx > (0.5 - 0.03) || mafx2 > (0.50 - 0.03)) {
						mod = "highfreqvariant";
						include = true;
					} else if (Math.abs(mafx - mafx2) < 0.03 || Math.abs(mafx - mafx2) > (0.5 - 0.03)) {
						// one allele is very rare, the other very common: probably one of the allleles is close to 0.5
						include = true;
						mod = "highorlowfreqfreqvariant";
					} else if ((minorallele.equals(otheralternate) && alternate.equals(otherMinorAllele))) {
						// classic allele flip
						include = true;
						mod = "flip allele";
					}
				}

				if (include) {
					out.writeln(snp
							+ "\t" + otherMinorAllele + "/" + otheralternate
							+ "\t" + otherMAF
							+ "\t" + minorallele + "/" + alternate
							+ "\t" + maf
							+ "\t" + othermissingnesssss
							+ "\t" + missingnessssss
							+ "\t" + Math.abs(mafx - mafx2)
							+ "\t" + (minorallele.equals(otherMinorAllele)
							+ "\t" + (minorallele.equals(otheralternate) && alternate.equals(otherMinorAllele)))
							+ "\t" + mod);
				} else {
					outlf.writeln(snp + "\t" + otherMAF + "\t" + maf + "\tfailedfilter " + mod);
				}

			} catch (NumberFormatException x) {
				outlf.writeln(snp + "\t" + otherMAF + "\t" + maf + "\tnullallele");
			}

			ln = tf2.readLine();
		}
		out.close();
		tf2.close();
		outlf.close();
	}

	public void rewriteMapFileChromosomeNames(String refMap, String inMap, String outMap) throws IOException {

		TextFile tf = new TextFile(refMap, TextFile.R);
		Set<String> refRS = tf.readAsSet(1, TextFile.tab);
		tf.close();

		TextFile tf2 = new TextFile(inMap, TextFile.R);
		Set<String> inRS = tf2.readAsSet(1, TextFile.tab);
		tf2.close();

		int counter = 0;
		for (String s : inRS) {
			if (!refRS.contains(s)) {
				System.out.println(s + " not found");
				counter++;
			}
		}
		System.out.println(counter);

		tf.open();
		HashMap<String, String> rsToChr = new HashMap<String, String>();
		HashMap<String, String> rsToChrPos = new HashMap<String, String>();

		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			rsToChr.put(elems[1], elems[0]);
			rsToChrPos.put(elems[1], elems[3]);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile outf = new TextFile(outMap, TextFile.W);
		tf2.open();
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[1];
			String chr = rsToChr.get(snp);
			String pos = rsToChrPos.get(snp);

			elems[0] = chr;
			elems[3] = pos;

			outf.writeln(Strings.concat(elems, Strings.tab));

			elems = tf2.readLineElems(TextFile.tab);
		}

		tf2.close();
		outf.close();


	}

	public void updateRSNames(String dbsnpvcf, String mapin, String mapout) throws IOException {
		HashMap<String, String> map = new HashMap<String, String>();
		TextFile tf = new TextFile(mapin, TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			map.put(elems[0] + "_" + elems[3], elems[1]);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();


		TextFile vcf = new TextFile(dbsnpvcf, TextFile.R);

		String[] lineElems = vcf.readLineElems(TextFile.tab);
		int ln = 0;
		while (lineElems != null) {

			if (lineElems[0].startsWith("#")) {
				// header
			} else {
				String query = lineElems[0] + "_" + lineElems[1];
				String rs = lineElems[2];
				if (map.containsKey(query)) {
					map.put(query, rs);
				}
			}

			ln++;
			if (ln % 1000000 == 0) {
				System.out.println(ln + " positions parsed...");
			}
			lineElems = vcf.readLineElems(TextFile.tab);
		}
		vcf.close();

		TextFile out = new TextFile(mapout, TextFile.W);
		tf.open();
		elems = tf.readLineElems(TextFile.tab);
		int nrReplaced = 0;
		while (elems != null) {
			String query = elems[0] + "_" + elems[3];
			String rs = map.get(query);
			if (!rs.equals(elems[1])) {
				System.out.println("Replacing " + elems[1] + " with rs: " + rs);
				nrReplaced++;
			}
			elems[1] = rs;

			out.writeln(Strings.concat(elems, Strings.tab));
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println(nrReplaced + " totally replaced...");


		out.close();
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


	private void rewriteRSNamesUsingDBSNP(String rsMergeFile, String mapfileIn, String mapfileOut) throws IOException {

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

	public void convertPostLiftOverMAP(String mapfile, String mapfileout, String liftoverbed) throws IOException {


		HashMap<String, Integer> variantToPos = new HashMap<String, Integer>();
		HashMap<String, String> variantToChr = new HashMap<String, String>();
		TextFile tf1 = new TextFile(liftoverbed, TextFile.R);
		String[] bedelems = tf1.readLineElems(TextFile.tab);
		while (bedelems != null) {
			String variant = bedelems[3];
			Integer pos = Integer.parseInt(bedelems[1]) - 1;
			variantToPos.put(variant, pos);
			variantToChr.put(variant, bedelems[0]);
			bedelems = tf1.readLineElems(TextFile.tab);
		}
		tf1.close();

		TextFile tf = new TextFile(mapfile, TextFile.R);
		TextFile out = new TextFile(mapfileout, TextFile.W);
		TextFile outf2 = new TextFile(mapfileout + "excludeThese.txt", TextFile.W);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String variant = elems[1];
			Integer newPosition = variantToPos.get(variant);
			if (newPosition == null) {
//				for (int i = 0; i < elems.length; i++) {
//					elems[i] = "hg18_" + elems[i];
//
//				}
				elems[1] = "hg18_" + elems[1];
				outf2.writeln(elems[1]);
			} else {
				elems[0] = "" + variantToChr.get(variant).replaceAll("chr", "");
				elems[3] = "" + newPosition;
			}

			out.writeln(Strings.concat(elems, Strings.tab));
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		outf2.close();
		out.close();
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

	public ArrayList<Pair<String, String>> parseFam(String famfile) throws IOException {
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

	public ArrayList<Feature> getVariantsFromMapFile(String mapfile) throws IOException {
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

	public Pair<byte[][], String[]> loadPedGenotypes(String ped,
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

	public void rewriteMapToBed(String mapIn, String bedOut) throws IOException {
		TextFile tf = new TextFile(mapIn, TextFile.R);
		TextFile out = new TextFile(bedOut, TextFile.W);
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {

			//chr rs cm/Mb pos
			String chrStr = elems[0];
			if (!elems[0].startsWith("chr")) {
				if (elems[0].equals("23")) {
					chrStr = "chrX";
				} else if (elems[0].equals("24")) {
					chrStr = "chrY";
				} else {
					chrStr = "chr" + elems[0];
				}

			}
			out.writeln(chrStr + "\t" + (Integer.parseInt(elems[3]) + 1) + "\t" + (Integer.parseInt(elems[3]) + 2) + "\t" + elems[1]);

			elems = tf.readLineElems(Strings.whitespace);
		}

		tf.close();
		out.close();


	}

	public ArrayList<String> readMAPFile(String mapfile) throws IOException {
		TextFile tf = new TextFile(mapfile, TextFile.R);


		HashSet<String> uniqueFeatures = new HashSet<String>();
		ArrayList<String> allFeats = new ArrayList<String>();
		HashSet<String> duplicates = new HashSet<String>();

		String[] lineElems = tf.readLineElems(TextFile.tab);
		while (lineElems != null) {

			String f = lineElems[0] + "_" + lineElems[3] + "@" + lineElems[1];
			if (uniqueFeatures.contains(f)) {
				duplicates.add(f);
			}
			uniqueFeatures.add(f);


			lineElems = tf.readLineElems(TextFile.tab);
		}
		tf.close();


		tf.open();
		lineElems = tf.readLineElems(TextFile.tab);
		while (lineElems != null) {
			String f = lineElems[0] + "_" + lineElems[3] + "@" + lineElems[1];
			if (duplicates.contains(f)) {
				f = lineElems[0] + "_" + lineElems[3] + "@" + lineElems[1];
			}
			allFeats.add(f);
			lineElems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		return allFeats;
	}


	public void updateMapFileRsIdsUsingMapFile(String refMap, String inMap, String outMap) throws IOException {

		TextFile tf = new TextFile(refMap, TextFile.R);
		Set<String> refRS = tf.readAsSet(1, TextFile.tab);
		tf.close();

		TextFile tf2 = new TextFile(inMap, TextFile.R);
		Set<String> inRS = tf2.readAsSet(1, TextFile.tab);
		tf2.close();

		int counter = 0;
		for (String s : inRS) {
			if (!refRS.contains(s)) {
				// System.out.println(s + " not found");
				counter++;
			}
		}
		System.out.println(counter + " variants not found using RS id");

		tf.open();
		HashMap<String, HashSet<String>> chrToRS = new HashMap<String, HashSet<String>>();


		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String query = elems[0] + "_" + elems[3];
			HashSet<String> rsssssss = chrToRS.get(query);
			if (rsssssss == null) {
				rsssssss = new HashSet<String>();
			}
			rsssssss.add(elems[1]);

			chrToRS.put(query, rsssssss);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();


		Set<String> dups = chrToRS.keySet();
		int dupctr = 0;
		for (String keuy : dups) {
			if (chrToRS.get(keuy).size() > 1) {
				dupctr++;
			}
		}
		System.out.println(dupctr + " dups in mapfile");

		int counter2 = 0;
		TextFile outf = new TextFile(outMap, TextFile.W);
		tf2.open();
		elems = tf2.readLineElems(TextFile.tab);

		// also ensure uniqueness of the rsIds
		HashSet<String> visitedRsIds = new HashSet<String>();
		while (elems != null) {
			String query = elems[0] + "_" + elems[3];

			HashSet<String> rsssssss = chrToRS.get(query);
			String currentVar = elems[1];
			if (rsssssss == null) {
				//System.out.println("Variant not found: " + query + "\t" + currentVar);
				counter2++;
			} else {
				if (rsssssss.contains(currentVar)) {
					// don't change
					//System.out.println("not replacing: " + currentVar);
					visitedRsIds.add(currentVar);
				} else {
					String[] rssses = rsssssss.toArray(new String[0]);
					if (rssses.length > 1) {
						System.out.println("replacing: " + currentVar + " with " + Strings.concat(rssses, Strings.comma) + " for query: " + query);

						// also don't change.
						visitedRsIds.add(currentVar);
					} else {
						int ctrrr = 0;
						String newVar = rssses[0];
						while (visitedRsIds.contains(newVar)) {
							newVar = rssses[0] + "_" + ctrrr;
							ctrrr++;
						}
						elems[1] = newVar;

					}

				}
			}


			outf.writeln(Strings.concat(elems, Strings.tab));

			elems = tf2.readLineElems(TextFile.tab);
		}

		System.out.println(counter2 + " variants not found using position");
		tf2.close();
		outf.close();


	}

	public void identifyDuplicatesInMap(String liftedmergedMapfile, String liftedmergedMapFileDedupped) throws IOException {
		TextFile tf = new TextFile(liftedmergedMapfile, TextFile.R);

		HashSet<String> set1 = new HashSet<String>();
		HashSet<String> set2 = new HashSet<String>();
		HashSet<String> dups1 = new HashSet<String>();
		HashSet<String> dups2 = new HashSet<String>();

		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 4) {
				String snp = elems[1];
				String pos = elems[0] + "_" + elems[3];
				if (set1.contains(snp)) {
					dups1.add(snp);
				}
				if (set2.contains(pos)) {
					dups2.add(pos);
				}
				set1.add(snp);
				set2.add(pos);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		System.out.println(set1.size());
		System.out.println(set2.size());

	}

	public void filterMap(String mapFile, String regionFile, String variantsToKeep) throws IOException {
		ArrayList<Feature> regions = readRegionFile(regionFile);

		TextFile in = new TextFile(mapFile, TextFile.R);

		String[] elems = in.readLineElems(TextFile.tab);
		TextFile tfout = new TextFile(variantsToKeep, TextFile.W);
		while (elems != null) {
			String chr = elems[0];
			Integer start = Integer.parseInt(elems[3]);
			Integer stop = Integer.parseInt(elems[3]);
			Feature f = new Feature();
			f.setChromosome(Chromosome.parseChr(chr));
			f.setStart(start);
			f.setStop(stop);
			boolean overlap = false;
			for (Feature r : regions) {
				if (r.overlaps(f)) {
					overlap = true;
				}
			}

			if (overlap) {
				// write SNP ids for PLINK
				tfout.writeln(elems[1]);
			}

			elems = in.readLineElems(TextFile.tab);
		}
		tfout.close();

		in.close();

	}

	public ArrayList<Feature> readRegionFile(String bed) throws IOException {
		TextFile tf = new TextFile(bed, TextFile.R);
		ArrayList<Feature> features = new ArrayList<Feature>();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String chr = elems[0];
			Integer start = Integer.parseInt(elems[1]);
			Integer stop = Integer.parseInt(elems[2]);
			Feature f = new Feature();
			f.setChromosome(Chromosome.parseChr(chr));
			f.setStart(start);
			f.setStop(stop);
			features.add(f);
			elems = tf.readLineElems(TextFile.tab);
		}
		return features;
	}

	public void filterOutIncompleteTrios(String famFileIn, String samplesToExclude) throws IOException {

		TextFile in = new TextFile(famFileIn, TextFile.R);


		HashSet<String> allSamples = new HashSet<String>();
		String[] elems = in.readLineElems(TextFile.tab);
		while (elems != null) {
			allSamples.add(elems[1]);
			elems = in.readLineElems(TextFile.tab);
		}

		in.close();

		in.open();
		elems = in.readLineElems(TextFile.tab);
		HashSet<String> samplesToInclude = new HashSet<String>();
		while (elems != null) {

			String father = elems[2];
			String mother = elems[3];

			if (allSamples.contains(father) && allSamples.contains(mother)) {
				samplesToInclude.add(father);
				samplesToInclude.add(mother);
				samplesToInclude.add(elems[1]);
			}
			elems = in.readLineElems(TextFile.tab);
		}
		in.close();


		TextFile out = new TextFile(samplesToExclude, TextFile.W);
		for (String s : allSamples) {
			if (!samplesToInclude.contains(s)) {
				out.writeln(s);
			}
		}
		out.close();
	}

	public void filterMapForVCFVariants(String mapIn, String vcfIn, String listOfVariantsToExclude) throws IOException {
		TextFile in = new TextFile(vcfIn, TextFile.R);
		String ln = in.readLine();
		HashSet<String> variantsToFilter = new HashSet<String>();
		while (ln != null) {
			if (ln.startsWith("##")) {
// metadata

			} else if (ln.startsWith("#CHROM")) {
// header

			} else {
				String[] elems = ln.split("\t");
				String variant = elems[0] + "_" + elems[1];
				variantsToFilter.add(variant);
			}

			ln = in.readLine();
		}
		in.close();


		TextFile out = new TextFile(listOfVariantsToExclude, TextFile.W);
		TextFile tf = new TextFile(mapIn, TextFile.R);

		String[] elems2 = tf.readLineElems(TextFile.tab);
		while (elems2 != null) {
			String variant = elems2[0] + "_" + elems2[3];
			if (variantsToFilter.contains(variant)) {
				out.writeln(elems2[1]);
			}
			elems2 = tf.readLineElems(TextFile.tab);
		}
		out.close();

	}

	public void replaceFamPhenotypes(String famin, String refFam, String famout) throws IOException {

		HashMap<String, String> sampletoPhenotype = new HashMap<String, String>();
		HashMap<String, String> sampletoFamid = new HashMap<String, String>();
		HashMap<String, String> sampletoFid = new HashMap<String, String>();
		HashMap<String, String> sampletoMid = new HashMap<String, String>();
		HashMap<String, String> sampletoSex = new HashMap<String, String>();

		TextFile tf = new TextFile(refFam, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {
			String sample = elems[1];
			String phenotype = elems[5];
			sampletoFamid.put(sample, elems[0]);
			sampletoFid.put(sample, elems[2]);
			sampletoMid.put(sample, elems[3]);
			sampletoSex.put(sample, elems[4]);


			sampletoPhenotype.put(sample, phenotype);
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();

		TextFile tf2 = new TextFile(famin, TextFile.R);
		TextFile tfout = new TextFile(famout, TextFile.W);
		elems = tf2.readLineElems(Strings.whitespace);
		while (elems != null) {
			String sample = elems[1];
			String phenotype = sampletoPhenotype.get(sample);
			if (phenotype == null) {
				System.err.println("ERROR: " + sample + " not found in reference.");
			}

			elems[2] = sampletoFid.get(sample);
			elems[3] = sampletoMid.get(sample);
			elems[4] = sampletoSex.get(sample);
			elems[5] = phenotype;

			tfout.writeln(Strings.concat(elems, Strings.tab));


			elems = tf2.readLineElems(Strings.whitespace);
		}
		tfout.close();
		tf2.close();

	}
}
