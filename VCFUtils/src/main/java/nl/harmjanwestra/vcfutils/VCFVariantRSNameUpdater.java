package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.io.trityper.util.BaseAnnot;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by Harm-Jan on 02/10/16.
 */
public class VCFVariantRSNameUpdater {

	public void updateRSNames(String dbsnpvcf, String input) throws IOException {

		String[] filenames = input.split(",");
		if (filenames.length == 1) {
			if (input.contains("CHR")) {
				filenames = new String[22];
				for (int i = 1; i < 23; i++) {
					filenames[i - 1] = input.replaceAll("CHR", "" + i);
				}
			}
		}

		System.out.println(filenames.length + " files to process");
		HashSet<String> mapOrig = new HashSet<String>();
		HashSet<String> mapOrigPos = new HashSet<String>();

		for (String vcfin : filenames) {
			System.out.println("Parsing: " + vcfin);
			TextFile tf = new TextFile(vcfin, TextFile.R);

			String ln = tf.readLine();
			int lnctr = 0;
			while (ln != null) {
				if (!ln.startsWith("#")) {

					String substr = ln;
					if (ln.length() > 200) {
						substr = ln.substring(0, 200);
					}
					String[] elems = substr.split("\t");

					// also add alleles
					mapOrigPos.add(Chromosome.parseChr(elems[0]).toString() + "_" + elems[1]);
					mapOrig.add(Chromosome.parseChr(elems[0]).toString() + "_" + elems[1] + "_" + elems[3] + "_" + elems[4]);
				}
				lnctr++;
				if (lnctr % 10000 == 0) {
					System.out.print(lnctr + " lines parsed\r");
				}

				ln = tf.readLine();
			}
			tf.close();
			System.out.println();
			System.out.println(mapOrig.size() + " variants in list");
		}

		TextFile vcf = new TextFile(dbsnpvcf, TextFile.R);

		System.out.println("parsing DBSNP VCF: " + dbsnpvcf);
		String ln = vcf.readLine();

		int lnctr = 0;
		HashMap<String, String> mapDbSNP = new HashMap<String, String>();
		System.out.println("Parsing reference data..");

		int nrMissingRSPos = 0;
		while (ln != null) {
			if (ln.charAt(0) == '#') {
				// header
			} else {
				if (ln.length() > 500) {
					ln = ln.substring(0, 500);
				}
				// elems[0] + "_" + elems[1] + "_" + elems[3] + "_" + elems[4]
				String[] lineElems = TextFile.tab.split(ln);
				String query = Chromosome.parseChr(lineElems[0]).toString() + "_" + lineElems[1];

				if (mapOrigPos.contains(query)) {
					String rspos = getRSPos(lineElems[7]);

					if (rspos == null) {
//						System.out.println("Missing RSPOS: " + query + " using chrom pos instead.");
						rspos = lineElems[1];
						nrMissingRSPos++;
					}

					query = Chromosome.parseChr(lineElems[0]).toString() + "_" + rspos;
					String[] elems = new String[]{
							lineElems[0],
							rspos,
							lineElems[2],
							lineElems[3],
							lineElems[4],
					};

					String snpId = Strings.concat(elems, Strings.tab);

//					String snpId = query + "_" + lineElems[3] + "_" + lineElems[4] + "_" + rs; // chr pos
					if (mapDbSNP.containsKey(query)) {
						String availableValues = mapDbSNP.get(query);
						availableValues += "#####" + snpId;
						mapDbSNP.put(query, availableValues);
					} else {
						mapDbSNP.put(query, snpId);
					}

				}
			}

			lnctr++;

			if (lnctr % 2500000 == 0) {
				System.out.print(ln + " positions parsed...\r");
			}
			ln = vcf.readLine();
		}
		vcf.close();
		System.out.println(ln + "total lines read.");
		if (nrMissingRSPos > 0) {
			System.out.println(nrMissingRSPos + " with missing RSPOS column. Is this a dbSNP VCF file?");
		}
		System.out.println(mapDbSNP.size() + " variant annotations total");

		for (String vcfin : filenames) {
			String vcfout = vcfin + "-updatedRSId.vcf.gz";

			TextFile tf = new TextFile(vcfin, TextFile.R);

			TextFile out = new TextFile(vcfout, TextFile.W);

			tf.open();
			int nrReplaced = 0;
			int nrEqual = 0;
			int numnull = 0;
			boolean headerwritten = false;
			String line = tf.readLine();
			lnctr = 0;
			int siteswithmultiplersids = 0;
			int siteswithmultiplersidsresolved = 0;
			int normalsites = 0;
			int sitesNotResolved = 0;
			while (line != null) {
				if (!line.startsWith("#")) {
					String[] elems = line.split("\t");
					String query = Chromosome.parseChr(elems[0]).toString() + "_" + elems[1];
					String replacement = mapDbSNP.get(query);
					String replacementRS = null;
					if (replacement == null) {
						replacementRS = elems[0] + ":" + elems[1] + ":" + elems[2];
						numnull++;
					} else {

						String[] replacementelems = replacement.split("#####");
						if (replacementelems.length == 1) {
							String[] replacementelems2 = replacement.split("\t");
							replacementRS = replacementelems2[2];
							normalsites++;
						} else {
							siteswithmultiplersids++;
							// multiple RS ids at this position...
							// try to match alleles
							String refAllele = elems[3];
							String altAllele = elems[4];
							String[] altAlleleElems = altAllele.split(",");
							String[] allAlleles = new String[1 + altAlleleElems.length];
							allAlleles[0] = refAllele;
							for (int q = 0; q < altAlleleElems.length; q++) {
								allAlleles[q + 1] = altAlleleElems[q];
							}

							HashMap<String, Integer> set = new HashMap<String, Integer>();
							HashMap<String, String> setAlleles = new HashMap<String, String>();

							// for each possible combination count the number of alleles that are shared.
							for (String rep : replacementelems) {
								String[] replacementelems2 = rep.split("\t");
								if (replacementelems2.length >= 4) {
									String chr = replacementelems2[0];
									String pos = replacementelems2[1];
									String allele1 = replacementelems2[3];
									String allele2 = replacementelems2[4];
									String rsid = replacementelems2[2];

									// try another combo
									// is the variant an indel? maybe we can only match one allele
									String[] allele2elems = allele2.split(",");
									String[] dbSNPAlleles = new String[1 + allele2elems.length];

									dbSNPAlleles[0] = allele1;
									for (int q = 0; q < allele2elems.length; q++) {
										dbSNPAlleles[q + 1] = allele2elems[q];
									}

									int equalCtr = 0;
									for (String a : allAlleles) {
										for (String b : dbSNPAlleles) {
											if (a.equals(b)) {
												equalCtr++;
											}
										}
									}

									if (equalCtr == 0) {
										System.out.println("Could not find match:\t" + elems[0] + ":" + elems[1] + ":" + elems[2] + "\t" + refAllele + "/" + altAllele + " with " + rep);
										// try reverse complement
										String complementAllele1 = "";
										for (int i = 0; i < allele1.length(); i++) {
											complementAllele1 += "" + BaseAnnot.getComplement("" + allele1.charAt(i));

										}
										System.out.println("Allele1: " + complementAllele1);
										String[] complementAllele2 = new String[allele2elems.length];
										for (int j = 0; j < complementAllele2.length; j++) {
											complementAllele2[j] = "";
											String a = allele2elems[j];
											for (int i = 0; i < a.length(); i++) {
												complementAllele2[j] += "" + BaseAnnot.getComplement("" + a.charAt(i));
											}
											System.out.println("Allele2: " + complementAllele2[j]);
										}


										dbSNPAlleles[0] = complementAllele1;
										for (int q = 0; q < allele2elems.length; q++) {
											dbSNPAlleles[q + 1] = complementAllele2[q];
										}

										equalCtr = 0;
										for (String a : allAlleles) {
											for (String b : dbSNPAlleles) {
												if (a.equals(b)) {
													equalCtr++;
												}
											}
										}
									}
									set.put(rep, equalCtr);
									setAlleles.put(rep, Strings.concat(dbSNPAlleles, Strings.comma));
								} else {
									System.err.println("Was expecting >4 elems for string, separated by _: " + rep);
									System.err.println("Input: " + query);
									System.err.println("Replacement: " + replacement);
									System.err.println("");
								}

							}

							// determine which variant has the highest number of equal alleles.
							if (replacementRS == null) {
								Set<String> keys = set.keySet();
								Integer max = 0;
								String mostlikelymatch = "";
								for (String key : keys) {
									Integer ctr = set.get(key);
									if (ctr > max) {
										max = ctr;
										mostlikelymatch = key;
									}
								}

								if (max > 0) {
									String[] replacementelems2 = mostlikelymatch.split("\t");
									replacementRS = replacementelems2[2];
									System.out.println(elems[0] + ":" + elems[1] + ":" + elems[2] + ":" + Strings.concat(allAlleles, Strings.comma) + "\t" + refAllele + "/" + altAllele + "\t-->\t" + replacementRS + ":" + setAlleles.get(mostlikelymatch) + "\t" + max);
									for (String s : keys) {
										if (!mostlikelymatch.equals(s)) {
											System.out.println(s);
										}
									}
									System.out.println();


								} else {

									replacementRS = elems[0] + ":" + elems[1] + ":" + elems[2];
//									System.out.println(elems[0] + ":" + elems[1] + ":" + elems[2] + ":" + Strings.concat(allAlleles, Strings.comma) + "\t" + refAllele + "/" + altAllele + "\t-->\t" + replacementRS + ":" + setAlleles.get(mostlikelymatch) + "\t" + max);
									System.out.println("Could not find match:\t" + elems[0] + ":" + elems[1] + ":" + elems[2] + "\t" + refAllele + "/" + altAllele + " --> " + replacementRS);
								}
							} else {
								siteswithmultiplersidsresolved++;
							}
						}

						// System.out.println("Replacing " + elems[1] + " with rs: " + rs);
						if (replacement.equals(elems[2])) {
							nrEqual++;
						}
						nrReplaced++;
					}

					if (replacementRS == null) {
						replacementRS = elems[0] + ":" + elems[1] + ":" + elems[2];
						sitesNotResolved++;
					}

					elems[2] = replacementRS;
					out.writeln(Strings.concat(elems, Strings.tab));
				} else {
					out.writeln(line);
					if (!headerwritten) {
						out.writeln("##RSIdsReplaced=" + dbsnpvcf);
						headerwritten = true;
					}
				}

				lnctr++;
				if (lnctr % 10000 == 0) {
					System.out.print(lnctr + " lines parsed\r");
				}
				line = tf.readLine();
			}
			tf.close();
			System.out.println(nrReplaced + " (" + nrEqual + " equal / " + numnull + " null) totally replaced in " + vcfin);
			System.out.println("Sites w/multiple rsids: " + siteswithmultiplersids + "\tMulti-sites resolved: " + siteswithmultiplersidsresolved + "\tSingle rsId sites: " + normalsites + "\tNot resolved: " + sitesNotResolved);
			out.close();
		}

	}

	private String getRSPos(String info) {
		String[] infoElems = info.split(";");

		for (String ielem : infoElems) {
			if (ielem.startsWith("RSPOS")) {
				String[] rsposelems = ielem.split("=");
				if (rsposelems.length == 2) {
					return rsposelems[1];
				}
			}
		}
		return null;
	}
}
