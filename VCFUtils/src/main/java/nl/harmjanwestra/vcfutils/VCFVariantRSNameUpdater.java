package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
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
					String substr = ln.substring(0, 200);
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
		String[] lineElems = vcf.readLineElems(TextFile.tab);
		int ln = 0;
		HashMap<String, String> mapDbSNP = new HashMap<String, String>();
		while (lineElems != null) {

			if (lineElems[0].startsWith("#")) {
				// header
			} else {
				// elems[0] + "_" + elems[1] + "_" + elems[3] + "_" + elems[4]

				String query = Chromosome.parseChr(lineElems[0]).toString() + "_" + lineElems[1];

				if (mapOrigPos.contains(query)) {
					String rspos = getRSPos(lineElems[7]);

					if (rspos == null) {
						System.out.println("Missing RSPOS: " + query);
						rspos = lineElems[1];
					}

					query = Chromosome.parseChr(lineElems[0]).toString() + "_" + rspos;
					String rs = lineElems[2];
					String snpId = query + "_" + lineElems[3] + "_" + lineElems[4] + "_" + rs;
					if (mapDbSNP.containsKey(query)) {
						String availableValues = mapDbSNP.get(query);
						availableValues += ";" + snpId;
						mapDbSNP.put(query, availableValues);
					} else {
						mapDbSNP.put(query, snpId);
					}

				}
			}

			ln++;

			if (ln % 2500000 == 0) {
				System.out.print(ln + " positions parsed...\r");
			}
			lineElems = vcf.readLineElems(TextFile.tab);
		}
		vcf.close();
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
			int lnctr = 0;
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

						String[] replacementelems = replacement.split(";");
						if (replacementelems.length == 1) {
							String[] replacementelems2 = replacement.split("_");
							replacementRS = replacementelems2[replacementelems2.length - 1];
						} else {
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

							for (String rep : replacementelems) {
								String[] replacementelems2 = rep.split("_");
								String chr = replacementelems2[0];
								String pos = replacementelems2[1];
								String allele1 = replacementelems2[2];
								String allele2 = replacementelems2[3];
								String rsid = replacementelems2[4];


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
								set.put(rep, equalCtr);
								setAlleles.put(rep, Strings.concat(dbSNPAlleles, Strings.comma));
							}

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
									String[] replacementelems2 = mostlikelymatch.split("_");
									replacementRS = replacementelems2[4];
									System.out.println(elems[0] + ":" + elems[1] + ":" + elems[2] + ":" + Strings.concat(allAlleles, Strings.comma) + "\t-->\t" + replacementRS + ":" + setAlleles.get(mostlikelymatch) + "\t" + max);
									for (String s : keys) {
										if (!mostlikelymatch.equals(s)) {
											System.out.println(s);
										}
									}
									System.out.println();


								} else {
									replacementRS = elems[0] + ":" + elems[1] + ":" + elems[2];
									System.out.println("Could not find match:\t" + elems[0] + ":" + elems[1] + ":" + elems[2] + " --> " + replacementRS);
								}
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
