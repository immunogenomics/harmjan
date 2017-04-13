import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by hwestra on 11/26/15.
 */
public class LiftOverTools {


	public static void main(String[] args) {
		if (args.length >= 1) {
			if (args[0].equals("post")) {


				if (args.length > 6) {
					String regionFile = args[1];
					String dsName = args[2];
					String hg18map = args[3];
					String liftedBed = args[4];
					String dbsnpvcf = args[5];
					String workdir = args[6];
					LiftOverTools t = new LiftOverTools();
					try {
						t.post(regionFile, dsName, hg18map, liftedBed, dbsnpvcf, workdir);
					} catch (IOException e) {
						e.printStackTrace();
					}
				} else {
					System.out.println("Usage post regionfile dsname hg18map liftedbed dbsnpvcf workdir");
				}

			} else if (args[0].equals("pre")) {
				if (args.length > 4) {
					String plinkDataset = args[1];
					String dsName = args[2];
					String refmap = args[3];
					String workdir = args[4];
					LiftOverTools t = new LiftOverTools();
					try {
						t.pre(plinkDataset, dsName, refmap, workdir);
					} catch (IOException e) {
						e.printStackTrace();
					}
				} else {
					System.out.println("Usage pre plinkdataset dsname refmap workdir");
				}
			} else {
				System.out.println("Usage: pre or post");
			}
		} else {
			System.out.println("Usage: pre or post");
		}
	}

	public void pre(String plinkDataset, String dsName, String refmap, String workdir) throws IOException {

		String map = plinkDataset + ".map";
		// String refmap = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC-RsIDsFromT1DStudy.map";
		String chrupdate = workdir + dsName + "-rewrittenChrNames.map";
		rewriteMapFileChromosomeNames(refmap, map, chrupdate);

//		String rsnameref = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC-RsIDsFromT1DStudy.map";
		String rsupdate = plinkDataset + "-rewrittenChrNames-updatedRS.map";
		updateMapFileRsIdsUsingMapFile(refmap, chrupdate, rsupdate);

		String rsupdatededup = workdir + dsName + "-rewrittenChrNames-updatedRS-dedup.map";
		deduplicateMAP(rsupdate, rsupdatededup);

		String bedout = workdir + dsName + "-rewrittenChrNames-dedup.bed";
		rewriteMapToBed(rsupdatededup, bedout);

	}

	public void post(String regionFile, String dsName, String hg18map, String liftedBed, String dbsnpvcf, String workdir) throws IOException {

		String hg19map = workdir + dsName + "-hg19.map";
		convertPostLiftOverMAP(hg18map, hg19map, liftedBed);

		String hg19mapupd = workdir + dsName + "-hg19-updRS.map";
		updateRSNames(dbsnpvcf, hg19map, hg19mapupd);

		// dedup
		String hg19mapdedup = workdir + dsName + "-hg19-updRS-dedup.map";
		deduplicateMAP(hg19mapupd, hg19mapdedup);

		String variantSelect = workdir + dsName + "-hg19-updRS-dedup-selectVariants.txt";
		filterMap(hg19mapdedup, regionFile, variantSelect);

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

	public void updateRSNames(String dbsnpvcf, String mapin, String mapout) throws IOException {
		System.out.println("Updating RS names: " + mapin + " to " + mapout);
		HashMap<String, String> map = new HashMap<String, String>();
		TextFile tf = new TextFile(mapin, TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			map.put(elems[0] + "_" + elems[3], elems[1]);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();


		TextFile vcf = new TextFile(dbsnpvcf, TextFile.R);

		System.out.println("parsing DBSNP VCF: " + dbsnpvcf);
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

			if (ln % 2500000 == 0) {
				System.out.println(ln + " positions parsed...");
			}
			lineElems = vcf.readLineElems(TextFile.tab);
		}
		vcf.close();
		System.out.println(map.size() + " variant annotations readAsTrack");

		TextFile out = new TextFile(mapout, TextFile.W);
		tf.open();
		elems = tf.readLineElems(TextFile.tab);
		int nrReplaced = 0;
		while (elems != null) {
			String query = elems[0] + "_" + elems[3];
			String rs = map.get(query);
			if (!rs.equals(elems[1])) {
				// System.out.println("Replacing " + elems[1] + " with rs: " + rs);
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

	public void rewriteMapToBed(String mapIn, String bedOut) throws IOException {
		TextFile tf = new TextFile(mapIn, TextFile.R);
		TextFile out = new TextFile(bedOut, TextFile.W);
		System.out.println("converting " + mapIn + " to " + bedOut);
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

	public void deduplicateMAP(String map1, String map2) throws IOException {

		System.out.println("Dedupping: " + map1 + " to " + map2);

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


	public void updateMapFileRsIdsUsingMapFile(String refMap, String inMap, String outMap) throws IOException {

		TextFile referenceMapFile = new TextFile(refMap, TextFile.R);
		Set<String> refRS = referenceMapFile.readAsSet(1, TextFile.tab);
		referenceMapFile.close();

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

		referenceMapFile.open();
		HashMap<String, HashSet<String>> chrToRS = new HashMap<String, HashSet<String>>();


		String[] elems = referenceMapFile.readLineElems(TextFile.tab);
		while (elems != null) {
			String query = elems[0] + "_" + elems[3];
			HashSet<String> rsssssss = chrToRS.get(query);
			if (rsssssss == null) {
				rsssssss = new HashSet<String>();
			}
			rsssssss.add(elems[1]);

			chrToRS.put(query, rsssssss);
			elems = referenceMapFile.readLineElems(TextFile.tab);
		}
		referenceMapFile.close();


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

	public void rewriteMapFileChromosomeNames(String refMap, String inMap, String outMap) throws IOException {

		TextFile tf = new TextFile(refMap, TextFile.R);
		Set<String> refRS = tf.readAsSet(1, TextFile.tab);
		tf.close();


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

		System.out.println("Replacing chromosome names");
		System.out.println("Loaded " + rsToChr.size() + " snps from " + refMap);

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

		System.out.println(inRS.size() + " SNPs in " + inMap + "\n" + counter + " not in ref");
		tf2.open();
		elems = tf2.readLineElems(TextFile.tab);

		int counternotupdated = 0;
		int counterupdated = 0;

		int diffpos = 0;

		while (elems != null) {
			String snp = elems[1];
			String chr = rsToChr.get(snp);
			if (chr == null) {
				counternotupdated++;
			} else {
				counterupdated++;
				String pos = rsToChrPos.get(snp);
				elems[0] = chr;

				if (!pos.equals(elems[3])) {
					diffpos++;
				}
				elems[3] = pos;
			}

			outf.writeln(Strings.concat(elems, Strings.tab));

			elems = tf2.readLineElems(TextFile.tab);
		}
		System.out.println(counterupdated + " updated / " + counternotupdated + " not updated");
		System.out.println(diffpos + " had actual different position");

		tf2.close();
		outf.close();


	}

}
