package nl.harmjanwestra.miscscripts;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class RsMergeRewriter {

	public static void main(String[] args) {
		try {
			String rsMergeFile = "/Data/dbSNP/b142/RsMergeArch.bcp";
			String rsRemoveFile = "";
			String mapfileIn = "/Data/Projects/2014-FR-Reseq/ImmunoChip/eur_T1D_hg19_lifted.map";
			String mapfileOut = "/Data/Projects/2014-FR-Reseq/ImmunoChip/eur_T1D_hg19_lifted_newRSIds.map";
			int dbsnpB = 138;

// do the actual stuffffff.....
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
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}