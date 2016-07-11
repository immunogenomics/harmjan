package nl.harmjanwestra.vcfutils;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by Harm-Jan on 02/10/16.
 */
public class VCFVariantRSNameUpdater {

	public void updateRSNames(String dbsnpvcf, String input) throws IOException {


		String[] filenames = input.split(",");
		System.out.println(filenames.length + " files to process");
		HashSet<String> mapOrig = new HashSet<String>();
		for (String vcfin : filenames) {
			System.out.println("Parsing: " + vcfin);
			TextFile tf = new TextFile(vcfin, TextFile.R);

			String ln = tf.readLine();
			int lnctr = 0;
			while (ln != null) {
				if (!ln.startsWith("#")) {
					String substr = ln.substring(0, 200);
					String[] elems = substr.split("\t");
					mapOrig.add(elems[0] + "_" + elems[1]);
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
				String query = lineElems[0] + "_" + lineElems[1];
				String rs = lineElems[2];
				if (mapOrig.contains(query)) {
					mapDbSNP.put(query, rs);
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
					String query = elems[0] + "_" + elems[1];
					String rs = mapDbSNP.get(query);
					if (rs == null) {
						rs = elems[0] + ":" + elems[1] + ":" + elems[2];
						numnull++;
					} else {
						// System.out.println("Replacing " + elems[1] + " with rs: " + rs);
						if (rs.equals(elems[2])) {
							nrEqual++;
						}
						nrReplaced++;
					}
					elems[2] = rs;
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
}
