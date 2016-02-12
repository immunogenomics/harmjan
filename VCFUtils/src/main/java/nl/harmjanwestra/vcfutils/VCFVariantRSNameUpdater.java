package nl.harmjanwestra.vcfutils;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

/**
 * Created by Harm-Jan on 02/10/16.
 */
public class VCFVariantRSNameUpdater {

	public void updateRSNames(String dbsnpvcf, String vcfin, String vcfout) throws IOException {
		System.out.println("Updating RS names: " + vcfin + " to " + vcfout);
		HashMap<String, String> map = new HashMap<String, String>();
		TextFile tf = new TextFile(vcfin, TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (!elems[0].startsWith("#")) {
				map.put(elems[0] + "_" + elems[1], elems[1]);
			}

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

		TextFile out = new TextFile(vcfout, TextFile.W);
		tf.open();
		elems = tf.readLineElems(TextFile.tab);
		int nrReplaced = 0;
		while (elems != null) {
			if (!elems[0].startsWith("#")) {
				String query = elems[0] + "_" + elems[1];
				String rs = map.get(query);
				if (!rs.equals(elems[2])) {
					// System.out.println("Replacing " + elems[1] + " with rs: " + rs);
					nrReplaced++;
				}
				elems[2] = rs;
			}
			out.writeln(Strings.concat(elems, Strings.tab));
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println(nrReplaced + " totally replaced...");
		out.close();
	}
}
