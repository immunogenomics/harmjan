package nl.harmjanwestra.vcfutils;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

/**
 * Created by Harm-Jan on 02/15/16.
 */
public class VCFSampleNameReplace {

	public void replace(String replacefile, String vcfin, String vcfout) throws IOException {

		HashMap<String, String> replacements = new HashMap<String, String>();
		TextFile tf = new TextFile(replacefile, TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			replacements.put(elems[0], elems[1]);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		System.out.println(replacements.size() + " replacements loaded from " + replacefile);

		TextFile tf2 = new TextFile(vcfin, TextFile.R);
		TextFile tf3 = new TextFile(vcfout, TextFile.W);

		int found = 0;
		String ln = tf2.readLine();
		while (ln != null) {
			if (ln.startsWith("#CHROM")) {
				String[] elemsStr = ln.split("\t");
				for (int i = 9; i < elemsStr.length; i++) {
					String replacement = replacements.get(elemsStr[i]);
					if (replacement != null) {
						elemsStr[i] = replacement;
						found++;
					} else {
						System.out.println("not found: "+elemsStr[i]);
					}
				}
				System.out.println(found + " replacements made out of " + (elemsStr.length - 9));
				tf3.writeln(Strings.concat(elemsStr, Strings.tab));
			} else {
				tf3.writeln(ln);
			}
			ln = tf2.readLine();
		}

		tf3.close();

		tf2.close();


	}
}
