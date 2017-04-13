package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

/**
 * Created by Harm-Jan on 02/15/16.
 */
public class VCFSampleNameReplace {

	public void replacePlinkDups(String famfile, String vcfIn, String vcfOut) throws IOException {

		HashMap<String, String> sampleToSample = new HashMap<String, String>();
		TextFile fam = new TextFile(famfile, TextFile.R);
		String[] elems = fam.readLineElems(Strings.whitespace);
		while (elems != null) {
			sampleToSample.put(elems[0] + "_" + elems[1], elems[1]);
			elems = fam.readLineElems(Strings.whitespace);
		}
		fam.close();

		TextFile in = new TextFile(vcfIn, TextFile.R);
		TextFile out = new TextFile(vcfOut, TextFile.W);

		String ln = in.readLine();
		while (ln != null) {
			if (ln.startsWith("#CHROM")) {
				elems = ln.split("\t");
				int replaced = 0;
				for (int i = 9; i < elems.length; i++) {
					String sample = elems[i];
					String replacement = sampleToSample.get(sample);
					if (replacement != null) {
						elems[i] = replacement;
						replaced++;
					}
				}
				out.writeln(Strings.concat(elems, Strings.tab));
				System.out.println(replaced + " sample names replaced out of " + (elems.length - 9));
			} else {
				out.writeln(ln);
			}
			ln = in.readLine();
		}
		out.close();
		in.close();

	}

	public void replace(String replacefile, String vcfin, String vcfout) throws IOException {

		HashMap<String, String> replacements = new HashMap<String, String>();
		TextFile tf = new TextFile(replacefile, TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length >= 2) {
				replacements.put(elems[0], elems[1]);
			}
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
						System.out.println("not found: " + elemsStr[i]);
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
