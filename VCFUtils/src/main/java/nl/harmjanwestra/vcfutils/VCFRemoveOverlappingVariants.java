package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashSet;

/**
 * Created by hwestra on 2/17/16.
 */
public class VCFRemoveOverlappingVariants {

	public void remove(String vcf1, String vcf2, String out) throws IOException {
		System.out.println("Reference: " + vcf1);
		System.out.println("To Filter: " + vcf2);
		System.out.println("Out: " + out);

		HashSet<String> variantsInVCF1 = new HashSet<String>();

		TextFile tf1 = new TextFile(vcf1, TextFile.R);
		String ln1 = tf1.readLine();

		while (ln1 != null) {
			if (!ln1.startsWith("#")) {
				VCFVariant var = new VCFVariant(ln1, VCFVariant.PARSE.HEADER);
				variantsInVCF1.add(var.toString());
			}
			ln1 = tf1.readLine();
		}
		tf1.close();

		System.out.println(variantsInVCF1.size() + " variants in reference");

		TextFile tf = new TextFile(vcf2, TextFile.R);
		TextFile tfout = new TextFile(out, TextFile.W);
		String ln = tf.readLine();
		int included = 0;
		int excluded = 0;
		int lnr = 0;
		while (ln != null) {
			if (ln.startsWith("#")) {
				tfout.writeln(ln);
			} else {
				VCFVariant var = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
				if (!variantsInVCF1.contains(var.toString())) {
					tfout.writeln(ln);
					included++;
				} else {
					excluded++;
				}
			}
			ln = tf.readLine();
			lnr++;
			if (lnr % 1000 == 0) {
				System.out.println(lnr + " variants processed");
			}
		}
		System.out.println(included + " included and " + excluded + " excluded");
		tfout.close();
		tf.close();
	}
}
