package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by Harm-Jan on 01/12/16.
 */
public class FrequencyFilter {

	public void filter(String in, String out) throws IOException {

		TextFile tf1 = new TextFile(in, TextFile.R);
		TextFile tf2 = new TextFile(out, TextFile.W);
		String ln = tf1.readLine();
		while (ln != null) {
			if (ln.startsWith("#")) {
				tf2.writeln(ln);
			} else {
				VCFVariant var = new VCFVariant(ln);
				if (var.getCallrate() > 0 && var.getMAF() > 0) {
					tf2.writeln(ln);
				}
			}
			ln = tf1.readLine();
		}
		tf1.close();
		tf2.close();


	}
}
