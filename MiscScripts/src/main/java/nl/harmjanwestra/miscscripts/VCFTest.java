package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by hwestra on 6/13/16.
 */
public class VCFTest {

	public static void main(String[] args) {
		String file = "/Data/tmp/2016-06-10/RA-Beagle1kg-regionfiltered-COSMO-ImpQualsReplaced-chr22.vcf.gz";

		try {
			TextFile tf = new TextFile(file, TextFile.R);
			String ln = tf.readLine();
			while (ln != null) {

				if (!ln.startsWith("#")) {
					VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.ALL);
				}
				ln = tf.readLine();
			}
			tf.close();
		} catch (IOException e) {

		}

	}

}
