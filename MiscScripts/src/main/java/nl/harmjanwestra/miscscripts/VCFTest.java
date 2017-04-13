package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 6/13/16.
 */
public class VCFTest {

	public static void main(String[] args) {
		String file = "d:\\tmp\\COSMO-chr22-ImpQualsUpdated.vcf.gz";

		try {
			TextFile tf = new TextFile(file, TextFile.R);
			String ln = tf.readLine();
			int ctr = 0;
			ArrayList<String> id = new ArrayList<String>();
			while (ln != null) {

				if (!ln.startsWith("#")) {
					VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.ALL);
					id.add(variant.getId());
					ctr++;
					if (ctr % 100 == 0) {
						System.out.println(ctr);

					}
				}
				ln = tf.readLine();
				if(ctr == 1000){
					break;
				}
			}
			tf.close();
			System.out.println(ctr);

			VCFGenotypeData d = new VCFGenotypeData(file);
			ctr = 0;
			while (d.hasNext()) {
				VCFVariant v = d.next();
				if (!id.get(ctr).equals(v.getId())) {
					System.out.println("found " + v.getId() + "\texpected " + id.get(ctr));
				}
				ctr++;
				if (ctr % 100 == 0) {
					System.out.println(ctr);

				}
			}

			System.out.println(ctr);
		} catch (IOException e) {

		}

	}

}
