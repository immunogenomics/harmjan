package nl.harmjanwestra.utilities.vcf;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by Harm-Jan on 06/02/16.
 */
public class VCFTest {

	public static void main(String[] args) {

		String file = "/Data/tmp/2016-06-06/rs139965430.txt";
		try {
			TextFile tf = new TextFile(file, TextFile.R);

			String ln = tf.readLine();
			while (ln != null) {
				if (!ln.startsWith("#")) {
					VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.ALL);
					VCFImputationQualScoreBeagle q = new VCFImputationQualScoreBeagle(variant, true);
					VCFImputationQualScoreImpute q2 = new VCFImputationQualScoreImpute();
					q2.computeAutosomal(variant);
					System.out.println(variant.getId() + "\t" + q.allelicR2() + "\t" + q.doseR2() + "\t" + q2.getImpinfo() + "\t" + q2.getInfo());
					System.out.println(q.toString());
					System.out.println(variant.getMAF() + "\t" + variant.getMinorAllele());
					System.out.println();
				} else {

				}
				ln = tf.readLine();
			}
			tf.close();
		} catch (IOException e) {

		}
	}

}
