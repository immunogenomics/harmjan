package nl.harmjanwestra.utilities.vcf;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by Harm-Jan on 06/02/16.
 */
public class VCFTest {

	public static void main(String[] args) {

		String file = "D:\\tmp\\2016-06-02\\problemvars\\rs35149334.txt";
		try {
			TextFile tf = new TextFile(file, TextFile.R);

			String ln = tf.readLine();
			while (ln != null) {
				if (!ln.startsWith("#")) {
					VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.ALL);
					VCFImputationQualScore q = new VCFImputationQualScore(variant, true);
					System.out.println(variant.getId() + "\t" + q.allelicR2() + "\t" + q.doseR2());
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
