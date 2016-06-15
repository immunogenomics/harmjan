package nl.harmjanwestra.utilities.vcf;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by Harm-Jan on 06/02/16.
 */
public class VCFTest {

	public static void main(String[] args) {

//		String file = "/Data/tmp/2016-06-06/rs139965430.txt";
//		try {
//			TextFile tf = new TextFile(file, TextFile.R);
//
//			String ln = tf.readLine();
//			while (ln != null) {
//				if (!ln.startsWith("#")) {
//					VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.ALL);
//					VCFImputationQualScoreBeagle q = new VCFImputationQualScoreBeagle(variant, true);
//					VCFImputationQualScoreImpute q2 = new VCFImputationQualScoreImpute();
//					q2.computeAutosomal(variant);
//					System.out.println(variant.getId() + "\t" + q.allelicR2() + "\t" + q.doseR2() + "\t" + q2.getImpinfo() + "\t" + q2.getInfo());
//					System.out.println(q.toString());
//					System.out.println(variant.getMAF() + "\t" + variant.getMinorAllele());
//					System.out.println();
//				} else {
//
//				}
//				ln = tf.readLine();
//			}
//			tf.close();

		try {
			String f = "/Data/tmp/2016-06-10/test.vcf";

			TextFile tf = new TextFile(f, TextFile.R);
			String ln = tf.readLine();
			while (ln != null) {
				if (!ln.startsWith("#")) {
					VCFVariant var = new VCFVariant(ln, VCFVariant.PARSE.ALL);
					DoubleMatrix2D d1 = var.getGenotypeDosagesAsMatrix2D();
					DoubleMatrix2D d2 = var.getDosagesAsMatrix2D();

					for (int i = 0; i < d1.rows(); i++) {
						System.out.println(d1.get(i, 0) + "\t" + d2.get(i, 0));
					}

				}
				ln = tf.readLine();
			}
			tf.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}


}
