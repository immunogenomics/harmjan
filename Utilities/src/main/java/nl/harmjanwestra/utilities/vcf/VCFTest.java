package nl.harmjanwestra.utilities.vcf;

import org.apache.tools.ant.util.StringUtils;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.util.RunTimer;

import java.io.IOException;
import java.util.Vector;

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
			String f = "/Data/tmp/2016-06-10/RA-Beagle1kg-regionfiltered-COSMO-ImpQualsReplaced-chr22.vcf.gz";
			RunTimer t = new RunTimer();
			t.start();

			for (int i = 0; i < 3; i++) {
				TextFile tf = new TextFile(f, TextFile.R);
				String ln = tf.readLine();
				long sum = 0;
				while (ln != null) {
					if (!ln.startsWith("#")) {
						String[] elems = ln.split("\t");
						for (int j = 0; j < elems.length; j++) {
							sum++;
						}
					}
					ln = tf.readLine();
				}
				System.out.println(i + "/" + 1000 + "\t1\t" + sum);
			}
			String time1 = t.stop();

			t.start();
			for (int i = 0; i < 3; i++) {
				TextFile tf = new TextFile(f, TextFile.R);
				String ln = tf.readLine();
				long sum = 0;
				while (ln != null) {
					if (!ln.startsWith("#")) {
						Vector numbers = StringUtils.split(ln, '\t');
						for (int j = 0; j < numbers.size(); j++) {
							sum++;
						}


					}
					ln = tf.readLine();
				}
				System.out.println(i + "/" + 1000 + "\t2\t" + sum);
			}
			String time2 = t.stop();

			System.out.println(time1);
			System.out.println(time2);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}


}
