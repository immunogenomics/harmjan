package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by hwestra on 2/18/15.
 */
public class VCFFilter {

	public static void main(String[] args) {
		try {
			String file = "/Data/Projects/2014-FR-Reseq/2015-01-31-HaplotypeCallerVariants-762Samples/merged.vcf";
			String querySample = "t1d-run3a/40458101";
			String queryVariant = "6:167368741";


			TextFile tf = new TextFile(file, TextFile.R);

			String ln = tf.readLine();
			int sampleColumn = 0;
			while (ln != null) {

				if (ln.startsWith("##")) {

				} else if (ln.startsWith("#")) {
					String[] samples = ln.split("\t");
					for (int i = 0; i < samples.length; i++) {
						if (samples[i].equals(querySample)) {
							sampleColumn = i;
							break;
						}
					}
				} else {
					// GT:AD:DP:GQ:PL
					String[] elems = ln.split("\t");
					String variant = elems[0] + ":" + elems[1];
					if (variant.equals(queryVariant)) {
						System.out.println(variant + "\t" + querySample + "\t" + elems[3] + "\t" + elems[4] + "\t" + elems[sampleColumn]);
					}

				}

				ln = tf.readLine();


			}

			tf.close();

		} catch (IOException e) {

		}
	}
}
