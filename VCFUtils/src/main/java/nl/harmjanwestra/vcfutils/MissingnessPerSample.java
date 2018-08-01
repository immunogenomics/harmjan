package nl.harmjanwestra.vcfutils;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 6/30/17.
 */
public class MissingnessPerSample {

	public void run(String vcf, String out) throws IOException {

		VCFGenotypeData d = new VCFGenotypeData(vcf);
		ArrayList<String> samples = d.getSamples();

		TextFile tf = new TextFile(vcf, TextFile.R);
		String ln = tf.readLine();
		int[] nrCalled = new int[samples.size()];
		int nrVars = 0;
		while (ln != null) {
			if (!ln.startsWith("#")) {
				VCFVariant v = new VCFVariant(ln);
				DoubleMatrix2D d2d = v.getGenotypeAllelesAsMatrix2D();
				for (int i = 0; i < samples.size(); i++) {
					if (d2d.get(i, 0) > -1) {
						nrCalled[i]++;
					}
				}
				nrVars++;
			}
			ln = tf.readLine();
		}

		tf.close();

		TextFile outf = new TextFile(out, TextFile.W);
		for (int i = 0; i < samples.size(); i++) {
			outf.writeln(samples.get(i) + "\t" + nrCalled[i] + "\t" + ((double) nrCalled[i] / nrVars));
		}
		outf.close();

	}
}
