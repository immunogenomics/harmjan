package nl.harmjanwestra.miscscripts;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;

/**
 * Created by hwestra on 1/13/15.
 */
public class SplitTargetRegionsPerChr {
	public static void main(String[] args) {

		// not 8
		// not 3
		try {
			String inf = "/Data/Projects/2014-FR-Reseq/TargetRegions/2015-01-07-targetRegions.bed";

			for (int i = 1; i < 23; i++) {
				String outf = "/Data/Projects/2014-FR-Reseq/TargetRegions/perchr/" + i + ".bed";
				filter(inf, outf, "" + i);
			}

			String outf = "/Data/Projects/2014-FR-Reseq/TargetRegions/perchr/X.bed";
			filter(inf, outf, "X");

		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public static void filter(String inf, String outf, String istr) throws IOException {

		TextFile tf = new TextFile(inf, TextFile.R);
		TextFile out = new TextFile(outf, TextFile.W);
		String[] elems = tf.readLineElems(TextFile.tab);

		while (elems != null) {
			if (elems[0].equals(istr)) {
				out.writeln(Strings.concat(elems, Strings.tab));
			}


			elems = tf.readLineElems(TextFile.tab);
		}
		out.close();
		tf.close();
	}
}
