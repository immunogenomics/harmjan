package nl.harmjanwestra.finemapping.annotation;

import nl.harmjanwestra.utilities.enums.Chromosome;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by hwestra on 9/25/16.
 */
public class eQTLFilter {

	public static void main(String[] args) {
		String eqtl = "/Data/eQTLs/BiosEQTLs/eQTLsFDR0.05-ProbeLevel.txt.gz";
		Chromosome chr = Chromosome.TWO;
		int start = 204446380;
		int stop = 204816382;

		// Chr2_204446380-204816382

		String out = "/Data/tmp/2016-09-25/CD28locus-2.txt";

		eQTLFilter f = new eQTLFilter();
		try {
			f.run(eqtl,chr,start,stop,out);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void run(String eqtl, Chromosome chr, int start, int stop, String out) throws IOException {


		TextFile tf = new TextFile(eqtl, TextFile.R);
		TextFile tfout = new TextFile(out, TextFile.W);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 3) {
				Chromosome chre = Chromosome.parseChr(elems[2]);
				Integer snppos = Integer.parseInt(elems[3]);

				if (chre.equals(chr) && snppos >= start && snppos <= stop) {
					tfout.writeln(
							elems[0]
									+ "\t" + (-Math.log10(Double.parseDouble(elems[0])))
									+ "\t" + elems[1]
									+ "\t" + elems[2]
									+ "\t" + elems[3]
									+ "\t" + elems[16] // elems[elems.length-2]
					);
				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();
		tfout.close();


	}
}
