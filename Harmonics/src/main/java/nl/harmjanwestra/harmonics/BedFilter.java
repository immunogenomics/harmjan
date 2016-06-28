package nl.harmjanwestra.harmonics;

import nl.harmjanwestra.utilities.enums.Chromosome;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;

/**
 * Created by hwestra on 9/9/15.
 */
public class BedFilter {


	public static void main(String[] args) {

		try {
			if (args.length < 2) {
				System.out.println("Usage: bedin bedout");
			}
			BedFilter f = new BedFilter();
			f.run(args[0], args[1]);

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String filein, String fileout) throws IOException {

		TextFile tf = new TextFile(filein, TextFile.R);
		TextFile tfout = new TextFile(fileout, TextFile.W);

		System.out.println(filein + "\t" + fileout);

		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {
			Chromosome chr = Chromosome.parseChr(elems[0]);
			// {print $1,$2,$3,$7,$5,$6}
			if (!chr.equals(Chromosome.NA)) {
				elems[3] = elems[6];
				tfout.writeln(Strings.concat(elems, Strings.tab, 0, 6));
			} else {
				//System.out.println("cannot parse chr: " + chr + " in file: " + filein);
			}
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();
		tfout.close();
	}
}
