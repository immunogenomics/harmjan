package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by hwestra on 3/2/17.
 */
public class VCFListVariantsAsBed {
	public void run(String input, String output) throws IOException {
		TextFile in = new TextFile(input, TextFile.R);
		TextFile out = new TextFile(output, TextFile.W);
		String ln = in.readLine();
		while (ln != null) {
			if (!ln.startsWith("#")) {
				String[] elems = ln.substring(0, 200).split("\t");
				Integer pos = Integer.parseInt(elems[1]);
				out.writeln(elems[0] + "\t" + elems[1] + "\t" + pos + 1 + "\t" + elems[2]);
			}

			ln = in.readLine();
		}
		out.close();
		in.close();


	}
}
