package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;

/**
 * Created by hwestra on 6/4/16.
 */
public class StripInfo {

	public void strip(String file) throws IOException {

		TextFile in = new TextFile(file, TextFile.R);
		TextFile out = new TextFile(file + "_tmp.vcf.gz", TextFile.W);
		String ln = in.readLine();
		while (ln != null) {
			if (ln.startsWith("#")) {
				out.writeln(ln);
			} else {
				String[] elems = ln.split("\t");
				elems[7] = ".";
				out.writeln(Strings.concat(elems, Strings.tab));
			}
			ln = in.readLine();
		}
		in.close();
		out.close();

		Gpio.moveFile(out.getFileName(), in.getFileName());

	}
}
