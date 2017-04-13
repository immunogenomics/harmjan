package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 3/23/17.
 */
public class VCFRemoveSamplesWithMissingGT {

	public static void main(String[] args) {
		String file = "/Data/tmp/sh2b3fix/test.vcf";
		String fileout = "/Data/tmp/sh2b3fix/test-filter.vcf";
		try {
			TextFile tf = new TextFile(file, TextFile.R);
			TextFile tf2 = new TextFile(fileout, TextFile.W);
			String ln = tf.readLine();
			String[] header = null;
			boolean[] includesamples = null;
			ArrayList<String> output = new ArrayList<String>();
			while (ln != null) {
				if (ln.startsWith("##")) {
					tf2.writeln(ln);
				} else if (ln.startsWith("#")) {
					header = ln.split("\t");
					includesamples = new boolean[header.length];
					for (int i = 0; i < header.length; i++) {
						includesamples[i] = true;
					}
				} else {
					String[] elems = ln.split("\t");
					for (int i = 9; i < elems.length; i++) {
						if (elems[i].startsWith("./.")) {
							includesamples[i] = false;
						}
					}
					output.add(ln);
				}
				ln = tf.readLine();
			}

			String outln = header[0];
			System.out.println(header.length + " header elements");
			int wr = 1;
			for (int q = 1; q < header.length; q++) {
				if (includesamples[q]) {
					outln += "\t" + header[q];
					wr++;
				}
			}
			tf2.writeln(outln);
			System.out.println(wr + " written");

			for (String l : output) {
				String[] elems = l.split("\t");
				outln = elems[0] + "";
				int wr2 = 1;
				for (int i = 1; i < elems.length; i++) {
					if (includesamples[i]) {
						outln += "\t" + elems[i];
						wr2++;
					}
				}
				System.out.println(wr2 + " written");
				tf2.writeln(outln);
			}
			tf2.close();
			tf.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}


}
