package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashSet;

/**
 * Created by hwestra on 5/4/16.
 */
public class FilterSet {


	public static void main(String[] args) {
		String in1 = "/Data/tmp/2016-05-04/ICVariants/T1D-AllICVariants.txt";
		String in2 = "/Data/tmp/2016-05-04/ICVariants/T1D-Cosmo-merged-variants.txt";
		String out = "/Data/tmp/2016-05-04/ICVariants/AllVsCosmoComparison";

		FilterSet f = new FilterSet();
		try {
			f.filter(in1,in2,out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void filter(String in1, String in2, String out) throws IOException {

		HashSet<String> variants1 = getVariants(in1);
		HashSet<String> variants2 = getVariants(in2);

		TextFile outf1 = new TextFile(out + "-shared.txt", TextFile.W);
		TextFile outf2 = new TextFile(out + "-uniqueIn1.txt", TextFile.W);
		TextFile outf3 = new TextFile(out + "-uniqueIn2.txt", TextFile.W);

		for (String s : variants1) {
			if (variants2.contains(s)) {
				outf1.writeln(s);
			} else {
				outf2.writeln(s);
			}
		}

		for (String s : variants2) {
			if (!variants1.contains(s)) {
				outf3.writeln(s);
			}
		}

		outf1.close();
		outf2.close();
		outf3.close();

	}

	public HashSet<String> getVariants(String in1) throws IOException {
		HashSet<String> output = new HashSet<>();

		TextFile tf = new TextFile(in1, TextFile.R);

		String ln = tf.readLine();
		while (ln != null) {

			String[] elems = ln.split("\t");
			if (elems.length == 3) {
				output.add(ln);
			}

			ln = tf.readLine();
		}

		tf.close();

		return output;
	}

}
