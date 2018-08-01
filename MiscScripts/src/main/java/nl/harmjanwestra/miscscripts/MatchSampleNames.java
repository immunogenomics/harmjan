package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by Harm-Jan on 02/15/16.
 */
public class MatchSampleNames {

	public static void main(String[] args) {
		MatchSampleNames n = new MatchSampleNames();


		n.method2();

	}

	public void method1() {
		String f1 = "/Data/tmp/nixup/samples/ic.txt";
		String f2 = "/Data/tmp/nixup/samples/seq.txt";
		String out = "/Data/tmp/nixup/samples/mapout.txt";

		HashSet<String> samples1 = new HashSet<String>();
		try {
			TextFile tf1 = new TextFile(f1, TextFile.R);
			String ln = tf1.readLine();
			while (ln != null) {
				samples1.add(ln);
				ln = tf1.readLine();
			}
			tf1.close();


			TextFile tf2 = new TextFile(f2, TextFile.R);

			ln = tf2.readLine();
			TextFile tfout = new TextFile(out, TextFile.W);
			while (ln != null) {
				if (samples1.contains(ln)) {
					System.out.println(ln + "\t" + ln);
					tfout.writeln(ln + "\t" + ln);
				} else {
					String[] elems = ln.split("/");
					if (elems.length == 2) {
						String sample = elems[1];
						if (samples1.contains(sample)) {
							System.out.println(ln + "\t" + sample);
							tfout.writeln(ln + "\t" + sample);
						}
					}
				}
				ln = tf2.readLine();
			}
			tfout.close();
			tf2.close();
		} catch (IOException e) {

		}
	}

	public void method2() {
		String f1 = "/Sync/OneDrive/RASamples.txt";
		String f2 = "/Data/tmp/nixup/samples/seq.txt";
		String out = "/Data/tmp/nixup/samples/mapout-RA.txt";

		HashMap<String, String> samples1 = new HashMap<String, String>();
		try {
			TextFile tf1 = new TextFile(f1, TextFile.R);
			String[] elems = tf1.readLineElems(TextFile.tab);
			while (elems != null) {
				samples1.put(elems[0], elems[1]);
				elems = tf1.readLineElems(TextFile.tab);
			}
			tf1.close();

			TextFile tf2 = new TextFile(f2, TextFile.R);
			String ln = tf2.readLine();
			TextFile tfout = new TextFile(out, TextFile.W);
			while (ln != null) {
				if (samples1.containsKey(ln)) {
					System.out.println(ln + "\t" + samples1.get(ln));
					tfout.writeln(ln + "\t" + ln);
				} else {
					elems = ln.split("/");
					if (elems.length == 2) {
						String sample = elems[1];
						if (samples1.containsKey(sample)) {
							System.out.println(ln + "\t" + samples1.get(sample));
							tfout.writeln(ln + "\t" + samples1.get(sample));
						}
					}
				}
				ln = tf2.readLine();
			}
			tfout.close();
			tf2.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
