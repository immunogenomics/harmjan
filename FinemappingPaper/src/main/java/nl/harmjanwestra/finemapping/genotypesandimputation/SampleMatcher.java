package nl.harmjanwestra.finemapping.genotypesandimputation;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 2/10/17.
 */
public class SampleMatcher {

	public static void main(String[] args) {
		String a = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/SequencingPanel/SeqIdToRAImmunoChip.txt";
		String b = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/SequencingPanel/uniquelySequencedSamples.txt";
		String c = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/SequencingPanel/uniquelySequencedSamplesRASamplesReplaced.txt";
		SampleMatcher m = new SampleMatcher();
		try {
//			m.run(a, b, c);

			String d = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/ImmunoChipSampleLists/sampleListT1D-imputed.txt";
			String e = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/ImmunoChipSampleLists/sampleListRA-imputed.txt";
			m.countSamples(b, d);
			m.countSamples(c, e);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String f1, String f2, String out) throws IOException {
		TextFile tf = new TextFile(f1, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, String> strToStr = new HashMap<String, String>();
		while (elems != null) {
			if (elems.length >= 2) {
				String a = elems[0];
				String b = elems[1];
				strToStr.put(a, b);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile tf2 = new TextFile(f2, TextFile.R);
		TextFile outf = new TextFile(out, TextFile.W);
		elems = tf2.readLineElems(TextFile.tab);
		int replacements = 0;
		while (elems != null) {
			String a = elems[0];

			String b = strToStr.get(a);
			if (b != null) {
				outf.writeln(b);
				System.out.println("Replacing " + a + " with " + b);
				replacements++;
			} else {
				outf.writeln(a);
			}

			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		outf.close();

		System.out.println(replacements + " replacements made.");
	}

	public void countSamples(String a, String b) throws IOException {
		HashSet<String> samplesA = new HashSet<String>();
		HashSet<String> samplesB = new HashSet<String>();

		TextFile tf = new TextFile(a, TextFile.R);
		String ln = tf.readLine();
		while (ln != null) {
			samplesA.add(ln.trim());
			ln = tf.readLine();
		}
		tf.close();
		System.out.println(samplesA.size());

		TextFile tf2 = new TextFile(b, TextFile.R);
		String ln2 = tf2.readLine();
		while (ln2 != null) {
			samplesB.add(ln2.trim());
			ln2 = tf2.readLine();
		}
		tf2.close();

		int shared = 0;
		for (String sample : samplesB) {
			if (samplesA.contains(sample)) {
				shared++;
			}
		}
		System.out.println(shared + " out of " + samplesB.size() + " shared");

		for (String sample : samplesA) {
			if (!samplesB.contains(sample)) {
				System.out.println("Could not find:\t" + sample);
			}
		}




	}

}
