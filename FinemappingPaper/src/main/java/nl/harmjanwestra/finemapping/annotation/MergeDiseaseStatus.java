package nl.harmjanwestra.finemapping.annotation;

import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 7/18/17.
 */
public class MergeDiseaseStatus {

	public static void main(String[] args) {
		try {
			String in1 = "/Data/tmp/2017-07-18/2016-03-11-T1D-diseaseStatusWithPseudos.txt";
			String in2 = "/Data/tmp/2017-07-18/covarmerged.txtmergeddisease.txt";
			String out = "/Data/tmp/2017-07-18/diseaseStatusMulti.txt";

			run(in1, in2, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void run(String in1, String in2, String out) throws IOException {

		HashMap<String, DiseaseStatus> status1 = read(in1);
		HashMap<String, DiseaseStatus> status2 = read(in2);
		HashSet<String> sharedSamples = new HashSet<String>();
		HashSet<String> allSamples = new HashSet<String>();
		allSamples.addAll(status1.keySet());
		allSamples.addAll(status2.keySet());
		System.out.println(allSamples.size() + " unique samples");

		for (String sample : status1.keySet()) {
			if (status2.containsKey(sample)) {
				DiseaseStatus s1 = status1.get(sample);
				DiseaseStatus s2 = status2.get(sample);
				if (!s1.equals(s2)) {
					System.out.println("WARNING: different disease status for sample: " + sample + "\t" + status1 + "\t" + status2);
				}
				sharedSamples.add(sample);

			}

		}
		System.out.println(sharedSamples.size() + " shared samples");

		TextFile outf = new TextFile(out, TextFile.W);
		for (String sample : allSamples) {


			DiseaseStatus s1 = status1.get(sample);
			if (s1 == null) {
				s1 = DiseaseStatus.CONTROL;
			}
			DiseaseStatus s2 = status2.get(sample);
			if (s2 == null) {
				s2 = DiseaseStatus.CONTROL;
			}

			String ln = sample + "\t" + s1.toString() + "\t" + s2.toString();
			outf.writeln(ln);
		}

		outf.close();

	}

	public static HashMap<String, DiseaseStatus> read(String in) throws IOException {
		TextFile tf = new TextFile(in, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, DiseaseStatus> diseaseStatus1 = new HashMap<String, DiseaseStatus>();
		while (elems != null) {
			String sample = elems[0];
			DiseaseStatus s = DiseaseStatus.parseStatus(elems[1]);
			diseaseStatus1.put(sample, s);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return diseaseStatus1;
	}
}
