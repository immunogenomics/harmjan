package nl.harmjanwestra.goshifter;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Created by hwestra on 2/26/15.
 */
public class GoShifterResultSummarizer {

	public void run(String[] args) {
		String dir = null;
		String out = null;

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("--in")) {
				if (i + 1 < args.length) {
					dir = args[i + 1];
				}
			} else if (args[i].equals("--out")) {
				if (i + 1 < args.length) {
					out = args[i + 1];
				}
			}
		}

		if (dir != null && out != null) {
			try {
				if (!Gpio.exists(dir)) {
					System.err.println("Error: dir not found: " + dir);
				}
				System.out.println("Summarizing: " + dir + " output to" + out);
				summarize(dir, out);

			} catch (IOException e) {
				e.printStackTrace();
			}
		} else {
			printUsage();
		}

	}

	private void printUsage() {
		System.out.println("Usage: --mode summarize --in /dir/ --out output.txt");
	}

	private void summarize(String dir, String out) throws IOException {
		String[] filesInDir = Gpio.getListOfFiles(dir);
		ArrayList<String> enrichFiles = new ArrayList<String>();
		dir = Gpio.formatAsDirectory(dir);
		for (String s : filesInDir) {
			if (s.endsWith(".enrich")) {
				enrichFiles.add(s);
			}
		}

		System.out.println(enrichFiles + " enrichment files found.");

		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln("File\tPval");
		for (String file : enrichFiles) {
			TextFile tf = new TextFile(dir + file, TextFile.R);
// skip the header
			tf.readLine();

			String[] elems = tf.readLineElems(TextFile.tab);
			double realperc = 1;
			double realnroverlap = 0;
			double realnrloci = 0;
			ArrayList<Double> permVals = new ArrayList<Double>();
			while (elems != null) {
				if (elems.length > 2) {
					Integer d1 = Integer.parseInt(elems[0]);
					Integer d2 = Integer.parseInt(elems[1]);
					Integer d3 = Integer.parseInt(elems[2]);

					if (d1 > 0) {
						permVals.add((double) d2 / d3);
					} else {
						realperc = (double) d2 / d3;
						realnroverlap = d2;
						realnrloci = d3;
					}
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

			Collections.sort(permVals);
			int ct = 0;
			for (double d : permVals) {
				if (d >= realperc) {
					ct++;
				}
			}

			double pval = (double) ct / permVals.size();
			outf.writeln(file + "\t" + realnroverlap + "\t" + realnrloci + "\t" + realperc + "\t" + pval);

		}
		outf.close();
	}

}
