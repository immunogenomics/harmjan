package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by Harm-Jan on 05/31/17.
 */
public class AFDist {


	public static void main(String[] args) {

		String file = "D:\\tmp\\HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz";
		String outfile = "D:\\tmp\\dist.txt";
		AFDist d = new AFDist();
		try {
			d.run(file, outfile);
		} catch (IOException e) {
			e.printStackTrace();
		}


	}

	public void run(String file, String outfile) throws IOException {
		int nrbins = 100;

		// place most bins < 0.1

		int[] freq = new int[nrbins + 1];

		TextFile tf = new TextFile(file, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		int total = 0;
		double minfreq = 1;
		while (elems != null) {
			String af = elems[7];
			Double d = Double.parseDouble(af);

			if (d < minfreq) {
				minfreq = d;
			}
			int binno = 0;
			if (d == 1d) {
				binno = 0;
			} else if (d == 0d) {
				binno = nrbins;
			} else {
				d = -Math.log10(d) * 10;
				binno = (int) Math.floor(d);
			}

			if (binno < 0) {
				binno = 0;
			}
			if (binno > nrbins) {
				binno = nrbins;
			}

			freq[binno]++;
			total++;
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		System.out.println(minfreq + " min freq");
		TextFile out = new TextFile(outfile, TextFile.W);
		for (int i = freq.length - 1; i > -1; i--) {
			double bin = Math.pow(10, -((double) i / 10));
			String ln = "" + i + "\t" + bin + "\t" + freq[i] + "\t" + ((double) freq[i] / total);
			out.writeln(ln);
		}
		out.close();

	}
}
