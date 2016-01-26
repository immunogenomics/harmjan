package nl.harmjanwestra.miscscripts;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.StringTokenizer;

/**
 * Created by hwestra on 1/20/16.
 */
public class VCFToMap {

	public static void main(String[] args) {
		try {

			TextFile out = new TextFile("/Sync/Dropbox/PlotsWoICStudy/2016-01-19/4.VariantLists/1kg.map", TextFile.W);

			for (int i = 1; i < 23; i++) {
				String file = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-11-29-AllChr-EurOnly-WithX/5-1KG-SequencedRegions-HighFreq/chr" + i + "-filtered.vcf.gz";

				TextFile in = new TextFile(file, TextFile.R);

				if (Gpio.exists(file)) {
					String ln = in.readLine();

					while (ln != null) {

						if (ln.startsWith("#")) {
							// skip
						} else {
							StringTokenizer tokenizer = new StringTokenizer(ln);
							int ctr = 0;
							String[] header = new String[9];
							while (tokenizer.hasMoreTokens() && ctr < 9) {
								header[ctr] = tokenizer.nextToken();
								ctr++;
							}
							// chr snp 0 pos
							//X	154445759	rs5940453
							out.writeln(header[0] + "\t" + header[2] + "\t" + 0 + "\t" + header[1]);
						}
						ln = in.readLine();
					}
				}

				in.close();

			}

			out.close();

			String file2 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/MixupsFixed/merged-ICIds-MixupsFixed.vcf";
			out = new TextFile("/Sync/Dropbox/PlotsWoICStudy/2016-01-19/4.VariantLists/seq-unfiltered.map", TextFile.W);
			TextFile in = new TextFile(file2, TextFile.R);

			if (Gpio.exists(file2)) {
				String ln = in.readLine();

				while (ln != null) {

					if (ln.startsWith("#")) {
						// skip
					} else {
						StringTokenizer tokenizer = new StringTokenizer(ln);
						int ctr = 0;
						String[] header = new String[9];
						while (tokenizer.hasMoreTokens() && ctr < 9) {
							header[ctr] = tokenizer.nextToken();
							ctr++;
						}
						// chr snp 0 pos
						//X	154445759	rs5940453
						out.writeln(header[0] + "\t" + header[2] + "\t" + 0 + "\t" + header[1]);
					}
					ln = in.readLine();
				}
			}

			in.close();
			out.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
