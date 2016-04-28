package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.features.Chromosome;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashSet;

/**
 * Created by Harm-Jan on 04/27/16.
 */
public class CompareOverlappingVariants {

	public void run(String vcf1, String vcf2list) throws IOException {
		TextFile tf1 = new TextFile(vcf1, TextFile.R);

		String ln = tf1.readLine();
		HashSet<String> varStrs = new HashSet<String>();
		while (ln != null) {
			if (!ln.startsWith("#")) {
				int strlen = ln.length();
				int substrlen = 1000;
				if (strlen < substrlen) {
					substrlen = strlen;
				}
				String lnheader = ln.substring(0, substrlen);
				String[] lnheaderelems = lnheader.split("\t");

				Chromosome chr = Chromosome.parseChr(lnheaderelems[0]);
				if (!chr.equals(Chromosome.X)) {
					// return this.chr + "_" + this.pos + "_" + this.id;
					String varStr = lnheaderelems[0] + "_" + lnheaderelems[1] + "_" + lnheaderelems[2];
					varStrs.add(varStr);
				}
			}
			ln = tf1.readLine();
		}
		tf1.close();

		System.out.println(varStrs.size() + " variants loaded");

		int nrSeenKnownAtAll = 0;
		int nrSeenUnknownAtAll = 0;
		int nrSeenKnownImpQual = 0;
		int nrSeenUnkownImpQual = 0;
		String[] vcf2arr = vcf2list.split(",");
		int ctr = 0;
		for (String file : vcf2arr) {
			TextFile tf = new TextFile(file, TextFile.R);
			System.out.println(ctr + "/" + vcf2arr.length + "\tParsing: " + file);
			ln = tf.readLine();
			int lnctr = 0;
			while (ln != null) {
				if (!ln.startsWith("#")) {
					int strlen = ln.length();
					int substrlen = 1000;
					if (strlen < substrlen) {
						substrlen = strlen;
					}
					String lnheader = ln.substring(0, substrlen);
					String[] lnheaderelems = lnheader.split("\t");

					// return this.chr + "_" + this.pos + "_" + this.id;
					String varStr = lnheaderelems[0] + "_" + lnheaderelems[1] + "_" + lnheaderelems[2];

					boolean present = varStrs.contains(varStr);

					if (present) {

						boolean qualok = false;
						String[] info = lnheaderelems[7].split(";");
						for (int i = 0; i < info.length; i++) {
							if (info[i].startsWith("INFO")) {
								String[] scorelemes = info[i].split("=");
								Double q = Double.parseDouble(scorelemes[1]);
								if (q > 0.3) {
									qualok = true;
								}
							}
						}


						if (lnheaderelems[2].startsWith("rs")) {
							nrSeenKnownAtAll++;
							if (qualok) {
								nrSeenKnownImpQual++;
							}
						} else {
							if (qualok) {
								nrSeenUnkownImpQual++;
							}
							nrSeenUnknownAtAll++;
						}


					}

					lnctr++;
					if (lnctr % 10000 == 0) {
						System.out.println(lnctr + " lines parsed");
					}
				}
				ln = tf.readLine();

			}
			tf.close();
			ctr++;
		}

		int sum = nrSeenKnownAtAll + nrSeenUnknownAtAll;
		int sumQual = nrSeenUnkownImpQual + nrSeenKnownImpQual;


		System.out.println("Known seen at all: " + nrSeenKnownAtAll);
		System.out.println("Known seen rsq>0.3: " + nrSeenKnownImpQual);
		System.out.println("UnKnown seen at all: " + nrSeenUnknownAtAll);
		System.out.println("UnKnown seen rsq>0.3: " + nrSeenUnkownImpQual);
		System.out.println("Sum: " + sum);
		System.out.println("Sum rsq>0.3: " + sumQual);


	}
}
