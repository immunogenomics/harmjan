package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.enums.Chromosome;
import umcg.genetica.console.FileReadProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;

/**
 * Created by hwestra on 12/15/16.
 * rewrite blueprint QTL files to tab file
 */
public class BlueprintQTLRewrite {

	public static void main(String[] args) {


		if (args.length < 2) {
			System.out.println("usage: filein fileout");

		} else {
			BlueprintQTLRewrite f = new BlueprintQTLRewrite();
			try {
				f.run(args[0], args[1]);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	}

	public void run(String in, String out) throws IOException {

		/*
		column_number   column_label
1               chr:pos_ref_alt
2               rsid
3               phenotypeID
4               p.value
5               beta
6               Bonferroni.p.value
7               FDR
8               alt_allele_frequency
9               std.error_of_beta
		 */

		TextFile tf = new TextFile(in, TextFile.R);
		TextFile outf = new TextFile(out, TextFile.W);
		String[] elems = tf.readLineElems(Strings.whitespace);
		FileReadProgressBar pb = new FileReadProgressBar(1000000, true);
		while (elems != null) {
			// filter for FDR < 5%
			if (elems.length < 6) {

			} else {
				Double fdr = Double.parseDouble(elems[6]);
				if (fdr < 0.05) {
					String[] poselems = elems[0].split(":");
					Chromosome chr = Chromosome.parseChr(poselems[0]);
					String[] poselemselems = poselems[1].split("_");
					String pos = poselemselems[0];

					String rsid = elems[1];
					String pheno = elems[2];
					String pval = elems[3];

					String outln = chr.toString() + "\t" + pos + "\t" + rsid + "\t" + pheno + "\t" + pval;
					outf.writeln(outln);
					pb.iteratectr2();
				}
				pb.iterate();
			}

			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();
		outf.close();
		pb.close();


	}
}
