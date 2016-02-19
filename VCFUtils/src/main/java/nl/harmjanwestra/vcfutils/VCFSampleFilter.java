package nl.harmjanwestra.vcfutils;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

/**
 * Created by hwestra on 2/10/16.
 */
public class VCFSampleFilter {

	public void filter(String fileIn, String fileout, String sampleFile) throws IOException {


		System.out.println("Sample Filter");
		System.out.println("in: " + fileIn);
		System.out.println("out: " + fileout);
		System.out.println("file: " + sampleFile);

		TextFile tf1 = new TextFile(sampleFile, TextFile.R);
		String[] elems = tf1.readLineElems(TextFile.tab);
		HashSet<String> samples = new HashSet<String>();
		while (elems != null) {
			samples.add(elems[0]);
			elems = tf1.readLineElems(TextFile.tab);
		}
		tf1.close();

		System.out.println(samples.size() + " samples loaded from: " + sampleFile);

		TextFile vcfin = new TextFile(fileIn, TextFile.R);
		TextFile out = new TextFile(fileout, TextFile.W);
		String ln = vcfin.readLine();
		int lnCtr = 0;
		boolean[] includecol = null;
		while (ln != null) {
			if (ln.startsWith("##")) {
				out.writeln(ln);
			} else if (ln.startsWith("#CHROM")) {
				elems = ln.split("\t");
				includecol = new boolean[elems.length];
				int samppleCtr = 0;

				for (int i = 0; i < elems.length; i++) {
					if (i < 9) {
						includecol[i] = true;
					} else if (samples.contains(elems[i])) {
						includecol[i] = true;

						samppleCtr++;
					}
				}
				System.out.println("After parsing header: " + samppleCtr + " samples found");
				out.writeln("##VCFSampleFilter=" + samppleCtr + "/" + (elems.length - 9) + " samples selected using " + sampleFile);
				out.writeln(Strings.concat(elems, includecol, Strings.tab));
			} else {
				if (includecol != null) {
					elems = ln.split("\t");
					out.writeln(Strings.concat(elems, includecol, Strings.tab));
				}


				lnCtr++;

				if (lnCtr % 1000 == 0) {
					System.out.println(lnCtr + " variants parsed");

				}
			}
			ln = vcfin.readLine();
		}
		out.close();
		vcfin.close();


	}

}