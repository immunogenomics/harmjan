package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

/**
 * Created by hwestra on 2/10/16.
 */
public class VCFSampleFilter {

	public void filteroverlapping(String toFilter, String ref, String vcfout, boolean keep) throws IOException {

		System.out.println("Filtering: " + toFilter);
		System.out.println("for samples present in: " + ref);
		System.out.println("out: " + vcfout);


		HashSet<String> samples1Hash = new HashSet<String>();
		VCFGenotypeData data1 = new VCFGenotypeData(ref);
		samples1Hash.addAll(data1.getSamples());
		data1.close();

		filter(toFilter, samples1Hash, vcfout, keep);

	}

	public void filter(String fileIn, String fileout, String sampleFile, boolean keep) throws IOException {

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

		filter(fileIn, samples, fileout, keep);

	}

	public void filter(String toFilter, HashSet<String> sampleHash, String vcfout, boolean keep) throws IOException {
		TextFile out = new TextFile(vcfout, TextFile.W);
		TextFile in = new TextFile(toFilter, TextFile.R);

		String ln = in.readLine();

		boolean[] includecol = null;
		int excluded = 0;
		while (ln != null) {
			if (ln.startsWith("##")) {
				out.writeln(ln);
			} else if (ln.startsWith("#CHROM")) {
				String[] elems = ln.split("\t");

				for (int i = 0; i < 9; i++) {
					includecol[i] = true;
				}

				includecol = new boolean[elems.length];
				if (!keep) {
					for (int i = 0; i < includecol.length; i++) {
						includecol[i] = true;
					}
				}

				for (int i = 9; i < elems.length; i++) {
					if (sampleHash.contains(elems[i])) {
						if (keep) {
							includecol[i] = true;
						} else {
							includecol[i] = false;
						}
						excluded++;
					}
				}

				out.writeln(Strings.concat(elems, includecol, Strings.tab));
				System.out.println(excluded + " / " + (elems.length - 9) + " samples will be removed");
			} else {
				String[] elems = ln.split("\t");
				out.writeln(Strings.concat(elems, includecol, Strings.tab));
			}
			ln = in.readLine();
		}

		in.close();
		out.close();
	}


}