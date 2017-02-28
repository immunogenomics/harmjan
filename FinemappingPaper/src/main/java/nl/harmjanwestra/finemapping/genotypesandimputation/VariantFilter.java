package nl.harmjanwestra.finemapping.genotypesandimputation;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 2/13/17.
 */
public class VariantFilter {


	public static void main(String[] args) {
		String a = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/SequencingPanel/seqpanelfiltered-maf1percent.txt";
		String b = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/SequencingPanel/regions.bed";
		String c = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/SequencingPanel/seqpanelfiltered-maf1percent-regionsfiltered.txt";


		BedFileReader reader = new BedFileReader();
		try {
			ArrayList<Feature> feats = reader.readAsList(b);


			TextFile out = new TextFile(c, TextFile.W);
			TextFile tf = new TextFile(a, TextFile.R);

			out.writeln(tf.readLine());

			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {

				if (elems.length >= 2) {
					Chromosome chr = Chromosome.parseChr(elems[0]);
					Integer start = Integer.parseInt(elems[1]);
					Integer stop = start + 1;
					Feature f = new Feature(chr, start, stop);
					boolean overlap = false;
					for (Feature r : feats) {
						if (r.overlaps(f)) {
							overlap = true;
						}
					}

					if (overlap) {
						out.writeln(Strings.concat(elems, Strings.tab));
					}
				}

				elems = tf.readLineElems(TextFile.tab);
			}

			tf.close();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
