package nl.harmjanwestra.utilities.vcf;

import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 10/6/16.
 */
public class VCFTabix {


	public static boolean[] getSampleFilter(String tabixfile, String tabixsamplelimit) throws IOException {
		boolean[] samplesToInclude = null;


		VCFGenotypeData d = new VCFGenotypeData(tabixfile);
		ArrayList<String> tabixSamples = d.getSamples();
		samplesToInclude = new boolean[tabixSamples.size()];

		TextFile tf = new TextFile(tabixsamplelimit, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);
		HashSet<String> allSamplesInListSet = new HashSet<String>();
		while (elems != null) {
			allSamplesInListSet.add(elems[0]);
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();

		for (int s = 0; s < tabixSamples.size(); s++) {
			if (allSamplesInListSet.contains(tabixSamples.get(s))) {
				samplesToInclude[s] = true;
			}
		}

		return samplesToInclude;
	}

	public static TabixReader.Iterator query(String tabixfile, Feature region) throws IOException {
		TabixReader treader = new TabixReader(tabixfile);
		TabixReader.Iterator window = treader.query(region.getChromosome().getNumber() + ":" + (region.getStart() - 10) + "-" + (region.getStop() + 10));
		return window;
	}

}
