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

	private final String tabixfile;
	private final TabixReader treader;

	public VCFTabix(String filename) throws IOException {
		this.tabixfile = filename;
		treader = new TabixReader(tabixfile);

	}

	public boolean[] getSampleFilter(String tabixsamplelimit) throws IOException {
		if (tabixsamplelimit == null) {
			return null;
		}
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

	public TabixReader.Iterator query(Feature region) throws IOException {
		int start = region.getStart() - 10;
		if (start < 0) {
			start = 1;
		}
		int stop = region.getStop() + 10;
		if (stop > region.getChromosome().getLength()) {
			stop = region.getChromosome().getLength();
		}

		TabixReader.Iterator window = treader.query(region.getChromosome().getNumber() + ":" + start + "-" + stop);
		return window;
	}

	public void close() {
		treader.close();
	}

	public ArrayList<VCFVariant> getAllVariants(Feature f, boolean[] samplefilter) throws IOException {
		TabixReader.Iterator iterator = query(f);
		String next = iterator.next();
		ArrayList<VCFVariant> output = new ArrayList<>();
		while (next != null) {
			output.add(new VCFVariant(next, VCFVariant.PARSE.ALL, samplefilter));
			next = iterator.next();
		}
		return output;
	}
}
