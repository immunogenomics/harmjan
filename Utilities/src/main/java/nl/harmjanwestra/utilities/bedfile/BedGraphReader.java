package nl.harmjanwestra.utilities.bedfile;

import nl.harmjanwestra.utilities.features.BedGraphFeature;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by Harm-Jan on 04/24/16.
 */
public class BedGraphReader {

	public ArrayList<BedGraphFeature> read(String file, boolean fileHasStrandInfo, ArrayList<Feature> regionList) throws IOException {

		System.out.println("Reading bedgraph: " + file);
		TextFile tf = new TextFile(file, TextFile.R);

		ArrayList<BedGraphFeature> output = new ArrayList<>();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

			if (elems.length >= 4) {
				Chromosome chr = Chromosome.parseChr(elems[0]);

				Integer start = Integer.parseInt(elems[1]);
				Integer stop = Integer.parseInt(elems[2]);
				double value = Double.parseDouble(elems[3]);
				BedGraphFeature feat = new BedGraphFeature(chr, start, stop);
				feat.setValue(value);
				if (elems.length > 3) {
					if (fileHasStrandInfo) {
						double valuepos = Double.parseDouble(elems[4]);
						double valueneg = Double.parseDouble(elems[5]);
						feat.setValue(valuepos, valueneg);
					}
				}

				boolean overlap = false;
				if (regionList != null) {
					for (Feature f : regionList) {
						if (f.overlaps(feat)) {
							overlap = true;
						}
					}
				}
				if (overlap || regionList == null) {
					output.add(feat);
				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		System.out.println();

		return output;

	}
}
