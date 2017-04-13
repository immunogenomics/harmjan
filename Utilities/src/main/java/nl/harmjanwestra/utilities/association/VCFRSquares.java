package nl.harmjanwestra.utilities.association;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 10/27/15.
 */
public class VCFRSquares {
	public ArrayList<Feature> loadVariantInfo(String file, Feature region) throws IOException {
		ArrayList<Feature> output = new ArrayList<Feature>();
		TextFile tf = new TextFile(file, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {
			if (elems.length > 1 && !elems[0].startsWith("#")) {
				Chromosome chr = Chromosome.parseChr(elems[0]);
				Integer pos = Integer.parseInt(elems[1]);
//					System.out.println(chr.toString() + "\t" + pos);
				Feature f = new Feature();
				f.setChromosome(chr);
				f.setStart(pos);
				f.setStop(pos);
				if (region.overlaps(f)) {
					output.add(f);
				}
			}
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();
		return output;
	}

	public ArrayList<Pair<Feature, Double>> loadRSquareds(String file, Feature region) throws IOException {
		ArrayList<Pair<Feature, Double>> output = new ArrayList<Pair<Feature, Double>>();
		TextFile tf = new TextFile(file, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {
			if (elems.length > 1 && !elems[0].startsWith("#")) {
				Chromosome chr = Chromosome.parseChr(elems[0]);
				Integer pos = Integer.parseInt(elems[1]);
//					System.out.println(chr.toString() + "\t" + pos);
				Feature f = new Feature();
				f.setChromosome(chr);
				f.setStart(pos);
				f.setStop(pos);
				if (region.overlaps(f)) {
					String[] infoElems = elems[7].split(";");
					for (int e = 0; e < infoElems.length; e++) {
						if (infoElems[e].startsWith("AR2")) {
							String[] infoElemsElems = infoElems[e].split("=");
							Double ar2 = Double.parseDouble(infoElemsElems[1]);
							// plot
							Pair<Feature, Double> p = new Pair<Feature, Double>(f, ar2);
							output.add(p);
						}
					}

				}
			}
			elems = tf.readLineElems(Strings.whitespace);
		}
		return output;
	}
}
