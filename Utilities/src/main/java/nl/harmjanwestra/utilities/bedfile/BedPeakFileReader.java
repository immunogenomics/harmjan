package nl.harmjanwestra.utilities.bedfile;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.enums.Strand;

/**
 * Created by hwestra on 1/6/15.
 */
public class BedPeakFileReader extends BedFileReader {
	@Override
	protected Feature parseElems(String[] elems) {

		int len = 0;
		Chromosome featureChr = Chromosome.parseChr(elems[0]);

		Strand featureStrand = Strand.NA;

		featureStrand = Strand.parseStr(elems[elems.length - 1]);

		int featureStart = -1;
		int featureStop = -1;
		try {
			featureStart = Integer.parseInt(elems[1]);
		} catch (NumberFormatException e) {
			System.out.println("Could not parse chromosome start position: " + elems[1]);
		}

		try {
			featureStop = Integer.parseInt(elems[2]);
		} catch (NumberFormatException e) {
			System.out.println("Could not parse chromosome stop position: " + elems[2]);
		}


		len = featureStop - featureStart;
		featureLengthSum += len;
		nrFeatures++;
		Feature f = new Feature();
		f.setChromosome(featureChr);
		f.setStrand(featureStrand);
		f.setStart(featureStart);
		f.setStop(featureStop);


		int peakPos = 0;
		double peakQ = 0;
		double foldChange = 0;

		if (elems.length > 4) {
			try {
				peakPos = Integer.parseInt(elems[4]);
				peakQ = Double.parseDouble(elems[elems.length - 1]);
				foldChange = Double.parseDouble(elems[elems.length - 2]);

			} catch (NumberFormatException e) {
				System.out.println("Could not parse chromosome stop position: " + elems[elems.length - 1] + " or " + elems[elems.length - 2]);

			}
		}


		if (filter == null || filter.passesFilter(f)) {
			return f;
		} else {
			return null;
		}

	}
}
