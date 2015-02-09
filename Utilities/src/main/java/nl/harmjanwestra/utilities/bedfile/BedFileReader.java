/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.utilities.bedfile;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.Strand;
import nl.harmjanwestra.utilities.features.Track;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * @author Harm-Jan
 */
public class BedFileReader {


	protected int nrFeatures;
	protected int featureLengthSum;
	protected BedFileFeatureFilter filter = null;

	public void setFilter(BedFileFeatureFilter filter) {
		this.filter = filter;
	}

	public Track read(String file, String name) throws IOException {
		nrFeatures = 0;
		featureLengthSum = 0;
		TextFile tf = new TextFile(file, TextFile.R);

		System.out.println("Reading file: " + file);
		Feature filter = null;


		// chr1	8128340	8128539	C011PABXX110504:4:2203:14692:158380	0	-
		String[] elems = tf.readLineElems(TextFile.tab);
		Track track = new Track(name);

		while (elems != null) {
			Feature f = parseElems(elems);
			if (f != null) {
				track.addFeature(f);
			}
			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();

		System.out.println("Average feature featureLengthSum: " + ((double) featureLengthSum / nrFeatures) + "\tNumber of elements: " + nrFeatures);
		track.printNrFeatures();
		return track;
	}

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
		if (filter == null || filter.passesFilter(f)) {
			return f;
		} else {
			return null;
		}
	}

}
