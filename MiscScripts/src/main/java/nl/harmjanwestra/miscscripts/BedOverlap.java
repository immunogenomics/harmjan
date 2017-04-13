package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 8/14/16.
 */
public class BedOverlap {

	public static void main(String[] args) {
		String bed1 = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci-woMHC.txt";
		String bed2 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed";
		String out1 = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/OverlapWithSequencing/ICLociOverlappingWithSequencingRegions.bed";
		String out2 = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/OverlapWithSequencing/SequencingRegionsOverlappingWithIc.bed";
		BedOverlap b = new BedOverlap();
		try {
			b.run(bed1, bed2, out1, out2);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String bed1, String bed2, String out1, String out2) throws IOException {

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> list1 = reader.readAsList(bed1);
		ArrayList<Feature> list2 = reader.readAsList(bed2);

		TextFile outf1 = new TextFile(out1, TextFile.W);
		TextFile outf2 = new TextFile(out2, TextFile.W);
		HashSet<Feature> overlapIn2 = new HashSet<Feature>();
		for (Feature f : list1) {
			boolean overlap = false;
			for (Feature f2 : list2) {
				if (f.overlaps(f2)) {
					overlap = true;
					overlapIn2.add(f2);
				}
			}
			if (overlap) {
				outf1.writeln(f.toBedString());

			}
		}

		for (Feature f : overlapIn2) {
			outf2.writeln(f.toBedString());
		}
		outf2.close();
		outf1.close();


	}
}
