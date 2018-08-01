package nl.harmjanwestra.utilities.coverage;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.filter.AggregateFilter;
import nl.harmjanwestra.utilities.bamfile.BamFileReader;
import nl.harmjanwestra.utilities.features.Feature;

import java.util.ArrayList;

/**
 * Created by hwestra on 9/4/15.
 */
public class BasicStats {


	public BasicStats(BamFileReader reader, AggregateFilter filter, int distributionSize) {

		// iterate BAM path, collect:
		// average insert size
		// total mapped reads
		// total mapped fragments
		// reads mapped per chr
		// duplicates per chr

	}

	public BasicStats(BamFileReader reader, AggregateFilter filter, int distributionSize, ArrayList<Feature> regions) {

		// iterate BAM path, collect:
		// average insert size
		// total mapped reads
		// total mapped fragments
		// reads mapped per chr
		// duplicates per chr

	}

	int totalReadsMapped = 0;
	int totalReadsMappedPassFilter = 0;
	int totalFragmentsMapped = 0;
	int totalFragmentsMappedPassFilter = 0;
	double averageInsertSize = 0;


	private void stats(BamFileReader reader, AggregateFilter filter, ArrayList<Feature> regions) {
		// count total number of mapping reads..
		int totalReadsMapped = 0;
		SAMRecordIterator it = reader.iterator();
		int readcounter = 0;
		while (it.hasNext()) {
			SAMRecord record = it.next();
			readcounter++;
			if (readcounter % 100000 == 0) {
				System.out.println(readcounter + " reads processed.");
			}
			if (filter == null || !filter.filterOut(record)) {
				totalReadsMapped++;
			}
		}
		it.close();
	}


}
