package nl.harmjanwestra.utilities.bamfile.filters;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * Created by hwestra on 9/18/15.
 */
public class MappingQualityUnavailableFilter implements SamRecordFilter {

	@Override
	public boolean filterOut(SAMRecord samRecord) {
		int mapq = samRecord.getMappingQuality();
		return (mapq == samRecord.NO_MAPPING_QUALITY || mapq == samRecord.UNKNOWN_MAPPING_QUALITY);

	}

	@Override
	public boolean filterOut(SAMRecord samRecord, SAMRecord samRecord1) {
		return false;
	}
}
