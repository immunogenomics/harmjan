package nl.harmjanwestra.utilities.bamfile.filters;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * Created by hwestra on 9/18/15.
 */
public class UnmappedReadFilter implements SamRecordFilter {

	@Override
	public boolean filterOut(SAMRecord samRecord) {
		return samRecord.getReadUnmappedFlag() || samRecord.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START;

	}

	@Override
	public boolean filterOut(SAMRecord samRecord, SAMRecord samRecord1) {
		return false;
	}


}
