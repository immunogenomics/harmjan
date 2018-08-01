package nl.harmjanwestra.utilities.bamfile.filters;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * Created by hwestra on 9/18/15.
 */
public class FailsVendorQualityCheckFilter implements SamRecordFilter {

	@Override
	public boolean filterOut(SAMRecord samRecord) {
		return samRecord.getReadFailsVendorQualityCheckFlag();
	}

	@Override
	public boolean filterOut(SAMRecord samRecord, SAMRecord samRecord1) {
		return false;
	}
}
