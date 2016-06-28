package nl.harmjanwestra.utilities.bamfile.filters;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;
import nl.harmjanwestra.utilities.enums.Chromosome;

/**
 * Created by hwestra on 10/13/15.
 */
public class NotMappingToAutosomeFilter implements SamRecordFilter {
	@Override
	public boolean filterOut(SAMRecord samRecord) {
		String refName = samRecord.getReferenceName();
		Chromosome chr = Chromosome.parseChr(refName);
		return !chr.isAutosome();
	}

	@Override
	public boolean filterOut(SAMRecord samRecord, SAMRecord samRecord1) {
		return false;
	}
}
