package nl.harmjanwestra.harmonics;

import htsjdk.samtools.*;
import htsjdk.samtools.filter.*;
import nl.harmjanwestra.utilities.bamfile.filters.FailsVendorQualityCheckFilter;
import nl.harmjanwestra.utilities.bamfile.filters.MappingQualityUnavailableFilter;
import nl.harmjanwestra.utilities.bamfile.filters.NotMappingToAutosomeFilter;
import nl.harmjanwestra.utilities.bamfile.filters.UnmappedReadFilter;
import nl.harmjanwestra.utilities.features.Chromosome;
import umcg.genetica.io.Gpio;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

/**
 * Created by hwestra on 9/8/15.
 */
public class FilterAndShift {

	int[] shift = new int[3];


	public static void main(String[] args) {
		try {

			if (args.length < 3) {
				System.out.println("Usage: bamin bamout boolAtacShift");
			} else {
				FilterAndShift t = new FilterAndShift();
				t.run(args[0], args[1], Boolean.parseBoolean(args[2]));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String bamin, String bamout, boolean atacshift) throws IOException {


		shift[0] = 5;
		shift[1] = -4;

		SamReader reader = SamReaderFactory.makeDefault().open(new File(bamin));

		SAMSequenceDictionary dictionary = reader.getFileHeader().getSequenceDictionary();
		List<SAMSequenceRecord> sequences = dictionary.getSequences();
		HashSet<String> allowedSequences = new HashSet<String>();
		for (SAMSequenceRecord record : sequences) {
			Chromosome chr = Chromosome.parseChr(record.getSequenceName());
			if (chr.equals(Chromosome.NA) || chr.equals(Chromosome.MT)) {

			} else {
				allowedSequences.add(record.getSequenceName());
				System.out.println("Found sequence: " + record.getSequenceName());
			}
		}


		int maxInsert = 10000;

		ArrayList<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();
		AggregateFilter filter = new AggregateFilter(filters);

		filters.add(new NotPrimaryAlignmentFilter());
		filters.add(new FailsVendorQualityCheckFilter());
		filters.add(new DuplicateReadFilter());
		filters.add(new UnmappedReadFilter());
		filters.add(new MappingQualityUnavailableFilter());
		filters.add(new MappingQualityFilter(20));
		filters.add(new NotMappingToAutosomeFilter());

		filters.add(new InsertSizeFilter(0, maxInsert));

		SamRecordFilter f1 = new NotPrimaryAlignmentFilter();
		SamRecordFilter f2 = new FailsVendorQualityCheckFilter();
		SamRecordFilter f3 = new DuplicateReadFilter();
		SamRecordFilter f4 = new UnmappedReadFilter();
		SamRecordFilter f5 = new MappingQualityUnavailableFilter();
		SamRecordFilter f6 = new MappingQualityFilter(20);
		SamRecordFilter f7 = new NotMappingToAutosomeFilter();


		if (Gpio.exists(bamout)) {
			Gpio.delete(bamout);
		}

		SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), false, new File(bamout));
		SAMRecordIterator iterator = reader.iterator();
		int nrPrinted = 0;
		int nrRemoved = 0;
		int nrTotal = 0;

		int ctr1 = 0;
		int ctr2 = 0;
		int ctr3 = 0;
		int ctr4 = 0;
		int ctr5 = 0;
		int ctr6 = 0;
		int ctr7 = 0;


		while (iterator.hasNext()) {
			SAMRecord next = iterator.next();

			if (f1.filterOut(next)) {
				ctr1++;
			}

			if (f2.filterOut(next)) {
				ctr2++;
			}

			if (f3.filterOut(next)) {
				ctr3++;
			}

			if (f4.filterOut(next)) {
				ctr4++;
			}

			if (f5.filterOut(next)) {
				ctr5++;
			}

			if (f6.filterOut(next)) {
				ctr6++;
			}

			if (f7.filterOut(next)) {
				ctr7++;
			}


			if (!filter.filterOut(next)) {

				if (atacshift) {
					if (next.getReadNegativeStrandFlag()) {
						next.setAlignmentStart(next.getAlignmentStart() + shift[1]);
						next.setMateAlignmentStart(next.getAlignmentStart() + shift[0]);
					} else {
						next.setAlignmentStart(next.getAlignmentStart() + shift[0]);
						next.setMateAlignmentStart(next.getAlignmentStart() + shift[1]);
					}
				}
				outputSam.addAlignment(next);
				nrPrinted++;


			} else {
				nrRemoved++;
			}

			nrTotal++;
			if (nrTotal % 10000000 == 0) {
				System.out.println(nrTotal + " reads " + nrRemoved + " removed " + nrPrinted + " kept");
				System.out.println("NotPrimaryAlignmentFilter: " + ctr1);
				System.out.println("FailsVendorQualityCheckFilter: " + ctr2);
				System.out.println("DuplicateReadFilter: " + ctr3);
				System.out.println("UnmappedReadFilter: " + ctr4);
				System.out.println("MappingQualityUnavailableFilter: " + ctr5);
				System.out.println("MappingQualityFilter: " + ctr6);
				System.out.println("NotMappingToAutosomeFilter: " + ctr7);

			}
		}
		iterator.close();
		outputSam.close();
		reader.close();

	}
}
