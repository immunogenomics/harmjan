package nl.harmjanwestra.ngs;

import htsjdk.samtools.*;
import nl.harmjanwestra.utilities.bamfile.BamFileReader;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.FeatureMerger;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

/**
 * Created by hwestra on 6/22/15.
 */
public class SAMFilter {

	public static void main(String[] args) {
		if (args.length < 3) {
			System.out.println("Usage: bamin bamout bedfile");
		} else {
			SAMFilter filter = new SAMFilter();
			try {
				filter.filterRegions(args[0], args[1], args[2]);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	public void filterRegions(String samfilename, String samout, String bedfilename) throws IOException {

		BedFileReader bed = new BedFileReader();
		ArrayList<Feature> regions = bed.readAsList(bedfilename);

		FeatureMerger merger = new FeatureMerger();
		regions = merger.merge(regions, false);

		System.out.println(regions.size() + " regions after merging....");

		// get overlapping chromosomenames...


		BamFileReader reader = new BamFileReader(samfilename);


		SAMFileHeader header = reader.getHeader();

//		SAMFileHeader newHeader = new SAMFileHeader();
//		newHeader.setReadGroups(header.getReadGroups());
//		newHeader.setComments(header.getComments());
//		newHeader.setProgramRecords(header.getProgramRecords());
//		newHeader.setGroupOrder(header.getGroupOrder());
//
//		newHeader.setSortOrder(header.getSortOrder());
//		newHeader.setTextHeader(header.getTextHeader());
//		newHeader.setValidationErrors(header.getValidationErrors());
//
//		SAMSequenceDictionary newDictionary = new SAMSequenceDictionary();
//
//		newHeader.setSequenceDictionary(reader.);
//		CoverageTask t = new CoverageTask();
		HashMap<String, String> chromosomeToSequence = reader.matchChromosomeNames(regions);
//		Set<String> chroms = chromosomeToSequence.keySet();
//		for (String s : chroms) {
//			SAMSequenceRecord record = reader.getHeader().getSequence(s);
//			if (record != null) {
//				newDictionary.addSequence(record);
//			}
//		}


		File samoutfile = new File(samout);

		SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, samoutfile);

		Collections.sort(regions, new FeatureComparator(false));

		int nrWritten = 0;
		for (Feature f : regions) {
			Chromosome c = f.getChromosome();
			int start = f.getStart();
			int stop = f.getStop();
			String chrName = chromosomeToSequence.get(c.getName());

			if (chrName != null) {
				SAMRecordIterator it = reader.query(chrName, start - 1000, stop + 1000, true);
				while (it.hasNext()) {
					SAMRecord record = it.next();
					outputSam.addAlignment(record);
					nrWritten++;
				}
				it.close();
			}

		}

		reader.close();
		outputSam.close();

		System.out.println(nrWritten + " records written");
	}

}
