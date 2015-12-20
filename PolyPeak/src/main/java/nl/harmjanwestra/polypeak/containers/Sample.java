/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.polypeak.containers;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import nl.harmjanwestra.utilities.bamfile.BamFileReader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

/**
 * @author hwestra
 */
public class Sample {
	private final boolean hasReaders;
	private final File[] bamfiles;
	private boolean isControl = false;
	private BamFileReader[] readers;
	private final HashSet<String> readgroups;
	private final String samplename;
	private ArrayList<String> sequences;


	public Sample(String samplename, File[] bamfiles, boolean checkSampleNamesInBAMFiles, boolean isControl) throws IOException {
		this(samplename, bamfiles, checkSampleNamesInBAMFiles);
		this.isControl = isControl;
	}

	public Sample(String samplename, File[] bamfiles, boolean checkSampleNamesInBAMFiles) throws IOException {

		if (samplename == null) {
			throw new IllegalArgumentException("Sample name not set!");
		}
		this.samplename = samplename;

		if (bamfiles == null) {
			throw new IllegalArgumentException("BAM file location not set!");
		}

		this.bamfiles = bamfiles;
		this.readgroups = new HashSet<String>();

		this.checkBAMFilesAndInitialize(checkSampleNamesInBAMFiles);
		if (readers == null || readers.length == 0) {
			throw new IllegalArgumentException(bamfiles.length + " BAM file(s) defined for sample: " + samplename + " but sample not found in BAM file(s)");
		} else {
			this.hasReaders = true;
		}
	}

	private void checkBAMFilesAndInitialize(boolean checkSampleNamesInBAMFiles) throws IOException {

		ArrayList<BamFileReader> tmpReaders = new ArrayList<BamFileReader>();
		for (int i = 0; i < bamfiles.length; i++) {
			// check whether file exists
			File f = bamfiles[i];
			if (!f.exists()) {
				throw new IOException("File: " + bamfiles[i] + " does not exist.");
			}

			BamFileReader r = new BamFileReader(bamfiles[i]);
			if (checkSampleNamesInBAMFiles) {
				List<SAMReadGroupRecord> allBamFileReadGroups = r.getReadGroups();

				if (allBamFileReadGroups.isEmpty()) {
					System.err.println("WARNING: trying to check sample names in BAM file: " + bamfiles[i].getAbsolutePath() + " but no readAsTrack group information found.");
					System.err.println("Forcing: " + bamfiles[i].getAbsolutePath() + " as BAM file for sample: " + samplename);
					tmpReaders.add(r);
				} else {
					boolean fileContainsSampleInfo = false;
					for (SAMReadGroupRecord record : allBamFileReadGroups) {
						String sample = record.getSample();
						if (sample.equals(samplename)) {
							readgroups.add(record.getReadGroupId());
							tmpReaders.add(r);
							fileContainsSampleInfo = true;
						}
					}
					if (!fileContainsSampleInfo) {
						r.close();
						System.err.println("Warning: BAM file " + bamfiles[i].getAbsolutePath() + " defined for sample: " + samplename + " but samplename not found in BAM file.");
					}
				}
			} else {
				System.out.println("Adding BAM file: " + bamfiles[i] + " for sample: " + samplename);
				tmpReaders.add(r);
			}

		}

		if (tmpReaders.isEmpty()) {
			System.err.println("WARNING: no sample");
		}

		this.readers = tmpReaders.toArray(new BamFileReader[0]);

		System.out.println("Sample: " + samplename + " has " + readers.length + " BAM file readers");


		// check the sequence libraries: do they conflict?
		for (BamFileReader reader : readers) {

			SAMSequenceDictionary dictionary = reader.getHeader().getSequenceDictionary();
			List<SAMSequenceRecord> sequenceList = dictionary.getSequences();

			for (SAMSequenceRecord record : sequenceList) {
				boolean b = SequenceLibrary.checkSequence(record.getSequenceName(), record.getSequenceLength());

				if (!b) {
					Sequence s = new Sequence(record.getSequenceName(), record.getSequenceLength());
					Sequence otherSequence = SequenceLibrary.getSequence(s);
					System.err.println("Warning: sequence length differs for sample: " + samplename + ". Sequence: " + record.getSequenceName() + " differs from " + otherSequence + " in length: " + record.getSequenceLength() + " vs " + otherSequence.getLength());
					if (s.getLength() > otherSequence.getLength()) {
						SequenceLibrary.setSequence(s);
					}
				}
			}
		}
	}

	public String getName() {
		return this.samplename;
	}

	public boolean hasReaders() {
		return hasReaders;
	}

	public void close() throws IOException {
		for (BamFileReader reader : readers) {
			if (reader != null) {
				reader.close();
			}
		}
	}


	public BamFileReader[] getReaders() {
		return readers;
	}
}
