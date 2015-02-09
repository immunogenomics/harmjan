/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.utilities.bamfile;

import htsjdk.samtools.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * @author Harm-Jan
 */
public class BamFileReader {

	private final File bamFile;
	private final SamReader reader;

	private static final Logger logger = LogManager.getLogger();

	public BamFileReader(String bamFile) throws IOException {
		this(new File(bamFile), false);
	}

	public BamFileReader(File bamFile) throws IOException {
		this(bamFile, false);
	}

	public BamFileReader(File bamFile, boolean indexFileWhenNotIndexed) throws IOException {

		this.bamFile = bamFile;

		SamReader tmpreader = SamReaderFactory.make()
				.enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES)
				.validationStringency(ValidationStringency.DEFAULT_STRINGENCY)
				.open(bamFile);


		if (tmpreader.hasIndex()) {

			reader = tmpreader;
			SAMFileHeader header = reader.getFileHeader();
			if (!header.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
				logger.warn("Sort order is not based on coordinate, so access may be slow for: " + bamFile.getAbsolutePath());
			}
		} else {
			tmpreader.close();
			logger.info("BAM file does not have an index: " + bamFile.getAbsolutePath());
			if (indexFileWhenNotIndexed) {
				// index file
				logger.info("Indexing BAM file: " + bamFile.getAbsolutePath());
				SAMFileReader samfileReader = new SAMFileReader(bamFile);
				String outputPath = bamFile.getAbsolutePath();
//                outputPath = outputPath.substring(3);
				outputPath += ".bai";
				BAMIndexer.createIndex(samfileReader, new File(outputPath));
				samfileReader.close();
				logger.info("Done indexing: " + outputPath);
			}
			tmpreader = SamReaderFactory.make()
					.enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES)
					.validationStringency(ValidationStringency.DEFAULT_STRINGENCY)
					.open(bamFile);
			reader = tmpreader;

		}

		SAMFileHeader header = reader.getFileHeader();
		logger.info("Summary:");
		logger.info("File name: " + bamFile.getAbsolutePath());
		logger.info("Creator: " + header.getCreator());
		logger.info("");
		int nrAligned = 0;
		int nrNonAligned = 0;
		int nrTotal = 0;
		int nrSequences = 0;
		for (SAMSequenceRecord record : header.getSequenceDictionary().getSequences()) {
//
			nrSequences++;
//			logger.debug("Sequence: " + record.getSequenceName());
//			logger.debug("Sequence index: " + record.getSequenceIndex());
//			logger.debug("Sequence len: " + record.getSequenceLength());
			BAMIndexMetaData metaData = reader.indexing().getIndex().getMetaData(record.getSequenceIndex());
			int tmpnrAligned = metaData.getAlignedRecordCount();
			int tmpnrNonAligned = metaData.getAlignedRecordCount();
			int tmpnrTotal = tmpnrAligned + tmpnrNonAligned;

			nrAligned += tmpnrAligned;
			nrNonAligned += tmpnrNonAligned;
			nrTotal += tmpnrTotal;
//			logger.debug("Aligned: " + tmpnrAligned);
//			logger.debug("NonAligned: " + tmpnrNonAligned);
//			logger.debug("NrTotal: " + tmpnrTotal);

		}

		logger.info("Sequences: " + nrSequences);
		logger.info("Aligned reads: " + nrAligned);
		logger.info("NonAligned reads: " + nrNonAligned);
		logger.info("NrTotal reads: " + nrTotal);
		logger.info("Readgroups: " + header.getReadGroups().size());

	}

	public SAMRecordIterator iterator() {
		return reader.iterator();
	}

	public SAMRecordIterator query(String chr, int positionBegin, int positionEnd, boolean overlapRequired) {
		SAMRecordIterator iterator = reader.query(chr, positionBegin, positionEnd, overlapRequired);
		return iterator;
	}


	public BrowseableBAMIndex getIndexMetaData() {
		return reader.indexing().getBrowseableIndex();
	}

	public SAMFileHeader getHeader() {
		return reader.getFileHeader();
	}

	public File getFile() {
		return bamFile;
	}

	public void close() throws IOException {
		reader.close();
	}

	public List<SAMReadGroupRecord> getReadGroups() {
		return reader.getFileHeader().getReadGroups();
	}

}
