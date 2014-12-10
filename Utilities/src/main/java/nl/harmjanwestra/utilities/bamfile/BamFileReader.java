/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.utilities.bamfile;

import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.BAMIndexer;
import htsjdk.samtools.BrowseableBAMIndex;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import umcg.genetica.containers.Pair;

/**
 *
 * @author Harm-Jan
 */
public class BamFileReader {

    private final File bamFile;
    private final SamReader reader;
    private final AggregateFilter filter;
    private static final Logger logger = LogManager.getLogger();

    public BamFileReader(File bamFile, AggregateFilter filter, boolean indexFileWhenNotIndexed) throws IOException {
        this.filter = filter;
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
        logger.info("Sequences: " + header.getSequenceDictionary().getSequences().size());
        logger.info("");
        int nrAligned = 0;
        int nrNonAligned = 0;
        int nrTotal = 0;
        for (SAMSequenceRecord record : header.getSequenceDictionary().getSequences()) {
//
            logger.debug("Sequence: " + record.getSequenceName());
            logger.debug("Sequence index: " + record.getSequenceIndex());
            logger.debug("Sequence len: " + record.getSequenceLength());
            BAMIndexMetaData metaData = reader.indexing().getIndex().getMetaData(record.getSequenceIndex());
            int tmpnrAligned = metaData.getAlignedRecordCount();
            int tmpnrNonAligned = metaData.getAlignedRecordCount();
            int tmpnrTotal = tmpnrAligned + tmpnrNonAligned;
            
            nrAligned += tmpnrAligned;
            nrNonAligned += tmpnrNonAligned;
            nrTotal += tmpnrTotal;
            logger.debug("Aligned: " + tmpnrAligned);
            logger.debug("NonAligned: " + tmpnrNonAligned);
            logger.debug("NrTotal: " + tmpnrTotal);

            
        }

        logger.info("");

        logger.info("Total Aligned: " + nrAligned);
        logger.info("Total NonAligned: " + nrNonAligned);
        logger.info("Total NrTotal: " + nrTotal);
        logger.info("");

    }

    public SAMRecordIterator iterator() {
        return reader.iterator();
    }

    public SAMRecordIterator query(String chr, int positionBegin, int positionEnd, boolean overlapRequired) {
        SAMRecordIterator iterator = reader.query(chr, positionBegin, positionEnd, overlapRequired);
        return iterator;
    }

    public Pair<Integer, Integer> getNrOfAlignedReadsPassingFilter() {
        SAMRecordIterator iterator = reader.iterator();
        int nrPassingFilter = 0;
        int nrNotPassingFilter = 0;
        while (iterator.hasNext()) {
            SAMRecord record = iterator.next();
            if (filter == null || !filter.filterOut(record)) {
                nrPassingFilter++;
            } else {
                nrNotPassingFilter++;
            }
        }

        return new Pair<Integer, Integer>(nrPassingFilter, nrNotPassingFilter);
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
}
