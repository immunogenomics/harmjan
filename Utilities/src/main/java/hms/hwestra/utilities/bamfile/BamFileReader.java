/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.utilities.bamfile;

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
import java.util.logging.Level;

import umcg.genetica.containers.Pair;

/**
 *
 * @author Harm-Jan
 */
public class BamFileReader {

    private final File bamFile;
    private final SamReader reader;
    private final AggregateFilter filter;
    static final java.util.logging.Logger LOG = java.util.logging.Logger.getLogger(BamFileReader.class.getName());
//    private static final Logger LOG = Logger.getLogger(BamFileReader.class);

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
                LOG.log(Level.WARNING, "Sort order is not based on coordinate, so access may be slow for: {0}", bamFile.getAbsolutePath());
            }
        } else {
            tmpreader.close();
            LOG.log(Level.INFO, "BAM file does not have an index: {0}", bamFile.getAbsolutePath());
            if (indexFileWhenNotIndexed) {
                // index file
                LOG.log(Level.INFO, "Indexing BAM file: {0}", bamFile.getAbsolutePath());
                SAMFileReader samfileReader = new SAMFileReader(bamFile);
                String outputPath = bamFile.getAbsolutePath();
//                outputPath = outputPath.substring(3);
                outputPath += ".bai";
                BAMIndexer.createIndex(samfileReader, new File(outputPath));
                samfileReader.close();
                LOG.log(Level.INFO, "Done indexing: {0}", outputPath);
            }
            tmpreader = SamReaderFactory.make()
                    .enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES)
                    .validationStringency(ValidationStringency.DEFAULT_STRINGENCY)
                    .open(bamFile);
            reader = tmpreader;

        }

        SAMFileHeader header = reader.getFileHeader();
        LOG.info("Summary:");
        LOG.log(Level.INFO, "File name: "+bamFile.getAbsolutePath());
        LOG.log(Level.INFO, "Creator: {0}", header.getCreator());
        LOG.log(Level.INFO, "Sequences: {0}", header.getSequenceDictionary().getSequences().size());
        LOG.info("");
        int nrAligned = 0;
        int nrNonAligned = 0;
        int nrTotal = 0;
        for (SAMSequenceRecord record : header.getSequenceDictionary().getSequences()) {
//
            LOG.log(Level.INFO, "Sequence: {0}", record.getSequenceName());
            LOG.log(Level.INFO, "\tSequence index: {0}", record.getSequenceIndex());
            LOG.log(Level.INFO, "\tSequence len: {0}", record.getSequenceLength());
            BAMIndexMetaData metaData = reader.indexing().getIndex().getMetaData(record.getSequenceIndex());
            int tmpnrAligned = metaData.getAlignedRecordCount();
            int tmpnrNonAligned = metaData.getAlignedRecordCount();
            int tmpnrTotal = tmpnrAligned + tmpnrNonAligned;
            nrAligned += tmpnrAligned;
            nrNonAligned += tmpnrNonAligned;
            nrTotal += tmpnrTotal;
            LOG.log(Level.INFO, "\tAligned: {0}", tmpnrAligned);
            LOG.log(Level.INFO, "\tNonAligned: {0}", tmpnrNonAligned);
            LOG.log(Level.INFO, "\tNrTotal: {0}\n", tmpnrTotal);

        }

        System.out.println("");

        LOG.log(Level.INFO, "Total Aligned: {0}", nrAligned);
        LOG.log(Level.INFO, "Total NonAligned: {0}", nrNonAligned);
        LOG.log(Level.INFO, "Total NrTotal: {0}", nrTotal);
        LOG.info("");

    }

    public SAMRecordIterator iterator() {
        return reader.iterator();
    }

    public SAMRecord[] query(String chr, int positionBegin, int positionEnd, boolean overlapRequired) {
        SAMRecordIterator iterator = reader.query(chr, positionBegin, positionEnd, overlapRequired);

        ArrayList<SAMRecord> records = new ArrayList<SAMRecord>();
        while (iterator.hasNext()) {
            SAMRecord record = iterator.next();
            if (filter == null || !filter.filterOut(record)) {
                records.add(record);
            }
        }
        return records.toArray(new SAMRecord[0]);
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
