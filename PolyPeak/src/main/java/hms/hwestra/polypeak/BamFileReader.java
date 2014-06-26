/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.polypeak;

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
import org.apache.log4j.Logger;
import umcg.genetica.containers.Pair;

/**
 *
 * @author Harm-Jan
 */
public class BamFileReader {

    private final File bamFile;
    private final SamReader reader;
    private final AggregateFilter filter;

    private static final Logger LOG = Logger.getLogger(BamFileReader.class);

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
                LOG.warn("Sort order is not based on coordinate, so access may be slow for: " + bamFile.getAbsolutePath());
            }
        } else {
            tmpreader.close();
            LOG.info("BAM file does not have an index: " + bamFile.getAbsolutePath());
            if (indexFileWhenNotIndexed) {
                // index file
                LOG.info("Indexing BAM file: " + bamFile.getAbsolutePath());
                SAMFileReader samfileReader = new SAMFileReader(bamFile);
                String outputPath = bamFile.getAbsolutePath();
//                outputPath = outputPath.substring(3);
                outputPath += ".bai";
                BAMIndexer.createIndex(samfileReader, new File(outputPath));
                samfileReader.close();
                LOG.info("Done indexing: " + outputPath);
            }
            tmpreader = SamReaderFactory.make()
                    .enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES)
                    .validationStringency(ValidationStringency.DEFAULT_STRINGENCY)
                    .open(bamFile);
            reader = tmpreader;

            SAMFileHeader header = reader.getFileHeader();
            LOG.info("Summary:");
            LOG.info("Creator: " + header.getCreator());
            LOG.info("Sequences: " + header.getSequenceDictionary().getSequences().size());
            LOG.info("");
            int nrAligned = 0;
            int nrNonAligned = 0;
            int nrTotal = 0;
            for (SAMSequenceRecord record : header.getSequenceDictionary().getSequences()) {

                LOG.info("Sequence: " + record.getSequenceName());
                LOG.info("\tSequence index: " + record.getSequenceIndex());
                LOG.info("\tSequence len: " + record.getSequenceLength());
                BAMIndexMetaData metaData = reader.indexing().getIndex().getMetaData(record.getSequenceIndex());
                int tmpnrAligned = metaData.getAlignedRecordCount();
                int tmpnrNonAligned = metaData.getAlignedRecordCount();
                int tmpnrTotal = tmpnrAligned + tmpnrNonAligned;
                nrAligned += tmpnrAligned;
                nrNonAligned += tmpnrNonAligned;
                nrTotal += tmpnrTotal;
                LOG.info("\tAligned: " + tmpnrAligned);
                LOG.info("\tNonAligned: " + tmpnrNonAligned);
                LOG.info("\tNrTotal: " + tmpnrTotal + "\n");

            }
            
            System.out.println("");
            
            LOG.info("Total Aligned: " + nrAligned);
            LOG.info("Total NonAligned: " + nrNonAligned);
            LOG.info("Total NrTotal: " + nrTotal);
            LOG.info("");

        }

    }

    public SAMRecord[] query(String chr, int positionBegin, int positionEnd, boolean overlapRequired) {
        SAMRecordIterator iterator = reader.query(chr, positionBegin, positionEnd, overlapRequired);
        ArrayList<SAMRecord> records = new ArrayList<SAMRecord>();
        while (iterator.hasNext()) {
            SAMRecord record = iterator.next();
            if (!filter.filterOut(record)) {
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
            if (!filter.filterOut(record)) {
                nrPassingFilter++;
            } else {
                nrNotPassingFilter++;
            }
        }
        return new Pair<Integer, Integer>(nrPassingFilter, nrNotPassingFilter);
    }
    
    public BrowseableBAMIndex getIndexMetaData(){
        return reader.indexing().getBrowseableIndex();
    }
    
    public SAMFileHeader getHeader(){
        return reader.getFileHeader();
    }

    public File getFile() {
        return bamFile;
    }

    public void close() throws IOException {
        reader.close();
    }
}
