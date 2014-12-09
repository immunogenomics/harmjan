/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.polypeak;

import hms.hwestra.utilities.bamfile.BamFileReader;
import hms.hwestra.utilities.features.Chromosome;
import java.io.File;
import java.io.IOException;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import java.util.List;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import umcg.genetica.util.RunTimer;

/**
 *
 * @author Harm-Jan
 */
public class PolyPeak {

    private static final Logger LOG = LogManager.getLogger();
    

    String tmpFile = null;

    public static void main(String[] args) {
        String file = "/Data/ATAC-seq/GSE47753/SRR891269/SRR891269-sort-dedup-sort.bam";
        PolyPeak p = new PolyPeak(file);

        
    }

    int windowsize = 10000;

    public PolyPeak(String file) {
        this.tmpFile = file;
        try {
            run();
        } catch (IOException ex) {
            LOG.error(ex.getMessage());
            System.exit(-1);
        }
        System.exit(0);
    }

    private void run() throws IOException {

        // set up windows
        // we should know the sizes of each of the chromosomes beforehand
        // they are hardcoded right now
        // keep a running total coverage distribution (e.g. coverage vs count)
        //
        // for a given sample
        // -- open bam file reader
        BamFileReader reader = new BamFileReader(new File(tmpFile), null, false);
            // -- keep a distribution of d versus correlation 500x100

        // determine the reference chromosome names +  from the header
        // some constants
        int basequalthreshold = 10;
        int windowSize = 10000;
        int mapQThreshold = 20;

        boolean pairedReadData = true;
        int windowOverlap = windowSize / 10;
        int maxD = windowSize / 2;
        int[][] crossCorrelationDistribution = new int[maxD][100]; // [d][corr]

        int windowSizeWithOverlap = windowSize + (2 * windowOverlap);

        // iterate chromosomes
        Chromosome[] allChromosomes = Chromosome.values();
        long basesprocessed = 0;
        long fragmentsprocessed = 0;
        RunTimer timer = new RunTimer();
        int maxInsertSize = 5000;
        int maxReadLen = 200;
        for (Chromosome chr : allChromosomes) {
            LOG.info("Processing: " + chr.toString());
            int chromosomeLength = chr.getLength();
            int windowStart = 0;
            while (windowStart < chromosomeLength) {
                // get all reads in window
                // for a given window
                int windowStop = windowStart + windowSize;
                int nA = 0;
                int nT = 0;
                int nC = 0;
                int nG = 0;
                int nN = 0;

                SAMRecordIterator recordIterator = null;
                try {
                    recordIterator = reader.query("" + 1, windowStart - windowOverlap, windowStop + windowOverlap, true);
                } catch (IllegalStateException e) {
                    LOG.error(e.getMessage());
                }

                int nrReadsWindowTotal = 0;
                if (recordIterator != null) {
                    // -- determine read coverage per strand
                    int[][] coverage = new int[2][windowSizeWithOverlap]; // [strand][position] 
                    int[][] insertLengthDist = new int[2][maxInsertSize];
                    int[][] readLengthDistribution = new int[2][maxReadLen];

                    //    -- get read
                    int actualWindowStart = windowStart - windowOverlap;
                    while (recordIterator.hasNext()) {
                        //    -- get mapq > 20
                        int matchingBases = 0;
                        SAMRecord record = recordIterator.next();
                        int mapQ = record.getMappingQuality();
                        int insertSize = record.getInferredInsertSize();
                        if (mapQ > mapQThreshold && insertSize > 0) {
                            //    -- if paired end, check whether other pair is present
                            if (!pairedReadData || properPairFilter(record)) {
                                //    -- parse cigar
                                nrReadsWindowTotal++;
                                Cigar cigar = record.getCigar();
                                List<CigarElement> cigarElements = cigar.getCigarElements();

                                int mapPos = record.getAlignmentStart();
                                int windowRelativePosition = mapPos - actualWindowStart;

                                boolean negStrandFlag = record.getReadNegativeStrandFlag();
                                int strand = 0;
                                if (negStrandFlag) {
                                    strand = 1;
                                }
                                int readPosition = 0;

                                byte[] baseQual = record.getBaseQualities();
                                byte[] bases = record.getReadBases();
                                if (insertSize > maxInsertSize) {
                                    insertLengthDist[strand][maxInsertSize - 1]++;
                                } else {
                                    insertLengthDist[strand][insertSize]++;
                                }

                                fragmentsprocessed++;
                                for (CigarElement e : cigarElements) {
                                    int cigarElementLength = e.getLength();

                                    switch (e.getOperator()) {
                                        case H:
                                            break;
                                        case P:
                                            break;
                                        case S: // soft clip
                                            readPosition += cigarElementLength;
                                            break;
                                        case N: // ref skip
                                            windowRelativePosition += cigarElementLength;
                                            break;
                                        case D: // deletion
                                            windowRelativePosition += cigarElementLength;
                                            break;
                                        case I: // insertion
                                            windowRelativePosition += cigarElementLength;
                                            break;
                                        case M:
                                        case EQ:
                                        case X:
                                            int endPosition = readPosition + cigarElementLength;
                                            for (int pos = readPosition; pos < endPosition; pos++) {
                                                byte base = bases[pos];
                                                boolean properbase = false;
                                                if (windowRelativePosition >= 0) { // the read could overlap the leftmost edge of this window
                                                    if (base == 78 || base == 110) {
                                                        nN++;
                                                    }

                                                    //    -- for each base pos: check whether basequal > 50
                                                    if (baseQual[readPosition] > basequalthreshold) {
                                                        //    -- determine number of A/T/C/G/N bases
                                                        if (base == 64 || base == 97) {
                                                            nA++;
                                                            properbase = true;
                                                            matchingBases++;
                                                        } else if (base == 67 || base == 99) {
                                                            nC++;
                                                            properbase = true;
                                                            matchingBases++;
                                                        } else if (base == 71 || base == 103) {
                                                            nG++;
                                                            properbase = true;
                                                            matchingBases++;
                                                        } else if (base == 84 || base == 116) { // extend to capture U?
                                                            nT++;
                                                            properbase = true;
                                                            matchingBases++;
                                                        }
                                                    }

                                                }
                                                if (properbase) {
                                                    coverage[strand][windowRelativePosition]++;
                                                    basesprocessed++;
                                                }
                                            } // if pos < readposition
                                            readPosition += cigarElementLength;
                                            break;
                                        default:
                                            LOG.warn("Unknown CIGAR operator found: " + e.getOperator().toString());
                                            LOG.warn("In read: " + record.toString());
                                            break;
                                    } // switch operator
                                } // for each cigarelement
                                readLengthDistribution[strand][matchingBases]++;
                            } // if paired read or properpair 
                        } // if mapq > threshold

                    } // for each samrecord
                    if (fragmentsprocessed % 100000 == 0) {
                        long diff = timer.getTimeDiff();
                        long seconds = diff / 1000;
                        double basesPerSecond = fragmentsprocessed / seconds;
                        LOG.info(fragmentsprocessed + " fragments processed in " + timer.getTimeDesc(diff) + " - " + basesPerSecond + " bases/second");
                    }
                    recordIterator.close();
                    // determine mean insert length and variance
                    double meanInsertLengthPosStrand = JSci.maths.ArrayMath.mean(insertLengthDist[0]);
                    double varInsertPosStrand = JSci.maths.ArrayMath.variance(insertLengthDist[0]);

                    double meanInsertLengthNegStrand = JSci.maths.ArrayMath.mean(insertLengthDist[1]);
                    double varInsertNegStrand = JSci.maths.ArrayMath.variance(insertLengthDist[1]);

                    // determine mean read length and variance
                    double meanReadLengthPosStrand = JSci.maths.ArrayMath.mean(readLengthDistribution[0]);
                    double varReadPosStrand = JSci.maths.ArrayMath.variance(readLengthDistribution[0]);

                    double meanReadLengthNegStrand = JSci.maths.ArrayMath.mean(readLengthDistribution[1]);
                    double varReadNegStrand = JSci.maths.ArrayMath.variance(readLengthDistribution[1]);

                    // determine mean coverage and variance
                    double meanCoveragePosStrand = JSci.maths.ArrayMath.mean(coverage[0]);
                    double varCoveragePosStrand = JSci.maths.ArrayMath.variance(coverage[0]);

                    double meanCoverageNegStrand = JSci.maths.ArrayMath.mean(coverage[1]);
                    double varCoverageNegStrand = JSci.maths.ArrayMath.variance(coverage[1]);

                    LOG.info("Chromosome position: " + chr.getName() + ":" + windowStart + "-" + windowStop + "\tNr reads: " + nrReadsWindowTotal);
                    LOG.debug("---");
                    LOG.debug("QTY\t\tPositive Strand\tNegative Strand");
                    LOG.debug("Mean insert size:\t" + meanInsertLengthPosStrand + "\t" + meanInsertLengthNegStrand);
                    LOG.debug("Variance insert:\t" + varInsertPosStrand + "\t" + varInsertNegStrand);
                    LOG.debug("Mean read length:\t" + meanReadLengthPosStrand + "\t" + meanReadLengthNegStrand);
                    LOG.debug("Variance read len:\t" + varReadPosStrand + "\t" + varReadNegStrand);
                    LOG.debug("Mean coverage:\t" + meanCoveragePosStrand + "\t" + meanCoverageNegStrand);
                    LOG.debug("Variance coverage:\t" + varCoveragePosStrand + "\t" + varCoverageNegStrand);
                    LOG.debug("---");
                    long sum = nA + nC + nT + nG + nN;
                    LOG.debug("nrA: " + nA);
                    LOG.debug("nrC: " + nC);
                    LOG.debug("nrT: " + nT);
                    LOG.debug("nrG: " + nG);
                    LOG.debug("nrN: " + nN);
                    LOG.debug("sum: " + sum);
                    LOG.debug("");

                }

                windowStart += windowSize;
            } // while windowstart < chr size
        } // for all chromsoomes

        // -- over a random selection of windows:
        //    -- determine correlation between strands
        //    -- set d += half window size
        //    -- recalculate correlation
        //    -- determine d with highest correlation
        //  
        reader.close();

    }

    private boolean properPairFilter(SAMRecord record) {
        // -- if paired end, check whether other pair is present
        // -- check whether read pair maps to same chr;
        // -- check whether read pair maps to same strand;

        return true;
    }
}
