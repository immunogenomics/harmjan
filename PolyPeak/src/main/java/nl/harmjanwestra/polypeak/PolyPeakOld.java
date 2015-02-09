/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.polypeak;

import nl.harmjanwestra.utilities.bamfile.BamFileReader;
import nl.harmjanwestra.utilities.features.Chromosome;
import java.io.File;
import java.io.IOException;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
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
public class PolyPeakOld {

    private static final Logger LOG = LogManager.getLogger();
    

    String tmpFile = null;

    public static void main(String[] args) {
        String file = "/Data/ATAC-seq/GSE47753/SRR891269/SRR891269-sort-dedup-sort.bam";
        PolyPeakOld p = new PolyPeakOld(file);

        
    }

    int windowsize = 10000;

    public PolyPeakOld(String file) {
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
        BamFileReader reader = new BamFileReader(new File(tmpFile));
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
