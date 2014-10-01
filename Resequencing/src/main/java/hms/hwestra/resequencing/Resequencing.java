/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.resequencing;

import htsjdk.samtools.BAMIndexer;
import htsjdk.samtools.MergingSamRecordIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamFileHeaderMerger;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author hwestra
 */
public class Resequencing {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        Resequencing s = new Resequencing(args);
    }

    private final CLI cli;

    private Resequencing(String[] args) {
        cli = new CLI(args);

        // stuff for merger
        // get list of samples to use
        try {
            CoverageMatrix mat = new CoverageMatrix(cli.getCoverageMatrix());

            String indir = cli.getBamFileIndir();
            String outdir = cli.getBamFileOutdir();

            indir = Gpio.formatAsDirectory(indir);
            outdir = Gpio.formatAsDirectory(outdir);
            Gpio.createDir(outdir);
            String tmpdir = outdir + "tmp";
            tmpdir = Gpio.formatAsDirectory(tmpdir);
            Gpio.createDir(tmpdir);
            String suffix = cli.getBamSuffix();

            ArrayList<Pair<String, ArrayList<String>>> samples = mat.selectSamples(cli.getCoverageThresholdColumn(), cli.getCoverageThreshold(), outdir);
            System.out.println(samples.size() + " samples selected");
            TextFile outfile = new TextFile(outdir + "SamplesToUse.txt", TextFile.W);

            for (Pair<String, ArrayList<String>> p : samples) {

                String sampleName = p.getLeft();
                ArrayList<String> dups = p.getRight();

                if (dups.size() > 1) {
                    // check whether files are actually there..
                    int d = 0;
                    ArrayList<File> finalFiles = new ArrayList<File>();
                    System.out.println(sampleName + "\thas " + dups.size() + " duplicates");
                    for (String f : dups) {
                        if (Gpio.exists(indir + f + "." + suffix)) {
                            d++;
                            finalFiles.add(new File(indir + f + "." + suffix));
                        } else {
                            System.out.println(sampleName + "\tFile not found: " + indir + f + "." + suffix);
                        }
                    }
                    System.out.println(sampleName + "\tFinal duplication size: " + finalFiles.size());
                    if (finalFiles.size() > 1) {

                        // merge bam files
                        File outputFile = new File(outdir + sampleName + "_merged.bam");
                        System.out.println(sampleName + "\tMerging BAM files");
                        mergeBamFiles(sampleName, finalFiles, outputFile);
                        // dedup

                        File deduppedOut = new File(outdir + sampleName + "_merged.dedup.bam");
                        File deduppedMetricsOut = new File(outdir + sampleName + "_merged.dedupmetrics.txt");
                        System.out.println(sampleName + "\tDeduplicating alignments in BAM file" + outputFile.getAbsolutePath());
                        System.out.println(sampleName + "\tWriting to: " + deduppedOut.getAbsolutePath());
                        MarkDups md = new MarkDups(outputFile, deduppedOut, deduppedMetricsOut, new File(tmpdir));
// md.go();
                        // sort file
                        File deduppedSortOut = new File(outdir + sampleName + "_merged.dedup.sort.bam");
                        SortBAM sb = new SortBAM(deduppedOut, deduppedSortOut, new File(tmpdir));

                        // index
                        indexBAM(deduppedSortOut);
                        outfile.writeln(deduppedSortOut.getAbsolutePath());
                    } else {
                        // print filename to file
                        if (finalFiles.size() == 1) {
                            // check whether the file has been indexed..

                            if (!Gpio.exists(finalFiles.get(0).getAbsolutePath() + ".bai")) {
                                System.out.println(sampleName + "\tIndexing BAM: " + finalFiles.get(0).getAbsolutePath());
                                indexBAM(finalFiles.get(0));
                            }
                            outfile.writeln(finalFiles.get(0).getAbsolutePath());
                        }
                    }

                } else {
                    // check whether file is there,
                    String f = dups.get(0);
                    System.out.println(sampleName + "\thas no duplicates.");
                    if (Gpio.exists(indir + f + "." + suffix)) {
                        outfile.writeln(indir + f + "." + suffix);
                        // check whether file has been indexed...
                        if (!Gpio.exists(indir + f + "." + suffix + ".bai")) {
                            System.out.println(sampleName + "\tIndexing BAM: " + indir + f + "." + suffix);
                            indexBAM(new File(indir + f + "." + suffix));
                        }
                        System.out.println(sampleName + "\tFound file: " + indir + f + "." + suffix);
                    } else {
                        System.out.println(sampleName + "\tCould not find file: " + indir + f + "." + suffix);
                    }
                }

            }
            outfile.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private void mergeBamFiles(String sample, ArrayList<File> finalFiles, File outfile) {
        ArrayList<SAMFileReader> readers = new ArrayList<SAMFileReader>();
        ArrayList<SAMFileHeader> headers = new ArrayList<SAMFileHeader>();
        for (File f : finalFiles) {
            System.out.println(sample + "\tIncluding file: " + f.getAbsolutePath());
            SAMFileReader tmpreader = new SAMFileReader(f);
            headers.add(tmpreader.getFileHeader());
            readers.add(tmpreader);
        }

        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate, headers, true);
        MergingSamRecordIterator iterator = new MergingSamRecordIterator(headerMerger, readers, false);

        SAMFileHeader header = headerMerger.getMergedHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        List<SAMReadGroupRecord> var = header.getReadGroups();
        System.out.println(var.size() + " readgroups found");
        for (SAMReadGroupRecord r : var) {
            System.out.println(String.format(sample + "\tDetected read group ID=%s PL=%s LB=%s SM=%s%n", r.getId(), r.getPlatform(), r.getLibrary(), r.getSample()));
        }

        String RGID = sample;
        String RGLB = var.get(0).getLibrary();
        String RGPL = var.get(0).getPlatform();
        String RGSM = sample;
        String RGPU = var.get(0).getPlatformUnit();

        final SAMReadGroupRecord rg = new SAMReadGroupRecord(RGID);
        rg.setLibrary(RGLB);
        rg.setPlatform(RGPL);
        rg.setSample(RGSM);
        rg.setPlatformUnit(RGPU);
        System.out.println(String.format(sample + "\tCreated read group ID=%s PL=%s LB=%s SM=%s%n", rg.getId(), rg.getPlatform(), rg.getLibrary(), rg.getSample()));
        header.setReadGroups(Arrays.asList(rg));

        SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(headerMerger.getMergedHeader(), false, outfile);
        System.out.println(sample + "\tWriting to: " + outfile);
        while (iterator.hasNext()) {
            SAMRecord record = iterator.next();
            record.setAttribute(SAMTag.RG.name(), RGID);
            outputSam.addAlignment(record);

        }

        outputSam.close();

        for (SAMFileReader r : readers) {
            r.close();
        }

        indexBAM(outfile);

    }

    private void sortBAM(String in, String out) {

    }

    private void indexBAM(File in) {
        SAMFileReader reader = new SAMFileReader(in);
        BAMIndexer.createIndex(reader, new File(in.getAbsolutePath() + ".bai"));
    }

}
