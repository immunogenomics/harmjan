/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.utilities.bedfile;

import hms.hwestra.utilities.features.Chromosome;
import hms.hwestra.utilities.features.Strand;
import java.io.IOException;
import umcg.genetica.io.text.TextFile;



/**
 *
 * @author Harm-Jan
 */
public class BedFileReader {

    public Track read(String file, String name, Chromosome chr, int start, int stop, boolean peakFiles) throws IOException {
        TextFile tf = new TextFile(file, TextFile.R);

        System.out.println("Reading file: " + file);
        System.out.println("Filtering for Chr: " + chr.getName() + "\t " + start + " - " + stop);

        // chr1	8128340	8128539	C011PABXX110504:4:2203:14692:158380	0	-
        String[] elems = tf.readLineElems(TextFile.tab);
        Track track = new Track(name, start, stop);
        int nrReads = 0;
        long length = 0;
        while (elems != null) {

            int len = 0;
            Chromosome lineChr = Chromosome.parseChr(elems[0]);
            Strand lineStrand = Strand.parseStr(elems[5]);
            int lineStart = -1;
            int lineStop = -1;
            try {
                lineStart = Integer.parseInt(elems[1]);
            } catch (NumberFormatException e) {
                System.out.println("Could not parse chromosome start position: " + elems[1]);
            }

            try {
                lineStop = Integer.parseInt(elems[2]);
            } catch (NumberFormatException e) {
                System.out.println("Could not parse chromosome stop position: " + elems[2]);
            }

            int peakPos = 0;
            double peakQ = 0;
            double foldChange = 0;
            if (peakFiles) {
                try {
                    peakPos = Integer.parseInt(elems[4]);
                    peakQ = Double.parseDouble(elems[elems.length - 1]);
                    foldChange = Double.parseDouble(elems[elems.length - 2]);
                } catch (NumberFormatException e) {
                    System.out.println("Could not parse chromosome stop position: " + elems[elems.length - 1] + " or " + elems[elems.length - 2]);

                }
            }

            len = lineStop - lineStart;

            if (chr == lineChr) {

                if (lineStart <= start && lineStop <= stop && lineStop >= start) {
                    Feature f = new Feature(lineChr, lineStrand, track, lineStart, lineStop);
                    if(peakFiles){
                        f.setPeakValues(peakPos, foldChange, peakQ);
                    }
                    track.addFeature(f);
                }
                if (lineStart >= start && lineStop >= stop && lineStart <= stop) {
                    Feature f = new Feature(lineChr, lineStrand, track, lineStart, lineStop);
                    f.setPeakValues(peakPos, foldChange, peakQ);
                    track.addFeature(f);
                }
                if (lineStart >= start && lineStop <= stop) {
                    Feature f = new Feature(lineChr, lineStrand, track, lineStart, lineStop);
                    f.setPeakValues(peakPos, foldChange, peakQ);
                    track.addFeature(f);
                }
            }
            length += len;
            nrReads++;
            track.addReadLengthToDist(len, lineChr);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println("Average length: " + ((double) length / nrReads) + "\tNumber of reads: " + nrReads);
        track.printNrFeatures();
        return track;
    }

}
