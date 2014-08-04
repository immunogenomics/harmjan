/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.utilities.bedfile;

import hms.hwestra.utilities.features.Track;
import hms.hwestra.utilities.features.Chromosome;
import hms.hwestra.utilities.features.Feature;
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
        if (chr != null) {
            System.out.println("Filtering for Chr: " + chr.getName() + "\t " + start + " - " + stop);
        }

        // chr1	8128340	8128539	C011PABXX110504:4:2203:14692:158380	0	-
        String[] elems = tf.readLineElems(TextFile.tab);
        Track track = new Track(name, start, stop);
        int nrReads = 0;
        long length = 0;
        while (elems != null) {

            int len = 0;
            Chromosome featureChr = Chromosome.parseChr(elems[0]);

            Strand featureStrand = Strand.NA;

            featureStrand = Strand.parseStr(elems[elems.length - 1]);

            int featureStart = -1;
            int featureStop = -1;
            try {
                featureStart = Integer.parseInt(elems[1]);
            } catch (NumberFormatException e) {
                System.out.println("Could not parse chromosome start position: " + elems[1]);
            }

            try {
                featureStop = Integer.parseInt(elems[2]);
            } catch (NumberFormatException e) {
                System.out.println("Could not parse chromosome stop position: " + elems[2]);
            }

            int peakPos = 0;
            double peakQ = 0;
            double foldChange = 0;
            if (peakFiles && elems.length > 4) {
                try {
                    peakPos = Integer.parseInt(elems[4]);
                    peakQ = Double.parseDouble(elems[elems.length - 1]);
                    foldChange = Double.parseDouble(elems[elems.length - 2]);
                } catch (NumberFormatException e) {
                    System.out.println("Could not parse chromosome stop position: " + elems[elems.length - 1] + " or " + elems[elems.length - 2]);

                }
            }

            len = featureStop - featureStart;
            boolean add = false;
            if (chr == null) {
                add = true;
            } else if (chr == featureChr) {

                if (featureStart <= start && featureStop <= stop && featureStop >= start) {
                    add = true;
                }
                if (featureStart >= start && featureStop >= stop && featureStart <= stop) {
                    add = true;
                }
                if (featureStart >= start && featureStop <= stop) {
                    add = true;
                }
            }
            if (add) {
                Feature f = new Feature();
                f.setChromosome(featureChr);
                f.setStrand(featureStrand);
                f.setStart(featureStart);
                f.setStop(featureStop);
                if (peakFiles) {
//                        f.setPeakValues(peakPos, foldChange, peakQ);
                }
//                if (track.containsFeature(f)) {
//                    System.out.println("Duplicate feature: " + f.toString());
//                } else {
                track.addFeature(f);
//                }
                // System.out.println(f.toString());
            }

            length += len;
            nrReads++;
            elems = tf.readLineElems(TextFile.tab);
        }

        tf.close();

        System.out.println(
                "Average length: " + ((double) length / nrReads) + "\tNumber of reads: " + nrReads
        );
        track.printNrFeatures();
        return track;
    }

}
