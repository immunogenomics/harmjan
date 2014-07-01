/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.utilities.bedfile;

import hms.hwestra.utilities.features.Chromosome;
import hms.hwestra.utilities.features.Strand;
import java.util.NavigableSet;
import java.util.TreeSet;

/**
 *
 * @author Harm-Jan
 */
public class Track {

    private final TreeSet<BedFileFeature> features;
    private final String name;
    private final int start;
    private final int stop;

    private final int[][] readLengthDist = new int[25][500];

    public Track(String name, int start, int stop) {
        this.name = name;
        this.features = new TreeSet<>(new BedFileFeatureComparator());
        this.start = start;
        this.stop = stop;
    }

    public void addFeature(BedFileFeature f) {
        features.add(f);
    }

    public void addReadLengthToDist(int len) {

    }

    public void addReadLengthToDist(int len, Chromosome lineChr) {
        if (len >= 500) {
            len = 499;
        }
        readLengthDist[lineChr.getNumber()][len]++;
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    public int getStart() {
        return start;
    }

    public int getStop() {
        return stop;
    }

    public Iterable<BedFileFeature> getFeatures() {
        return features;
    }

    public int getNrReads() {
        return features.size();
    }

    public NavigableSet<BedFileFeature> getFeatureSet(Chromosome chr, Strand strand, int start, int end) {
        BedFileFeature left = new BedFileFeature(chr, strand, null, start, start);
        BedFileFeature right = new BedFileFeature(chr, strand, null, end, end);
        return features.subSet(left, true, right, true);

    }

    public void printNrFeatures() {
        System.out.println(features.size() + " features in track.");
    }

}
