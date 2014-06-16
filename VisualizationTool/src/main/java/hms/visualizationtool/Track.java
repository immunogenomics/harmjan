/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.visualizationtool;

import java.util.NavigableSet;
import java.util.TreeSet;

/**
 *
 * @author Harm-Jan
 */
public class Track {

    private final TreeSet<Feature> features;
    private final String name;
    private final int start;
    private final int stop;

    private final int[][] readLengthDist = new int[25][500];

    public Track(String name, int start, int stop) {
        this.name = name;
        this.features = new TreeSet<>(new FeatureComparator());
        this.start = start;
        this.stop = stop;
    }

    public void addFeature(Feature f) {
        features.add(f);
    }

    void addReadLengthToDist(int len) {

    }

    void addReadLengthToDist(int len, Chromosome lineChr) {
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

    Iterable<Feature> getFeatures() {
        return features;
    }

    int getNrReads() {
        return features.size();
    }

    public NavigableSet<Feature> getFeatureSet(Chromosome chr, Strand strand, int start, int end) {
        Feature left = new Feature(chr, strand, null, start, start);
        Feature right = new Feature(chr, strand, null, end, end);
        return features.subSet(left, true, right, true);

    }

    void printNrFeatures() {
        System.out.println(features.size() + " features in track.");
    }

}
