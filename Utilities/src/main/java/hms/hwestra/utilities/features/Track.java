/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.utilities.features;

import java.util.NavigableSet;
import java.util.TreeSet;

/**
 *
 * @author Harm-Jan
 */
public class Track extends Feature {

    private TreeSet<Feature> features;

    public Track(String name, int start, int stop) {
        this.name = name;
        this.features = new TreeSet<Feature>(new FeatureComparator(false));
        this.start = start;
        this.stop = stop;
    }

    public void addFeature(Feature f) {
        features.add(f);
    }

    public Iterable<Feature> getFeatures() {
        return features;
    }

    public int getNrReads() {
        return features.size();
    }

    public NavigableSet<Feature> getFeatureSet(Chromosome chr, int start, int end) {
        Feature left = new Feature();
        Feature right = new Feature();
        left.setChromosome(chr);
        left.setStart(start);
        left.setStop(start);
        left.setStrand(Strand.POS);
        right.setChromosome(chr);
        right.setStart(end);
        right.setStop(end);
        right.setStrand(Strand.NEG);
        NavigableSet<Feature> set = features.subSet(left, true, right, true);
        System.out.println(set.size() + "\t" + left.toString() + "\t" + right.toString());
        return set;
    }

    public void printNrFeatures() {
        System.out.println(features.size() + " features in track.");
    }

    public void setFeatures(NavigableSet<Feature> f) {
        this.features = new TreeSet<Feature>(new FeatureComparator(false));
        features.addAll(f);
    }

    public Track getSubset(Chromosome chr, int start, int stop) {
        Track t = new Track(this.name, start, stop);
        NavigableSet<Feature> set = getFeatureSet(chr, start, stop);

        t.setFeatures(set);
        return t;
    }

    public int getNrFeatures() {
        return features.size();
    }

    public boolean containsFeature(Feature f) {
        return features.contains(f);
    }
    
    public void addFeatures(Track t){
        features.addAll(t.features);
    }

}
