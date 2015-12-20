/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.utilities.features;

import java.util.ArrayList;
import java.util.NavigableSet;
import java.util.TreeSet;

/**
 * @author Harm-Jan
 */
public class Track extends Feature {

	private TreeSet<Feature> features;
	private ArrayList<Feature> allFeatures;

	public Track(String name) {
		this(name, 0, Integer.MAX_VALUE);
	}

	public Track(String name, int start, int stop) {
		this.name = name;
		this.features = new TreeSet<Feature>(new FeatureComparator(false));
		this.allFeatures = new ArrayList<Feature>();
		this.start = start;
		this.stop = stop;
	}

	public void addFeature(Feature f) {
		features.add(f);
	}

	public Iterable<Feature> getFeatures() {
		return features;
	}

	public ArrayList<Feature> getAllFeatures() {
		return allFeatures;
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
//		System.out.println(set.size() + "\t" + left.toString() + "\t" + right.toString());
		return set;
	}

	public void printNrFeatures() {
		System.out.println(features.size() + " features in track.");
	}

	public void setFeatures(NavigableSet<Feature> f) {
		this.features = new TreeSet<Feature>(new FeatureComparator(false));
		this.allFeatures = new ArrayList<Feature>();
		features.addAll(f);
		allFeatures.addAll(f);
	}

	public Track getSubset(Chromosome chr, int start, int stop) {
		Track t = new Track(this.name, start, stop);
		NavigableSet<Feature> set = getFeatureSet(chr, start, stop);

		t.setFeatures(set);
		return t;
	}

	public int getNrUniqueFeatures() {
		return features.size();
	}

	public int getNrFeatures() {
		return allFeatures.size();
	}

	public boolean containsFeature(Feature f) {
		return features.contains(f);
	}

	public void addFeatures(Track t) {
		features.addAll(t.features);
		allFeatures.addAll(t.features);
	}

	public void addFeatures(ArrayList<Feature> features) {
		for (Feature f : features) {
			this.features.add(f);
		}
		allFeatures.addAll(features);
	}

	public NavigableSet<Feature> getFeatureSet(Feature feat1) {
		return getFeatureSet(feat1.getChromosome(), feat1.getStart(), feat1.getStop());
	}
}
