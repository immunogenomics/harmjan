/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.visualizationtool;

import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.Gene;
import java.util.ArrayList;

/**
 *
 * @author hwestra
 */
public class Locus extends Feature {

    ArrayList<Feature> proxies = new ArrayList<Feature>();
    ArrayList<Boolean> proxiesOverlap = new ArrayList<Boolean>();
    ArrayList<Gene> genes = new ArrayList<Gene>();

    double score = 0;
    boolean overlap = false;
    private int leftBound;
    private int rightBound;

    public ArrayList<Feature> getProxies() {
        return proxies;
    }

    public void setProxies(ArrayList<Feature> proxies) {
        this.proxies = proxies;
    }

    public double getScore() {
        return score;
    }

    public void setScore(double score) {
        this.score = score;
    }

    public boolean isOverlap() {
        return overlap;
    }

    public void setOverlap(boolean overlap) {
        this.overlap = overlap;
    }

    public void addProxy(String snp, boolean overlap) {
        Feature f = new Feature();
        f.setName(snp);

        proxies.add(f);
        proxiesOverlap.add(overlap);
    }

    public int getNrOfGenes() {
        return genes.size();
    }

    public void addGene(Gene g) {
        genes.add(g);
    }

    public ArrayList<Gene> getGenes() {
        return genes;
    }

    public ArrayList<ArrayList<Gene>> determineNonOverlappingGenes() {

        ArrayList<ArrayList<Gene>> nonOverLappingGenes = new ArrayList<ArrayList<Gene>>();
        ArrayList<Gene> currentLevel = new ArrayList<Gene>();
        boolean[] hasBeenPlaced = new boolean[genes.size()];
        int maxX = 0;
        int total = 0;
        System.out.println("Determining overlap for: " + genes.size() + " genes... ");
        for (int g = 0; g < genes.size(); g++) {

            if (!hasBeenPlaced[g]) {
                currentLevel = new ArrayList<Gene>();
                currentLevel.add(genes.get(g));
                hasBeenPlaced[g] = true;

                for (int g2 = g + 1; g2 < genes.size(); g2++) {
                    if (!hasBeenPlaced[g2]) {
                        boolean b = determineGeneOverlap(genes.get(g), genes.get(g2));
                        if (!b) {
                            currentLevel.add(genes.get(g2));
                            hasBeenPlaced[g2] = true;
                        }
                    }
                }
                nonOverLappingGenes.add(currentLevel);
                if (currentLevel.size() > maxX) {
                    maxX = currentLevel.size();

                }
                total += currentLevel.size();
            }
        }

        System.out.println("Level of genes: x: " + maxX + " y: " + nonOverLappingGenes.size() + " total: " + total + " original: " + genes.size());
        return nonOverLappingGenes;
    }

    private boolean determineGeneOverlap(Gene g1, Gene g2) {

        return g1.overlaps(g2);

    }

    public void setLeftBound(int maxStart) {
        this.leftBound = maxStart;
    }

    public void setRightBound(int maxStop) {
        this.rightBound = maxStop;
    }

    public ArrayList<Boolean> getProxiesOverlap() {
        return proxiesOverlap;
    }

    public int getLeftBound() {
        return leftBound;
    }

    public int getRightBound() {
        return rightBound;
    }

}
