/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.polypeak.containers;

import hms.hwestra.utilities.features.Chromosome;

/**
 *
 * @author hwestra
 */
public class Window {

    private int[][] coverage;
    private final int start;
    private final int stop;
    private final Chromosome chr;
    private int maxCoverage = 0;
    private int minCoverage = 0;

    public Window(int start, int stop, Chromosome chr) {
        this.start = start;
        this.stop = stop;
        this.chr = chr;
    }

    public int[][] getCoverage() {
        return coverage;
    }

    public int getStart() {
        return start;
    }

    public int getStop() {
        return stop;
    }

    public Chromosome getChr() {
        return chr;
    }

    public int getMaxCoverage() {
        return maxCoverage;
    }

    public int getMinCoverage() {
        return minCoverage;
    }

}
