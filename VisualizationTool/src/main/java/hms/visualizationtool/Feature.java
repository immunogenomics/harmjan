/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.visualizationtool;

import java.util.Objects;

/**
 *
 * @author Harm-Jan
 */
public class Feature {

    private final Chromosome chr;
    private final Strand str;
    private final Track track;
    private final int start;
    private final int stop;
    private int peakPos;
    private double foldChange;
    private double peakQ;

    public Feature(Chromosome chr, Strand str, Track track, int start, int stop) {
        this.chr = chr;
        this.str = str;
        this.track = track;
        this.start = start;
        this.stop = stop;
    }

    public int getStop() {
        return stop;
    }

    public Strand getStrand() {
        return str;
    }

    public int getStart() {
        return start;
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 17 * hash + Objects.hashCode(this.chr);
        hash = 17 * hash + Objects.hashCode(this.str);
        hash = 17 * hash + Objects.hashCode(this.track);
        hash = 17 * hash + this.start;
        hash = 17 * hash + this.stop;
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Feature other = (Feature) obj;
        if (this.chr != other.chr) {
            return false;
        }
        if (this.str != other.str) {
            return false;
        }
        if (!Objects.equals(this.track, other.track)) {
            return false;
        }
        if (this.start != other.start) {
            return false;
        }
        return this.stop == other.stop;
    }

    public int compareTo(Object obj) {
        if (this.equals(obj)) {
            return 0;
        }
        if (getClass() != obj.getClass()) {
            return -1;
        }

        final Feature other = (Feature) obj;
        if (!Objects.equals(this.track, other.track)) {
            return -1;
        }
        if (this.chr.getNumber() > other.chr.getNumber()) {
            return 1;
        } else if (this.chr.getNumber() < other.chr.getNumber()) {
            return -1;
        } else {
            if (this.str == other.str) {
                if (this.start > other.start) {
                    return 1;
                } else if (this.start == other.start) {
                    if (this.stop == other.stop) {
                        return 0;
                    } else {
                        if (this.stop < other.stop) {
                            return -1;
                        } else {
                            return 1;
                        }
                    }
                } else {
                    return -1;
                }
            } else {
                if (this.str == Strand.POS && other.str == Strand.NEG) {
                    return 1;
                } else {
                    return -1;
                }
            }
        }

    }

    public boolean overlaps(Feature that) {
        if (this.chr != that.chr) {
            return false;
        }

        if (this.start >= that.start && this.stop <= that.stop) {
            return true;
        }
        if (this.start >= that.start && this.start < that.stop) {
            return true;
        }
        return this.stop > that.start && this.stop < that.stop;
    }

    public boolean overlapsStrand(Feature that) {
        if (this.str != that.str) {
            return false;
        }
        return overlaps(that);
    }

    public int compare(Feature otherFeature) {
        FeatureComparator comp = new FeatureComparator();
        return comp.compare(this, otherFeature);
    }

    void setPeakValues(int peakPos, double foldChange, double peakQ) {
        this.peakPos = peakPos;
        this.foldChange = foldChange;
        this.peakQ = peakQ;
    }

    double getPeakQ(){
        return peakQ;
    }
    
    double getFoldChange(){
        return foldChange;
    }
    
    int getPeakPos(){
        return peakPos;
    }
}
