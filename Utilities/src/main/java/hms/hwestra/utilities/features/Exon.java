/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.utilities.features;

import java.util.ArrayList;
import java.util.Objects;

/**
 *
 * @author Harm-Jan
 */
public class Exon {

    private final Chromosome chromosome;
    private final String name;

    private final Strand strand;

    private  ArrayList<Transcript> transcripts;
    private final Gene gene;

    int start = Integer.MAX_VALUE;
    int stop = -Integer.MAX_VALUE;

    public Exon(String name, Chromosome chr, Strand strand, Gene gene, int start, int stop) {
        this.name = name;
        this.chromosome = chr;
        this.strand = strand;
        
        this.gene = gene;
        this.start = start;
        this.stop = stop;
    }

    public Chromosome getChromosome() {
        return chromosome;
    }

    public String getName() {
        return name;
    }

    public Strand getStrand() {
        return strand;
    }

    public ArrayList<Transcript> getTranscripts() {
        return transcripts;
    }

    public Gene getGene() {
        return gene;
    }

    public int getStart() {
        return start;
    }

    public int getStop() {
        return stop;
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 37 * hash + Objects.hashCode(this.chromosome);
        hash = 37 * hash + Objects.hashCode(this.name);
        hash = 37 * hash + Objects.hashCode(this.strand);
        hash = 37 * hash + Objects.hashCode(this.gene);
        hash = 37 * hash + this.start;
        hash = 37 * hash + this.stop;
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
        final Exon other = (Exon) obj;
        if (this.chromosome != other.chromosome) {
            return false;
        }
        if (!Objects.equals(this.name, other.name)) {
            return false;
        }
        if (this.strand != other.strand) {
            return false;
        }
        if (!Objects.equals(this.gene, other.gene)) {
            return false;
        }
        if (this.start != other.start) {
            return false;
        }
        if (this.stop != other.stop) {
            return false;
        }
        return true;
    }

    public void addTranscript(Transcript t) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    
}
