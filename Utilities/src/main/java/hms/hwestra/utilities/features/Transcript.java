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
public class Transcript {

    private final String name;
    private final Chromosome chromosome;
    private final Strand strand;
    private final Gene gene;
    private ArrayList<Exon> exons;

    int start = Integer.MAX_VALUE;
    int stop = -Integer.MAX_VALUE;

    public Transcript(String name, Chromosome chromosome, Strand strand, Gene gene) {
        this.name = name;
        this.chromosome = chromosome;
        this.strand = strand;
        this.gene = gene;

    }

    public String getName() {
        return name;
    }

    public Chromosome getChromosome() {
        return chromosome;
    }

    public Strand getStrand() {
        return strand;
    }

    public Gene getGene() {
        return gene;
    }

    public ArrayList<Exon> getExons() {
        return exons;
    }

    public int getStart() {
        return start;
    }

    public int getStop() {
        return stop;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 47 * hash + Objects.hashCode(this.name);
        hash = 47 * hash + Objects.hashCode(this.chromosome);
        hash = 47 * hash + Objects.hashCode(this.strand);
        hash = 47 * hash + Objects.hashCode(this.gene);
        hash = 47 * hash + this.start;
        hash = 47 * hash + this.stop;
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
        final Transcript other = (Transcript) obj;
        if (!Objects.equals(this.name, other.name)) {
            return false;
        }
        if (this.chromosome != other.chromosome) {
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

    public void addExon(Exon e) {
        int estart = e.getStart();
        int estop = e.getStop();
        if(estart < start){
            start = estart;
        }
        if(estop > stop){
            stop = estop;
        }
        exons.add(e);
    }

}
