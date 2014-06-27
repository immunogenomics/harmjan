/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.utilities.features;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Objects;

/**
 *
 * @author Harm-Jan
 */
public class Gene {

    Chromosome chromosome;
    String name;
    Strand strand;

    ArrayList<Transcript> transcripts;

    int start = Integer.MAX_VALUE;
    int stop = -Integer.MAX_VALUE;

    public Gene(String name, Chromosome chromosome, Strand strand) {
        this.name = name;
        this.chromosome = chromosome;
        this.strand = strand;
    }

    public void addTranscript(Transcript t) {
        if (transcripts == null) {
            transcripts = new ArrayList<Transcript>();
        }
        transcripts.add(t);
        if(t.getStart() < start){
            start = t.getStart();
        }
        if(t.getStop() > stop){
            stop = t.getStop();
        }
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

    public int getStart() {
        return start;
    }

    public int getStop() {
        return stop;
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 67 * hash + Objects.hashCode(this.chromosome);
        hash = 67 * hash + Objects.hashCode(this.name);
        hash = 67 * hash + Objects.hashCode(this.strand);
        hash = 67 * hash + this.start;
        hash = 67 * hash + this.stop;
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
        final Gene other = (Gene) obj;
        if (this.chromosome != other.chromosome) {
            return false;
        }
        if (!Objects.equals(this.name, other.name)) {
            return false;
        }
        if (this.strand != other.strand) {
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

    public boolean overlaps(Gene that) {
        if (chromosome != that.chromosome) {
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

}
