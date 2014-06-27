/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.rnaaccess.probeannotation;

import hms.hwestra.utilities.features.Chromosome;
import java.util.Objects;

/**
 *
 * @author Harm-Jan
 */
public class Probe {

    private final Chromosome chromosome;
    private final int start;
    private final int stop;
    private final String name;

    public Probe(Chromosome chromosome, int start, int stop, String name) {
        this.chromosome = chromosome;
        this.start = start;
        this.stop = stop;
        this.name = name;
    }

    public Chromosome getChromosome() {
        return chromosome;
    }

    public int getStart() {
        return start;
    }

    public int getStop() {
        return stop;
    }

    public String getName() {
        return name;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 19 * hash + Objects.hashCode(this.chromosome);
        hash = 19 * hash + this.start;
        hash = 19 * hash + this.stop;
        hash = 19 * hash + Objects.hashCode(this.name);
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
        final Probe other = (Probe) obj;
        if (this.chromosome != other.chromosome) {
            return false;
        }
        if (this.start != other.start) {
            return false;
        }
        if (this.stop != other.stop) {
            return false;
        }
        if (!Objects.equals(this.name, other.name)) {
            return false;
        }
        return true;
    }
    
    

    public boolean overlaps(Probe that) {
        if (chromosome != that.getChromosome()) {
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
