/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.utilities.features;

import nl.harmjanwestra.utilities.enums.Chromosome;

/**
 *
 * @author Harm-Jan
 */
public class Probe extends Feature {

    public Probe(Chromosome chromosome, int start, int stop, String name) {
        this.chromosome = chromosome;
        this.start = start;
        this.stop = stop;
        this.name = name;
    }

    @Override
    public String toString() {
        return "Probe{" + "chromosome=" + chromosome + ", start=" + start + ", stop=" + stop + ", name=" + name + '}';
    }

}
