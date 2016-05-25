/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.utilities.features;

import java.util.ArrayList;

/**
 *
 * @author Harm-Jan
 */
public class Transcript extends Feature {

    private final Gene gene;
    private ArrayList<Exon> exons;

    public Transcript(String name, Chromosome chromosome, Strand strand, Gene gene) {
        this.name = name;
        this.chromosome = chromosome;
        this.strand = strand;
        this.gene = gene;
    }

    public Transcript(String name, Chromosome chromosome, Strand strand, Gene gene, int start, int stop) {
        this.name = name;
        this.chromosome = chromosome;
        this.strand = strand;
        this.gene = gene;
        this.start = start;
        this.stop = stop;
    }

    public Gene getGene() {
        return gene;
    }

    public ArrayList<Exon> getExons() {
        return exons;
    }

    public void addExon(Exon e) {
        if(exons == null){
            exons = new ArrayList<Exon>();
        }
        int estart = e.getStart();
        int estop = e.getStop();
        if (estart < start) {
            start = estart;
        }
        if (estop > stop) {
            stop = estop;
        }
        exons.add(e);
    }

    @Override
    public String toString() {
        return "Transcript{" + "name=" + name + ", chromosome=" + chromosome + ", strand=" + strand + ", gene=" + gene.getName() + ", start=" + start + ", stop=" + stop + '}';
    }


    

  

}
