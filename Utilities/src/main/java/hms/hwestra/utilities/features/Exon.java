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
public class Exon extends Feature {

    private ArrayList<Transcript> transcripts;
    private final Gene gene;

    public Exon(String name, Chromosome chr, Strand strand, Gene gene, int start, int stop) {
        this.name = name;
        this.chromosome = chr;
        this.strand = strand;

        this.gene = gene;
        this.start = start;
        this.stop = stop;

    }

    public ArrayList<Transcript> getTranscripts() {
        return transcripts;
    }

    public Gene getGene() {
        return gene;
    }

    public void addTranscript(Transcript t) {
        if (this.transcripts == null) {
            this.transcripts = new ArrayList<Transcript>();
        }
        this.transcripts.add(t);
    }

    @Override
    public String toString() {
        return "Exon{" + "chromosome=" + chromosome + ", name=" + name + ", strand=" + strand + ", gene=" + gene.getName() + ", start=" + start + ", stop=" + stop + '}';
    }

}
