/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.utilities.features;

import java.util.ArrayList;
import java.util.Objects;

/**
 *
 * @author Harm-Jan
 */
public class Gene extends Feature {

    ArrayList<Transcript> transcripts;
    private String geneId;

    public Gene(String name, Chromosome chromosome, Strand strand) {
        this.name = name;
        this.chromosome = chromosome;
        this.strand = strand;
    }

    public Gene(String name, Chromosome chromosome, Strand strand, int start, int stop) {
        this.name = name;
        this.chromosome = chromosome;
        this.strand = strand;
        this.start = start;
        this.stop = stop;

    }

    public void addTranscript(Transcript t) {
        if (transcripts == null) {
            transcripts = new ArrayList<Transcript>();
        }
        transcripts.add(t);
    }

    public ArrayList<Transcript> getTranscripts() {
        return transcripts;
    }

    @Override
    public String toString() {
        return "Gene{" + "chromosome=" + chromosome + ", geneId=" + geneId + ", name=" + name + ", strand=" + strand + ", start=" + start + ", stop=" + stop + '}';
    }

    public void getBounds() {
        for (Transcript t : transcripts) {

            if (t.getStart() < start) {
                start = t.getStart();
            }
            if (t.getStop() > stop) {
                stop = t.getStop();
            }

        }
    }

    public void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    public String getGeneId(){
        return geneId;
    }
}
