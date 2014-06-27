/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.utilities.gtf;

import hms.hwestra.utilities.features.Exon;
import hms.hwestra.utilities.features.Gene;
import hms.hwestra.utilities.features.Transcript;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.TreeSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author Harm-Jan
 */
public class GTFAnnotation {

    private String annotationLocation;
    private Collection<Gene> genes;

    public GTFAnnotation(String annotation) throws IOException {
        this.annotationLocation = annotation;
        readAnnotation();
    }

    private void readAnnotation() throws IOException {

        TextFile tf = new TextFile(annotationLocation, TextFile.R);
        String ln = tf.readLine();

        HashMap<String, Gene> strToGene = new HashMap<String, Gene>();
        HashMap<String, Transcript> strToTranscript = new HashMap<String, Transcript>();
        HashMap<String, Exon> strToExon = new HashMap<String, Exon>();

// this all assumes the file is sorted on genomic coordinates..
        while (ln != null) {
            GTFLine lineObj = new GTFLine(ln);
            Gene currentGene = null;
            Transcript currentTranscript = null;

            String geneName = lineObj.getGeneName();
            String transcriptName = lineObj.getTranscriptId();
            Exon e = new Exon(lineObj.getType(), lineObj.getChr(), lineObj.getStr(), currentGene, lineObj.getStart(), lineObj.getStop());
            if (strToGene.containsKey(geneName)) {
                // continue annotating transcripts and exons with this gene
                currentGene = strToGene.get(geneName);

            } else {
                currentGene = new Gene(geneName, lineObj.getChr(), lineObj.getStr());
                strToGene.put(geneName, currentGene);
            }

            if (strToTranscript.containsKey(transcriptName)) {
                currentTranscript = strToTranscript.get(transcriptName);
                if (strToExon.containsKey(e.toString())) {
                    // no work to be done
                    System.out.println("Duplicate entry in file: " + lineObj.toString());
                } else {
                    currentTranscript.addExon(e);
                    strToExon.put(e.toString(), e);
                }
            } else {
                currentTranscript = new Transcript(lineObj.getTranscriptId(), lineObj.getChr(), lineObj.getStr(), currentGene);
                currentGene.addTranscript(currentTranscript);
                strToTranscript.put(lineObj.getTranscriptId(), currentTranscript);
                if (strToExon.containsKey(e.toString())) {
                    e = strToExon.get(e.toString());
                    currentTranscript.addExon(e);
                } else {
                    currentTranscript.addExon(e);
                    strToExon.put(e.toString(), e);
                }
            }

            ln = tf.readLine();
        }
        tf.close();
        System.out.println(annotationLocation + ", Genes: " + strToGene.size() + "\tTranscripts: " + strToTranscript.size() + "\tExons: " + strToExon.size());

        // set the relative start and end positions of each gene and transcript
        genes = strToGene.values();

    }

    public TreeSet<Gene> buildGeneTree() {
        TreeSet<Gene> geneTree = new TreeSet<Gene>(new GeneComparator());
        
        return geneTree;
    }

    public void getGenes() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
