/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.rnaaccess.probeannotation;

import hms.hwestra.utilities.features.Chromosome;
import hms.hwestra.utilities.gtf.GTFAnnotation;
import java.io.IOException;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author Harm-Jan
 */
public class ProbeAnnotation {

    public static void main(String[] args) {
        String probeFile = args[0];
        String genomeFasta = args[1];
        String geneGTF = args[2];

        // annotate RNA Access probes
        // --------------------------
        // load gene/transcript/exon annotation
        ProbeAnnotation p = new ProbeAnnotation();
        try {
            p.readGTF(geneGTF);

            // load RNA Access probes
            p.readProbes(probeFile);

            // for every probe, try to find overlapping genes, transcripts, exons
            p.overlapWithGenes();
            
            // unload gene/transcript/exon annotation
            // sort probes by chromosome position
            // read fasta line by line in chunks of 10Mb with 2mb overlap
            // convert to character[]
            // get all probes in sequence window
            // get specific subsequences for probes
            // store sequence in probe fasta file
            // calculate G/C content for specific probe
            // 
        } catch (IOException ex) {
            Logger.getLogger(ProbeAnnotation.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    GTFAnnotation gtf;
    private TreeSet<Probe> probes;

    public void readGTF(String gtfFile) throws IOException {
        gtf = new GTFAnnotation(gtfFile);
    }

    private void readProbes(String probeFile) throws IOException {
        TextFile tf = new TextFile(probeFile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);

        probes = new TreeSet<Probe>(new ProbeComparator());
        while (elems != null) {
            Chromosome chr = Chromosome.parseChr(elems[0]);
            int start = Integer.parseInt(elems[1]);
            int stop = Integer.parseInt(elems[2]);
            String probeName = elems[3];
            Probe probe = new Probe(chr, start, stop, probeName);
            probes.add(probe);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

    }

    private void overlapWithGenes() {
        gtf.getGenes();
    }

}
