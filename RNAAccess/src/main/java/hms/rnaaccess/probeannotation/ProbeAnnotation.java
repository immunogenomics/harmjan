/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.rnaaccess.probeannotation;

import hms.hwestra.utilities.features.Probe;
import hms.hwestra.utilities.features.Chromosome;
import hms.hwestra.utilities.features.Exon;
import hms.hwestra.utilities.features.FeatureComparator;
import hms.hwestra.utilities.features.Gene;
import hms.hwestra.utilities.features.Strand;
import hms.hwestra.utilities.features.Transcript;
import hms.hwestra.utilities.gtf.GTFAnnotation;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author Harm-Jan
 */
public class ProbeAnnotation {

    public static void main(String[] args) {
        String probeFile = "/Data/Projects/RNAAccess/Smita's data/nexterarapidcapture_exome_probes_v1.1.bed";
        String genomeFasta = "/Data/Annotation/UCSC/genome.fa";
        String geneGTF = "/Data/Annotation/UCSC/genes.gtf";

        // annotate RNA Access probes
        // --------------------------
        // load gene/transcript/exon annotation
        ProbeAnnotation p = new ProbeAnnotation();
        try {
            p.readGTF(geneGTF);

            // load RNA Access probes
            p.readProbes(probeFile);

            // for every probe, try to find overlapping genes, transcripts, exons
            String outfile = probeFile + "-AnnotatedWithGenes.txt";
            p.overlapWithGenes(outfile);

            // determine GC content for transcripts and probes
            p.determineGCContentForTranscriptsAndProbes(genomeFasta, outfile + "-GCContent.txt", geneGTF + "-GCContent.txt");

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

        FeatureComparator comp = new FeatureComparator();
        probes = new TreeSet<Probe>(new FeatureComparator());
        ArrayList<Probe> probeList = new ArrayList<Probe>();
        while (elems != null) {
            Chromosome chr = Chromosome.parseChr(elems[0]);
            int start = Integer.parseInt(elems[1]);
            int stop = Integer.parseInt(elems[2]);
            String probeName = elems[3];
            Probe probe = new Probe(chr, start, stop, probeName);
            if (probes.contains(probe)) {

                System.out.println("Duplicate: " + probe);
                for (Probe p : probes) {
                    if (comp.compare(probe, p) == 0) {
                        System.out.println(p);
                        System.out.println(probe);
                        System.out.println();
                        p.equals(probe);
                    }
                }

                System.exit(0);
            } else {
                probes.add(probe);
            }
            probeList.add(probe);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(probeList.size() + " probes total");
        System.out.println(probes.size() + "\tprobes loaded");
        int nrInThere = 0;
        int nrMissing = 0;
        for (Probe p : probeList) {
            if (probes.contains(p)) {
                nrInThere++;
            } else {
                // System.out.println("Probe not found: " + p.toString());
                nrMissing++;
            }
        }
        System.out.println(nrInThere + "\t" + nrMissing);

    }

    private void overlapWithGenes(String outputFile) throws IOException {
        TreeSet<Gene> genes = gtf.getGeneTree();

        TextFile outfile = new TextFile(outputFile, TextFile.W);
        outfile.writeln("Probe\tProbeChr\tProbeStart\tProbeStop\tGenes\tGeneChr\tTranscript\tTranscriptChr\tExon\tExonChr");
        for (Probe p : probes) {
            Gene geneStart = new Gene("", p.getChromosome(), Strand.POS, p.getStart(), p.getStart());
            Gene geneStop = new Gene("", p.getChromosome(), Strand.POS, p.getStop(), p.getStop());

            SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);
            // now we have a set of genes, overlap the transcripts as well..

            String overlappingGeneStr = null;
            String overlappingGeneAnnotStr = null;
            String overlappingTranscriptStr = null;
            String overlappingTranscriptAnnotStr = null;
            String overlappingExonStr = null;
            String overlappingExonAnnotStr = null;
            if (overlappingGenes.isEmpty()) {
                overlappingGeneStr = "-";
                overlappingGeneAnnotStr = "-";
                overlappingTranscriptStr = "-";
                overlappingTranscriptAnnotStr = "-";
                overlappingExonStr = "-";
                overlappingExonAnnotStr = "-";
            } else {
                for (Gene g : overlappingGenes) {
                    ArrayList<Transcript> transcripts = g.getTranscripts();
                    ArrayList<Transcript> overlappingTranscripts = new ArrayList<Transcript>();

                    if (overlappingGeneStr == null) {
                        overlappingGeneStr = g.getName();
                        overlappingGeneAnnotStr = g.getChromosome().getName() + ":" + g.getStrand() + ":" + g.getStart() + "-" + g.getStop();
                    } else {
                        overlappingGeneStr += ";" + g.getName();
                        overlappingGeneAnnotStr += ";" + g.getChromosome().getName() + ":" + g.getStrand() + ":" + g.getStart() + "-" + g.getStop();
                    }

                    for (Transcript t : transcripts) {
                        if (p.getStart() >= t.getStart() && p.getStop() <= t.getStop()) {
                            // perfect overlap
                            overlappingTranscripts.add(t);

                        } else if (p.getStart() < t.getStart() && p.getStop() > t.getStart()) {
                            // partial overlap
                            overlappingTranscripts.add(t);
                        } else if (p.getStart() < t.getStop() && p.getStop() > t.getStop()) {
                            // more partial overlap
                            overlappingTranscripts.add(t);
                        }
                    }

                    // now find the overlapping exons
                    for (Transcript t : overlappingTranscripts) {
                        ArrayList<Exon> exons = t.getExons();
                        ArrayList<Exon> overlappingExons = new ArrayList<Exon>();
                        if (overlappingTranscriptStr == null) {
                            overlappingTranscriptStr = t.getName();
                            overlappingTranscriptAnnotStr = t.getChromosome().getName() + ":" + t.getStrand() + ":" + t.getStart() + "-" + t.getStop();
                        } else {
                            overlappingTranscriptStr += ";" + t.getName();
                            overlappingTranscriptAnnotStr += ";" + t.getChromosome().getName() + ":" + t.getStrand() + ":" + t.getStart() + "-" + t.getStop();
                        }
                        for (Exon e : exons) {
                            if (p.getStart() >= e.getStart() && p.getStop() <= e.getStop()) {
                                // perfect overlap
                                overlappingExons.add(e);
                            } else if (p.getStart() < e.getStart() && p.getStop() > e.getStart()) {
                                // partial overlap
                                overlappingExons.add(e);
                            } else if (p.getStart() < e.getStop() && p.getStop() > e.getStop()) {
                                // more partial overlap
                                overlappingExons.add(e);
                            }
                        }

                        for (Exon e : overlappingExons) {

                            if (overlappingExonStr == null) {
                                overlappingExonStr = e.getName();
                                overlappingExonAnnotStr = e.getChromosome().getName() + ":" + e.getStrand() + ":" + e.getStart() + "-" + e.getStop();
                            } else {
                                overlappingExonStr += ";" + e.getName();
                                overlappingExonAnnotStr += ";" + e.getChromosome().getName() + ":" + e.getStrand() + ":" + e.getStart() + "-" + e.getStop();
                            }

                        }

                    }
                }
            }

            outfile.writeln(p.getName() + "\t" + p.getChromosome().getName() + "\t" + p.getStart() + "\t" + p.getStop()
                    + "\t" + overlappingGeneStr + "\t" + overlappingGeneAnnotStr
                    + "\t" + overlappingTranscriptStr + "\t" + overlappingTranscriptAnnotStr
                    + "\t" + overlappingExonStr + "\t" + overlappingExonAnnotStr);

        }
        outfile.close();

    }

    private void determineGCContentForTranscriptsAndProbes(String genomeFasta, String probeOut, String transcriptOut) throws IOException {
        TreeSet<Transcript> transcripts = gtf.getTranscriptTree();

        // open fasta
        FastaSequenceFile fastaFile = new FastaSequenceFile(new File(genomeFasta), false);
        ReferenceSequence seq = fastaFile.nextSequence();
        int windowSize = 10000000;
        TextFile tfProbeAnnotOut = new TextFile(probeOut, TextFile.W);
        tfProbeAnnotOut.writeln("Probe\tChr\tChrStart\tChrStop\tnrGC\tnrAT\tnrN\t%GC");
        TextFile tfProbeSequenceOut = new TextFile(probeOut + ".fa.gz", TextFile.W);
        TextFile tfTranscriptAnnotOut = new TextFile(transcriptOut, TextFile.W);
        tfTranscriptAnnotOut.writeln("Gene\tGeneChr\tGeneChrStart\tGeneChrStop\t"
                + "Transcript\tTranscriptChr\tTranscriptChrStart\tTranscriptChrStop\tnrExons\t"
                + "totalNrBasesInExons\tnrGC\tnrAT\tnrN\t%GC");
        while (seq != null) {
            System.out.println(seq.getName() + "\t" + seq.length());
            byte[] bases = seq.getBases();
            // get all transcript on this sequence

            Set<Transcript> transcriptsOnSequence = transcripts.subSet(new Transcript("", Chromosome.parseChr(seq.getName()), Strand.POS, null, 0, 0), true,
                    new Transcript("", Chromosome.parseChr(seq.getName()), Strand.POS, null, bases.length + 1, bases.length + 1), true);
            System.out.println("Nr Transcripts on this sequence: " + transcriptsOnSequence.size());

            // are positions one based or 0-based? does it matter?
            // bases[] == 0 based
            // GTF (genes) == 1 based
            // BED (probes)== 0 based
            // annotate probes
            Set<Probe> probesOnSequence = probes.subSet(new Probe(Chromosome.parseChr(seq.getName()), 0, 0, ""), true, new Probe(Chromosome.parseChr(seq.getName()), bases.length + 1, bases.length + 1, ""), true);
            System.out.println("Probes on sequence: " + probesOnSequence.size());

            for (Probe p : probesOnSequence) {
                int start = p.getStart();
                int stop = p.getStop();

                byte[] subsequence = new byte[stop - start];
                // System.out.println(p.getName() + "\t" + p.getChromosome().getName() + "\t" + start + "\t" + stop + "\t" + bases.length + "\t" + subsequence.length);
                System.arraycopy(bases, start, subsequence, 0, subsequence.length);
                Pair<Triple<Integer, Integer, Integer>, Double> gcContent = determineGCContent(subsequence);
                tfProbeAnnotOut.writeln(p.getName() + "\t" + p.getChromosome().getName() + "\t" + p.getStart() + "\t" + p.getStop()
                        + "\t" + gcContent.getLeft().getLeft()
                        + "\t" + gcContent.getLeft().getMiddle()
                        + "\t" + gcContent.getLeft().getRight()
                        + "\t" + gcContent.getRight()
                );
                writeSequence(tfProbeSequenceOut, p.getName(), subsequence);
            }
            //System.exit(-1);

            for (Transcript t : transcriptsOnSequence) {
                // get exons
                ArrayList<Exon> exons = t.getExons();
                int totalLength = 0;
                for (Exon e : exons) {
                    totalLength += e.getStop() - e.getStart();
                }
                byte[] subsequence = new byte[totalLength];

                int ctr = 0;
                for (Exon e : exons) {
                    int len = e.getStop() - e.getStart();
                    System.arraycopy(bases, e.getStart() - 1, subsequence, ctr, len);
                    ctr += len;
                }

                Pair<Triple<Integer, Integer, Integer>, Double> gcContent = determineGCContent(subsequence);
                Gene parentGene = t.getGene();
                tfTranscriptAnnotOut.writeln(parentGene.getName() + "\t" + parentGene.getChromosome().getName() + "\t" + parentGene.getStart() + "\t" + parentGene.getStop()
                        + "\t" + t.getName() + "\t" + t.getChromosome().getName() + "\t" + t.getStart() + "\t" + t.getStop() + "\t" + exons.size() + "\t" + totalLength
                        + "\t" + gcContent.getLeft().getLeft()
                        + "\t" + gcContent.getLeft().getMiddle()
                        + "\t" + gcContent.getLeft().getRight()
                        + "\t" + gcContent.getRight()
                );
            }
            seq = fastaFile.nextSequence();
        }
        tfTranscriptAnnotOut.close();
        tfProbeAnnotOut.close();
        tfProbeSequenceOut.close();
    }

    private Pair<Triple<Integer, Integer, Integer>, Double> determineGCContent(byte[] sequence) {

        int nrGC = 0;
        int nrAT = 0;
        int nrN = 0;
        for (int i = 0; i < sequence.length; i++) {

            byte c = sequence[i];
            if (c == 'c' || c == 'C' || c == 'g' || c == 'G') {
                nrGC++;
            } else if (c == 'a' || c == 'A' || c == 't' || c == 'T') {
                nrAT++;
            } else {
                nrN++;
            }

        }

        double gcContent = (double) nrGC / (nrGC + nrAT + nrN);

        return new Pair<Triple<Integer, Integer, Integer>, Double>(new Triple<Integer, Integer, Integer>(nrGC, nrAT, nrN), gcContent);
    }

    private void writeSequence(TextFile tfProbeSequenceOut, String name, byte[] sequence) throws IOException {
        tfProbeSequenceOut.writeln(">" + name);
        tfProbeSequenceOut.writeln(new String(sequence, "UTF-8"));
    }

}
