/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.miscscripts;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

/**
 *
 * @author hwestra
 */
public class NewMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        NewMain.correlationBonferroni(0.05, 48000, 860);
//        String f1 = "/Data/Projects/BrandonPierce/mediated_trans_eQTLs_HJW.txt";
//        String f2 = "/Data/ProbeAnnotation/2013-07-18-HT12v3.txt";
//
//        NewMain m = new NewMain();
//        try {
//            m.remapProbesBrandonPierce(f2, f1);
//        } catch (IOException e) {
//
//        }
    }

    private static void correlationBonferroni(double pthresh, int tests, int sampleSize) {
        double bonferroni = pthresh / tests;
        double stepsize = 0.001;
        Correlation.correlationToZScore(sampleSize);
        int iterations = 0;
        double corr = 0;
        double z = 0;
        double p = 0;
        System.out.println(bonferroni);

        while (corr < 1) {
            z = Correlation.convertCorrelationToZScore(sampleSize, corr);
            p = ZScores.zToP(z);
//            if (p < bonferroni) {

                System.out.println((p<bonferroni)+"\t"+iterations + "\t" + stepsize + "\t" + corr + "\t" + z + "\t" + p);

//            }

            corr += stepsize;
        }
    }

    public void remapProbesBrandonPierce(String annot, String file) throws IOException {
        TextFile tf1 = new TextFile(annot, TextFile.R);
        String[] elems = tf1.readLineElems(TextFile.tab);
        HashMap<String, String> seqToProbe = new HashMap<String, String>();
        while (elems != null) {
            String seq = elems[elems.length - 1];
            String probe = elems[1];
            seqToProbe.put(seq, probe);
            elems = tf1.readLineElems(TextFile.tab);
        }
        tf1.close();

        TextFile tfout = new TextFile(file + "-wcombos.txt", TextFile.W);
        TextFile tf3 = new TextFile(file, TextFile.R);

        elems = tf3.readLineElems(TextFile.tab);
        while (elems != null) {

            String seq1 = elems[3];
            String seq2 = elems[4];
            tfout.writeln(Strings.concat(elems, Strings.tab) + "\t" + seqToProbe.get(seq1) + "\t" + seqToProbe.get(seq2));
            elems = tf3.readLineElems(TextFile.tab);
        }
        tf3.close();
        tfout.close();
    }

}
