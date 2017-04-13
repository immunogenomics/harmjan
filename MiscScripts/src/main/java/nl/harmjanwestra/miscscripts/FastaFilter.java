/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.miscscripts;

import java.io.IOException;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

/**
 *
 * @author hwestra
 */
public class FastaFilter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String file = "/Data/Ref/human_g1k_v37.fasta";
        String fileOut = "/Data/Ref/human_g1k_v37-AutosomalSexAndMT.fasta";
        FastaFilter f = new FastaFilter();
        try {
            f.filter(file, fileOut);
        } catch (IOException ex) {
            Logger.getLogger(FastaFilter.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void filter(String filein, String fileout) throws IOException {

        TextFile f = new TextFile(filein, TextFile.R);
        TextFile fout = new TextFile(fileout, TextFile.W);
        String ln = f.readLine();
        HashSet<Character> set = new HashSet<Character>();
        for (int i = 1; i < 23; i++) {
            String s = "" + i;
            set.add(s.charAt(0));
        }
        set.add('X');
        set.add('Y');
        set.add('M');
        boolean printBlock = false;
        while (ln != null) {

            if (ln.startsWith(">")) {
                // is header
                char c = ln.charAt(1);
                if (set.contains(c)) {
                    fout.writeln(ln);
                    printBlock = true;
                    System.out.println("Keeping: " + ln);
                } else {
                    printBlock = false;
                }

            } else if (printBlock) {
                fout.writeln(ln);
            }

            ln = f.readLine();
        }
        f.close();
        fout.close();
    }

}
