/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.visualizationtool;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author Harm-Jan
 */
public class Filter {

    public static void main(String[] args) {
        try {
            String file = "C:\\Work\\eQTLs.txt";
            HashSet<String> s = new HashSet<String>();
            s.add("rs354033");
            s.add("rs874628");
            s.add("rs1790100");
            s.add("rs1062158");
            s.add("rs1132200");

            String outfile = "c:\\work\\selection2.txt";

            TextFile tf = new TextFile(file, TextFile.R);
            TextFile tfOut = new TextFile(outfile, TextFile.W);
            String[] elems = tf.readLineElems(TextFile.tab);
            tfOut.writeln(Strings.concat(elems, Strings.tab));
            elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {

                if(elems[16].contains("ICAM3")){
                    tfOut.writeln(Strings.concat(elems, Strings.tab));
                }

                elems = tf.readLineElems(TextFile.tab);
            }
            tfOut.close();
            tf.close();
        } catch (IOException e) {

        }
    }
}
