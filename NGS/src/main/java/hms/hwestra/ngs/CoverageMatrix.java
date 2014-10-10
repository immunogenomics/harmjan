/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.ngs;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.text.Strings;

/**
 *
 * @author hwestra
 */
public class CoverageMatrix {

    private final DoubleMatrixDataset<String, String> matrix;

    public CoverageMatrix(String coverageMatrix) throws IOException {
        matrix = DoubleMatrixDataset.loadDoubleData(coverageMatrix);
    }

    public ArrayList<Pair<String, ArrayList<String>>> selectSamples(int column, double threshold, String outdir) throws IOException {

        System.out.println("Selecting samples from column: " + matrix.getColObjects().get(column) + " with threshold: " + threshold);

        // filter samples on 20x coverage
        ArrayList<String> selectedSamples = new ArrayList<String>();
        for (int s = 0; s < matrix.rows(); s++) {
            if (matrix.getElement(s, column) >= threshold) {
                String sample = matrix.getRowObjects().get(s);
                selectedSamples.add(sample.trim());
            }
        }
        System.out.println(selectedSamples.size() + " samples with coverage above threshold");

// detect duplicates
        HashMap<String, ArrayList<String>> duplicates = new HashMap<String, ArrayList<String>>();
        for (String sampleName : selectedSamples) {
            String[] sampleElems = sampleName.split("/");
            String realSample = sampleElems[1].trim();
            realSample = realSample.split("_")[0];

            ArrayList<String> dups = duplicates.get(realSample);
            if (dups == null) {
                dups = new ArrayList<String>();
            }
            dups.add(sampleName);
            duplicates.put(realSample, dups);
        }

        System.out.println("");
        System.out.println("Duplicates:");
        TextFile outfile = new TextFile(outdir + "sampleSelection.txt", TextFile.W);
        outfile.writeln("Sample\tDups");
        Set<String> keys = duplicates.keySet();
        ArrayList<Pair<String, ArrayList<String>>> selectedSamplesWDups = new ArrayList<Pair<String, ArrayList<String>>>();
        for (String key : keys) {

            ArrayList<String> dups = duplicates.get(key);

            outfile.writeln(key + "\t" + Strings.concat(dups.toArray(new String[0]), Strings.semicolon));
            if (dups.size() > 1) {
                System.out.println(key + "\t" + Strings.concat(dups.toArray(new String[0]), Strings.semicolon));
                Pair<String, ArrayList<String>> sampleWDups = new Pair<String, ArrayList<String>>(key, dups);
                selectedSamplesWDups.add(sampleWDups);
            } else {
                Pair<String, ArrayList<String>> sampleWDups = new Pair<String, ArrayList<String>>(key, dups);
                selectedSamplesWDups.add(sampleWDups);
            }

        }
        outfile.close();

        return selectedSamplesWDups;
    }
}
