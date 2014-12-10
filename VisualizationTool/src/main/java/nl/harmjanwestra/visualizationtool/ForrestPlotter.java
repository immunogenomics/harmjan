/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.visualizationtool;

import umcg.genetica.graphics.ForestPlot;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author Harm-Jan
 */
public class ForrestPlotter {

    public static void main(String[] args) {
        ForestPlot f = new ForestPlot();
        String[] snps = new String[]{"rs6265", "rs10767664"};

        String infile = "C:\\Work\\2012-12-21-CisAssociationsProbeLevelFDR0.5.txt";
        String outdir = "C:\\Work\\";
        int[] datasetSizes = new int[]{891, 963, 1240, 229, 762, 509, 611, 43, 63};
        //EGCUT, SHIP-TREND, Fehrmann-HT12v3, Fehrmann-H8v2, Rotterdam Study, DILGOM, InChianti, HVH-HT12v3, HVH-HT12v4
        try {
            for (String snp : snps) {
                TextFile tf = new TextFile(infile, TextFile.R);
                tf.readLine();
                String[] elems = tf.readLineElems(TextFile.tab);
                while (elems != null) {
                    if (elems.length > 1) {
                        if (elems[1].equals(snp)) {
                            String gene = elems[13];
                            String probe = elems[4];
                            double fdr = Double.parseDouble(elems[14]);
                            double pval = Double.parseDouble(elems[0]);
                            String[] datasetNames = elems[11].split(",");
                            String[] datasetZScoresStr = elems[12].split(",");

                            Double[] datasetZscores = new Double[datasetNames.length + 1];
                            String[] rowNames = new String[datasetZscores.length];
                            int[] weights = new int[datasetNames.length + 1];
                            for (int i = 0; i < datasetNames.length; i++) {

                                weights[i] = 0;
                                try {
                                    datasetZscores[i] = Double.parseDouble(datasetZScoresStr[i]);
                                    weights[i] = datasetSizes[i];
                                    weights[weights.length - 1] += datasetSizes[i];
                                    rowNames[i] = datasetNames[i] + " (n = " + datasetSizes[i] + ")";
                                } catch (NumberFormatException e) {
                                    datasetZscores[i] = null;//Double.NaN;
                                    rowNames[i] = "-";
                                }
                            }
                            double metaZ = Double.parseDouble(elems[10]);
                            rowNames[rowNames.length - 1] = "Meta-analysis (n = " + weights[weights.length - 1] + ")";
                            datasetZscores[datasetZscores.length - 1] = metaZ;
                            byte[] chr = new byte[1];
                            chr[0] = ChrAnnotation.parseChr(elems[5]);
                            int[] chrPos = new int[1];
                            chrPos[0] = Integer.parseInt(elems[6]);

                            String outfileName = outdir + snp + "-" + probe + "-" + gene + ".pdf";
                            double[] thresholds = new double[]{3.824, 2.918};

                            f.setGeneNames(new String[]{gene});
                            f.drawForrestPlot("ZScore", rowNames, datasetZscores, outfileName, ForestPlot.Output.PDF, thresholds, -5d, 5d, weights, chr, chrPos, datasetZscores.length - 1);

                            //drawForrestPlot(String xAxisName, String[] yAxisNames, Double[] xValues, String filename, Output output, double[] significanceThresholds, Double minX, Double maxX, int[] weights, byte[] chr, int[] chrpos, int metaRow) throws IOException, DocumentException {
                        }
                    }
                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}
