/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.rnaaccess;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.containers.Pair;
import umcg.genetica.graphics.ViolinBoxPlot;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.WilcoxonMannWhitney;

/**
 *
 * @author Harm-Jan
 */
public class SaturationAnalysis {

    HashMap<String, Integer> geneIndex = new HashMap<String, Integer>();
    ArrayList<String> geneNames = new ArrayList<String>();

    public static void main(String[] args) {

        if (args.length < 3) {
            System.out.println("Usage: inputfile.txt nriter outputdir");
            System.exit(-1);
        }

        String inputfile = args[0];
        int nrIter = Integer.parseInt(args[1]);
        String outdir = args[2];

        // indexAllGenes
        SaturationAnalysis s = new SaturationAnalysis();

        ArrayList<String> dirs = new ArrayList<String>();
        ArrayList<String> sampleNames = new ArrayList<String>();

        try {
            TextFile tf = new TextFile(inputfile, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String sample = elems[0];
                String path = elems[1];
                dirs.add(path);
                sampleNames.add(sample);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            s.run(dirs.toArray(new String[0]), sampleNames.toArray(new String[0]), nrIter, outdir);
        } catch (IOException ex) {
            Logger.getLogger(SaturationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(SaturationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void run(String[] sampleLocations, String[] sampleNames, int nrIter, String outdir) throws IOException, Exception {

        ArrayList<HashMap<Integer, Integer>> readCounts = new ArrayList<HashMap<Integer, Integer>>();
        ArrayList<double[][]> fpkms = new ArrayList<double[][]>();
        for (String dir : sampleLocations) {
            indexAllGenes(dir, nrIter);
            readCounts.add(getReadCounts(dir, nrIter));
            fpkms.add(getfpkms(dir, nrIter));
        }

        // everything is loaded, now do some statistics...
        // every dir is 1 sample. 
        // compare iterations within sample
        SpearmansCorrelation corr = new SpearmansCorrelation();
        WilcoxonMannWhitney mwm = new WilcoxonMannWhitney();
        for (int d = 0; d < sampleLocations.length; d++) {
            TextFile tfOut = new TextFile(outdir + "PairWiseCorrelations-" + sampleNames[d] + ".txt", TextFile.W);
            TextFile tfOut2 = new TextFile(outdir + "PairWiseWilcoxonTests-" + sampleNames[d] + ".txt", TextFile.W);
            String s = "-";
            for (int i = 0; i < nrIter; i++) {
                s += "\tIter" + i;
            }
            tfOut.writeln(s);

            double[][] fpkm = fpkms.get(d);
            for (int i = 0; i < nrIter; i++) {
                double[] x = fpkm[i];
                StringBuilder sb = new StringBuilder();
                sb.append("Iter");
                sb.append(i);

                StringBuilder sb2 = new StringBuilder();
                sb2.append("Iter");
                sb2.append(i);

                for (int j = i; j < nrIter; j++) {
                    double[] y = fpkm[j];
                    Pair<double[], double[]> p = trim(x, y);

                    x = p.getLeft();
                    y = p.getRight();
                    sb.append("\t");
                    sb2.append("\t");
                    if (x.length > 3 && y.length > 3) {

                        double spearman = corr.correlation(x, y);
                        double mwmp = mwm.returnWilcoxonMannWhitneyPValue(x, y);

                        double auc = mwm.getAUC();

                        sb2.append(mwmp);
                        sb2.append("(").append(auc).append(")");

                        sb.append(spearman).append("(").append(x.length).append(")");
                    } else {
                        sb.append(Double.NaN);
                        sb2.append(Double.NaN);
                    }
                }
                tfOut2.writeln(sb2.toString());
                tfOut.writeln(sb.toString());
            }
            tfOut.close();
            tfOut2.close();

            // plot forestplot
            ViolinBoxPlot plot = new ViolinBoxPlot();

            double[][][] data = new double[1][nrIter][0];
            String[][] xLabels = new String[1][nrIter];
            for (int i = 0; i < nrIter; i++) {
                double[] x = fpkm[i];
                xLabels[0][i] = "" + ((i + 1) * 10);
                int ctr = 0;
                for (int q = 0; q < x.length; q++) {
                    if (x[q] != 0) {
                        ctr++;
                    }
                }
                double[] vals = new double[x.length - ctr];
                ctr = 0;
                for (int q = 0; q < x.length; q++) {
                    if (x[q] != 0) {
                        vals[ctr] = Math.log(x[q]) / Math.log(2);
                        if (vals[ctr] > 250) {
                            vals[ctr] = 250;
                        }
                        ctr++;
                    }
                }
                data[0][i] = vals;
            }
            plot.draw(data, new String[]{sampleNames[d]}, xLabels, "Log2 FPKM", ViolinBoxPlot.Output.PDF, outdir + "Distributions-" + sampleNames[d] + ".pdf");

            // plot histograms for each iteration...
            CombineTrackingFiles f = new CombineTrackingFiles();

            DoubleMatrixDataset ds = new DoubleMatrixDataset();
            ds.setMatrix(data[0]);
            ds.setRowObjects(Arrays.asList(xLabels[0]));
            ds.setColObjects(geneNames);
            /// drawHistograms(DoubleMatrixDataset<String, String> ds, String histogramOut, boolean binZeros, double maxVal, boolean log10, Output output) throws DocumentException, FileNotFoundException, IOException {
            f.drawHistograms(ds, outdir + "Histograms-" + sampleNames[d] + ".pdf", true, 250, true, CombineTrackingFiles.Output.PDF);

            double[][] foldchanges = new double[nrIter - 1][];
            double[][] averages = new double[nrIter - 1][];
            double[] maxVals = fpkm[fpkm.length - 1];
            String[] plotTitles = new String[nrIter - 1];
            for (int iter = 0; iter < nrIter - 1; iter++) {
                plotTitles[iter] = "" + iter;
                double[] vals = fpkm[iter];
                int ctr = 0;
                for (int q = 0; q < maxVals.length; q++) {
                    if (maxVals[q] != 0 && vals[q] != 0) {
                        ctr++;
                    }
                }

                foldchanges[iter] = new double[ctr];
                averages[iter] = new double[ctr];
                ctr = 0;
                for (int q = 0; q < maxVals.length; q++) {
                    if (maxVals[q] != 0 && vals[q] != 0) {
                        foldchanges[iter][ctr] = log2(vals[q] / maxVals[q]);
                        averages[iter][ctr] = (log2(vals[q]) + log2(maxVals[q])) / 2;
                        ctr++;
                    }
                }
            }

            ScatterPlot2 scatter = new ScatterPlot2(averages, foldchanges, false, plotTitles, "Log2 Average FPKM", "Log2 FoldChange FPKM",
                    sampleNames[d], 250, 250, 50, 100, 5, "MAPlots-" + sampleNames[d] + ".pdf");

            /*
             double[][] x, double[][] y, boolean interceptAtZero, String[] plotTitles, String xAxisTitle, String yAxisTitle,
             String figureTitle, int plotWidth, int plotHeight, int plotMargin, int margin, int nrColumns, String outputFileName
             */
        }

//            TextFile tfOut = new TextFile(outdir + "ComparisonTopVsIterations.txt", TextFile.W);
//            for (int d = 0; d < dirs.length; d++) {
//                // get the max iter
//                // compare against iterations from other samples
//            }
//            tfOut.close();
        // correlate genes' fpkm with readcount over all iterations
//            for (int d = 0; d < sampleLocations.length; d++) {
//                TextFile tfOut = new TextFile(outdir + "PairWiseCorrelations-" + sampleNames[d] + ".txt", TextFile.W);
//
//                tfOut.close();
//            }
        // plot and test read distributions per iteration.
        // test the difference between genes between samples:
        // pick pair of samples, calculate fold-change per gene for full set of reads
        // plot FPKM vs fold-change
    }

    private double log2(double d) {
        return Math.log(d) / Math.log(2);
    }

    private void indexAllGenes(String dir, int nrIter) throws IOException {
        for (int i = 1; i < nrIter + 1; i++) {
            double d = (double) i / nrIter;
            String dStr = "" + d;
            if (d < 1) {
                dStr = dStr.substring(1, 3);
            } else {
                dStr = dStr.substring(0, 3);
            }
            System.out.println(dStr);
            String readFile1 = dir + "Iteration" + dStr + "\\genes.fpkm_tracking";

            TextFile tf = new TextFile(readFile1, TextFile.R);
            tf.readLine(); // skip header
            String[] elems = tf.readLineElems(TextFile.tab);

            while (elems != null) {
                String gene = elems[0] + "-" + elems[4] + "-" + elems[5];
                if (!geneIndex.containsKey(gene)) {
                    geneIndex.put(gene, geneIndex.size());
                    geneNames.add(gene);
                }

                elems = tf.readLineElems(TextFile.tab);
            }

            tf.close();
        }
    }

    public HashMap<Integer, Integer> getReadCounts(String dir, int nrIter) throws IOException {
        // read the number of reads
        HashMap<Integer, Integer> nrReads = new HashMap<>();
        for (int i = 1; i < nrIter + 1; i++) {
            double d = (double) i / nrIter;
            String dStr = "" + d;
            if (d < 1) {
                dStr = dStr.substring(1, 3);
            } else {
                dStr = dStr.substring(0, 3);
            }
            String readFile1 = dir + "aligned_reads-counts." + dStr + ".BAM";

            TextFile tf = new TextFile(readFile1, TextFile.R);

            Integer it = Integer.parseInt(tf.readLine());

            tf.close();

            nrReads.put(i - 1, it);
        }
        return nrReads;
    }

    private double[][] getfpkms(String dir, int nrIter) throws IOException {
        double[][] fpkms = new double[nrIter][geneIndex.size()];
        for (int i = 1; i < nrIter + 1; i++) {
            double d = (double) i / nrIter;
            String dStr = "" + d;
            if (d < 1) {
                dStr = dStr.substring(1, 3);
            } else {
                dStr = dStr.substring(0, 3);
            }
            String readFile1 = dir + "Iteration" + dStr + "\\genes.fpkm_tracking";

            TextFile tf = new TextFile(readFile1, TextFile.R);
            tf.readLine(); // skip header
            String[] elems = tf.readLineElems(TextFile.tab);

            while (elems != null) {
                String gene = elems[0] + "-" + elems[4] + "-" + elems[5];
                int geneIndexNr = geneIndex.get(gene);

                double fpkm = Double.parseDouble(elems[10]);
                fpkms[i - 1][geneIndexNr] = fpkm;

                elems = tf.readLineElems(TextFile.tab);
            }

            tf.close();
        }
        return fpkms;
    }

    private Pair<double[], double[]> trim(double[] x, double[] y) {
        int ctr = 0;
        for (int q = 0; q < x.length; q++) {
            if (x[q] == 0 && y[q] == 0) {
                ctr++;
            }
        }
        System.out.println(ctr + " nulls");
        double[] a = new double[x.length - ctr];
        double[] b = new double[x.length - ctr];
        ctr = 0;
        for (int q = 0; q < x.length; q++) {
            if (x[q] == 0 && y[q] == 0) {

            } else {
                a[ctr] = x[q];
                b[ctr] = y[q];
                ctr++;
            }
        }
        return new Pair<>(a, b);
    }
}
