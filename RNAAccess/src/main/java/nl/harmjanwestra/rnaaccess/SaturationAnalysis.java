/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.rnaaccess;

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
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.Regression;
import umcg.genetica.math.stats.WilcoxonMannWhitney;
import umcg.genetica.util.Primitives;

/**
 *
 * @author Harm-Jan
 */
public class SaturationAnalysis {

    HashMap<String, Integer> geneIndex = new HashMap<>();
    ArrayList<String> geneNames = new ArrayList<>();

    public static void main(String[] args) {

        if (args.length < 3) {
            System.out.println("Usage: inputfile.txt nriter outputdir");
            System.out.println("Inputfileformat: samplename path");
            System.exit(-1);
        }

        String inputfile = args[0];
        int nrIter = Integer.parseInt(args[1]);
        String outdir = args[2];

        // indexAllGenes
        SaturationAnalysis s = new SaturationAnalysis();

        ArrayList<String> dirs = new ArrayList<>();
        ArrayList<String> sampleNames = new ArrayList<>();

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

        ArrayList<int[]> readCounts = new ArrayList<>();
        ArrayList<double[][]> fpkms = new ArrayList<>();
        for (String dir : sampleLocations) {
            indexAllGenes(dir, nrIter);
            readCounts.add(getReadCounts(dir, nrIter));
            fpkms.add(getfpkms(dir, nrIter));
        }

        // everything is loaded, now do some statistics...
        // every dir is 1 sample. 
        // compare iterations within sample
        for (int d = 0; d < sampleLocations.length; d++) {

            ArrayList<String> colNames = new ArrayList<>();
            for (int i = 0; i < nrIter; i++) {
                String s = sampleNames[d] + "-" + ((i + 1) * 10);
                colNames.add(s);
            }

            double[][] fpkm = fpkms.get(d);
            System.out.println(fpkm.length + "x" + fpkm[fpkm.length - 1].length + " size of matrix");
            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<>();
            ds.setMatrix(fpkm);
            ds.setMatrix(ds.getMatrix().viewDice());

            ds.setColObjects(colNames);
            ds.setRowObjects(geneNames);

            ds.save(outdir + sampleNames[d] + "-FPKMs.txt");

//            saveReadCounts(sampleNames[d], outdir);
            correlateIterations(fpkm, nrIter, sampleNames[d], outdir);

            // plot forestplot
            createViolinBoxPlots(fpkm, nrIter, sampleNames[d], outdir);

            // plot histograms for each iteration...
            createHistograms(fpkm, nrIter, sampleNames[d], outdir);

            // create foldchange plots
            createMAPlots(fpkm, nrIter, sampleNames[d], outdir);

            correlateWithReadCounts(fpkm, readCounts.get(d), nrIter, sampleNames[d], outdir);

        }

        // save full table
        double[][] fullTable = new double[geneNames.size()][fpkms.size() * nrIter];
        int[] fullReadCounts = new int[fpkms.size() * nrIter];
        int ctr = 0;
        ArrayList<String> colNames = new ArrayList<String>();
        TextFile tfOut = new TextFile(outdir + "AllSamples-ReadCounts.txt", TextFile.W);
        tfOut.writeln("Sample\tIter\tReadCount");
        for (int d = 0; d < sampleNames.length; d++) {
            double[][] fpkm = fpkms.get(d);
            for (int iter = 0; iter < nrIter; iter++) {
                String s = sampleNames[d] + "-" + ((iter + 1) * 10);
                colNames.add(s);
                for (int gene = 0; gene < geneNames.size(); gene++) {
                    fullTable[gene][ctr] = fpkm[iter][gene];

                }

                fullReadCounts[ctr] = readCounts.get(d)[iter];
                tfOut.writeln(s + "\t" + iter + "\t" + fullReadCounts[ctr]);
                ctr++;
            }
        }
        tfOut.close();

        DoubleMatrixDataset<String, String> dsOut = new DoubleMatrixDataset<String, String>();
        dsOut.setColObjects(colNames);
        dsOut.setRowObjects(geneNames);
        dsOut.setMatrix(fullTable);
        dsOut.save(outdir + "AllSamples.txt");

        correlateFullTableSamples(dsOut, fullReadCounts, outdir + "AllSamples-CorrelationsWithReadCounts.txt");

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

            String readFile1 = dir + "/cufflinks/Iteration" + dStr + "/genes.fpkm_tracking";
            System.out.println("Reading file: " + readFile1);
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

    public int[] getReadCounts(String dir, int nrIter) throws IOException {
        // read the number of reads
        int[] readCounts = new int[nrIter];
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

            readCounts[i - 1] = it;
        }
        return readCounts;
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
            String readFile1 = dir + "/cufflinks/Iteration" + dStr + "/genes.fpkm_tracking";

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

    private void correlateIterations(double[][] fpkm, int nrIter, String sampleName, String outputdir) throws Exception {
        SpearmansCorrelation corr = new SpearmansCorrelation();
        WilcoxonMannWhitney mwm = new WilcoxonMannWhitney();

        double[][] mwmArr = new double[nrIter][nrIter];
        double[][] corArr = new double[nrIter][nrIter];

        ArrayList<String> rowNames = new ArrayList<String>();

        for (int iter1 = 0; iter1 < nrIter; iter1++) {
            rowNames.add(sampleName + "-" + ((iter1 + 1) * 10));
        }

        for (int iter1 = 0; iter1 < nrIter; iter1++) {
            double[] x = fpkm[iter1];

            for (int iter2 = iter1; iter2 < nrIter; iter2++) {
                double[] y = fpkm[iter2];
                Pair<double[], double[]> p = trim(x, y);

                x = p.getLeft();
                y = p.getRight();

                if (x.length > 3 && y.length > 3) {
                    corArr[iter1][iter2] = corr.correlation(x, y);
                    mwmArr[iter1][iter2] = mwm.returnWilcoxonMannWhitneyPValue(x, y);

                } else {

                    corArr[iter1][iter2] = Double.NaN;
                    mwmArr[iter1][iter2] = Double.NaN;
                }
            }
        }

        DoubleMatrixDataset ds = new DoubleMatrixDataset();
        ds.setRowObjects(rowNames);
        ds.setColObjects(rowNames);
        ds.setMatrix(corArr);
        ds.save(outputdir + "CorrelationMatrix-" + sampleName + ".txt");
        ds.setMatrix(mwmArr);
        ds.save(outputdir + "WilcoxonMatrix-" + sampleName + ".txt");

    }

    private void createViolinBoxPlots(double[][] fpkm, int nrIter, String sampleName, String outdir) throws IOException {
        ViolinBoxPlot plot = new ViolinBoxPlot();

        double[][][] data = new double[1][nrIter][0];
        String[][] xLabels = new String[1][nrIter];

        // remove zeroes
        for (int iteration = 0; iteration < nrIter; iteration++) {
            double[] x = fpkm[iteration];
            xLabels[0][iteration] = sampleName + "-" + ((iteration + 1) * 10);
            int ctr = 0;
            for (int q = 0; q < x.length; q++) {
                if (x[q] != 0) {
                    ctr++;
                }
            }
            double[] vals = new double[ctr];
            ctr = 0;
            for (int q = 0; q < x.length; q++) {
                if (x[q] != 0) {
                    vals[ctr] = Math.log(x[q]) / Math.log(2);

                    ctr++;
                }
            }
            data[0][iteration] = vals;
        }
        plot.draw(data, new String[]{sampleName}, xLabels, "Log2 FPKM", ViolinBoxPlot.Output.PDF, outdir + "Distributions-" + sampleName + ".pdf");

        // dont remove zeroes.
        for (int iteration = 0; iteration < nrIter; iteration++) {
            double[] x = fpkm[iteration];
            double[] vals = new double[x.length];
            for (int q = 0; q < x.length; q++) {

                vals[q] = Math.log(x[q]) / Math.log(2);
                if (Double.isNaN(vals[q])) {
                    vals[q] = 0;
                }
            }

            data[0][iteration] = vals;
        }
        plot.draw(data, new String[]{sampleName}, xLabels, "Log2 FPKM", ViolinBoxPlot.Output.PDF, outdir + "Distributions-WithNulls-" + sampleName + ".pdf");
    }

    private void createHistograms(double[][] fpkm, int nrIter, String sampleName, String outdir) throws IOException, Exception {
        CombineTrackingFiles f = new CombineTrackingFiles();

        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset();
        ds.setMatrix(fpkm);
        String[] xLabels = new String[nrIter];
        for (int iteration = 0; iteration < nrIter; iteration++) {

            xLabels[iteration] = sampleName + "-" + ((iteration + 1) * 10);
        }
        ds.setMatrix(ds.getMatrix().viewDice());
        ds.setColObjects(Arrays.asList(xLabels));
        ds.setRowObjects(geneNames);
        /// drawHistograms(DoubleMatrixDataset<String, String> ds, String histogramOut, boolean binZeros, double maxVal, boolean log10, Output output) throws DocumentException, FileNotFoundException, IOException {
        f.drawHistograms(ds, outdir + "Histograms-" + sampleName + ".pdf", false, 250, true, CombineTrackingFiles.Output.PDF);

    }

    private void createMAPlots(double[][] fpkm, int nrIter, String sampleName, String outdir) throws IOException {
        double[][] foldchanges = new double[nrIter - 1][];
        double[][] averages = new double[nrIter - 1][];
        double[] maxVals = fpkm[fpkm.length - 1];
        String[] plotTitles = new String[nrIter - 1];

        for (int iter = 0; iter < nrIter - 1; iter++) {
            String iterSampleStr = sampleName + "-" + ((iter + 1) * 10);
            plotTitles[iter] = "" + iterSampleStr;
            double[] vals = fpkm[iter];

            System.out.println("FoldChange comp: " + iter + " - " + vals.length + " genes");
            foldchanges[iter] = new double[maxVals.length];
            averages[iter] = new double[maxVals.length];

            for (int gene = 0; gene < maxVals.length; gene++) {

                if (maxVals[gene] == 0 || vals[gene] == 0) {
                    foldchanges[iter][gene] = Double.NaN;
                    averages[iter][gene] = Double.NaN;
                } else {
                    foldchanges[iter][gene] = log2(vals[gene] / maxVals[gene]);
                    averages[iter][gene] = (log2(vals[gene]) + log2(maxVals[gene])) / 2;
                }

            }
        }

        TextFile tfout = new TextFile(outdir + "MA-" + sampleName + ".txt", TextFile.W);

        String header = "Gene";
        for (int iter = 0; iter < nrIter - 1; iter++) {
            String iterSampleStr = sampleName + "-" + ((iter + 1) * 10);
            header += "\tFC-" + iterSampleStr + "\tAvg-" + iterSampleStr;
        }
        tfout.writeln(header);
        for (int gene = 0; gene < geneNames.size(); gene++) {
            String ln = geneNames.get(gene);
            for (int iter = 0; iter < nrIter - 1; iter++) {

                ln += "\t" + averages[iter][gene] + "\t" + foldchanges[iter][gene];

            }
            tfout.writeln(ln);
        }

        tfout.close();

        System.out.println("Plotting MA plot");
        ScatterPlot2 scatter = new ScatterPlot2(averages, foldchanges, false, plotTitles, "Log2 average FPKM", "Log2 FoldChange (iteration/100%) FPKM",
                sampleName, 250, 250, 100, 100, 5, outdir + "MAPlots-" + sampleName + ".pdf");

    }

    private void correlateWithReadCounts(double[][] fpkm, int[] readCounts, int nrIter, String sampleName, String outdir) throws IOException {
        int nrGenes = fpkm[fpkm.length - 1].length;
        SpearmansCorrelation corr = new SpearmansCorrelation();

        double[] x = new double[nrIter];
        for (int iter = 0; iter < nrIter; iter++) {
            x[iter] = readCounts[iter];
        }

        TextFile out = new TextFile(outdir + "Stats-" + sampleName + ".txt", TextFile.W);
        out.writeln("Gene\tAvgReads\tAvgExp\tVarExp\tCorr");
        ArrayList<Double> xvals = new ArrayList<Double>();
        ArrayList<Double> yvals = new ArrayList<Double>();

        for (int gene = 0; gene < nrGenes; gene++) {
            double[] fpkmForGene = new double[nrIter];
            for (int iter = 0; iter < nrIter; iter++) {
                fpkmForGene[iter] = fpkm[iter][gene];
            }

            double correlation = corr.correlation(x, fpkmForGene);
            double avgExp = Descriptives.mean(fpkmForGene);
            double avgReads = Descriptives.mean(x);
            double sdExp = Math.sqrt(Descriptives.variance(fpkmForGene));

            if (!Double.isNaN(correlation)) {
                xvals.add(log2(avgExp));
                yvals.add(correlation);
            }

            out.writeln(geneNames.get(gene) + "\t" + avgReads + "\t" + avgExp + "\t" + sdExp + "\t" + correlation);

        }

        out.close();
        String outfilename = outdir + "ComparisonExpVsReadCountCorrelation-" + sampleName + ".pdf";
        String[] plotTitles = new String[]{sampleName};
        double[][] xd = new double[1][];
        double[][] yd = new double[1][];
        xd[0] = Primitives.toPrimitiveArr(xvals.toArray(new Double[0]));
        yd[0] = Primitives.toPrimitiveArr(yvals.toArray(new Double[0]));
        ScatterPlot2 plot = new ScatterPlot2(xd, yd, true, plotTitles, "Log2 Avg FPKM", "Correlation", "Expression vs Correlation with readcount", 500, 500, 100, 100, 1, outfilename);
    }

    private void correlateFullTableSamples(DoubleMatrixDataset<String, String> dsOut, int[] fullReadCounts, String outfile) throws Exception {
        int nrSamples = dsOut.columns();
        double[][] correlationMatrixWithNulls = new double[nrSamples][nrSamples];
        double[][] correlationMatrixWithoutNulls = new double[nrSamples][nrSamples];

        SpearmansCorrelation corr = new SpearmansCorrelation();
        for (int sample1 = 0; sample1 < nrSamples; sample1++) {
            double[] x = new double[dsOut.rows()];
            for (int row = 0; row < dsOut.rows(); row++) {
                x[row] = dsOut.getElement(row, sample1);
            }
            for (int sample2 = sample1 + 1; sample2 < nrSamples; sample2++) {
                double[] y = new double[dsOut.rows()];

                int nrNonNull = 0;
                for (int row = 0; row < dsOut.rows(); row++) {
                    y[row] = dsOut.getElement(row, sample2);
                    if (x[row] == 0 && x[row] == 0) {

                    } else {
                        nrNonNull++;
                    }
                }

                double[] xnonnull = new double[nrNonNull];
                double[] ynonnull = new double[nrNonNull];
                nrNonNull = 0;
                for (int row = 0; row < dsOut.rows(); row++) {

                    if (x[row] == 0 && x[row] == 0) {

                    } else {
                        xnonnull[nrNonNull] = x[row];
                        ynonnull[nrNonNull] = y[row];
                        nrNonNull++;
                    }
                }

                correlationMatrixWithNulls[sample1][sample2] = corr.correlation(x, y);
                correlationMatrixWithNulls[sample2][sample1] = correlationMatrixWithNulls[sample1][sample2];
                correlationMatrixWithoutNulls[sample1][sample2] = corr.correlation(xnonnull, ynonnull);
                correlationMatrixWithoutNulls[sample2][sample1] = correlationMatrixWithoutNulls[sample1][sample2];
            }
            correlationMatrixWithNulls[sample1][sample1] = 1;
            correlationMatrixWithoutNulls[sample1][sample1] = 1;
        }

        DoubleMatrixDataset ds = new DoubleMatrixDataset();
        ds.setMatrix(correlationMatrixWithNulls);
        ds.setColObjects(dsOut.getColObjects());
        ds.setRowObjects(dsOut.getColObjects());

        ds.save(outfile);

        ds.setMatrix(correlationMatrixWithoutNulls);
        ds.save(outfile + "-WithoutNulls.txt");

        TextFile out = new TextFile(outfile + "-GeneCorrelationWithReadCounts.txt", TextFile.W);
        out.writeln("Gene\tAverageExp\tMedianExp\tBeta\tSE\tCorrelation");

        double[] ctsDbl = new double[fullReadCounts.length];
        for (int q = 0; q < ctsDbl.length; q++) {
            ctsDbl[q] = fullReadCounts[q];
        }

        for (int row = 0; row < dsOut.rows(); row++) {
            String gene = dsOut.getRowObjects().get(row);
            double[] geneData = dsOut.getMatrix().viewRow(row).toArray();
            double[] coeff = Regression.getLinearRegressionCoefficients(geneData, ctsDbl);
            try {
                out.writeln(gene
                        + "\t" + JSci.maths.ArrayMath.mean(geneData)
                        + "\t" + JSci.maths.ArrayMath.median(geneData)
                        + "\t" + coeff[0]
                        + "\t" + coeff[2]
                        + "\t" + corr.correlation(geneData, ctsDbl));
            } catch (org.apache.commons.math3.exception.NotANumberException e) {

                System.out.println(gene + "\t" + row);
                for (double d : geneData) {
                    System.out.println(d);
                }

                e.printStackTrace();

            }

        }
        out.close();

        out = new TextFile(outfile + "-SamplesCorrelationWithReadCounts.txt", TextFile.W);
        out.writeln("Sample\tAverageExp\tMedianExp\tNrNonNull\tAverageExpNonNull\tMedianExpNonNull\tnrReads");
        for (int col = 0; col < dsOut.columns(); col++) {
            double[] sampleVals = dsOut.getMatrix().viewColumn(col).toArray();
            int nrNonNulls = 0;
            for (int row = 0; row < dsOut.rows(); row++) {

                if (sampleVals[row] != 0) {
                    nrNonNulls++;
                }
            }
            double[] sampleValsWoZeroes = new double[nrNonNulls];
            nrNonNulls = 0;
            for (int row = 0; row < dsOut.rows(); row++) {

                if (sampleVals[row] != 0) {
                    sampleValsWoZeroes[nrNonNulls] = sampleVals[row];
                    nrNonNulls++;
                }
            }

            out.writeln(dsOut.getColObjects().get(col)
                    + "\t" + JSci.maths.ArrayMath.mean(sampleVals)
                    + "\t" + JSci.maths.ArrayMath.median(sampleVals)
                    + "\t" + nrNonNulls
                    + "\t" + JSci.maths.ArrayMath.mean(sampleValsWoZeroes)
                    + "\t" + JSci.maths.ArrayMath.median(sampleValsWoZeroes)
                    + "\t" + fullReadCounts[col]
            );
        }
        out.close();
    }
}
