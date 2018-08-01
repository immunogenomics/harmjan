/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.rnaaccess;

import com.itextpdf.text.DocumentException;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import nl.harmjanwestra.utilities.legacy.genetica.graphics.ViolinBoxPlot;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.math.matrix2.DoubleMatrixDataset;
import nl.harmjanwestra.utilities.legacy.genetica.math.matrix2.MatrixHandling;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.Descriptives;
import nl.harmjanwestra.utilities.legacy.genetica.util.Primitives;

/**
 *
 * @author Harm-Jan
 */
public class CombineTrackingFiles {

    public static void main(String[] args) {
//        String inDir = "/Data/Projects/2014-Epipilot/rna-seq/tracking/samples/";
//        String selection = "genes.fpkm_tracking";
//        String outdir = "/Data/Projects/2014-Epipilot/rna-seq/tracking/";
//        String outFile = outdir + "genes.fpkm_tracking_combined";
        CombineTrackingFiles f = new CombineTrackingFiles();
        try {
//            File indir = new File(inDir);
//
//            File[] files = indir.listFiles();
//            f.combineTrackingFiles(files, selection, outdir, outFile);
//            double maxPercMissing = 1;
//            DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(outFile + ".txt");
//            HashSet<String> rowsToInclude = new HashSet<String>();
//            for (int r = 0; r < ds.rows(); r++) {
//                int missing = 0;
//                for (int c = 0; c < ds.columns(); c++) {
//                    if (ds.getElement(r, c) == 0d || Double.isNaN(ds.getElement(r, c))) {
//                        missing++;
//                    }
//                }
//
//                if (missing == ds.columns()) {
//                    System.out.println(ds.getRowObjects().get(r) + "\t" + r);
//                } else {
//                    rowsToInclude.add(ds.getRowObjects().get(r));
//                }
//            }
//            ds = DoubleMatrixDataset.loadSubsetOfTextDoubleData(outFile + ".txt", "\t", rowsToInclude, null);
//            ds.save(outFile + "-withoutZeros.txt");
        String normalizedDataset = "/Data/Projects/2014-Epipilot/rna-seq/tracking/genes.fpkm_tracking_combined-withoutZeros.ProbesWithZeroVarianceRemoved.CovariatesRemoved.txt.gz";
        
////            
//
//            // C:\Work\RNAAccess\CufflinksSamples\QNNormalized
//            String normalizedDataset = inDir + "\\genes.fpkm_tracking_combined.txt";
////            String normalizedDataset = inDir + "\\QNNormalized\\genes.fpkm_tracking_combined_RemovedGenesWithMaxPercZero-0.95.ProbesWithZeroVarianceRemoved.QuantileNormalized.txt.gz";
////            String normalizedDataset = inDir + "\\QNNormalized\\genes.fpkm_tracking_combined_RemovedGenesWithMaxPercZero-0.95.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz";
////            
            DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(normalizedDataset);
//            boolean trimZeros = false;
//            if (trimZeros) {
//                String trimmedFile = outFile + "_RemovedGenesWithMaxPercZero-" + maxPercMissing + ".txt";
//                ds = f.removeGenesWithNullExpression(0.95, ds, trimmedFile);
//                normalizedDataset = trimmedFile;
//            }
//            
            String vbpOutFileName = normalizedDataset + "_SampleDistributions.pdf";
            f.drawViolinPlotsPerSample(ds, vbpOutFileName, false);
            String correlationOutput = normalizedDataset + "_CorrelationMatrix";
            f.correlateSamples(ds, correlationOutput);
            String histogramOut = normalizedDataset + "_Histograms.pdf";
            boolean binZeros = true;
            double max = 500;
            boolean log10 = true;
            f.drawHistograms(ds, histogramOut, binZeros, max, log10, Output.PDF);
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void combineTrackingFiles(File[] inFiles, String selection, String outdir, String outFile) throws IOException, Exception {

        HashMap<String, Integer> geneIndex = new HashMap<String, Integer>();

        ArrayList<String> genes = new ArrayList<String>();
        ArrayList<String> samples = new ArrayList<String>();
        HashMap<String, Integer> sampleIndex = new HashMap<String, Integer>();

        for (File f : inFiles) {

            if (sampleIndex.containsKey(f.getName())) {
                System.out.println("Sample is already processed.");

            } else {
                samples.add(f.getName());
                sampleIndex.put(f.getName(), sampleIndex.size());
            }

            String fileName = f.getName();
            if (fileName.contains(selection)) {
                TextFile tf = new TextFile(f, TextFile.R);
                String[] lineElems = tf.readLineElems(TextFile.tab);
                int geneCol = 0;
                int trackingCol = 0;
                int tssIdCol = 0;
                int fpkmCol = 0;
                int fpkmHiCol = 0;
                int fpkmLoCol = 0;
                System.out.println("Header length: " + lineElems.length);
                for (int i = 0; i < lineElems.length; i++) {
                    String e = lineElems[i];
                    if (e.equals("gene_id")) {
                        geneCol = i;
                    }
                    //tss_id
                    if (e.equals("locus")) {
                        tssIdCol = i;
                    }
                    if (e.equals("tracking_id")) {
                        trackingCol = i;
                    }
                    if (e.equals("FPKM")) {
                        fpkmCol = i;
                    }
                    if (e.equals("FPKM_conf_hi")) {
                        fpkmHiCol = i;
                    }
                    if (e.equals("FPKM_conf_lo")) {
                        fpkmLoCol = i;
                    }
                }
                System.out.println(f.getName());
                System.out.println("geneCol\t" + geneCol);
                System.out.println("trackingCol\t" + trackingCol);
                System.out.println("fpkmCol\t" + fpkmCol);
                System.out.println("fpkmHiCol\t" + fpkmHiCol);
                System.out.println("fpkmLoCol\t" + fpkmLoCol);
                lineElems = tf.readLineElems(TextFile.tab);

                while (lineElems != null) {

                    String gene = lineElems[geneCol];
                    String tssId = lineElems[tssIdCol];
                    String track = lineElems[trackingCol];
                    String fpkm = lineElems[fpkmCol];
                    String fpkmHi = lineElems[fpkmHiCol];
                    String fpkmLo = lineElems[fpkmLoCol];

                    double dFpkm = Double.parseDouble(fpkm);
                    double dFpkmHi = Double.parseDouble(fpkmHi);
                    double dFpkmLo = Double.parseDouble(fpkmLo);

                    String combinedId = track + "_" + tssId;
                    if (geneIndex.containsKey(combinedId)) {
//                                System.out.println("Gene already processed: " + combinedId);
                    } else {
                        genes.add(combinedId);
                        geneIndex.put(combinedId, geneIndex.size());
                    }
                    lineElems = tf.readLineElems(TextFile.tab);
                }
                tf.close();
            }

        }

        System.out.println(geneIndex.size() + " genes detected. ");
        System.out.println(sampleIndex.size() + " samples detected");

        double[][] fpkms = new double[geneIndex.size()][sampleIndex.size()];
        double[][] fpkmsLo = new double[geneIndex.size()][sampleIndex.size()];
        double[][] fpkmsHi = new double[geneIndex.size()][sampleIndex.size()];

        for (File f : inFiles) {

            Integer sampleId = sampleIndex.get(f.getName());

            String fileName = f.getName();
            if (fileName.contains(selection)) {
                TextFile tf = new TextFile(f, TextFile.R);
                String[] lineElems = tf.readLineElems(TextFile.tab);
                int geneCol = 0;
                int trackingCol = 0;
                int tssIdCol = 0;
                int fpkmCol = 0;
                int fpkmHiCol = 0;
                int fpkmLoCol = 0;
                System.out.println("Header length: " + lineElems.length);
                for (int i = 0; i < lineElems.length; i++) {
                    String e = lineElems[i];
                    if (e.equals("gene_id")) {
                        geneCol = i;
                    }
                    //tss_id
                    if (e.equals("locus")) {
                        tssIdCol = i;
                    }
                    if (e.equals("tracking_id")) {
                        trackingCol = i;
                    }
                    if (e.equals("FPKM")) {
                        fpkmCol = i;
                    }
                    if (e.equals("FPKM_conf_hi")) {
                        fpkmHiCol = i;
                    }
                    if (e.equals("FPKM_conf_lo")) {
                        fpkmLoCol = i;
                    }
                }
                System.out.println(f.getName());
                System.out.println("geneCol\t" + geneCol);
                System.out.println("trackingCol\t" + trackingCol);
                System.out.println("fpkmCol\t" + fpkmCol);
                System.out.println("fpkmHiCol\t" + fpkmHiCol);
                System.out.println("fpkmLoCol\t" + fpkmLoCol);
                lineElems = tf.readLineElems(TextFile.tab);

                while (lineElems != null) {

                    String gene = lineElems[geneCol];
                    String tssId = lineElems[tssIdCol];
                    String track = lineElems[trackingCol];
                    String fpkm = lineElems[fpkmCol];
                    String fpkmHi = lineElems[fpkmHiCol];
                    String fpkmLo = lineElems[fpkmLoCol];

                    double dFpkm = Double.parseDouble(fpkm);
                    double dFpkmHi = Double.parseDouble(fpkmHi);
                    double dFpkmLo = Double.parseDouble(fpkmLo);

                    String combinedId = track + "_" + tssId;

                    Integer geneId = geneIndex.get(combinedId);

                    fpkms[geneId][sampleId] = dFpkm;
                    fpkmsLo[geneId][sampleId] = dFpkmLo;
                    fpkmsHi[geneId][sampleId] = dFpkmHi;

                    lineElems = tf.readLineElems(TextFile.tab);
                }
                tf.close();
            }

        }

        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>();
        ds.setMatrix(fpkms);
        ds.setRowObjects(genes);
        ds.setColObjects(samples);
        ds.save(outFile + ".txt");
        System.out.println(outFile + ".txt");
        ds = new DoubleMatrixDataset<String, String>();
        ds.setMatrix(fpkmsLo);
        ds.setRowObjects(genes);
        ds.setColObjects(samples);
        ds.save(outFile + "_lo.txt");
        System.out.println(outFile + "_lo.txt");
        ds = new DoubleMatrixDataset<String, String>();
        ds.setMatrix(fpkmsHi);
        ds.setRowObjects(genes);
        ds.setColObjects(samples);
        ds.save(outFile + "_hi.txt");
        System.out.println(outFile + "_hi.txt");
    }

    private DoubleMatrixDataset<String, String> removeGenesWithNullExpression(double maxPercentageMissing, DoubleMatrixDataset<String, String> ds, String outfile) throws IOException, Exception {
        MatrixHandling.RemoveRowsWithToManyMissingValues(ds, (int) Math.ceil(maxPercentageMissing * ds.getColObjects().size()), 0);
        ds.save(outfile);
        return ds;
    }

    private void drawViolinPlotsPerSample(DoubleMatrixDataset<String, String> ds, String outfileName, boolean log10transform) throws IOException {
        ViolinBoxPlot vbp = new ViolinBoxPlot();

        double[][][] vals = new double[1][ds.columns()][ds.rows()];
        String[][] xLabels = new String[1][ds.columns()];
        for (int c = 0; c < ds.columns(); c++) {
            xLabels[0][c] = ds.getColObjects().get(c);
        }
        String[] datasetNames = new String[]{"RNA-Access"};
        for (int c = 0; c < ds.columns(); c++) {
            for (int r = 0; r < ds.rows(); r++) {
                if (log10transform) {
                    vals[0][c][r] = Math.log10(ds.getElement(r, c));
                } else {
                    double elem = ds.getElement(r, c);
                    if (elem > 250) {
                        elem = 250;
                    }
                    vals[0][c][r] = elem;
                }
            }
        }

        vbp.draw(vals, datasetNames, xLabels, "FPKM", ViolinBoxPlot.Output.PDF, outfileName);
    }

    public void correlateSamples(DoubleMatrixDataset<String, String> ds, String outputFileName) throws Exception {

        double[][] correlationMatrix1 = new double[ds.columns()][ds.columns()];
        double[][] correlationMatrix2 = new double[ds.columns()][ds.columns()];

        for (int s0 = 0; s0 < ds.columns(); s0++) {
            for (int s1 = s0; s1 < ds.columns(); s1++) {

                double[] sample0 = ds.getMatrix().viewColumn(s0).toArray();
                double[] sample1 = ds.getMatrix().viewColumn(s1).toArray();
                int ctr = 0;
                int ctr1 = 0;
                int ctr2 = 0;
                for (int i = 0; i < sample0.length; i++) {
                    if (sample0[i] == 0d) {
                        ctr++;
                    }
                    if (sample1[i] == 0d) {
                        ctr1++;
                    }

                    if (sample1[i] == 0d && sample0[i] == 0d) {
                        ctr2++;

                    }
                }

                double[] xval0 = new double[sample1.length - ctr2];
                double[] xval1 = new double[sample1.length - ctr2];
                int q = 0;
                for (int i = 0; i < sample0.length; i++) {

                    if (sample1[i] == 0d && sample0[i] == 0d) {

                    } else {
                        xval0[q] = sample0[i];
                        xval1[q] = sample1[i];
                        q++;
                    }
                }

//                System.out.println(ctr + "\t" + ctr1 + "\t" + ctr2 + "\t" + sample0.length);
//                System.out.println(Correlation.correlate(sample0, sample1));
                SpearmansCorrelation corr = new SpearmansCorrelation();

                double correlation = corr.correlation(sample0, sample1);
                System.out.println(correlation);

                double correlationWithoutZeroes = corr.correlation(xval0, xval1);

                System.out.println(ds.getColObjects().get(s0) + " - " + ds.getColObjects().get(s1) + "\t" + correlation + "\t" + correlationWithoutZeroes);

                correlationMatrix1[s0][s1] = correlation;
                correlationMatrix2[s0][s1] = correlationWithoutZeroes;

                correlationMatrix1[s1][s0] = correlation;
                correlationMatrix2[s1][s0] = correlationWithoutZeroes;

            }
        }

        DoubleMatrixDataset<String, String> dsOut1 = new DoubleMatrixDataset<>(ds.columns(), ds.columns());
        dsOut1.setMatrix(correlationMatrix1);
        dsOut1.setColObjects(ds.getColObjects());
        dsOut1.setRowObjects(ds.getColObjects());
        dsOut1.save(outputFileName + ".txt");
        DoubleMatrixDataset<String, String> dsOut2 = new DoubleMatrixDataset<>(ds.columns(), ds.columns());
        dsOut2.setMatrix(correlationMatrix2);
        dsOut2.setColObjects(ds.getColObjects());
        dsOut2.setRowObjects(ds.getColObjects());
        dsOut2.save(outputFileName + "_removedZeros.txt");
    }

    public enum Output {

        PDF, PNG
    };
    private static final Font LARGE_FONT = new Font("Verdana", Font.PLAIN, 14);
    private static final Font LARGE_FONT_BOLD = new Font("Verdana", Font.BOLD, 14);
    private static final Font SMALL_FONT = new Font("Verdana", Font.PLAIN, 10);
    private static final int upMaxR = 108;
    private static final int upMaxG = 189;
    private static final int upMaxB = 69;

    private static final int downMaxR = 0;
    private static final int downMaxG = 174;
    private static final int downMaxB = 239;

    private static final float downMaxH = 196f / 360f;
    private static final float downMaxS = 1f;
    private static final float downMaxBr = 0.93f;

    private static final float upMaxH = 100f / 360f;
    private static final float upMaxS = 0.63f;
    private static final float upMaxBr = 0.74f;

    private static final Color red1Color = new Color(242, 101, 94, 50);
    private static final Color red2Color = new Color(242, 101, 94, 25);
    private static final Color gray1Color = new Color(237, 237, 237);
    private static final Color gray2Color = new Color(247, 247, 247);
    private static final Logger LOGGER = Logger.getLogger(CombineTrackingFiles.class.getName());
    private static final Stroke dashed = new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0, new float[]{4}, 0);
    private static final Stroke line2pt = new BasicStroke(2, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
    private static final Stroke line = new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);

    public void drawHistograms(DoubleMatrixDataset<String, String> ds, String histogramOut, boolean binZeros, double maxVal, boolean log10, Output output) throws DocumentException, FileNotFoundException, IOException {
        int nrHistograms = ds.columns();
        int nrCols = 5;
        int nrRows = nrHistograms / nrCols;

        int marginBetween = 100;
        int margin = 100;
        int plotSize = 250;
        int width = (margin * 2) + (nrCols * plotSize) + ((nrCols - 1) * marginBetween);
        int height = (margin * 2) + nrRows * plotSize + ((nrRows - 1) * marginBetween);

        Locale defaultLocale = Locale.getDefault();
        Locale.setDefault(Locale.US);
        // set up Graphics2D depending on required format using iText in case PDF
        Graphics2D g2d = null;
        com.itextpdf.text.Document document = null;
        com.itextpdf.text.pdf.PdfWriter writer = null;
        BufferedImage bi = null;

        bi = new BufferedImage(1, 1, BufferedImage.TYPE_INT_RGB);
        g2d = bi.createGraphics();

        g2d.setFont(LARGE_FONT);
        com.itextpdf.text.pdf.PdfContentByte cb = null;
        if (output == Output.PDF) {
            com.itextpdf.text.Rectangle rectangle = new com.itextpdf.text.Rectangle(width, height);
            document = new com.itextpdf.text.Document(rectangle);
            writer = com.itextpdf.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(histogramOut));

            document.open();
            cb = writer.getDirectContent();
            cb.saveState();

            g2d = cb.createGraphics(width, height);
        } else {
            bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
            g2d = bi.createGraphics();
        }

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(Color.white);
        g2d.fillRect(0, 0, width, height);

        // bin data
        int minVal = 0;
        int nrBins = plotSize;
        double maxCtr = 0;
        int[][] allBins = new int[ds.columns()][nrBins];
        double[] maxVals = new double[ds.columns()];
        double[] minVals = new double[ds.columns()];
        double[] medianVals = new double[ds.columns()];
        double[] meanVals = new double[ds.columns()];

        DescriptiveStatistics s = new DescriptiveStatistics();

        for (int sample = 0; sample < ds.columns(); sample++) {
            double[] sampleVals = ds.getMatrix().viewColumn(sample).toArray();
            maxVals[sample] = Primitives.max(sampleVals);
            minVals[sample] = Primitives.min(sampleVals);
            meanVals[sample] = Descriptives.mean(sampleVals);
            medianVals[sample] = JSci.maths.ArrayMath.percentile(sampleVals, 0.50d);
            allBins[sample] = binData(sampleVals, nrBins, minVal, maxVal, binZeros);
            for (int bin = 0; bin < nrBins; bin++) {

                if (allBins[sample][bin] > maxCtr) {
                    maxCtr = allBins[sample][bin];
                }

            }
        }

        System.out.println(width + "\t" + height);
        // plot histograms
        for (int sample = 0; sample < ds.columns(); sample++) {
            g2d.setColor(new Color(0, 128, 255));
            int row = sample / nrCols;
            int col = sample % nrCols;

            int startX = margin + (col * marginBetween) + (col * plotSize);
            int startY = margin + (row * marginBetween) + (row * plotSize);

            // System.out.println(sample + "\t" + row + "\t" + col + "\t" + x0 + "\t" + y0);
            // draw bars
            int barWidth = (plotSize / nrBins);
            for (int bin = 0; bin < nrBins; bin++) {

                int xPos = startX + (barWidth * bin);
                double val = allBins[sample][bin];
                double perc = val / maxCtr;
                if (log10) {
                    perc = Math.log10(val) / Math.log10(maxCtr);

                }

                int barHeigth = (int) Math.ceil(perc * plotSize);
                int yPos = startY + (plotSize - barHeigth);
                g2d.fillRect(xPos, yPos, barWidth, barHeigth);
            }
            FontMetrics fontmetrics = g2d.getFontMetrics();

            // draw Axis values
            DecimalFormat formatter = new java.text.DecimalFormat("##.###;-##.###", new java.text.DecimalFormatSymbols(java.util.Locale.US));
            String maxYValString = formatter.format(maxCtr);
            String minYValStr = formatter.format(0);
            // y
            g2d.drawString(maxYValString, startX - fontmetrics.stringWidth(maxYValString) - 15, startY + 5);
            g2d.drawString(minYValStr, startX - fontmetrics.stringWidth(minYValStr) - 15, startY + plotSize);

//X
            String maxXValString = formatter.format(maxVal);
            g2d.drawString(minYValStr, startX, startY + plotSize + 25);
            g2d.drawString(maxXValString, startX + plotSize - fontmetrics.stringWidth(maxXValString) + 5, startY + plotSize + 25);

            g2d.setColor(new Color(0, 0, 0));
            g2d.drawString(ds.getColObjects().get(sample), startX, startY - 30);

            // draw statistics and sample name
            String maxString = formatter.format(maxVals[sample]);
            String minString = formatter.format(minVals[sample]);
            String medianString = formatter.format(medianVals[sample]);
            String meanString = formatter.format(meanVals[sample]);
            String subtitle = "Min: " + minString + " Max: " + maxString + " Mean: " + meanString + " Median: " + medianString;
            g2d.drawString(subtitle, startX, startY - 15);

// draw box
            g2d.setColor(new Color(0, 0, 0, 64));
            g2d.drawRect(startX - 5, startY - 5, plotSize + 10, plotSize + 10);
        }
        g2d.dispose();
        if (output == Output.PDF) {
            g2d.dispose();
            cb.restoreState();
            document.close();
            writer.close();
        } else {
            bi.flush();
            ImageIO.write(bi, output.toString().toLowerCase(), new File(histogramOut));
        }

        Locale.setDefault(defaultLocale);
    }

    private int[] binData(double[] sampleVals, int nrBins, double minVal, double maxVal, boolean binZeros) {
        double range = maxVal - minVal;
        int[] bins = new int[nrBins];
        for (int i = 0; i < sampleVals.length; i++) {
            double v = sampleVals[i];
            if (v > maxVal) {
                v = maxVal;
            }
            if (v < minVal) {
                v = minVal;
            }
            int bin = (int) Math.ceil(((double) (v - minVal) / range) * nrBins);
            if (bin >= nrBins) {
                bin = nrBins - 1;
            }
            if (binZeros && v == 0) {

            } else {
                bins[bin]++;
            }
        }
        return bins;
    }
}
