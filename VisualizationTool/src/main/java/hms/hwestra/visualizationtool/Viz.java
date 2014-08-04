/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.visualizationtool;

import hms.hwestra.utilities.features.Track;
import hms.hwestra.utilities.features.BedFileFeature;
import hms.hwestra.utilities.features.Chromosome;
import hms.hwestra.utilities.features.Strand;
import hms.hwestra.utilities.bedfile.BedFileReader;
import com.lowagie.text.DocumentException;
import hms.hwestra.utilities.features.Feature;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Locale;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import umcg.genetica.util.Primitives;

/**
 *
 * @author Harm-Jan
 */
public class Viz {

    private String[] geneNames;

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
    private static final Logger LOGGER = Logger.getLogger(Viz.class.getName());
    private static final Stroke dashed = new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0, new float[]{4}, 0);
    private static final Stroke line2pt = new BasicStroke(2, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
    private static final Stroke line = new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);

    private final int maxReadLength = 1000;

    public static void main(String[] args) {
        try {
            String file = "/Data/Gosia/bedfiles/UCSF-UBC.Breast_vHMEC.H3K4me1.RM035.HS1994.bed.gz";
            BedFileReader f = new BedFileReader();

            int width = 1000;
            int margin = 100;
            
            Chromosome chr = Chromosome.FIVE;
            int start = 56025316; // 56025316-56041464
            int stop = 56041464;

            // read all data
            Track t = f.read(file, "NaivePrimaryCD8", chr, start, stop, false);

//            start = 117294211;
//            stop = 117294211 + 5000000;
            ArrayList<Track> tracks = new ArrayList<Track>();
            tracks.add(t);

            Viz v = new Viz();
            String outputDir = "/Data/Gosia/viz2.png";
            try {
                v.viz(tracks, chr, start, stop, outputDir, Output.PNG);
            } catch (DocumentException ex) {
                Logger.getLogger(Viz.class.getName()).log(Level.SEVERE, null, ex);
            }
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    private void viz(ArrayList<Track> tracks, Chromosome chr, int start, int stop, String outputFileName, Output output) throws IOException, DocumentException {
        Locale defaultLocale = Locale.getDefault();
        Locale.setDefault(Locale.US);
        // set up Graphics2D depending on required format using iText in case PDF
        Graphics2D g2d = null;
        com.lowagie.text.Document document = null;
        com.lowagie.text.pdf.PdfWriter writer = null;
        BufferedImage bi = null;
        int width = 1;
        int height = 1;
        bi = new BufferedImage(1, 1, BufferedImage.TYPE_INT_RGB);
        g2d = bi.createGraphics();

        g2d.setFont(LARGE_FONT);
        FontMetrics fontmetrics = g2d.getFontMetrics();

        int betweenplotMargin = 10;
        int borderMargin = 100;
        int windowPlotSize = 1000;
        int distPlotSize = 200;

        int plotHeight = distPlotSize;
        int withinPlotMargin = 10;

        // width = plotwidth + max string size + 3 margins
        width = (3 * distPlotSize) + (3 * betweenplotMargin) + (2 * borderMargin) + windowPlotSize;
        height = (2 * withinPlotMargin) + plotHeight + (2 * borderMargin) + 500;

        System.out.println(width);
        // height is topmargin * 2 + yAxisNames * fontsize

        int fontheight = fontmetrics.getHeight();

        // initialize plot
        com.lowagie.text.pdf.PdfContentByte cb = null;
        if (output == Output.PDF) {
            com.lowagie.text.Rectangle rectangle = new com.lowagie.text.Rectangle(width, height);
            document = new com.lowagie.text.Document(rectangle);
            writer = com.lowagie.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(outputFileName));

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
        g2d.setStroke(line);
        // draw datasetnames
        g2d.setColor(Color.gray);

        System.out.println(height + " x " + width);
        int plot1XStart = borderMargin + withinPlotMargin;
        int plot1XStop = plot1XStart + windowPlotSize;
        int plot1YStart = borderMargin + withinPlotMargin;
        int plot1YStop = plot1YStart + plotHeight;
        System.out.println(plot1XStart);
        System.out.println(plot1XStop);
        System.out.println(plot1YStart);
        System.out.println(plot1YStop);

        g2d.drawRect(plot1XStart,
                plot1YStart,
                windowPlotSize, plotHeight);
        g2d.setColor(Color.darkGray);
//        g2d.drawLine(plot1XStart, plot1YStart + (plotHeight / 2), plot1XStop, plot1YStart + (plotHeight / 2));
        g2d.setColor(Color.gray);

        int nrBpInWindow = stop - start;
        int nrBins = 1000;
        if (nrBpInWindow < 1000) {
            nrBins = nrBpInWindow;
        }

        double[][] coverage = binReadDepth(tracks.get(0), chr, start, stop, nrBins);
        double maxReads = java.lang.Math.max(Primitives.max(coverage[0]), Primitives.max(coverage[1]));

        int bpPerBin = nrBpInWindow / nrBins;
        g2d.drawString("Coverage per bin. Showing " + start + " - " + stop
                + " (" + (nrBpInWindow) + "bp), " + bpPerBin + " bp/bin. 100% == " + maxReads + " reads for given base"
                + " Total reads in track: " + tracks.get(0).getNrReads(),
                plot1XStart, plot1YStart - withinPlotMargin);

        g2d.drawString("100%", plot1XStart - fontmetrics.stringWidth("100%"), plot1YStart + (fontmetrics.getHeight() / 2));
        g2d.drawString("0%", plot1XStart - fontmetrics.stringWidth("0%"), plot1YStart + (plotHeight / 2) + (fontmetrics.getHeight() / 4));
        g2d.drawString("100%", plot1XStart - fontmetrics.stringWidth("100%"), plot1YStart + plotHeight);

        int pixelWidthPerBin = windowPlotSize / nrBins;

        for (int bin = 0; bin < nrBins; bin += 100) {

            int tickX = start + (bpPerBin * bin);

            g2d.drawString("" + tickX, plot1XStart + (bin * pixelWidthPerBin), plot1YStop + 15);
        }
        for (int bin = 0; bin < nrBins; bin++) {
            double nrPos = coverage[1][bin];
            double nrNeg = coverage[0][bin];

            int heightNeg = (int) Math.ceil(((double) nrNeg / maxReads) * (plotHeight / 2));
            int heightPos = (int) Math.ceil(((double) nrPos / maxReads) * (plotHeight / 2));
            int diff = heightNeg - heightPos;

            System.out.println(bin + "\t" + nrPos + "\t" + nrNeg + "\t" + heightPos + "\t" + heightNeg);

            int alpha = 128 + (int) Math.ceil(127 * (double) nrPos / maxReads);
            g2d.setColor(new Color(255, 128, 0, alpha));
            g2d.fillRect(plot1XStart + (bin * pixelWidthPerBin), plot1YStart + (plotHeight / 2) - heightPos, pixelWidthPerBin, heightPos);
            alpha = 128 + (int) Math.ceil(127 * (double) nrNeg / maxReads);
            g2d.setColor(new Color(0, 128, 255, alpha));
            g2d.fillRect(plot1XStart + (bin * pixelWidthPerBin), plot1YStart + (plotHeight / 2), pixelWidthPerBin, heightNeg);

            if (diff > 0) {
                g2d.setColor(new Color(0, 0, 0, 255));
                g2d.fillRect(plot1XStart + (bin * pixelWidthPerBin), plot1YStart + (plotHeight / 2), pixelWidthPerBin, diff);
            } else {
                g2d.setColor(new Color(0, 0, 0, 255));
                diff = diff * -1;
                g2d.fillRect(plot1XStart + (bin * pixelWidthPerBin), plot1YStart + (plotHeight / 2) - diff, pixelWidthPerBin, diff);
            }
        }

        // visualize reads
        if (nrBpInWindow >= 10000000) {
            System.out.println("Not plotting reads");
        } else {
            System.out.println("Plotting Reads:");
            visualizeReads(g2d, tracks, plot1XStart, plot1YStop + borderMargin, chr, start, stop, windowPlotSize);
        }

        g2d.setColor(Color.gray);
        int plot2XStart = plot1XStop + withinPlotMargin;

        g2d.drawRect(plot2XStart,
                plot1YStart,
                distPlotSize, plotHeight);

        // read length distribution
        int nrReadLengthBins = 200;
        g2d.drawString("Read length distribution within window",
                plot2XStart, plot1YStart - withinPlotMargin);

        int[][] bins = binReadLength(tracks.get(0), tracks.get(0).getStart(), tracks.get(0).getStop(), nrReadLengthBins);
        int maxFrequency = java.lang.Math.max(Primitives.max(bins[0]), Primitives.max(bins[1]));
        pixelWidthPerBin = plotHeight / nrReadLengthBins;
        for (int bin = 0; bin < nrReadLengthBins; bin++) {
            int nrPos = bins[1][bin];
            int nrNeg = bins[0][bin];

            int heightNeg = (int) Math.ceil(((double) nrNeg / maxFrequency) * plotHeight);
            int heightPos = (int) Math.ceil(((double) nrPos / maxFrequency) * plotHeight);

            g2d.setColor(new Color(255, 128, 0, 128));
            g2d.drawRect(plot2XStart + (bin * pixelWidthPerBin), plot1YStart + plotHeight - heightPos, pixelWidthPerBin, heightPos);
            g2d.setColor(new Color(0, 128, 255, 128));
            g2d.drawRect(plot2XStart + (bin * pixelWidthPerBin), plot1YStart + plotHeight - heightNeg, pixelWidthPerBin, heightNeg);
        }
        g2d.setColor(Color.gray);
        bpPerBin = maxReadLength / nrReadLengthBins;
        for (int bin = 0; bin < nrReadLengthBins; bin += nrReadLengthBins / 4) {
            int tickX = 0 + (bpPerBin * bin);

            g2d.drawString("" + tickX, plot2XStart + (bin * pixelWidthPerBin), plot1YStop + 15);
        }

        // coverage distribution
        // first get the max average 
        int plot3XStart = plot2XStart + distPlotSize + withinPlotMargin;
        int maxReadDepth = plotHeight;
        int[][] frequencyDist = binCoverage(tracks.get(0), chr, start, stop, maxReadDepth);
        maxFrequency = java.lang.Math.max(Primitives.max(frequencyDist[0]), Primitives.max(frequencyDist[1]));
        g2d.drawString("Coverage in window.",
                plot3XStart, plot1YStart - withinPlotMargin);

        
        
        
        
        g2d.drawRect(plot3XStart,
                plot1YStart,
                distPlotSize, plotHeight);

        int plot4XStart = plot3XStart + distPlotSize + withinPlotMargin;

        g2d.drawRect(plot4XStart,
                plot1YStart,
                distPlotSize, plotHeight);

        g2d.dispose();
        if (output == Output.PDF) {
            g2d.dispose();
            cb.restoreState();
            document.close();
            writer.close();
        } else {
            bi.flush();
            ImageIO.write(bi, output.toString().toLowerCase(), new File(outputFileName));
        }

        Locale.setDefault(defaultLocale);
    }

    private double determineUnit(double range) {

        double divisor = Math.log10(range);
        divisor = Math.floor(divisor);
        divisor = Math.pow(10, divisor);
        return divisor;
    }

    private double[][] binReadDepth(Track track, Chromosome chr, int start, int stop, int nrBins) {
        int bpPerBin = (stop - start) / nrBins;
        double[][] bins = new double[2][nrBins];

        Set<Feature> features = track.getFeatureSet(chr, start, stop);
        for (Feature f : features) {
            int fstart = f.getStart();
            int fstop = f.getStop();
            Strand s = f.getStrand();
            if (fstart < start) {
                fstart = start;
            }

            if (fstart > start && fstop > stop && fstart < stop) {
                fstop = stop;
            }
            if (fstart >= start && fstop <= stop) {
                for (int i = fstart; i < fstop; i++) {
                    int binNo = (int) Math.ceil((i - start) / bpPerBin);
                    if (binNo >= nrBins) {
                        binNo = nrBins - 1;
                    }
                    if (binNo < 0) {
                        binNo = 0;
                    }
                    

                    bins[s.getNumber()][binNo]++;

                }
            }
        }
        for (int bin = 0; bin < nrBins; bin++) {
            System.out.println(bin + "\t" + bins[0][bin] + "\t" + bins[1][bin] + "\t" + (bins[0][bin] / bpPerBin) + "\t" + (bins[1][bin] / bpPerBin));
            bins[0][bin] /= bpPerBin;
            bins[1][bin] /= bpPerBin;
        }
        return bins;
    }

    private int[][] binReadLength(Track track, int start, int stop, int nrBins) {

        int[][] bins = new int[2][nrBins];
        int maxLength = maxReadLength;

        for (Feature f : track.getFeatures()) {
            int fstart = f.getStart();
            int fstop = f.getStop();

            Strand s = f.getStrand();
            if (fstart < start && fstop < stop && fstop > start) {
                fstart = start;
            }

            if (fstart > start && fstop > stop && fstart < stop) {
                fstop = stop;
            }
            if (fstart >= start && fstop <= stop) {
                fstart = f.getStart();
                fstop = f.getStop();
                int len = fstop - fstart;
                int binNo = (int) Math.ceil(((double) len / maxLength) * nrBins);

                if (binNo >= nrBins) {
                    binNo = nrBins - 1;
                }

                bins[s.getNumber()][binNo]++;

            }
        }

        return bins;
    }

    private void visualizeReads(Graphics2D g2d, ArrayList<Track> tracks, int plotXStart, int plotYStart, Chromosome chr, int windowStart, int windowEnd, int windowPixelSize) {
        Track track = tracks.get(0);
        // + strand features first

        int y = 0;
        Set<Feature> features = null;
        features = track.getFeatureSet(chr, windowStart, windowEnd);

        System.out.println(features.size() + " features found within window...");
        HashSet<Feature> visitedFeatures = new HashSet<Feature>();
        int bpInWindow = windowEnd - windowStart;
        int featHeight = 5;
        int featMargin = 3;
//        boolean plot = false;
        while (visitedFeatures.size() != features.size()) {

            Set<Feature> otherFeaturesAtY = new HashSet<Feature>();
            int yStart = plotYStart + (y * featHeight) + (y * featMargin);
            for (Feature f : features) {

                Color boxcolor = null;
                Color linecolor = null;
                if (f.getStrand() == Strand.NEG) {
                    boxcolor = new Color(0, 128, 255, 128);
                    linecolor = new Color(0, 128, 255);
//                    plot = false;
//                    visitedFeatures.add(f);
                } else {
                    boxcolor = new Color(255, 128, 0, 128);
                    linecolor = new Color(255, 128, 0);
//                    plot = true;
                }

                if (!visitedFeatures.contains(f)) {

                    if (otherFeaturesAtY.isEmpty()) {
                        // plot 
                        int pixelEnd = 0;
                        int pixelStart = 0;
                        if (f.getStart() <= windowStart) {
                            pixelStart = 0;
                        } else {
                            pixelStart = (int) Math.ceil(((double) (f.getStart() - windowStart) / bpInWindow) * windowPixelSize);
                        }

                        if (f.getStop() >= windowEnd) {
                            pixelEnd = windowPixelSize;
                        } else {
                            pixelEnd = (int) Math.ceil(((double) (f.getStop() - windowStart) / bpInWindow) * windowPixelSize);
                        }

                        g2d.setColor(boxcolor);
                        g2d.fillRect(plotXStart + pixelStart, yStart, pixelEnd - pixelStart, featHeight);
                        g2d.setColor(linecolor);
                        g2d.drawRect(plotXStart + pixelStart, yStart, pixelEnd - pixelStart, featHeight);
                        //System.out.println((plotXStart + pixelStart) + "\t" + yStart + "\t" + (pixelEnd - pixelStart) + "\t" + visitedFeatures.size() + "/" + features.size());
                        visitedFeatures.add(f);
                        otherFeaturesAtY.add(f);

                    } else {
                        boolean overlaps = false;
                        for (Feature otherFeature : otherFeaturesAtY) {
                            if (f.overlaps(otherFeature)) {
                                overlaps = true;
                            }
                        }
                        if (!overlaps) {
                            // plot 
                            int pixelEnd = 0;
                            int pixelStart = 0;
                            if (f.getStart() <= windowStart) {
                                pixelStart = 0;
                            } else {
                                pixelStart = (int) Math.ceil(((double) (f.getStart() - windowStart) / bpInWindow) * windowPixelSize);
                            }

                            if (f.getStop() >= windowEnd) {
                                pixelEnd = windowPixelSize;
                            } else {
                                pixelEnd = (int) Math.ceil(((double) (f.getStop() - windowStart) / bpInWindow) * windowPixelSize);
                            }

                            g2d.setColor(boxcolor);
                            g2d.fillRect(plotXStart + pixelStart, yStart, pixelEnd - pixelStart, featHeight);
                            g2d.setColor(linecolor);
                            g2d.drawRect(plotXStart + pixelStart, yStart, pixelEnd - pixelStart, featHeight);
//                            System.out.println((plotXStart + pixelStart) + "\t" + yStart + "\t" + (pixelEnd - pixelStart) + "\t" + visitedFeatures.size() + "/" + features.size());
                            visitedFeatures.add(f);
                            otherFeaturesAtY.add(f);

                        }

                    }

                }

            }

            System.out.println("Last feature unset. Y ++ - level " + y);
            if (y % 10 == 0) {
                g2d.setColor(new Color(125, 125, 125, 125));
                g2d.drawLine(plotXStart - 10, yStart - 2, plotXStart + windowPixelSize + 10, yStart - 2);
                g2d.drawString("" + y, plotXStart + windowPixelSize + 15, yStart - 2);
            }
            y++;
        }

    }

    private int[][] binCoverage(Track get, Chromosome chr, int start, int stop, int maxReadDepth) {
        int[][] bins = new int[2][maxReadDepth];

        int windowsize = 10000;
        if (start + windowsize > stop) {
            windowsize = stop - start;
        }
        for (int i = start; i < stop; i += windowsize) {
            // get feats
            int[] negBaseSeen = new int[windowsize];
            int[] posBaseSeen = new int[windowsize];

            Set<Feature> features = get.getFeatureSet(chr, i, i + windowsize);
            for (Feature f : features) {
                int fstart = f.getStart();
                int fstop = f.getStop();
                Strand s = f.getStrand();
                int windowStart = fstart - i;
                if (windowStart < 0) {
                    windowStart = 0;
                }
                int windowStop = fstop - i;
                if (windowStop > windowsize - 1) {
                    windowStop = windowsize - 1;
                }

                for (int q = windowStart; q < windowStop; q++) {
                    if (s == Strand.NEG) {
                        negBaseSeen[q]++;
                    } else {
                        posBaseSeen[q]++;
                    }
                }

            }

            for (int j = 0; j < windowsize; j++) {
                int num1 = negBaseSeen[j];
                int num2 = posBaseSeen[j];
                if (num1 > maxReadDepth) {
                    num1 = maxReadDepth;
                }
                if (num2 > maxReadDepth) {
                    num2 = maxReadDepth;
                }

                bins[Strand.NEG.getNumber()][num1]++;
                bins[Strand.POS.getNumber()][num2]++;
            }
        }

        return bins;
    }
}
