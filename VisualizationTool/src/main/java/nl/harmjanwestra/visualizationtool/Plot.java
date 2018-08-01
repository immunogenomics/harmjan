/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.visualizationtool;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Exon;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.enums.Strand;
import nl.harmjanwestra.utilities.features.Track;
import nl.harmjanwestra.utilities.features.Transcript;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author hwestra
 */
public class Plot extends DefaultGraphics {

    private ArrayList<Locus> loci;
    private Track[][] bedFileTracks;
    int maxFigureWidth = 1000;
    int margin = 100;
    int betweenlocusMargin = 100;

    double bpPerPixel = 0;
    private int maxLocusNrPixels = 0;

    int betweenlocusElementMargin = 20;
    int locusGeneHeight = 10;
    int betweenLocusGeneMargin = 10;
    int locusAnnotationHeight = 200;
    int maxReadDepth = 20;
    Color grey = new Color(128, 128, 128);
    Color lightgrey = new Color(196, 196, 196);

    public void plot(ArrayList<Locus> loci, Track[][] bedFileTracks, Track[] peakTracks, String outputFileName) throws DocumentException, FileNotFoundException, IOException {
        this.loci = loci;
        this.bedFileTracks = bedFileTracks;

        determineFigureSize();

        initializePlot(outputFileName, figureWidth, figureHeight);

        // plot the actual things.
        int nextYStartPosition = margin;

        g2d.drawString(bpPerPixel + " bp/pixel. Max nr bp: " + maxLocusNrPixels, margin, 20);

        for (int locusId = 0; locusId < loci.size(); locusId++) {
            g2d.setColor(grey);
            int initialYStart = nextYStartPosition;
            Locus locus = loci.get(locusId);

            Track[] bedfileTracksForLocus = bedFileTracks[locusId];
            Iterable<Feature> peaks = peakTracks[locusId].getFeatures();

            int locusStart = locus.getLeftBound();
            int locusStop = locus.getRightBound();
            int locusWidth = locusStop - locusStart;

            int locusPixelWidth = (int) Math.ceil((double) locusWidth / bpPerPixel); // nrbins

            // draw line
            int xStart = margin;
            int xStop = margin + locusPixelWidth;

            // plot the variants
            int x1 = locus.getStart();
            int y = initialYStart - 20;
            int x = xStart + (int) Math.ceil((double) (x1 - locusStart) / bpPerPixel);

            int annotationstartY = initialYStart;
            int annotationstopY = (bedfileTracksForLocus.length * locusAnnotationHeight) + (bedfileTracksForLocus.length * betweenlocusElementMargin);

            g2d.setColor(new Color(0, 128, 255));
            for (Feature f : peaks) {
                if (f.overlaps(locus)) {
                    g2d.setColor(new Color(255, 0, 128));
                }
            }
            drawRotatedText(x, y - 10, locus.getName(), -45);
            g2d.drawLine(x, y, x, annotationstartY + annotationstopY);

            ArrayList<Feature> proxies = locus.getProxies();

            for (Feature proxy : proxies) {
                g2d.setColor(new Color(0, 0, 0, 64));
                x1 = proxy.getStart();
                y = initialYStart - 20;
                x = xStart + (int) Math.ceil((double) (x1 - locusStart) / bpPerPixel);

                for (Feature f : peaks) {
                    if (f.overlaps(proxy)) {
                        System.out.println("OVERLAPPPPP!!");
                        g2d.setColor(new Color(255, 0, 128));
                    }
                    System.out.println(proxy.toString());
                    System.out.println(f.toString());
                }

                drawRotatedText(x, y - 10, proxy.getName(), -45);
                g2d.drawLine(x, y, x, annotationstartY + annotationstopY);
            }

            g2d.drawString(locus.getName()
                    + " Width: " + locusWidth
                    + " Pixels: " + locusPixelWidth
                    + " Pos: " + locus.getChromosome().getName() + ":" + locus.getLeftBound() + "-" + locus.getRightBound()
                    + " Proxies: " + locus.getProxies().size()
                    + " Peaks: " + peakTracks[locusId].getNrFeatures()
                    + " Score: " + locus.getScore(), margin, 5);

            // plot chr positions
            for (int bin = 0; bin < locusPixelWidth; bin += 100) {

                int tickX = xStart + bin;
                int chrpos = locusStart + (int) Math.ceil(bin * bpPerPixel);

                g2d.drawString(String.format("%,d", chrpos), tickX, nextYStartPosition - 15);
                g2d.drawLine(tickX, nextYStartPosition, tickX, annotationstartY + annotationstopY);
            }

            //            g2d.drawLine(xStart, nextYStartPosition - 10, xStop, nextYStartPosition - 10);
            // plot genes
            ArrayList<ArrayList<Gene>> geneOverlapMatrix = locus.determineNonOverlappingGenes();
            int geneYLevel = 0;

            if (geneOverlapMatrix.size() > 0) {
                for (ArrayList<Gene> genes : geneOverlapMatrix) {
                    for (Gene g : genes) {
                        // plot the exons..
                        y = nextYStartPosition;
                        ArrayList<Transcript> transcripts = g.getTranscripts();
                        g.getBounds();
                        System.out.println("Plotting gene: "+g.getName()+", chrPos: "+g.getChromosome().getName()+":"+g.getStart()+"-"+g.getStop());
                        int geneStartX = g.getStart();
                        if (geneStartX < locusStart) {
                            geneStartX = locusStart;
                        }

                        int geneStopX = g.getStop();
                        if (geneStopX > locusStop) {
                            geneStopX = locusStop;
                        }

                        int genePixelWidth = (int) Math.ceil((double) (geneStopX - geneStartX) / bpPerPixel);
                        int genePixelStart = xStart + (int) Math.ceil((double) (geneStartX - locusStart) / bpPerPixel);
                        Strand geneStrand = g.getStrand();

                        g2d.setColor(lightgrey);
                        if (geneStrand.equals(Strand.NEG)) {
                            int arrowStart = genePixelStart + genePixelWidth + 10;
                            drawArrow(arrowStart, y + (locusGeneHeight / 2), genePixelStart - 10, y + (locusGeneHeight / 2));
                        } else {
                            drawArrow(genePixelStart - 10, y + (locusGeneHeight / 2), genePixelStart + genePixelWidth + 10, y + (locusGeneHeight / 2));
                        }
                        g2d.setColor(grey);
                        int lastX = xStart;

//                        int y = nextYStartPosition + (geneYLevel * locusGeneHeight) + (geneYLevel * betweenLocusGeneMargin);
                        for (Transcript t : transcripts) {
                            ArrayList<Exon> exons = t.getExons();
                            for (Exon e : exons) {
                                x1 = e.getStart();
                                int x2 = e.getStop();

                                if (x1 < locusStart) {
                                    x1 = locusStart;
                                }

                                if (x2 > locusStop) {
                                    x2 = locusStop;
                                }

                                if (x1 < locusStop) {
                                    // plot
                                    int width = (int) Math.ceil((double) (x2 - x1) / bpPerPixel);

                                    x = xStart + (int) Math.ceil((double) (x1 - locusStart) / bpPerPixel);
                                    g2d.fillRect(x, y, width, locusGeneHeight);
                                    if (x > lastX) {
                                        lastX = x;
                                    }
                                }

                            }
                        }
                        g2d.drawString(g.getName(), lastX + 25, y + locusGeneHeight);
                    }
                    geneYLevel++;
                    nextYStartPosition += locusGeneHeight + betweenLocusGeneMargin;
                }
            } else {
                System.out.println("No genes for locus: " + locus.getName() + "\t" + locus.getGenes().size());
                g2d.drawString("No genes in locus", margin, nextYStartPosition);
                g2d.setColor(new Color(255, 255, 255, 230));
                g2d.fillRect(xStart, nextYStartPosition, locusPixelWidth, locusGeneHeight);
                g2d.setColor(grey);
                nextYStartPosition += locusGeneHeight + betweenLocusGeneMargin;
            }

            // plot the peaks
            // plot the peaks on top of the annotations
            g2d.setColor(new Color(255, 0, 128));
            for (Feature e : peaks) {

                x1 = e.getStart();
                int x2 = e.getStop();

                if (x1 < locusStart) {
                    x1 = locusStart;
                }

                if (x2 > locusStop) {
                    x2 = locusStop;
                }

                if (x1 < locusStop) {
                    // plot
                    int width = (int) Math.ceil((double) (x2 - x1) / bpPerPixel);
                    x = xStart + (int) Math.ceil((double) (x1 - locusStart) / bpPerPixel);
                    g2d.fillRect(x, annotationstartY - 10, width, annotationstopY);
                }

            }

            g2d.setColor(grey);
            // plot annotations
            double maxPos = 0;
            double maxNeg = 0;
            double maxSum = 0;
//            for (int b = 0; b < bedfileTracksForLocus.length; b++) {
//
//                g2d.setColor(new Color(255, 255, 255, 230));
//                g2d.fillRect(xStart, nextYStartPosition, locusPixelWidth, locusAnnotationHeight);
//                g2d.setColor(grey);
//                g2d.drawRect(xStart, nextYStartPosition, locusPixelWidth, locusAnnotationHeight);
//                g2d.setColor(grey);
//                Track t = bedfileTracksForLocus[b];
//                if (locus.getRightBound() - locus.getLeftBound() <= 0) {
//                    System.err.println("Locus has zero size: " + locus.getName() + "\t" + locus.getLeftBound() + "-" + locus.getRightBound());
//                }
//
//                g2d.drawString(t.getName(), xStart + 10, nextYStartPosition + 15);
//
//                double[][] coverage = binReadDepthPerStrand(t, locus.getChromosome(), locus.getLeftBound(), locus.getRightBound(), locusPixelWidth);
//
//                g2d.drawLine(xStart, b, xStart + locusPixelWidth, b);
//
//                for (int bin = 0; bin < locusPixelWidth; bin++) {
//
//                    double nrreadsPos = coverage[0][bin];
//                    double nrreadsNeg = coverage[1][bin];
//                    double sum = nrreadsNeg + nrreadsPos;
//                    if (nrreadsPos > maxPos) {
//                        maxPos = nrreadsPos;
//                    }
//                    if (nrreadsNeg > maxNeg) {
//                        maxNeg = nrreadsNeg;
//                    }
//                    if (sum > maxSum) {
//                        maxSum = sum;
//                    }
//                    System.out.println(bin + "\t" + nrreadsPos + "\t" + nrreadsNeg);
//                    if (nrreadsPos > maxReadDepth) {
//                        nrreadsPos = maxReadDepth;
//                    }
//                    if (nrreadsNeg > maxReadDepth) {
//                        nrreadsNeg = maxReadDepth;
//                    }
//
//                    int pileupPixelsPos = (int) Math.ceil(((double) nrreadsPos / maxReadDepth) * (locusAnnotationHeight / 2));
//                    int pileupPixelsNeg = (int) Math.ceil(((double) nrreadsNeg / maxReadDepth) * (locusAnnotationHeight / 2));
//
////                    System.out.println(bin + "\t" + heightNeg + "\t" + nrNeg);
////                    System.out.println(bin + "\t" + heightPos + "\t" + nrPos);
//                    g2d.setColor(grey);
//                    g2d.fillRect(xStart + (bin * 1), nextYStartPosition + (locusAnnotationHeight / 2) - pileupPixelsPos, 1, pileupPixelsPos);
//
//                    g2d.setColor(lightgrey);
//                    g2d.fillRect(xStart + (bin * 1), nextYStartPosition + (locusAnnotationHeight / 2), 1, pileupPixelsNeg);
//
////                    System.out.println(bin + "\t" + nrPos + "\t" + nrNeg + "\t" + heightPos + "\t" + heightNeg);
////                    int alpha = 128 + (int) Math.ceil(127 * (double) pileupPixels / maxReadDepth);
////                    g2d.setColor(grey);
////                    g2d.fillRect(xStart + (bin * 1), nextYStartPosition + (locusAnnotationHeight) - pileupPixels, 1, pileupPixels);
//                }
//// bins[0]: bins for pos strand
//                // bins[1]: bins for neg strand
//
//                nextYStartPosition += betweenlocusElementMargin + locusAnnotationHeight;
//
//            }

            nextYStartPosition += locusGeneHeight + betweenlocusElementMargin;

            if (bedfileTracksForLocus.length == 1) {
                visualizeReads(g2d, bedfileTracksForLocus[0], xStart, nextYStartPosition, locus.getChromosome(), locusStart, locusStop, locusPixelWidth);
            }
            nextYStartPosition += betweenlocusMargin;
            System.out.println("Max pileup: pos: " + maxPos + " neg: " + maxNeg + " sum: " + maxSum);
        }

        close();
    }

    private void visualizeReads(Graphics2D g2d, Track t, int plotXStart, int plotYStart, Chromosome chr, int windowStart, int windowEnd, int windowPixelSize) {
        Track track = t;
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
            for (int i = 0; i < 2; i++) {
                for (Feature f : features) {

                    Color boxcolor = null;
                    Color linecolor = null;
                    if (f.getStrand() == Strand.NEG) {
                        boxcolor = lightgrey;
                        linecolor = new Color(96, 96, 96);
//                    plot = false;
//                    visitedFeatures.add(f);
                    } else {
                        boxcolor = grey;
                        linecolor = new Color(64, 64, 64);
//                    plot = true;
                    }

                    if (f.getStrand().getNumber() == i && !visitedFeatures.contains(f)) {

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
//                            g2d.setColor(linecolor);
//                            g2d.drawRect(plotXStart + pixelStart, yStart, pixelEnd - pixelStart, featHeight);
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
//                                g2d.setColor(linecolor);
//                                g2d.drawRect(plotXStart + pixelStart, yStart, pixelEnd - pixelStart, featHeight);
//                            System.out.println((plotXStart + pixelStart) + "\t" + yStart + "\t" + (pixelEnd - pixelStart) + "\t" + visitedFeatures.size() + "/" + features.size());
                                visitedFeatures.add(f);
                                otherFeaturesAtY.add(f);

                            }

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

    private final int ARR_SIZE = 4;

    void drawArrow(int x1, int y1, int x2, int y2) {
        Graphics2D g = g2d;

        AffineTransform at1 = g2d.getTransform();

        double dx = x2 - x1, dy = y2 - y1;
        double angle = Math.atan2(dy, dx);
        int len = (int) Math.sqrt(dx * dx + dy * dy);
        AffineTransform at = AffineTransform.getTranslateInstance(x1, y1);
        at.concatenate(AffineTransform.getRotateInstance(angle));
        g.transform(at);

        // Draw horizontal arrow starting in (0, 0)
        g.drawLine(0, 0, len, 0);
        g.fillPolygon(new int[]{len, len - ARR_SIZE, len - ARR_SIZE, len},
                new int[]{0, -ARR_SIZE, ARR_SIZE, 0}, 4);
        g2d.setTransform(at1);
    }

    void drawRotatedText(int x1, int y1, String str, double angle) {
        Graphics2D g = g2d;

        double rad = radian(angle);
        AffineTransform at1 = g2d.getTransform();

        AffineTransform at = AffineTransform.getTranslateInstance(x1, y1);
        at.concatenate(AffineTransform.getRotateInstance(rad));
        g.transform(at);

        // Draw horizontal arrow starting in (0, 0)
        g.drawString(str, 0, 0);

        g2d.setTransform(at1);
    }

    public double radian(double deg) {
        double rad = (deg * (Math.PI / 180));
        return rad;
    }

    private void determineFigureSize() {
        determineMaxWindowWidth(); // gets the max nr of BP over all loci
        System.out.println("Max locus size: " + maxLocusNrPixels);
        System.out.println("Max figure width: " + maxFigureWidth);
        bpPerPixel = (double) maxLocusNrPixels / maxFigureWidth;
        System.out.println("Bp/pixel: " + bpPerPixel);

        int ygenesTotal = 0;
        for (Locus locus : loci) {
            ArrayList<ArrayList<Gene>> geneOverlapMatrix = locus.determineNonOverlappingGenes();
            ygenesTotal += geneOverlapMatrix.size();
        }

        figureWidth = maxFigureWidth + margin * 2;
        figureHeight = (ygenesTotal * locusGeneHeight) + ((ygenesTotal - 1) * betweenlocusElementMargin) // genes
                + (loci.size() * bedFileTracks[loci.size() - 1].length * locusAnnotationHeight) + ((loci.size() - 1) * bedFileTracks[loci.size() - 1].length * betweenlocusElementMargin) // annotations
                + (loci.size() * locusGeneHeight) + betweenlocusElementMargin // peaks
                + (loci.size() * 300)
                + (2 * margin);

        System.out.println("Figure size: " + figureWidth + "x" + figureHeight);
    }

    private void determineMaxWindowWidth() {
        maxLocusNrPixels = 0;

        for (Locus l : loci) {
            int maxStart = l.getLeftBound();
            int maxStop = l.getRightBound();

            int windowSize = maxStop - maxStart;
            if (windowSize > maxLocusNrPixels) {
                maxLocusNrPixels = windowSize;
            }

        }
    }

    private double[] binReadDepth(Track track, Chromosome chr, int start, int stop, int nrBins) {
        double[] bins = new double[nrBins];
        if (stop - start <= 0) {
            return bins;
        }
        int bpPerBin = (stop - start) / nrBins;

        Set<Feature> features = track.getFeatureSet(chr, start, stop);
        System.out.println(features.size());
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

                    bins[binNo]++;

                }
            }
        }
        for (int bin = 0; bin < nrBins; bin++) {
//            System.out.println(bin + "\t" + bins[0][bin] + "\t" + bins[1][bin] + "\t" + (bins[0][bin] / bpPerBin) + "\t" + (bins[1][bin] / bpPerBin));
            bins[bin] /= bpPerBin;

        }
        return bins;
    }

    private double[][] binReadDepthPerStrand(Track track, Chromosome chr, int start, int stop, int nrBins) {
        double[][] bins = new double[2][nrBins];
        if (stop - start <= 0) {
            return bins;
        }
        int bpPerBin = (stop - start) / nrBins;

        Set<Feature> features = track.getFeatureSet(chr, start, stop);
        System.out.println(features.size());
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
//            System.out.println(bin + "\t" + bins[0][bin] + "\t" + bins[1][bin] + "\t" + (bins[0][bin] / bpPerBin) + "\t" + (bins[1][bin] / bpPerBin));
            bins[0][bin] /= bpPerBin;
            bins[1][bin] /= bpPerBin;

        }
        return bins;
    }

}
