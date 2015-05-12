package nl.harmjanwestra.ngs.graphics;

import com.lowagie.text.DocumentException;
import nl.harmjanwestra.ngs.GenotypeFormats.VCFFunctions;
import nl.harmjanwestra.ngs.GenotypeFormats.VCFVariant;
import nl.harmjanwestra.utilities.features.*;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.gtf.GTFAnnotation;
import umcg.genetica.io.text.TextFile;

import java.awt.*;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * Created by hwestra on 4/23/15.
 */
public class VariantPlot extends DefaultGraphics {

    private int margin;

    public VariantPlot(String name, int pageWidth, int pageHeight) throws FileNotFoundException, DocumentException {
        super(name, pageWidth, pageHeight);
        figureWidth = pageWidth;
    }

    public void setMargin(int margin) {
        this.margin = margin;
    }

    private double determineUnit(double range) {

        double divisor = Math.log10(range);
        divisor = Math.floor(divisor);
        divisor = Math.pow(10, divisor);
        return divisor;
    }

    public void plot(String[] variantFiles, String[] variantFileNames, String sequencedRegionFile, Feature region) throws IOException {
        TextFile tf1 = new TextFile(sequencedRegionFile, TextFile.R);

        HashSet<Feature> sequencedRegions = new HashSet<Feature>();

        String[] elems = tf1.readLineElems(TextFile.tab);
        while (elems != null) {
            Feature f = new Feature();
            if (elems.length > 2) {
                f.setChromosome(Chromosome.parseChr(elems[0]));
                f.setStart(Integer.parseInt(elems[1]));
                f.setStop(Integer.parseInt(elems[2]));
                sequencedRegions.add(f);
            }
            elems = tf1.readLineElems(TextFile.tab);
        }
        tf1.close();

        int regionSize = region.getStop() - region.getStart();
        int nrPixels = figureWidth - (2 * margin);

        System.out.println(nrPixels);

        GTFAnnotation annot = new GTFAnnotation("/Data/Annotation/UCSC/genes.gtf");
        TreeSet<Gene> genes = annot.getGeneTree();
        Gene geneStart = new Gene("", region.getChromosome(), Strand.POS, region.getStart(), region.getStart());
        Gene geneStop = new Gene("", region.getChromosome(), Strand.POS, region.getStop(), region.getStop());
        SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);

        g2d.setColor(new Color(64, 64, 64));
        for (Gene g : overlappingGenes) {
            ArrayList<Transcript> transcripts = g.getTranscripts();
            HashSet<Exon> exons = new HashSet<Exon>();
            for (Transcript t : transcripts) {
                exons.addAll(t.getExons());
            }

            for (Exon f : exons) {
                if (region.overlaps(f)) {
                    int start = f.getStart();
                    int featurewidth = f.getStop() - start;
                    int relativeStart = start - region.getStart();
                    if (relativeStart < 0) {
                        featurewidth -= Math.abs(relativeStart);
                        relativeStart = 0;
                    }

                    int relativeStop = relativeStart + featurewidth;
                    if (relativeStop > region.getStop()) {
                        relativeStop = regionSize;
                    }

                    double percStart = (double) relativeStart / regionSize;
                    double percStop = (double) relativeStop / regionSize;

                    int pixelStart = (int) Math.ceil(percStart * nrPixels);
                    int pixelStop = (int) Math.ceil(percStop * nrPixels);

                    System.out.println(f.toString());
                    System.out.println(pixelStart + "\t" + pixelStop);
                    int y1 = margin - 20;
                    int y2 = margin + 10;
                    g2d.fillRect(margin + pixelStart, y1, pixelStop - pixelStart, 10);
                }
            }

            int start = g.getStart();
            int featurewidth = g.getStop() - start;
            int relativeStart = start - region.getStart();
            if (relativeStart < 0) {
                featurewidth -= Math.abs(relativeStart);
                relativeStart = 0;
            }

            int relativeStop = relativeStart + featurewidth;
            if (relativeStop > region.getStop()) {
                relativeStop = regionSize;
            }

            double percStart = (double) relativeStart / regionSize;
            double percStop = (double) relativeStop / regionSize;

            int pixelStart = (int) Math.ceil(percStart * nrPixels);
            int pixelStop = (int) Math.ceil(percStop * nrPixels);

            System.out.println(g.toString());
            System.out.println(pixelStart + "\t" + pixelStop);
            int y1 = margin - 20;

            g2d.drawLine(margin + pixelStart, y1 + 5, margin + pixelStop, y1 + 5);
            g2d.drawString(g.getGeneId(), margin + pixelStop + 5, y1 + 5);

        }

        g2d.setColor(new Color(255, 1, 0, 202));
        for (Feature f : sequencedRegions) {
            if (f.overlaps(region)) {
                int start = f.getStart();
                int featurewidth = f.getStop() - start;
                int relativeStart = start - region.getStart();
                if (relativeStart < 0) {
                    featurewidth -= Math.abs(relativeStart);
                    relativeStart = 0;
                }

                int relativeStop = relativeStart + featurewidth;
                if (relativeStop > region.getStop()) {
                    relativeStop = regionSize;
                }

                double percStart = (double) relativeStart / regionSize;
                double percStop = (double) relativeStop / regionSize;

                int pixelStart = (int) Math.ceil(percStart * nrPixels);
                int pixelStop = (int) Math.ceil(percStop * nrPixels);

                System.out.println(f.toString());
                System.out.println(pixelStart + "\t" + pixelStop);
                int y1 = margin;
                int y2 = margin + 10;
                g2d.fillRect(margin + pixelStart, y1, pixelStop - pixelStart, 10);
            }
        }

        int betweenmargin = 25;
        int starty = margin + betweenmargin + 10;
        int dotsize = 50;
        int minDotSize = 2;

        VCFFunctions func = new VCFFunctions();
        HashMap<Feature, Integer> ctr = new HashMap<Feature, Integer>();
        for (int i = 0; i < variantFiles.length; i++) {
            String variantFile = variantFiles[i];
            TextFile tf = new TextFile(variantFile, TextFile.R);
            String ln = tf.readLine();
            while (ln != null) {

                if (!ln.startsWith("#")) {
                    VCFVariant v = new VCFVariant(ln, 10, 30, true);

                    Feature f = new Feature();

                    f.setChromosome(Chromosome.parseChr(v.getChr()));
                    f.setStart(v.getPos());
                    f.setStop(v.getPos());
                    Integer ct = ctr.get(f);
                    if (ct == null) {
                        ct = 0;
                    }
                    ct++;
                    ctr.put(f, ct);
                }
                ln = tf.readLine();
            }
            tf.close();

        }
        System.out.println(ctr.size() + " variants loaded");

        for (int i = 0; i < variantFiles.length; i++) {
            starty += betweenmargin + dotsize;
            String variantFile = variantFiles[i];
            String variantFileName = variantFileNames[i];

            TextFile tf = new TextFile(variantFile, TextFile.R);

            String ln = tf.readLine();

            Color grey = new Color(0, 0, 0, 28);
            g2d.setColor(grey);

            while (ln != null) {

                if (!ln.startsWith("#")) {
                    VCFVariant v = new VCFVariant(ln, 10, 30,true);

                    Feature f = new Feature();

                    f.setChromosome(Chromosome.parseChr(v.getChr()));
                    f.setStart(v.getPos());
                    f.setStop(v.getPos());

                    if (region.overlaps(f)) {
                        double maf = v.getMAF();
                        if (maf > 0 && maf < 1) {
                            int size = (int) Math.ceil((maf * 2) * dotsize);
                            System.out.println(maf + "\t" + size);
                            int relativeStart = f.getStart() - region.getStart();
                            if (relativeStart < 0) {
                                relativeStart = 0;
                            }

                            double percStart = (double) relativeStart / regionSize;
                            int pixelStart = (int) Math.ceil(percStart * nrPixels);

                            if (size < minDotSize) {
                                size = minDotSize;
                            }
                            int halfdsize = size / 2;

                            if (ctr.get(f) == null || ctr.get(f) == 1) {
                                g2d.setColor(new Color(255, 0, 0, 128));
                            } else {
                                g2d.setColor(grey);
                            }

                            g2d.fillOval(margin + pixelStart - halfdsize, starty - halfdsize, size, size);
                        }
                    }
                }
                ln = tf.readLine();
            }
            tf.close();
        }

        // plot coordinates
        starty += betweenmargin + 10;

        g2d.setColor(Color.BLACK);

        g2d.drawLine(margin, starty, margin + nrPixels, starty);

        int unit = (int) Math.ceil(determineUnit(regionSize));

        for (int i = region.getStart(); i < region.getStop(); i++) {
            if (i % unit == 0) {
                int relativeStart = i - region.getStart();
                double percStart = (double) relativeStart / regionSize;
                int pixelStart = (int) Math.ceil(percStart * nrPixels);
                g2d.drawLine(margin + pixelStart, starty - 5, margin + pixelStart, starty + 5);
                g2d.drawString("" + i, margin + pixelStart, starty + 10);
            }
        }

        starty += betweenmargin + 10;

        int startx = margin;
        int between = 5;
        int c = 0;
        int nextX = startx;
        for (double maf = 0; maf < 1; maf += 0.1) {
            int size = (int) Math.ceil((maf) * dotsize);

            g2d.setColor(new Color(255, 0, 0, 128));
            int halfdsize = size / 2;
            g2d.fillOval(margin + nextX - halfdsize, starty - halfdsize, size, size);

            nextX += (between + size);

            c++;
        }

        close();

    }

    public void plot2(String[] variantFiles, String[] variantFileNames, String sequencedRegionFile, String regionFile) throws IOException {
        TextFile tf1 = new TextFile(sequencedRegionFile, TextFile.R);

        HashSet<Feature> sequencedRegions = new HashSet<Feature>();

        String[] elems = tf1.readLineElems(TextFile.tab);
        while (elems != null) {
            Feature f = new Feature();
            if (elems.length > 2) {
                f.setChromosome(Chromosome.parseChr(elems[0]));
                f.setStart(Integer.parseInt(elems[1]));
                f.setStop(Integer.parseInt(elems[2]));
                sequencedRegions.add(f);
            }
            elems = tf1.readLineElems(TextFile.tab);
        }
        tf1.close();

        TextFile tf2 = new TextFile(regionFile, TextFile.R);
        elems = tf2.readLineElems(TextFile.tab);
        Feature impregion = new Feature();
        impregion.setChromosome(Chromosome.parseChr(elems[0]));
        impregion.setStart(Integer.parseInt(elems[1]));
        impregion.setStop(Integer.parseInt(elems[2]));
        tf2.close();

        ArrayList<Feature> overlappingRegions = new ArrayList<Feature>();

        for (Feature f : sequencedRegions) {
            if (f.overlaps(impregion)) {
                overlappingRegions.add(f);
            }
        }

        int nrPixels = figureWidth - (2 * margin);

        System.out.println(nrPixels);
        System.out.println(overlappingRegions.size() + " overlapping regions");
        int pixelsBetween = 5;
        int pixelsPerRegion = (nrPixels - (overlappingRegions.size() - 1) * pixelsBetween) / overlappingRegions.size();
        System.out.println(pixelsPerRegion + " pixels per region");

        int betweenmargin = 25;
        int origStarty = margin + betweenmargin + 10;
        int dotsize = 24;
        int minDotSize = 2;
        int origstartx = margin;

        HashMap<Feature, Integer> ctr = new HashMap<Feature, Integer>();
        for (int i = 0; i < variantFiles.length; i++) {
            String variantFile = variantFiles[i];
            TextFile tf = new TextFile(variantFile, TextFile.R);
            String ln = tf.readLine();
            while (ln != null) {

                if (!ln.startsWith("#")) {
                    VCFVariant v = new VCFVariant(ln, 10, 30, true);

                    Feature f = new Feature();

                    f.setChromosome(Chromosome.parseChr(v.getChr()));
                    f.setStart(v.getPos());
                    f.setStop(v.getPos());
                    Integer ct = ctr.get(f);
                    if (ct == null) {
                        ct = 0;
                    }
                    ct++;
                    ctr.put(f, ct);
                }
                ln = tf.readLine();
            }
            tf.close();

        }
        System.out.println(ctr.size() + " variants loaded");

        int fctr = 0;
        for (Feature region : overlappingRegions) {
            int starty = origStarty;

            int startx = origstartx + (fctr * pixelsBetween) + (fctr * pixelsPerRegion);

            int regionSize = region.getStop() - region.getStart();

            for (int i = 0; i < variantFiles.length; i++) {
                starty += betweenmargin + dotsize;
                String variantFile = variantFiles[i];
                String variantFileName = variantFileNames[i];

                TextFile tf = new TextFile(variantFile, TextFile.R);

                String ln = tf.readLine();

                Color grey = new Color(0, 0, 0, 64);
                g2d.setColor(grey);

                while (ln != null) {

                    if (!ln.startsWith("#")) {
                        VCFVariant v = new VCFVariant(ln, 10, 30, true);

                        Feature f = new Feature();

                        f.setChromosome(Chromosome.parseChr(v.getChr()));
                        f.setStart(v.getPos());
                        f.setStop(v.getPos());

                        if (region.overlaps(f)) {
                            double maf = v.getMAF();
                            if (maf > 0 && maf < 1) {
                                int size = (int) Math.ceil((maf * 2) * 10);
                                int relativeStart = f.getStart() - region.getStart();
                                if (relativeStart < 0) {
                                    relativeStart = 0;
                                }

                                double percStart = (double) relativeStart / regionSize;
                                int pixelStart = (int) Math.ceil(percStart * nrPixels);

                                int halfdotsize = dotsize / 2;
                                if (size < minDotSize) {
                                    size = minDotSize;
                                }
                                int halfdsize = size / 2;

                                if (ctr.get(f) == null || ctr.get(f) == 1) {
                                    g2d.setColor(new Color(255, 0, 0, 128));
                                } else {
                                    g2d.setColor(grey);
                                }
                                g2d.fillOval(startx + pixelStart - halfdsize, starty - halfdsize, size, size);
                            }
                        }
                    }
                    ln = tf.readLine();
                }
                tf.close();
            }
        }

        close();
    }
}
