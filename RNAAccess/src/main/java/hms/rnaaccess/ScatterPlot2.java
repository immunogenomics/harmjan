/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.rnaaccess;

import com.lowagie.text.Document;
import com.lowagie.text.DocumentException;
import com.lowagie.text.Rectangle;
import com.lowagie.text.pdf.PdfContentByte;
import com.lowagie.text.pdf.PdfWriter;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.text.DecimalFormat;
import umcg.genetica.util.Primitives;

/**
 *
 * @author Harm-Jan
 */
public class ScatterPlot2 {

    int plotWidth;
    int plotHeight;
    int plotMargin;
    int margin;
    int figureWidth;
    int figureHeight;
    double[][] x;
    double[][] y;
    private String[] plotTitles;
    private String figureTitle;
    private int nrColumns;
    private int nrRows;
    private ScatterPlot.OUTPUTFORMAT format;
    private Document document;
    private String outfilename;
    private PdfWriter writer;
    private PdfContentByte cb;
    private BufferedImage bi;
    private Graphics2D g2d;
    private Color color;
    private int fontheight;
    private double maxX;
    private double maxY;
    private double minX;
    private double minY;
    private boolean interceptAtZero;
    private double rangeX;
    private double rangeY;
    private double unitX = Double.NaN;
    private double unitY = Double.NaN;
    private final String xAxisTitle;
    private final String yAxisTitle;

    public ScatterPlot2(double[][] x, double[][] y, boolean interceptAtZero, String[] plotTitles, String xAxisTitle, String yAxisTitle,
            String figureTitle, int plotWidth, int plotHeight, int plotMargin, int margin, int nrColumns, String outputFileName) {
        this.x = x;
        this.y = y;
        this.plotHeight = plotHeight;
        this.plotWidth = plotWidth;
        this.plotMargin = plotMargin;
        this.margin = margin;
        this.plotTitles = plotTitles;
        this.figureTitle = figureTitle;
        this.nrColumns = nrColumns;
        this.outfilename = outputFileName;
        this.interceptAtZero = interceptAtZero;
        this.xAxisTitle = xAxisTitle;
        this.yAxisTitle = yAxisTitle;
        if (outputFileName.toLowerCase().endsWith(".jpg")) {
            this.format = ScatterPlot.OUTPUTFORMAT.JPG;
        } else if (outputFileName.toLowerCase().endsWith(".png")) {
            this.format = ScatterPlot.OUTPUTFORMAT.PNG;
        } else if (outputFileName.toLowerCase().endsWith(".pdf")) {
            this.format = ScatterPlot.OUTPUTFORMAT.PDF;
        } else {
            throw new IllegalArgumentException("File format not supported");
        }

        init();
        plot();
        write();
    }

    private void init() {

        if (nrColumns == 1) {
            figureWidth = (2 * margin) + plotWidth;
            figureHeight = (2 * margin) + (x.length * plotHeight) + ((x.length - 1) * plotMargin);
            nrRows = x.length;
        } else {
            int nrPlots = x.length;
            figureWidth = (2 * margin) + (plotWidth * nrColumns) + ((nrColumns - 1) * plotMargin);
            nrRows = (int) Math.ceil((double) nrPlots / nrColumns);
            figureHeight = (2 * margin) + (nrRows * plotHeight) + ((nrRows - 1) * plotMargin);
        }

        System.out.println("Scatterplot will be " + nrColumns + " columns and " + nrRows + " rows");

        if (format == ScatterPlot.OUTPUTFORMAT.PDF) {
            Rectangle rectangle = new Rectangle(figureWidth, figureHeight);
            document = new com.lowagie.text.Document(rectangle);

            if (!outfilename.toLowerCase().endsWith(".pdf")) {
                outfilename += ".pdf";
            }
            try {
                writer = com.lowagie.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(outfilename));

            } catch (DocumentException e) {
                e.printStackTrace();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
            document.open();
            cb = writer.getDirectContent();
            cb.saveState();

//            com.lowagie.text.pdf.DefaultFontMapper fontMap = new com.lowagie.text.pdf.DefaultFontMapper();
            g2d = cb.createGraphics(figureWidth, figureHeight);
        } else {
            bi = new java.awt.image.BufferedImage(figureWidth, figureHeight, java.awt.image.BufferedImage.TYPE_INT_RGB);
            g2d = bi.createGraphics();
        }

        color = new Color(255, 255, 255);

        g2d.setColor(color);
        g2d.setFont(new Font("Verdana", Font.PLAIN, 10));
        FontMetrics fontmetrics = g2d.getFontMetrics();
        fontheight = fontmetrics.getHeight();
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

        g2d.fillRect(0, 0, figureWidth, figureHeight);

        // draw axis
        color = new Color(0, 0, 0);
        g2d.setColor(color);
//        g2d.drawRect(margin, margin, graphHeight - (2 * margin), graphWidth - (2 * margin));

        g2d.drawString(figureTitle, margin, margin/2);
        // 
        determineRange();

        // draw axes
        drawAxis();

    }

    private void determineRange() {
        maxX = Primitives.max(x);
        maxY = Primitives.max(y);
        minX = Primitives.min(x);
        minY = Primitives.min(y);

        if (interceptAtZero) {
            if (minY > 0) {
                minY = 0;
            }
            if (minX > 0) {
                minX = 0;
            }

            if (maxY < 0) {
                maxY = 0;
            }
            if (maxX < 0) {
                maxX = 0;
            }
        }

        if (Double.isInfinite(maxX)) {
            maxX = Double.MAX_VALUE;
        }

        if (Double.isInfinite(minX)) {
            minX = -Double.MAX_VALUE;
        }

        if (Double.isInfinite(maxY)) {
            maxY = Double.MAX_VALUE;
        }

        if (Double.isInfinite(minY)) {
            minY = -Double.MAX_VALUE;
        }

        rangeX = Math.abs(minX - maxX);
        rangeY = Math.abs(minY - maxY);
//

//        System.out.println("MinX: " + minX + "\nMaxX: " + maxX + "\nMinY: " + minY + "\nMaxY: " + maxY + "\nRangeX: " + rangeX + "\nRangeY: " + rangeY + "\nUnitX: " + unitX + "\nUnitY: " + unitY);
        if (Double.isNaN(unitX)) {

            unitX = determineUnit(rangeX);
            if (unitX >= Math.abs(minX) && unitX >= Math.abs(maxX)) {
                unitX /= 10;
            }
            System.out.println("Determined unitX: " + unitX);
        }
        if (Double.isNaN(unitY)) {
            unitY = determineUnit(rangeY);
            // prevent excessive tickmarking..
            if (unitY >= Math.abs(minY) && unitY >= Math.abs(maxY)) {
                unitY /= 10;
            }
            System.out.println("Determined unitY: " + unitY);
        }

        // round off the limits towards the unit
        double remainder = Math.abs(maxX) % unitX;
        if (remainder > 0d && maxX != 0) {
            double diff = unitX - remainder;
            maxX += diff;
        }

        remainder = Math.abs(minX) % unitX;
//        System.out.println(unitX + " min: " + minX + " rem: " + remainder);
        if (remainder > 0d && minX != 0) {
            if (minX < 0) {
                minX -= (unitX - remainder);
            } else {
                minX -= remainder;
            }

        }

        rangeX = Math.abs(minX - maxX);
        if (rangeX == 0) {
            maxX = unitX;
        }

        // ensure the max Y is rounded off by to the next unitY
        remainder = Math.abs(maxY) % unitY;
        if (remainder > 0d && maxY != 0d) {
            double diff = unitY - remainder;
            maxY += diff;
        }

        remainder = Math.abs(minY) % unitY; // -8 , unit == 10, remainder = 8 // diff = 2
//        System.out.println(minY + "\t" + unitY + "\t" + remainder);
        if (remainder > 0d && minY != 0) {
            if (minY < 0) {
                minY -= (unitY - remainder);
            } else {
                minY -= remainder;
            }
        }
//        System.out.println(minY + "\t" + unitY + "\t" + remainder);
        rangeY = Math.abs(minY - maxY);
        if (rangeY == 0) {
            maxY = unitY;
        }

        if (Double.isInfinite(maxY)) {
            maxY = Double.MAX_VALUE;
        }

        if (Double.isInfinite(minY)) {
            minY = -Double.MAX_VALUE;
        }

        if (Double.isInfinite(maxX)) {
            maxX = Double.MAX_VALUE;
        }

        if (Double.isInfinite(minX)) {
            minX = -Double.MAX_VALUE;
        }

        rangeX = Math.abs(minX - maxX);
        rangeY = Math.abs(minY - maxY);

        System.out.println("MinX: " + minX + "\nMaxX: " + maxX + "\nMinY: " + minY + "\nMaxY: " + maxY + "\nRangeX: " + rangeX + "\nRangeY: " + rangeY + "\nUnitX: " + unitX + "\nUnitY: " + unitY);
    }

    private void drawAxis() {
        Color originalColor = g2d.getColor();
        Font oriFont = g2d.getFont();

        g2d.setFont(new Font("SansSerif", Font.PLAIN, 9));
        g2d.setColor(Color.gray);

        FontMetrics fontmetrics = g2d.getFontMetrics();
        int tickFontHeight = fontmetrics.getHeight();

        AffineTransform cc = new AffineTransform();
        cc.setToRotation( -(Math.PI / 2.0));
        AffineTransform norm = new AffineTransform();

        // x-axis
        for (int sample = 0; sample < x.length; sample++) {
            int r = sample / nrColumns;
            int c = sample % nrColumns;
            int columnXStart = margin + (c * plotWidth) + (c * plotMargin);
            int rowYStart = margin + (r * plotHeight) + (r * plotMargin);

            g2d.drawString(plotTitles[sample], columnXStart, rowYStart - 20);

            // draw box ?
            g2d.drawRect(columnXStart - 10, rowYStart - 10, plotWidth + 20, plotHeight + 20);

            g2d.drawString(xAxisTitle, columnXStart, rowYStart + plotHeight + 20);
            g2d.setTransform(cc);

            g2d.drawString(yAxisTitle, columnXStart - 40, rowYStart + plotHeight);

            g2d.setTransform(norm);

            int xposYAxis = 0;
            int yposXAxis = 0;
            if (minY >= 0) {
                // all Y values above 0, X axis starts at bottom left: |_
                g2d.drawLine(columnXStart, rowYStart + plotHeight, columnXStart + plotWidth, rowYStart + plotHeight);
                yposXAxis = rowYStart + plotHeight;
            } else if (maxY <= 0) {
                // all Y values below 0, X axis starts at top left 
                g2d.drawLine(columnXStart, rowYStart, columnXStart + plotWidth, rowYStart);
                yposXAxis = rowYStart;
            } else {
                // X-axis crosses the y-axis somewhere, at DIFF
                // this method does show some rounding errors at the moment..
                int diff = (int) Math.ceil((Math.abs(minY) / rangeY) * plotHeight);
                yposXAxis = rowYStart + plotHeight - diff;
                g2d.drawLine(columnXStart, yposXAxis, columnXStart + plotWidth, yposXAxis);

            }

            // y-axis
            if (minX > 0) {
                // all X values above 0, Y-axis starts at top left
                g2d.drawLine(columnXStart, rowYStart, columnXStart, rowYStart + plotHeight);
                xposYAxis = columnXStart;
            } else if (maxX <= 0) {
                // all X values below 0, axis starts at top right 
                g2d.drawLine(columnXStart + plotWidth, rowYStart, columnXStart + plotWidth, rowYStart + plotHeight);
                xposYAxis = columnXStart + plotWidth;
            } else {
                // Y-axis crosses the X-axis somewhere, at DIFF
                // this method does show some rounding errors at the moment..
                int diff = (int) Math.ceil((Math.abs(minX) / rangeX) * plotWidth);
                xposYAxis = columnXStart + diff;
                g2d.drawLine(xposYAxis, rowYStart, xposYAxis, rowYStart + plotHeight);

            }

            // X-ticks
            if (!Double.isNaN(unitX)) {
                // start drawing from minX, add one unitX at a time..
                // first round minX using the unitX
                double tickX = minX - (minX % unitX);

                while (tickX <= maxX) {

                    double nonroundedPerc = Math.abs(minX - tickX) / rangeX;
//                double roundedPerc = roundToDecimals((nonroundedPerc), 2);
                    int diff = (int) Math.floor(nonroundedPerc * plotWidth); // for the position relative to min and maxX

                    String tickLabelFormatted = null;
                    if (unitX > 10000 || unitX < 0.001) {
                        tickLabelFormatted = new DecimalFormat("0.#E0").format(tickX);
                    } else {
                        tickLabelFormatted = new DecimalFormat("#.###").format(tickX);
                    }

                    if (tickX == 0) {
                        g2d.drawLine(xposYAxis, yposXAxis - 3, xposYAxis, yposXAxis + 3);
                        int nrPixelsString = g2d.getFontMetrics().stringWidth(tickLabelFormatted) / 2;
                        g2d.drawString(tickLabelFormatted, columnXStart + diff - nrPixelsString, yposXAxis + tickFontHeight + 3);

                    } else {
                        g2d.drawLine(columnXStart + diff, yposXAxis - 3, columnXStart + diff, yposXAxis + 3);
                        int nrPixelsString = g2d.getFontMetrics().stringWidth(tickLabelFormatted) / 2;
                        g2d.drawString(tickLabelFormatted, columnXStart + diff - nrPixelsString, yposXAxis + tickFontHeight + 3);

                    }

                    tickX += unitX;

                }
            }

            // Y-ticks
            if (!Double.isNaN(unitY)) {
                double tickY = minY - (minY % unitY);
                while (tickY <= maxY) {

                    double nonroundedPerc = Math.abs(minY - tickY) / rangeY;

//                double perc = roundToDecimals((Math.abs(minY - tickY) / rangeY), 2);
                    int diff = (int) Math.floor(nonroundedPerc * plotHeight);
                    String tickLabelFormatted = null;
                    if (unitY > 100000 || unitY < 0.0001) {
                        tickLabelFormatted = new DecimalFormat("0.#E0").format(tickY);
                    } else {
                        tickLabelFormatted = new DecimalFormat("#.###").format(tickY);
                    }

                    if (tickY == 0) {
                        g2d.drawLine(xposYAxis - 3, yposXAxis, xposYAxis + 3, yposXAxis);
                        g2d.drawString(tickLabelFormatted, xposYAxis - g2d.getFontMetrics().stringWidth(tickLabelFormatted) - 5, rowYStart + plotHeight - diff + (tickFontHeight / 2) - 3);
                    } else {
                        g2d.drawLine(xposYAxis - 3, rowYStart + plotHeight - diff, xposYAxis + 3, rowYStart + plotHeight - diff);
                        g2d.drawString(tickLabelFormatted, xposYAxis - g2d.getFontMetrics().stringWidth(tickLabelFormatted) - 5, rowYStart + plotHeight - diff + (tickFontHeight / 2) - 3);
                    }

                    tickY += unitY;
                }
            }
        }

        g2d.setFont(oriFont);
        g2d.setColor(originalColor);

    }

    private void plot() {

        for (int sample = 0; sample < x.length; sample++) {
            double[] xtmp = x[sample];
            double[] ytmp = y[sample];

            int r = sample / nrColumns;
            int c = sample % nrColumns;

            int columnXStart = margin + (c * plotWidth) + (c * plotMargin);
            int rowYStart = margin + (r * plotHeight) + (r * plotMargin);
            g2d.setColor(new Color(0, 0, 0));
            for (int i = 0; i < xtmp.length; i++) {

//                if (category != null && colors != null) {
//                    g2d.setColor(colors[category[i]]);
//                } else {
//                    g2d.setColor(new Color(0, 128, 255, 64));
//                }
                double xval = xtmp[i];
                double yval = ytmp[i];

                if (Double.isInfinite(xval)) {
                    if (xval < 0) {
                        xval = -Double.MAX_VALUE;
                    } else {
                        xval = -Double.MAX_VALUE;
                    }
                }

                if (Double.isInfinite(yval)) {
                    if (yval < 0) {
                        yval = -Double.MAX_VALUE;
                    } else {
                        yval = -Double.MAX_VALUE;
                    }
                }

                int posX = columnXStart + (int) Math.ceil((Math.abs(minX - xval) / rangeX) * plotWidth);
                int posY = (rowYStart + plotHeight) - (int) Math.ceil((Math.abs(minY - yval) / rangeY) * plotHeight);

                g2d.fillOval(posX - 2, posY - 2, 4, 4);
            }

            c++;
        }

    }

    private void write() {
        try {
            g2d.dispose();
            if (format == ScatterPlot.OUTPUTFORMAT.JPG) {
                if (!outfilename.toLowerCase().endsWith(".jpg")) {
                    outfilename += ".jpg";
                }
                javax.imageio.ImageIO.write(bi, "jpg", new File(outfilename));
            } else if (format == ScatterPlot.OUTPUTFORMAT.PNG) {
                if (!outfilename.toLowerCase().endsWith(".png")) {
                    outfilename += ".png";
                }
                javax.imageio.ImageIO.write(bi, "png", new File(outfilename));
            } else {
                cb.restoreState();
                document.close();
                writer.close();
            }
        } catch (Exception e) {
            System.out.println(e.getMessage());
            System.out.println(e.getStackTrace());
        }

    }

    private double determineUnit(double range) {

        double divisor = Math.log10(range);
        divisor = Math.floor(divisor);
        divisor = Math.pow(10, divisor);
        return divisor;
    }
}
