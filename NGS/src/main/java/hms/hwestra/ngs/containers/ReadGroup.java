/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.ngs.containers;

import hms.hwestra.utilities.features.Chromosome;
import java.io.IOException;
import umcg.genetica.graphics.ViolinBoxPlot;

/**
 *
 * @author hwestra
 */
public class ReadGroup {

    public String name;
    public int[][] readsPerChr = new int[2][26];
    public int[][] readsPairedPerChr = new int[2][26];

    public int[][] mapqPerChrPosStrand = new int[26][61]; // [chr][mapq]
    public int[][] mapqPerChrNegStrand = new int[26][61]; // [chr][mapq]
    public int[][] dupsPerChr = new int[2][26];
    public int[][] baseQualPerPosPosStrand = new int[50][61]; // [pos][qual]
    public int[][] baseQualPerPosNegStrand = new int[50][61]; // [pos][qual]
    int nrReads = 0;

    public void plotMapQPerChr(String outputdir) throws IOException {
        ViolinBoxPlot plot = new ViolinBoxPlot(false);

        // vals[dataset][category][value]
        double[][][] data = new double[2][26][];
        String[] datasetLabel = new String[2];
        datasetLabel[0] = name+" - pos strand";
        datasetLabel[1] = name+" - neg strand";
        String[][] categoryLabel = new String[2][26];
        int sample = 100000;
        Chromosome[] allChr = Chromosome.values();

        for (Chromosome chr : allChr) {
            int c = chr.getNumber();

            long sum = 0;

            for (int q = 0; q < mapqPerChrPosStrand[chr.getNumber()].length; q++) {
                sum += mapqPerChrPosStrand[c][q];
            }

            double[] d = new double[sample];
            int y = 0;
            for (int q = 0; q < mapqPerChrPosStrand[c].length; q++) {
                int fraction = (int) Math.ceil(sample * ((double) mapqPerChrPosStrand[c][q] / sum));
//                System.out.println(c + "\t" + fraction + "\t" + sample + "\t" + mapqPerChr[c][q] + "\t" + sum);

                for (int z = 0; z < fraction; z++) {
                    if (y < sample) {
                        d[y] = q;
                        y++;
                    }
                }
            }
            categoryLabel[0][c] = chr.getName() + " - " + sum;

            data[0][c] = d;
        }
        
        for (Chromosome chr : allChr) {
            int c = chr.getNumber();

            long sum = 0;

            for (int q = 0; q < mapqPerChrNegStrand[chr.getNumber()].length; q++) {
                sum += mapqPerChrNegStrand[c][q];
            }

            double[] d = new double[sample];
            int y = 0;
            for (int q = 0; q < mapqPerChrNegStrand[c].length; q++) {
                int fraction = (int) Math.ceil(sample * ((double) mapqPerChrNegStrand[c][q] / sum));
//                System.out.println(c + "\t" + fraction + "\t" + sample + "\t" + mapqPerChr[c][q] + "\t" + sum);

                for (int z = 0; z < fraction; z++) {
                    if (y < sample) {
                        d[y] = q;
                        y++;
                    }
                }
            }
            categoryLabel[1][c] = chr.getName() + " - " + sum;

            data[1][c] = d;
        }
        System.out.println("output will go to: " + outputdir + "-MapQPerChr.pdf");
        plot.draw(data, datasetLabel, categoryLabel, "Ct", ViolinBoxPlot.Output.PDF, outputdir + "-MapQPerChr.pdf");

    }

    public void plotQualPerPos(String outputdir) throws IOException {
        ViolinBoxPlot plot = new ViolinBoxPlot(false);

        // vals[dataset][category][value]
        double[][][] data = new double[2][50][];
        String[] datasetLabel = new String[2];
        datasetLabel[0] = name+"pos strand";
        datasetLabel[1] = name+"neg strand";
        String[][] categoryLabel = new String[2][50];
        int sample = 10000;

        for (int pos = 0; pos < baseQualPerPosPosStrand.length; pos++) {
            long sum = 0; // sum of qualities for this position
            for (int qual = 0; qual < baseQualPerPosPosStrand[pos].length; qual++) {
                sum += baseQualPerPosPosStrand[pos][qual];
            }

            double[] d = new double[sample];
            int y = 0;
            for (int qual = 0; qual < baseQualPerPosPosStrand[pos].length; qual++) {
                int fraction = (int) Math.ceil(sample * ((double) baseQualPerPosPosStrand[pos][qual] / sum));
                System.out.println(pos + "\t" + fraction + "\t" + sample + "\t" + baseQualPerPosPosStrand[pos][qual] + "\t" + sum);

                for (int z = 0; z < fraction; z++) {
                    if (y < sample) {
                        d[y] = qual;
                        y++;
                    }
                }
            }
            categoryLabel[0][pos] = (pos + 1) + " - " + sum;

            data[0][pos] = d;

        }
        
        for (int pos = 0; pos < baseQualPerPosNegStrand.length; pos++) {
            long sum = 0; // sum of qualities for this position
            for (int qual = 0; qual < baseQualPerPosNegStrand[pos].length; qual++) {
                sum += baseQualPerPosNegStrand[pos][qual];
            }

            double[] d = new double[sample];
            int y = 0;
            for (int qual = 0; qual < baseQualPerPosNegStrand[pos].length; qual++) {
                int fraction = (int) Math.ceil(sample * ((double) baseQualPerPosNegStrand[pos][qual] / sum));
                System.out.println(pos + "\t" + fraction + "\t" + sample + "\t" + baseQualPerPosNegStrand[pos][qual] + "\t" + sum);

                for (int z = 0; z < fraction; z++) {
                    if (y < sample) {
                        d[y] = qual;
                        y++;
                    }
                }
            }
            categoryLabel[1][pos] = (pos + 1) + " - " + sum;

            data[1][pos] = d;

        }

        
        System.out.println("output will go to: " + outputdir + "-QualPerPos.pdf");
        plot.draw(data, datasetLabel, categoryLabel, "Ct", ViolinBoxPlot.Output.PDF, outputdir + "-QualPerPos.pdf");

    }
    
    

}
