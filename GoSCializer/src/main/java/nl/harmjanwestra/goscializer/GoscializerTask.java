/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.goscializer;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.SortableSNP;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author hwestra
 */
public class GoscializerTask implements Runnable {

    String snpfile;
    String tabixdir;
    double r2threshold;
    int ldWindowSize;
    int windowmargin;
    String[] listOfAnnotationFiles;
    String outputDir;
    int taskctr;
    private final boolean add1;

    public GoscializerTask(int taskctr, String snpfile, String tabixdir, double r2threshold,
            int ldWindowSize, int windowmargin, String[] listOfAnnotationFiles, String outputDir, boolean add1) {
        this.add1 = add1;
        this.taskctr = taskctr;
        this.snpfile = snpfile;
        this.tabixdir = tabixdir;
        this.r2threshold = r2threshold;
        this.ldWindowSize = ldWindowSize;
        this.windowmargin = windowmargin;
        this.listOfAnnotationFiles = listOfAnnotationFiles;
        this.outputDir = outputDir;
    }

    @Override
    public void run() {
        // load variants
        try {

            ArrayList<SortableSNP> snps = loadSNPs(snpfile);

            Collections.sort(snps);
//           .println(taskctr + "- Loaded SNPs: ");
            HashSet<Byte> chromosomes = new HashSet<Byte>();
            for (SortableSNP snp : snps) {

                chromosomes.add(snp.chr);
            }

            System.out.println(taskctr + "- Total: " + snps.size());

            // list chromosomes
            int total = 0;
            int locusCounter = 0;
            TextFile snpBedOutput = new TextFile(outputDir + "SNPs.bed", TextFile.W);
            TextFile regionBedOutput = new TextFile(outputDir + "Regions.bed", TextFile.W);
            ArrayList<Feature> windows = new ArrayList<Feature>();
            tabixdir = Gpio.formatAsDirectory(tabixdir);

            for (Byte chr : chromosomes) {
                if (chr > 0) {
//                    System.out.println(taskctr + "- Finding windows for chromosome: " + chr);
                    ArrayList<SortableSNP> snpsForChr = new ArrayList<SortableSNP>();
                    for (SortableSNP snp : snps) {
                        if (chr.byteValue() == snp.chr) {
                            snpsForChr.add(snp);
                        }
                    }

                    Chromosome chrE = Chromosome.parseChr("" + chr);

//                    System.out.println(taskctr + "- " + snpsForChr.size() + " SNPs found.");
                    total += snpsForChr.size();

                    if (snpsForChr.size() > 0) {

                        HashMap<SortableSNP, ArrayList<SortableSNP>> proxies = getProxiesForSNPGosiaFiles(tabixdir, chr, ldWindowSize, r2threshold, snpsForChr);// getProxiesForSNPs(tabixdir, chr, ldWindowSize, r2threshold, windowmargin, snpsForChr);
                        for (SortableSNP snp : snpsForChr) {
                            ArrayList<SortableSNP> ldSNPs = proxies.get(snp);
                            String locusName = "locus-" + locusCounter + "-chr" + chr;
                            int minPos = snp.chrpos;
                            int maxPos = snp.chrpos;
                            if (ldSNPs.size() == 1) {
                                if (add1) {
                                    snpBedOutput.writeln(locusName + "\t" + (snp.chrpos) + "\t" + (snp.chrpos+1) + "\t" + snp.name + "\t" + snp.effect);
                                } else {
                                    snpBedOutput.writeln(locusName + "\t" + (snp.chrpos) + "\t" + (snp.chrpos) + "\t" + snp.name + "\t" + snp.effect);
                                }
                                minPos = snp.chrpos;
                                maxPos = snp.chrpos;
                            } else {

                                for (SortableSNP s : ldSNPs) {
                                    if(add1){
                                        snpBedOutput.writeln(locusName + "\t" + s.chrpos + "\t" + (s.chrpos+1) + "\t" + s.name + "\t" + s.effect);
                                    } else {
                                        snpBedOutput.writeln(locusName + "\t" + s.chrpos + "\t" + (s.chrpos) + "\t" + s.name + "\t" + s.effect);
                                    }
                                    
                                    int pos = s.chrpos;
                                    if (pos > maxPos) {
                                        maxPos = pos;
                                    }
                                    if (pos < minPos) {
                                        minPos = pos;
                                    }
                                }
                            }

                            // define the window/region
                            // extend the region by certain bp
                            minPos -= windowmargin;
                            maxPos += windowmargin;
                            if (minPos < 0) {
                                minPos = 0;
                            }
                            
                            if(maxPos == 266 && minPos > 266){
                                System.out.println(taskctr+ " Warning: there's something wrong here: "+snp.chr+" - "+ snp.chrpos  +  "\t"+ldSNPs.size());
                            }

                            Feature f = new Feature();
                            f.setChromosome(chrE);
                            f.setName(locusName);
                            f.setStart(minPos);
                            f.setStop(maxPos);
                            windows.add(f);
                            regionBedOutput.writeln(locusName + "\t" + minPos + "\t" + maxPos);
                            locusCounter++;
                        }
                        /*
                         // we have our LD snps... create the window...
                         // write the SNPs as single bp features
                    
                         */
                    }
                }
            }

            regionBedOutput.close();
            snpBedOutput.close();

//            System.out.println(taskctr + "- " + listOfAnnotationFiles.length + " annotation files found...");
            // filter the annotations
            for (int i = 0; i < listOfAnnotationFiles.length; i++) {
                TextFile annotationInput = new TextFile(listOfAnnotationFiles[i], TextFile.R);

                File f = new File(listOfAnnotationFiles[i]);
                String n = f.getName();
//                System.out.println(taskctr + "- Filtering " + n + "\noutput: " + outputDir + n);

                String[] annelems = annotationInput.readLineElems(TextFile.tab);
                Chromosome prevChr = null;
                ArrayList<Feature> bedfileFeatures = new ArrayList<Feature>();
                while (annelems != null) {

                    Chromosome chrE = Chromosome.parseChr(annelems[0]);
                    int start = Integer.parseInt(annelems[1]);
                    int stop = Integer.parseInt(annelems[2]);
                    if (prevChr == null || !chrE.equals(prevChr)) {
                        prevChr = chrE;
                    }

                    Feature feat = new Feature();
                    feat.setChromosome(chrE);
                    feat.setStart(start);
                    feat.setStop(stop);
                    bedfileFeatures.add(feat);

                    annelems = annotationInput.readLineElems(TextFile.tab);
                }
                annotationInput.close();
//                System.out.println(taskctr + "- Found: " + bedfileFeatures.size() + " features");
//
//                System.out.println(taskctr + "- Writing file: " + outputDir + n);
                if (n.endsWith(".gz")) {
                    n = n.substring(0, n.length() - 3);
                }

                TextFile annotationOutput = new TextFile(outputDir + n, TextFile.W);
                for (int q = 0; q < bedfileFeatures.size(); q++) {
                    Feature f1 = bedfileFeatures.get(q);
                    for (int w = 0; w < windows.size(); w++) {
                        Feature f2 = windows.get(w);
                        if (f2.overlaps(f1)) {
                            annotationOutput.writeln(f2.getName() + "\t" + f1.getStart() + "\t" + f1.getStop());
                        }
                    }
                }
                annotationOutput.close();
            }

            System.out.println(taskctr + "- Processed total of " + total + " SNPs");
        } catch (IOException e) {

            System.out.println("Error for set: " + snpfile + " - " + tabixdir + " - " + r2threshold + " - "
                    + ldWindowSize + " - "
                    + windowmargin + " - "
                    + outputDir + " - "
                    + taskctr + "\n" + e.getMessage());
            System.exit(-1);
        }
    }

    private ArrayList<SortableSNP> loadSNPs(String inputfile) throws IOException {
        TextFile tf = new TextFile(inputfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        ArrayList<SortableSNP> snps = new ArrayList<SortableSNP>();
        int id = 0;
        while (elems != null) {
            String snp = elems[0];
            String chr = elems[1];
            Integer position = null;
            String pos = elems[2];
            try {
                position = Integer.parseInt(pos);
            } catch (NumberFormatException e) {
                System.out.println(pos + " is not a number..");
            }
            if (position != null) {
                byte chrb = (byte) Chromosome.parseChr(chr).getNumber();
                SortableSNP snpObj = new SortableSNP(snp, id, chrb, position, SortableSNP.SORTBY.CHRANDCHRPOS);
                snps.add(snpObj);
                id++;
            }

            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        return snps;
    }

    private HashMap<SortableSNP, ArrayList<SortableSNP>> getProxiesForSNPGosiaFiles(String dir, Byte chr, int ldWindowSize, double r2threshold, ArrayList<SortableSNP> snpsForChr) throws IOException {
        // open tabix file

        // iterate the variants
        HashMap<SortableSNP, ArrayList<SortableSNP>> proxies = new HashMap<SortableSNP, ArrayList<SortableSNP>>();

        for (SortableSNP snp : snpsForChr) {
            // get the variants within the windowsize of snp

            // open file 
            ArrayList<SortableSNP> ldSNPs = new ArrayList<SortableSNP>();
            int proxyId = 0;

            String snpldfile = dir + "chr" + chr + "/results_ld_" + snp.name + ".txt";
            if (!Gpio.exists(snpldfile)) {
                System.out.println("No LD information for: " + snp.name);
                ldSNPs.add(snp);
                proxies.put(snp, ldSNPs);
            } else {
                TextFile tf = new TextFile(snpldfile, TextFile.R);
                String proxyDataLn = tf.readLine();
                HashSet<String> visitedSNPs = new HashSet<String>();
                while (proxyDataLn != null) {

                    String[] lnElems = Strings.whitespace.split(proxyDataLn);

                    String snp1Str = lnElems[1];
                    String snp2Str = lnElems[3];
                    if (snp1Str.equals(snp.name) || snp2Str.equals(snp.name)) {

                        String r2Str = lnElems[5];
                        Double r2 = Double.parseDouble(r2Str);
                        if (r2 >= r2threshold) {

                            if (!visitedSNPs.contains(snp1Str)) {

                                String pos1Str = lnElems[2];

                                Integer pos1 = Integer.parseInt(pos1Str);
                                int distance = Math.abs(pos1 - snp.chrpos);

                                if (distance <= ldWindowSize) {
                                    // check LD
                                    // include variant
                                    ldSNPs.add(new SortableSNP(snp1Str, proxyId, chr, pos1, r2, SortableSNP.SORTBY.CHRANDCHRPOS));
                                    proxyId++;
                                }
                                visitedSNPs.add(snp1Str);

                            }

                            if (!visitedSNPs.contains(snp2Str)) {

                                String pos2Str = lnElems[4];
                                Integer pos2 = Integer.parseInt(pos2Str);
                                int distance = Math.abs(pos2 - snp.chrpos);

                                if (distance <= ldWindowSize) {
                                    // check LD
                                    // include variant
                                    ldSNPs.add(new SortableSNP(snp2Str, proxyId, chr, pos2, r2, SortableSNP.SORTBY.CHRANDCHRPOS));
                                    proxyId++;
                                }
                                visitedSNPs.add(snp2Str);

                            }

                        }
                    }

                    proxyDataLn = tf.readLine();
                }
                tf.close();

//            System.out.println(chr + "\t" + snp.name + "\tproxies: " + ldSNPs.size());
                proxies.put(snp, ldSNPs);

            }

        }

        return proxies;
    }

}
