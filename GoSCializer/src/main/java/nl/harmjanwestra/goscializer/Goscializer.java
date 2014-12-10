/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.goscializer;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author hwestra
 */
public class Goscializer {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try{
            TextFile tf1 = new TextFile(args[0], TextFile.R);
            TextFile tf2 = new TextFile(args[1], TextFile.W);
            tf2.writeln(tf1.readLine());
            String[] elems = tf1.readLineElems(TextFile.tab);
            while(elems!=null){
                
                if(Math.random() > 0.5){
                    elems[elems.length-1] = ""+0;
                } else {
                    elems[elems.length-1] = ""+1;
                }
                tf2.writeln(Strings.concat(elems, Strings.tab));
                
                elems = tf1.readLineElems(TextFile.tab);
            }
            tf1.close();
            tf2.close();
                    
            
            
        } catch (IOException e){
            e.printStackTrace();
        }
        
//        if (args.length < 13) {
//            System.out.println("usage:");
//            // System.out.println("snpfile tabixdir r2 ldwindow ldwindowextend annotationfile outputdir");
//            System.out.println("indir snpfilename numsets tabixdir r2 ldwindow ldwindowextend annotationfile outputdir nrbootstraps bootstrapout add1 nrthreads");
//        } else {
////            String snpfile = args[0];
////            String tabixdir = args[1];
////            double r2thresholdStr = Double.parseDouble(args[2]);
////            int ldWindowSizeStr = Integer.parseInt(args[3]);
////            int ldWindowMargin = Integer.parseInt(args[4]);
////            String outputdir = args[6];
////
////            try {
////                Goscializer gsc = new Goscializer();
////                String[] annotionFileList = gsc.loadAnnotationList(args[5]);
////                gsc.run(snpfile, tabixdir, r2thresholdStr, ldWindowSizeStr, ldWindowMargin, annotionFileList, outputdir);
////            } catch (IOException ex) {
////                Logger.getLogger(Goscializer.class.getName()).log(Level.SEVERE, null, ex);
////            }
//
//            String indir = args[0];
//            String snpfilename = args[1];
//            int numsets = Integer.parseInt(args[2]);
//
//            String tabixdir = args[3];
//            double r2thresholdStr = Double.parseDouble(args[4]);
//            int ldWindowSizeStr = Integer.parseInt(args[5]);
//            int ldWindowMargin = Integer.parseInt(args[6]);
//
//            String annotationfile = args[7];
//            String outputdir = args[8];
//            int nrbootstraps = Integer.parseInt(args[9]);
//            String bootstrapout = args[10];
//            boolean add1 = Boolean.parseBoolean(args[11]);
//            Integer nrThreads = Integer.parseInt(args[12]);
//            
//            try {
//                Goscializer gsc = new Goscializer();
//
//                gsc.runCreateRegionsForFigure2(tabixdir, ldWindowSizeStr, r2thresholdStr, ldWindowMargin, indir, snpfilename, numsets, annotationfile, outputdir, nrbootstraps, bootstrapout, add1, nrThreads);
//            } catch (IOException ex) {
//                Logger.getLogger(Goscializer.class.getName()).log(Level.SEVERE, null, ex);
//            }
//
//        }

        System.exit(0);
    }

    public String[] loadAnnotationList(String annotationFileList) throws IOException {
        TextFile annotF = new TextFile(annotationFileList, TextFile.R);
        String[] listOfAnnotationFiles = annotF.readAsArray();
        annotF.close();
        return listOfAnnotationFiles;
    }

    public void runCreateRegionsForFigure2(String tabixdir, int windowSize, double ldthreshold, int extendWindow,
            String inDir, String snpfilename, int numSets, String annotationfile,
            String outdir, int nrBootstraps, String bootstrapoutdir, boolean add1, Integer nrThreads) throws IOException {
        inDir = Gpio.formatAsDirectory(inDir);

        bootstrapoutdir = Gpio.formatAsDirectory(bootstrapoutdir);
        Gpio.createDir(bootstrapoutdir);
        TextFile commandsOut = new TextFile(bootstrapoutdir + "commands.sh", TextFile.W);

        ExecutorService executor = Executors.newFixedThreadPool(nrThreads);

        for (int set = 1; set < numSets + 1; set++) {
            String snpInputFile = inDir + "set" + set + "/" + snpfilename;
            if (Gpio.exists(snpInputFile)) {
                String[] annotfile = new String[]{annotationfile};

                File f = new File(annotationfile);

                String setOutputDir = Gpio.formatAsDirectory(outdir + "set" + set);
                Gpio.createDir(setOutputDir);

                GoscializerTask task = new GoscializerTask(set, snpInputFile, tabixdir,
                        ldthreshold, windowSize, extendWindow, annotfile, setOutputDir, add1);
                executor.execute(task);

                // run(snpInputFile, tabixdir, ldthreshold, windowSize, extendWindow, annotfile, setOutputDir);
                String annotationfilename = f.getName();
                if (annotationfilename.endsWith(".gz")) {
                    annotationfilename = annotationfilename.substring(0, annotationfilename.length() - 3);
                }
                String command = "python /home/unix/hwestra/Scripts/block_bootstrap-0.8.1/block_bootstrap.py -1 " + setOutputDir + "SNPs.bed "
                        + "-2 " + setOutputDir + annotationfilename + " "
                        + "-d " + setOutputDir + "Regions.bed "
                        + "-n " + nrBootstraps + " "
                        + "-r 0.01 -s 0.01 -t bc -B -o " + bootstrapoutdir + "null-set" + set + ".txt &> " + bootstrapoutdir + "set" + set + ".txt";

                commandsOut.writeln(command);

                // for each SNP
                // get all proxies
            }

        }

        commandsOut.close();

        executor.shutdown();
        while (!executor.isTerminated()) {
        }
        System.out.println("Finished all threads");
    }

}
