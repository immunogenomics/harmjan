/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.rnaaccess;

import java.io.File;
import java.io.IOException;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;

/**
 *
 * @author Harm-Jan
 */
public class CollectTrackingFiles {

    public static void main(String[] args) {

        try {
            String startdir = "/Data/Projects/2014-Epipilot/rna-seq/coreunit/";
            File f = new File(startdir);
            String outdir = "/Data/Projects/2014-Epipilot/rna-seq/tracking/";
            CollectTrackingFiles c = new CollectTrackingFiles();
            c.run(f, outdir);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void run(File dir, String outdir) throws IOException {

        File[] lsof = dir.listFiles();
        if (lsof == null) {
            System.err.println("ERROR: " + dir.getName() + " gives null results");
        } else {
            for (File f : lsof) {
                if (f.isDirectory()) {
                    run(f, outdir);
                } else {
                    if (f.getName().contains("fpkm_tracking")) {
//                        System.out.println(f+" in "+outdir);
                        if (f.getAbsolutePath().contains("genes")) {

                            String dirName = dir.getAbsolutePath();
//                            sampleName = sampleName.substring(i);

                            String[] elems = dirName.split("/");

                            String sampleName = elems[elems.length - 1];

                            System.out.println(sampleName);

                            System.out.println("Found path: " + dir.getAbsolutePath() + " / " + f.getName());
//                            Gpio.createDir(outdir + "/" + sampleName);
                            Gpio.copyFile(dir.getAbsolutePath() + "/" + f.getName(), outdir + "/" + sampleName + "-" + f.getName());
                        }

                    }
                }
            }
        }

    }

}
