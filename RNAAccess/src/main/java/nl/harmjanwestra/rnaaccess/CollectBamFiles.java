package nl.harmjanwestra.rnaaccess;

import java.io.File;
import java.io.IOException;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template path, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Harm-Jan
 */
public class CollectBamFiles {

    public static void main(String[] args) {
        if (args.length != 2) {
            System.out.println("Usage: collect indir outdir");
        } else {
            try {

                String startdir = args[0];
                File f = new File(startdir);
                String outdir = args[1];
                CollectBamFiles c = new CollectBamFiles();
                c.run(f, outdir);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        System.exit(0);
    }
    private TextFile pathsOut;

    public void run(File dir, String outdir) throws IOException {
        pathsOut = new TextFile(outdir + "BamFilePaths.txt", TextFile.W);
        iterate(dir);
        pathsOut.close();
    }

    private void iterate(File dir) throws IOException {
        File[] lsof = dir.listFiles();
        if (lsof == null) {
            System.err.println("ERROR: " + dir.getName() + " gives null results");
        } else {
            for (File f : lsof) {
                if (f.isDirectory()) {
                    System.out.println("Entering directory: " + f.getAbsolutePath());
                    iterate(f);
                } else {
                    if (f.getName().endsWith(".bam")) {
                        System.out.println(f.getAbsolutePath());
                        pathsOut.writeln(f.getAbsolutePath());
                    }
                }
            }
        }
    }
}
