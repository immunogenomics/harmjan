/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.visualizationtool;

/**
 *
 * @author hwestra
 */
public class Main {

    public static void main(String[] args) {

//            try {
//                String dir = args[0];
//                String outdir = args[1];
//
//                String[] files = Gpio.getListOfFiles(dir, "markers");
//                TextFile out = new TextFile(outdir, TextFile.W);
//                for (String f : files) {
//                    TextFile tf = new TextFile(f, TextFile.R);
//                    String[] elems = tf.readLineElems(TextFile.tab);
//
//                    String[] felems = f.split("\\.");
//                    String chr = felems[1];
//                    
//                    while (elems != null) {
//                        out.writeln(elems[0] + "\t" + chr + "\t" + elems[1]);
//                        
//                        elems = tf.readLineElems(TextFile.tab);
//                    }
//
//                    tf.close();
//
//                }
//                out.close();
//
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//        }
        GosiaViz2 viz = new GosiaViz2();

        if (args.length < 8) {
            System.out.println("Usage: Viz.jar locusscore snpoverlap gtf bed peak snpannot outputfilename extrawindowsize");
        } else {

            String locusScore = args[0];
            String snpOverlap = args[1];
            String gtf = args[2];
            String bed = args[3];
            String peakfile = args[4];
            String snpannot = args[5];
            String outputfilename = args[6];

            Integer extraWindowSize = 100;

            try {
                extraWindowSize = Integer.parseInt(args[7]);
            } catch (NumberFormatException e) {
                System.out.println(args[7] + " is not a number. Will default to extrawindowsize == 100");
            }

            try {
                viz.run(locusScore, snpOverlap, gtf, bed, peakfile, snpannot, outputfilename, extraWindowSize);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        System.exit(0);
    }

}
