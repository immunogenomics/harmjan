/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.resequencing;

import java.io.File;
import picard.sam.MarkDuplicates;

/**
 *
 * @author hwestra
 */
public class MarkDups extends MarkDuplicates {

    public MarkDups(File in, File out, File metricsFile, File tmp) {
        String[] args = new String[4];
        args[0] = "INPUT="+in.getAbsolutePath();
        args[1] = "OUTPUT="+out.getAbsolutePath();
        args[2] = "TMP_DIR="+tmp.getAbsolutePath();
        args[3] = "METRICS_FILE="+metricsFile.getAbsolutePath();
         
        new MarkDuplicates().instanceMain(args);
    }

//    void go() {
//        doWork();
//    }
}
