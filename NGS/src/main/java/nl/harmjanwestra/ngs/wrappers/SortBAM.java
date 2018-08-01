/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.ngs.wrappers;

import java.io.File;
import picard.sam.SortSam;

/**
 *
 * @author hwestra
 */
public class SortBAM extends SortSam {

    public SortBAM(File in, File out, File tmp) {
        /* 
         /tools/picard-tools-1.32/SortSam.jar \ 
         INPUT=/output/filename.b37_1kg.bam \ 
         OUTPUT=/output/filename.b37_1kg.sorted.bam \ 
         SORT_ORDER=coordinate \ 
         VALIDATION_STRINGENCY=LENIENT \ 
         TMP_DIR=/local
                
                
         */
        String[] args = new String[5];
        args[0] = "INPUT=" + in.getAbsolutePath();
        args[1] = "OUTPUT=" + out.getAbsolutePath();
        args[2] = "TMP_DIR=" + tmp.getAbsolutePath();
        args[3] = "SORT_ORDER=coordinate";
        args[4] = "VALIDATION_STRINGENCY=LENIENT";
        
        new SortSam().instanceMain(args);
    }
}
