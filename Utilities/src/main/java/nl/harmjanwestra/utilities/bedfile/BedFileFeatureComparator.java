/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.utilities.bedfile;

import nl.harmjanwestra.utilities.features.BedFileFeature;
import nl.harmjanwestra.utilities.enums.Strand;
import java.util.Comparator;

/**
 *
 * @author Harm-Jan
 */
public class BedFileFeatureComparator implements Comparator<BedFileFeature> {

    @Override
    public int compare(BedFileFeature t, BedFileFeature t1) {
        if(t.getChromosome().getNumber() > t1.getChromosome().getNumber()){
            return 1;
        } else if(t.getChromosome().getNumber() < t1.getChromosome().getNumber()){
            return -1;
        }
        
        if (t.equals(t1)) {
            return 0;
        }

        if (t.overlaps(t1)) {

            if (t.getStart() >= t1.getStart() && t.getStop() <= t1.getStop()) {
                if (t.getStrand() == t1.getStrand()) {
                    return 1;
                }
                if (t.getStrand() == Strand.NEG){
                    return 1;
                } else {
                    return -1;
                }
                        
                
            } else if (t.getStart() <= t1.getStart()) {
                return -1;
            } else {
                return 1;
            }
        } else {
            if (t.getStart() > t1.getStop()) {
                return 1;
            } else {
                return -1;
            }
        }
    }

}
