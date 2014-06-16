/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.visualizationtool;

import java.util.Comparator;

/**
 *
 * @author Harm-Jan
 */
public class FeatureComparator implements Comparator<Feature> {

    @Override
    public int compare(Feature t, Feature t1) {
        
        if (t.equals(t1)) {
            return 0;
        }

        if (t.overlaps(t1)) {

            if (t.getStart() >= t1.getStart() && t.getStop() <= t1.getStop()) {
                if (t.getStrand() == t1.getStrand()) {
                    return 0;
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
