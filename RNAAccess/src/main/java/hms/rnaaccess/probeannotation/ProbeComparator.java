/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.rnaaccess.probeannotation;

import hms.hwestra.utilities.features.Strand;
import java.util.Comparator;

/**
 *
 * @author Harm-Jan
 */
public class ProbeComparator implements Comparator<Probe> {

    @Override
    public int compare(Probe probe1, Probe probe2) {
        if(probe1.equals(probe2)){
            return 0;
        }
        if(probe1.getChromosome().getNumber() > probe2.getChromosome().getNumber()){
            return 1;
        } else if (probe1.getChromosome().getNumber() < probe2.getChromosome().getNumber()){
            return -1;
        } 
        if(probe1.overlaps(probe2)){
            return 0;
        } else {
            if (probe1.getStart() > probe2.getStop()) {
                return 1;
            } else {
                return -1;
            }
        }
    }

}
