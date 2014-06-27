/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.utilities.gtf;

import hms.hwestra.utilities.features.Gene;
import java.util.Comparator;

/**
 *
 * @author Harm-Jan
 */
class GeneComparator implements Comparator<Gene> {

    @Override
    public int compare(Gene gene1, Gene gene2) {
        if (gene1.equals(gene2)) {
            return 0;
        }
        if (gene1.getChromosome().getNumber() > gene2.getChromosome().getNumber()) {
            return 1;
        } else if (gene1.getChromosome().getNumber() > gene2.getChromosome().getNumber()) {
            return -1;
        }
        if (gene1.overlaps(gene2)) {
            return 0;
        } else {
            if (gene1.getStart() > gene2.getStop()) {
                return 1;
            } else {
                return -1;
            }
        }

    }

}
