/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.features.Chromosome;

/**
 *
 * @author hwestra
 */
public class ChrSize {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        Chromosome[] values = Chromosome.values();
        for (Chromosome v : values) {
            System.out.println(v.getName().toLowerCase() + "\t0\t" + v.getLength());

        }
    }

}
