/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.visualizationtool;

/**
 *
 * @author Harm-Jan
 */
public enum Strand {

    POS("+",1), NEG("-",2), NA("NA",0);

    static Strand parseStr(String string) {
        if (string.equals("+")) {
            return Strand.POS;
        }

        if (string.equals("-")) {
            return Strand.NEG;
        }

        return Strand.NA;
    }
    
    int num;
    String strand;

    private Strand(String s, int num) {
        this.strand = s;
        this.num = num;
    }

    @Override
    public String toString() {
        return strand;
    }

    int getNumber() {
        return num;
    }
}
