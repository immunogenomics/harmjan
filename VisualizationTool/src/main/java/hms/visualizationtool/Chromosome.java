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
public enum Chromosome {

    ONE(1, "Chr1", 248956422),
    TWO(2, "Chr2", 242193529),
    THREE(3, "Chr3", 198295559),
    FOUR(4, "Chr4", 19021455),
    FIVE(5, "Chr5", 181538259),
    SIX(6, "Chr6", 170805979),
    SEVEN(7, "Chr7", 159345973),
    EIGHT(8, "Chr8", 145138636),
    NINE(9, "Chr9", 138394717),
    TEN(10, "Chr10", 133797422),
    ELEVEN(11, "Chr11", 135086622),
    TWELVE(12, "Chr12", 133275309),
    THIRTEEN(13, "Chr13", 114364328),
    FOURTEEN(14, "Chr14", 107043718),
    FIFTEEN(15, "Chr15", 101991189),
    SIXTEEN(16, "Chr16", 90338345),
    SEVENTEEN(17, "Chr17", 83257441),
    EIGHTEEN(18, "Chr18", 80373285),
    NINETEEN(19, "Chr19", 58617616),
    TWENTY(20, "Chr20", 64444167),
    TWENTYONE(21, "Chr21", 46709983),
    TWENTYTWO(22, "Chr22", 50818468),
    X(23, "ChrX", 156040895),
    Y(24, "ChrY", 57227415),
    MT(25, "ChrMT", 1),
    NA(0, "N/A", 1);

    private final int number;
    private final String name;
    private final int length;

    private Chromosome(int num, String name, int length) {
        this.number = num;
        this.name = name;
        this.length = length;
    }

    public String getName() {
        return name;
    }

    public int getLength() {
        return length;
    }
    
     public int getNumber() {
        return number;
    }

    public static Chromosome parseChr(String chrStr) {
        chrStr = chrStr.toLowerCase().trim();
        if (chrStr.equals("chr1")) {
            return Chromosome.ONE;
        }
        if (chrStr.equals("chr2")) {
            return Chromosome.TWO;
        }
        if (chrStr.equals("chr3")) {
            return Chromosome.THREE;
        }
        if (chrStr.equals("chr4")) {
            return Chromosome.FOUR;
        }
        if (chrStr.equals("chr5")) {
            return Chromosome.FIVE;
        }
        if (chrStr.equals("chr6")) {
            return Chromosome.SIX;
        }
        if (chrStr.equals("chr7")) {
            return Chromosome.SEVEN;
        }
        if (chrStr.equals("chr8")) {
            return Chromosome.EIGHT;
        }
        if (chrStr.equals("chr9")) {
            return Chromosome.NINE;
        }
        if (chrStr.equals("chr10")) {
            return Chromosome.TEN;
        }
        if (chrStr.equals("chr11")) {
            return Chromosome.ELEVEN;
        }
        if (chrStr.equals("chr12")) {
            return Chromosome.TWELVE;
        }
        if (chrStr.equals("chr13")) {
            return Chromosome.THIRTEEN;
        }
        if (chrStr.equals("chr14")) {
            return Chromosome.FOURTEEN;
        }
        if (chrStr.equals("chr15")) {
            return Chromosome.FIFTEEN;
        }
        if (chrStr.equals("chr16")) {
            return Chromosome.SIXTEEN;
        }
        if (chrStr.equals("chr17")) {
            return Chromosome.SEVENTEEN;
        }
        if (chrStr.equals("chr18")) {
            return Chromosome.EIGHTEEN;
        }
        if (chrStr.equals("chr19")) {
            return Chromosome.NINETEEN;
        }
        if (chrStr.equals("chr20")) {
            return Chromosome.TWENTY;
        }
        if (chrStr.equals("chr21")) {
            return Chromosome.TWENTYONE;
        }
        if (chrStr.equals("chr22")) {
            return Chromosome.TWENTYTWO;
        }
        if (chrStr.equals("chry")) {
            return Chromosome.Y;
        }
        if (chrStr.equals("chrx")) {
            return Chromosome.X;
        }

        return Chromosome.NA;
    }
}
