/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.utilities.features;

/**
 * @author Harm-Jan
 */
public enum Chromosome {

	ONE(1, "Chr1", 248956422),
	TWO(2, "Chr2", 242193529),
	THREE(3, "Chr3", 198295559),
	FOUR(4, "Chr4", 190214555),
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
	NA(26, "N/A", 1);

	private final int number;
	private final String name;
	private final int length;

	private Chromosome(int num, String name, int length) {
		this.number = num;
		this.name = name;
		this.length = length;
	}

	public static Chromosome parseChr(String chrStr) {
		chrStr = chrStr.toLowerCase().trim();
		if (chrStr.equals("chr1") || chrStr.equals("1")) {
			return Chromosome.ONE;
		}
		if (chrStr.equals("chr2") || chrStr.equals("2")) {
			return Chromosome.TWO;
		}
		if (chrStr.equals("chr3") || chrStr.equals("3")) {
			return Chromosome.THREE;
		}
		if (chrStr.equals("chr4") || chrStr.equals("4")) {
			return Chromosome.FOUR;
		}
		if (chrStr.equals("chr5") || chrStr.equals("5")) {
			return Chromosome.FIVE;
		}
		if (chrStr.equals("chr6") || chrStr.equals("6")) {
			return Chromosome.SIX;
		}
		if (chrStr.equals("chr7") || chrStr.equals("7")) {
			return Chromosome.SEVEN;
		}
		if (chrStr.equals("chr8") || chrStr.equals("8")) {
			return Chromosome.EIGHT;
		}
		if (chrStr.equals("chr9") || chrStr.equals("9")) {
			return Chromosome.NINE;
		}
		if (chrStr.equals("chr10") || chrStr.equals("10")) {
			return Chromosome.TEN;
		}
		if (chrStr.equals("chr11") || chrStr.equals("11")) {
			return Chromosome.ELEVEN;
		}
		if (chrStr.equals("chr12") || chrStr.equals("12")) {
			return Chromosome.TWELVE;
		}
		if (chrStr.equals("chr13") || chrStr.equals("13")) {
			return Chromosome.THIRTEEN;
		}
		if (chrStr.equals("chr14") || chrStr.equals("14")) {
			return Chromosome.FOURTEEN;
		}
		if (chrStr.equals("chr15") || chrStr.equals("15")) {
			return Chromosome.FIFTEEN;
		}
		if (chrStr.equals("chr16") || chrStr.equals("16")) {
			return Chromosome.SIXTEEN;
		}
		if (chrStr.equals("chr17") || chrStr.equals("17")) {
			return Chromosome.SEVENTEEN;
		}
		if (chrStr.equals("chr18") || chrStr.equals("18")) {
			return Chromosome.EIGHTEEN;
		}
		if (chrStr.equals("chr19") || chrStr.equals("19")) {
			return Chromosome.NINETEEN;
		}
		if (chrStr.equals("chr20") || chrStr.equals("20")) {
			return Chromosome.TWENTY;
		}
		if (chrStr.equals("chr21") || chrStr.equals("21")) {
			return Chromosome.TWENTYONE;
		}
		if (chrStr.equals("chr22") || chrStr.equals("22")) {
			return Chromosome.TWENTYTWO;
		}
		if (chrStr.equals("chr23") || chrStr.equals("23")) {
			return Chromosome.X;
		}
		if (chrStr.equals("chry") || chrStr.equals("y")) {
			return Chromosome.Y;
		}
		if (chrStr.equals("chrx") || chrStr.equals("x")) {
			return Chromosome.X;
		}
		if (chrStr.equals("chrm") || chrStr.equals("chrmt") || chrStr.equals("m") || chrStr.equals("mt")) {
			return Chromosome.MT;
		}

		return Chromosome.NA;
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

	@Override
	public String toString() {
		return this.getName();
	}


	public int compare(Chromosome other) {
		if (this.equals(other)) {
			return 0;
		} else if (this.number > other.number) {
			return 1;
		} else {
			return 0;
		}
	}

	public boolean equals(Chromosome other) {
		return this.number == other.number;
	}

	public boolean isAutosome() {
		if (this.getNumber() > 0 && this.getNumber() < 23) {
			return true;
		} else {
			return false;
		}
	}

}
