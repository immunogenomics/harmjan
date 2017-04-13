package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;

/**
 * Created by hwestra on 1/16/15.
 */
public class DetermineFileSizeForGenome {

	public static void main(String[] args) {
		Chromosome[] chromosomes = Chromosome.values();
		long sumOfBases = 0;
		for (Chromosome chr : chromosomes) {
			int nr = chr.getNumber();

			if (nr < 23) {
				int len = chr.getLength();
				sumOfBases += len;
				System.out.println(chr.getName() + "\t" + nr + "\t" + len + "\t" + sumOfBases+"\t"+Gpio.humanizeFileSize(len*8)+"\t"+len/1000);
			}
		}
		sumOfBases = sumOfBases * 8;
		String fs = Gpio.humanizeFileSize(sumOfBases);
		System.out.println(fs);
		System.out.println(sumOfBases/1000);
	}
}
