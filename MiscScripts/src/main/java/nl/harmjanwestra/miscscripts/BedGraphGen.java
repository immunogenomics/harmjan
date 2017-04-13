package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by hwestra on 4/26/16.
 */
public class BedGraphGen {

	public static void main(String[] args) {

		Feature region = new Feature(Chromosome.TWO, 204446380, 204816382);
		try {

			BedGraphGen g = new BedGraphGen();
			Gpio.createDir("/Data/tmp/2016-04-26/");
//			g.gen("/Data/tmp/2016-04-26/sinus.bed.gz", region);
			g.gen2("/Data/tmp/2016-04-26/sinusSum.bed.gz", region);


		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void gen(String out, Feature region) throws IOException {

		TextFile outf = new TextFile(out, TextFile.W);

		double steps = (2 * Math.PI) / 10000;

		int ctr = 0;
		for (int i = region.getStart(); i < region.getStop(); i++) {
			double q = ctr * steps;
			if (q > (2 * Math.PI)) {
				ctr = 0;
			} else {
				ctr++;
			}

			double v = 50 * Math.sin(q); // max 50, min -50;
			double pos = 0;
			double neg = 0;
			if (v > 0) {
				pos = v;
			} else {
				neg = Math.abs(v);
			}
			v = Math.abs(v);

			outf.writeln(region.getChromosome().toString() + "\t" + i + "\t" + i + "\t" + v + "\t" + pos + "\t" + neg);
		}
		outf.close();
	}

	public void gen2(String out, Feature region) throws IOException {

		TextFile outf = new TextFile(out, TextFile.W);

		double steps = (2 * Math.PI) / 10000;

		int ctr = 0;
		for (int i = region.getStart(); i < region.getStop(); i++) {
			double q = ctr * steps;
			if (q > (2 * Math.PI)) {
				ctr = 0;
			} else {
				ctr++;
			}

			double v = 100 * Math.sin(q); // max 50, min -50;
			v = Math.abs(v);
			double pos = v / 2;
			double neg = v / 2;


			outf.writeln(region.getChromosome().toString() + "\t" + i + "\t" + i + "\t" + v + "\t" + pos + "\t" + neg);
		}
		outf.close();
	}
}
