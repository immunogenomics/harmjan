package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;

import java.io.IOException;

/**
 * Created by hwestra on 5/17/16.
 */
public class GTTest {

	public static void main(String[] args) {
		try {

			GTTest t = new GTTest();
			t.run();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run() throws IOException {

		String vcf1 = "/Data/ImmunoChip/RA/2015-09-23-IncX/ES/genotypes-filtered-sorted-Chr1.vcf.gz";
		VCFGenotypeData data1 = new VCFGenotypeData(vcf1);

		while (data1.hasNext()) {
			VCFVariant var = data1.next();
			Chromosome chr = Chromosome.parseChr(var.getChr());

			if (!chr.equals(Chromosome.X)) {
				String varln = data1.getNextLn();

				String[] lnelems = varln.split("\t");

				String id1 = var.getId();
				String id2 = lnelems[2];

				if (!id1.equals(id2)) {
					System.out.println("found " + id1 + " expected " + id2);

				}

//				var.parseGenotypes(varln, VCFVariant.PARSE.ALL);

			}

		}
		data1.close();
	}
}
