package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by Harm-Jan on 04/20/16.
 */
public class VCFListSamples {

	public void run(String vcf) throws IOException {
		VCFGenotypeData data = new VCFGenotypeData(vcf);
		ArrayList<String> samples = data.getSamples();
		for (int i = 0; i < samples.size(); i++) {
			System.out.println(samples.get(i));
		}
		data.close();
	}
}
