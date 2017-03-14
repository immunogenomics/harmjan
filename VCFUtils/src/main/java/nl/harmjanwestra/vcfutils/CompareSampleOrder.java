package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 3/13/17.
 */
public class CompareSampleOrder {

	public void run(String vcf1, String vcf2) throws IOException {
		System.out.println("VCF1: " + vcf1);
		System.out.println("VCF2: " + vcf2);
		VCFGenotypeData d1 = new VCFGenotypeData(vcf1);
		VCFGenotypeData d2 = new VCFGenotypeData(vcf2);


		ArrayList<String> l1 = d1.getSamples();
		ArrayList<String> l2 = d2.getSamples();
		System.out.println(l1.size() + " samples in vcf1");
		System.out.println(l2.size() + " samples in vcf2");

		int unequal = 0;

		HashSet<String> s1 = new HashSet<String>();
		s1.addAll(l1);

		HashSet<String> s2 = new HashSet<String>();
		s2.addAll(l2);

		int shared = 0;
		for (String sample : s1) {
			if (s2.contains(sample)) {
				shared++;
			}
		}
		System.out.println(shared + " samples in vcf1 also in vcf2");

		shared = 0;
		for (String sample : s2) {
			if (s1.contains(sample)) {
				shared++;
			}
		}
		System.out.println(shared + " samples in vcf2 also in vcf1");

		if (l1.size() != l2.size()) {
			System.out.println("different sizes: " + l1.size() + " vs " + l2.size());
		} else {
			for (int i = 0; i < l1.size(); i++) {
				String sa1 = l1.get(i);
				String sa2 = l2.get(i);
				if (!sa1.equals(sa2)) {
					System.out.println(i + " expecting\t" + sa1 + "\tfound\t" + sa2);
					unequal++;
				}
			}
		}
		System.out.println(unequal + " samples with different ordering");


	}
}
