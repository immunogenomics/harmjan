package nl.harmjanwestra.gwas.tasks;

import cern.colt.matrix.tbit.BitVector;
import nl.harmjanwestra.utilities.vcf.VCFVariant;

import java.util.ArrayList;
import java.util.concurrent.Callable;

/**
 * Created by Harm-Jan on 06/23/16.
 */
public class LRTestHaploTask implements Callable<BitVector[]> {

	ArrayList<VCFVariant> variants;
	int i = -1;

	public LRTestHaploTask(int i, ArrayList<VCFVariant> variants) {
		this.i = i;
		this.variants = variants;
	}

	@Override
	public BitVector[] call() throws Exception {
		BitVector haplotype1 = new BitVector(variants.size());
		BitVector haplotype2 = new BitVector(variants.size());

		for (int v = 0; v < variants.size(); v++) {
			VCFVariant variant = variants.get(v);
			double[][] alleles = variant.getGenotypeAlleles();
			if (alleles[i][0] == 1) {
				haplotype1.set(v);
			}
			if (alleles[i][1] == 1) {
				haplotype2.set(v);
			}
		}
		return new BitVector[]{haplotype1, haplotype2};
	}
}
