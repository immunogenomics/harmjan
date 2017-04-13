package nl.harmjanwestra.finemappingtools.gwas.tasks;

import cern.colt.matrix.tbit.BitVector;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Triple;

import java.util.ArrayList;
import java.util.concurrent.Callable;

/**
 * Created by Harm-Jan on 06/23/16.
 */
public class LRTestHaploTask implements Callable<Triple<BitVector[], Integer, Boolean>> {

	private final double genotypeProbThreshold;
	private final ArrayList<VCFVariant> variants;
	private final int i;


	public LRTestHaploTask(int i, ArrayList<VCFVariant> variants, double genotypeProbThreshold) {
		this.i = i;
		this.variants = variants;
		this.genotypeProbThreshold = genotypeProbThreshold;
	}

	@Override
	public Triple<BitVector[], Integer, Boolean> call() throws Exception {
		BitVector haplotype1 = new BitVector(variants.size());
		BitVector haplotype2 = new BitVector(variants.size());
		boolean[] callUncertain = new boolean[variants.size()];
		boolean hasUncertainCalls = false;

		for (int v = 0; v < variants.size(); v++) {
			VCFVariant variant = variants.get(v);
			DoubleMatrix2D probs = variant.getGenotypeProbabilies();
			double[][] alleles = variant.getGenotypeAlleles();

			int sum = (int) alleles[i][0] + (int) alleles[i][1];

			double prob = probs.get(i, sum);
			if (prob < genotypeProbThreshold) {
				callUncertain[v] = true;
				hasUncertainCalls = true;
			}
			if (alleles[i][0] == 1) {
				haplotype1.set(v);
			}
			if (alleles[i][1] == 1) {
				haplotype2.set(v);
			}
		}
		return new Triple<>(new BitVector[]{haplotype1, haplotype2}, i, hasUncertainCalls);
	}
}
