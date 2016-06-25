package nl.harmjanwestra.gwas.tasks;

import cern.colt.matrix.tbit.BitVector;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;

import java.util.ArrayList;
import java.util.concurrent.Callable;

/**
 * Created by Harm-Jan on 06/23/16.
 */
public class LRTestHaploTask implements Callable<Pair<BitVector[], boolean[]>> {

	private final double genotypecallqualthreshold;
	private final ArrayList<VCFVariant> variants;
	private final int i;

	public LRTestHaploTask(int i, ArrayList<VCFVariant> variants, double genotypecallqualthreshold) {
		this.i = i;
		this.variants = variants;
		this.genotypecallqualthreshold = genotypecallqualthreshold;
	}

	@Override
	public Pair<BitVector[], boolean[]> call() throws Exception {
		BitVector haplotype1 = new BitVector(variants.size());
		BitVector haplotype2 = new BitVector(variants.size());
		boolean[] callUncertain = new boolean[variants.size()];

		for (int v = 0; v < variants.size(); v++) {
			VCFVariant variant = variants.get(v);
			DoubleMatrix2D probs = variant.getGenotypeProbabilies();
			double[][] alleles = variant.getGenotypeAlleles();

			int sum = (int) alleles[i][0] + (int) alleles[i][1];

			double prob = probs.get(i, sum);
			if (prob < genotypecallqualthreshold) {
				callUncertain[v] = true;
			}
			if (alleles[i][0] == 1) {
				haplotype1.set(v);
			}
			if (alleles[i][1] == 1) {
				haplotype2.set(v);
			}
		}
		return new Pair<>(new BitVector[]{haplotype1, haplotype2}, callUncertain);
	}
}
