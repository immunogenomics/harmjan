package nl.harmjanwestra.utilities.vcf.filter;

import nl.harmjanwestra.utilities.vcf.VCFVariant;

/**
 * Created by hwestra on 2/8/16.
 */
public class AllelicDepthFilter implements VCFGenotypeFilter {

	private int minimumReadDepth;
	private double abCutoff = 0.15;

	public AllelicDepthFilter() {

	}

	public AllelicDepthFilter(double allelicBalanceCutoff, int minimumReadDepth) {
		this.abCutoff = allelicBalanceCutoff;
		this.minimumReadDepth = minimumReadDepth;
	}

	public void filter(VCFVariant variant) {

		short[][] allelicDepth = variant.getAllelicDepth();
		if (allelicDepth != null) {
			byte[][] alleles = variant.getGenotypeAlleles();
			for (int i = 0; i < alleles[0].length; i++) {
				if (alleles[0][i] != -1) {
					int sum = 0;
					for (int j = 0; j < allelicDepth.length; j++) {
						sum += allelicDepth[j][i];
					}
					if (alleles[0][i] != alleles[1][i]) {
						double ab1 = (double) allelicDepth[alleles[0][i]][i] / sum;
						double ab2 = (double) allelicDepth[alleles[1][i]][i] / sum;
						boolean ok = true;
						if (ab1 < abCutoff || ab1 > (1 - abCutoff)) {
							ok = false;
						}
						if (ab2 < abCutoff || ab2 > (1 - abCutoff)) {
							ok = false;
						}

						if (!ok || sum < minimumReadDepth) {
							alleles[0][i] = -1;
							alleles[1][i] = -1;
						}
					} else {

						if (allelicDepth[alleles[0][i]][i] < minimumReadDepth) {
							alleles[0][i] = -1;
							alleles[1][i] = -1;
						}

					}
				}
			}
		}
	}
}
