package nl.harmjanwestra.utilities.vcf.filter;

import nl.harmjanwestra.utilities.vcf.VCFVariant;

/**
 * Created by hwestra on 2/8/16.
 */
public class ReadDepthFilter implements VCFGenotypeFilter {


	private int minimalReadDepth = 10;

	public ReadDepthFilter() {

	}

	public ReadDepthFilter(int minimalReadDepth) {
		this.minimalReadDepth = minimalReadDepth;
	}

	public void filter(VCFVariant variant) {
		byte[][] alleles = variant.getGenotypeAlleles();
		short[][] allelicDepth = variant.getAllelicDepth();
		short[] approximateDepth = variant.getApproximateDepth();

		if (approximateDepth != null) {
			if (allelicDepth != null) {
				// account for a bug in newer GATK output
				for (int i = 0; i < alleles[0].length; i++) {
					short indSum = 0;
					for (int j = 0; j < allelicDepth.length; j++) {
						indSum += allelicDepth[j][i];
					}
					if (approximateDepth[i] != indSum) {
						approximateDepth[i] = indSum;
					}
				}
			}

			for (int i = 0; i < alleles[0].length; i++) {
				int depth = approximateDepth[i];
				if (depth < minimalReadDepth) {
					alleles[0][i] = -1;
					alleles[1][i] = -1;
				}
			}
		} else {
			// should everything be -1 now???
		}
	}
}
