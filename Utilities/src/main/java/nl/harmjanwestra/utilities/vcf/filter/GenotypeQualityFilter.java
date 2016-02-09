package nl.harmjanwestra.utilities.vcf.filter;

import nl.harmjanwestra.utilities.vcf.VCFVariant;

/**
 * Created by hwestra on 2/8/16.
 */
public class GenotypeQualityFilter implements VCFGenotypeFilter {

	private int minimalGenotypeQual = 30;

	public GenotypeQualityFilter() {
	}

	public GenotypeQualityFilter(int minimalGenotypeQual) {
		this.minimalGenotypeQual = minimalGenotypeQual;
	}

	public void filter(VCFVariant variant) {
		byte[][] alleles = variant.getGenotypeAlleles();
		short[] genotypeQuals = variant.getGenotypeQuals();
		if (genotypeQuals != null) {
			for (int i = 0; i < genotypeQuals.length; i++) {
				int qual = genotypeQuals[i];
				if (qual < minimalGenotypeQual) {
					alleles[0][i] = -1;
					alleles[1][i] = -1;
				}
			}
		} else {
			// should everything be -1?
		}
	}
}
