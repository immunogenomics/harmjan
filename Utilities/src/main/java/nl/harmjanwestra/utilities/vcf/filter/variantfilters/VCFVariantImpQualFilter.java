package nl.harmjanwestra.utilities.vcf.filter.variantfilters;

import nl.harmjanwestra.utilities.vcf.VCFVariant;

/**
 * Created by hwestra on 3/21/17.
 */
public class VCFVariantImpQualFilter implements VCFVariantFilter {
	private double threshold = 0.3;

	public VCFVariantImpQualFilter(double threshold) {
		this.threshold = threshold;
	}

	@Override
	public boolean passesThreshold(VCFVariant variant) {
		Double impqual = variant.getImputationQualityScore();
		if (impqual != null) {
			return impqual > threshold;
		} else {
			return true;
		}

	}

	@Override
	public String toString() {
		return "VCFVariantImpQualFilter{" +
				"threshold=" + threshold +
				'}';
	}
}
