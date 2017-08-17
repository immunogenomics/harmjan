package nl.harmjanwestra.utilities.vcf.filter.variantfilters;

import nl.harmjanwestra.utilities.vcf.VCFVariant;

/**
 * Created by hwestra on 6/21/17.
 */
public class VCFVariantMissingnessPFilter implements VCFVariantFilter {

	private double thresh = 5E-8;

	public VCFVariantMissingnessPFilter(double threshold) {
		this.thresh = threshold;
	}
	
	@Override
	public boolean passesThreshold(VCFVariant variant) {

		if(variant.getDiffMissingnessP() < thresh){
			return false;
		} else {
			return true;
		}
	}
	
	@Override
	public String toString() {
		return "VCFVariantMissingnessPFilter{" +
				"threshold=" + thresh +
				'}';
	}
}
