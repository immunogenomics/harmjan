package nl.harmjanwestra.utilities.vcf.filter.variantfilters;

/**
 * Created by hwestra on 6/21/17.
 */
public class VCFVariantMissingnessPFilter {

	private double thresh = 5E-8;

	public VCFVariantMissingnessPFilter(double threshold) {
		this.thresh = threshold;
	}
}
