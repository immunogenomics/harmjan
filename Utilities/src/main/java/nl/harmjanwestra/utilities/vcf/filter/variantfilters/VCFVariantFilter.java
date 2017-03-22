package nl.harmjanwestra.utilities.vcf.filter.variantfilters;

import nl.harmjanwestra.utilities.vcf.VCFVariant;

/**
 * Created by hwestra on 3/21/17.
 */
public interface VCFVariantFilter {
	public boolean passesThreshold(VCFVariant variant);
}
