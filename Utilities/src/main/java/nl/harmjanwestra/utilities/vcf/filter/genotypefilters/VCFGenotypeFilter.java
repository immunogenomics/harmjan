package nl.harmjanwestra.utilities.vcf.filter.genotypefilters;

import nl.harmjanwestra.utilities.vcf.VCFVariant;

/**
 * Created by hwestra on 2/8/16.
 */
public interface VCFGenotypeFilter {
	void filter(VCFVariant variant);
}
