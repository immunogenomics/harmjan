package nl.harmjanwestra.utilities.vcf.filter.variantfilters;

import nl.harmjanwestra.utilities.vcf.VCFVariant;

import java.util.ArrayList;

/**
 * Created by hwestra on 3/21/17.
 */
public class VCFVariantFilters {

	ArrayList<VCFVariantFilter> filters;

	public VCFVariantFilters() {
		filters = new ArrayList<VCFVariantFilter>();
	}

	public void addFilter(VCFVariantFilter filter) {
		filters.add(filter);
	}

	public boolean passesFilters(VCFVariant v) {
		for (VCFVariantFilter f : filters) {
			if (!f.passesThreshold(v)) {
				return false;
			}
		}
		return true;
	}

	public String toString() {
		String output = "VCF Variant filters:\n";
		for (VCFVariantFilter filter : filters) {
			output += filter.toString() + "\n";
		}
		return output;
	}

	public int size() {
		return filters.size();
	}
}

