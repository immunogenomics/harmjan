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
		checkRegionFilter();
	}

	public void add(VCFVariantFilter filter) {
		addFilter(filter);
	}

	public void addFilter(VCFVariantFilter filter) {
		filters.add(filter);
		checkRegionFilter();
	}

	private void checkRegionFilter() {
		if (hasRegionOrSetFilter == null) {
			for (VCFVariantFilter f : filters) {
				if (f instanceof VCFVariantSetFilter || f instanceof VCFVariantRegionFilter) {
					hasRegionOrSetFilter = true;
				}
			}
		}
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

	public boolean passesRegionOrVariantFilter(VCFVariant v) {
		boolean passesfilter = true;
		for (VCFVariantFilter f : filters) {
			if (f instanceof VCFVariantSetFilter) {
				passesfilter = f.passesThreshold(v);
			}
		}
		return passesfilter;

	}

	public Boolean hasRegionOrSetFilter = null;

	public boolean hasRegionOrVariantSetFilter() {

		return hasRegionOrSetFilter;
	}
}

