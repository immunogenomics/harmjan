package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.vcf.VCFVariantLoader;
import nl.harmjanwestra.utilities.vcf.filter.variantfilters.VCFVariantFilter;

import java.util.ArrayList;

/**
 * Created by hwestra on 5/1/17.
 */
public class FindIndividualsWithMissingData {

	public static void main(String[] args) {
		String vcf = "";


		ArrayList<VCFVariantFilter> filters = new ArrayList<VCFVariantFilter>();


		VCFVariantLoader loader = new VCFVariantLoader();
//		loader.run(vcf, filters);

	}

}
