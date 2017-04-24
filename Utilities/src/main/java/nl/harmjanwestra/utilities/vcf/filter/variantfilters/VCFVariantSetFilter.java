package nl.harmjanwestra.utilities.vcf.filter.variantfilters;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.vcf.VCFVariant;

import java.io.IOException;
import java.util.HashSet;

/**
 * Created by hwestra on 4/18/17.
 */
public class VCFVariantSetFilter implements VCFVariantFilter {

	private HashSet<String> toFilter = null;

	public VCFVariantSetFilter(HashSet<String> s) {
		this.toFilter = s;
	}

	public VCFVariantSetFilter(String file) throws IOException {
		TextFile tf = new TextFile(file, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		toFilter = new HashSet<>();
		while (elems != null) {
			String snp = elems[0];
			String[] snpelems = snp.split("_");
			if (snpelems.length > 1) {
				Chromosome chr = Chromosome.parseChr(snpelems[0]);
				String poselems = snpelems[1];
				String str = chr.toString() + "_" + poselems;
				toFilter.add(str);
			} else {
				System.err.println("Error loading: " + file + " expected format: Chrx_pos_rsid");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
	}

	@Override
	public boolean passesThreshold(VCFVariant variant) {
		String str = variant.getChrObj().toString() + "_" + variant.getPos();
		return toFilter.contains(str);
	}

	@Override
	public String toString() {
		return "VCFVariantSetFilter{" +
				"toFilter=" + toFilter.size() +
				'}';
	}
}
