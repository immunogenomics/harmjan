package nl.harmjanwestra.utilities.vcf;

import nl.harmjanwestra.utilities.features.Chromosome;

import java.util.Comparator;

/**
 * Created by Harm-Jan on 03/20/16.
 */
public class VCFVariantComparator implements Comparator<VCFVariant> {


	@Override
	public int compare(VCFVariant var1, VCFVariant var2) {
		Chromosome chr1 = Chromosome.parseChr(var1.getChr());
		Chromosome chr2 = Chromosome.parseChr(var2.getChr());

		if (chr1.equals(chr2)) {
			// compare positions
			if (var1.getPos() > var2.getPos()) {
				return 1;
			} else if (var1.getPos() < var2.getPos()) {
				return -1;
			} else {
				return 0;
			}
		} else {
			return chr1.compare(chr2);
		}
	}
}
