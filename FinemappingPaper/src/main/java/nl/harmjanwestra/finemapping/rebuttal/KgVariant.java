package nl.harmjanwestra.finemapping.rebuttal;

import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.SNPFeature;

public class KgVariant implements Comparable<KgVariant> {
	public SNPFeature f;
	public int nproxies;
	public double maf;
	public double hwep;
	public double score;
	
	@Override
	public int compareTo(KgVariant o) {
		FeatureComparator c = new FeatureComparator();
		return c.compare(this.f, o.f);
	}
	
	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;
		
		KgVariant kgVariant = (KgVariant) o;
		
		return f != null ? f.equals(kgVariant.f) : kgVariant.f == null;
	}
	
	@Override
	public int hashCode() {
		return f != null ? f.hashCode() : 0;
	}
	
	public int getMafBin() {
		return (int) Math.floor(100 * maf * 2);
	}
}


