package nl.harmjanwestra.finemapping.rebuttal;

public class KgVariantPair implements Comparable<KgVariantPair> {
	
	KgVariant v1;
	KgVariant v2;
	
	public KgVariantPair(KgVariant v1, KgVariant v2) {
		this.v1 = v1;
		this.v2 = v2;
	}
	
	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;
		
		KgVariantPair that = (KgVariantPair) o;
		if (this.v1.equals(v1) && this.v2.equals(that.v2) ||
				this.v2.equals(v1) && this.v1.equals(that.v2)) {
			return true;
		} else {
			return false;
		}
	}
	
	@Override
	public int hashCode() {
		int result = v1 != null ? v1.hashCode() : 0;
		result = 31 * result + (v2 != null ? v2.hashCode() : 0);
		return result;
	}
	
	@Override
	public int compareTo(KgVariantPair o) {
		if (this.equals(o)) {
			return 0;
		} else {
			return this.v1.compareTo(o.v1);
		}
	}
}
