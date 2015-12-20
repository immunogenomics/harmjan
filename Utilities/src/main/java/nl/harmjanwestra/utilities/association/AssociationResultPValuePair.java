package nl.harmjanwestra.utilities.association;

import java.util.Objects;

/**
 * Created by hwestra on 11/11/15.
 */
public class AssociationResultPValuePair implements Comparable<AssociationResultPValuePair> {

	double p = 0;
	AssociationResult r;
	boolean reverseSort;

	public AssociationResultPValuePair(AssociationResult r, double p, boolean reverseSort) {
		this.p = p;
		this.r = r;
		this.reverseSort = reverseSort;
	}

	@Override
	public int compareTo(AssociationResultPValuePair o) {
		if (this.equals(o)) {
			return 0;
		}
		if (reverseSort) {
			if (this.p > o.p) {
				return 1;
			} else {
				return -1;
			}
		} else {
			if (this.p > o.p) {
				return -1;
			} else {
				return 1;
			}
		}

	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof AssociationResultPValuePair)) return false;
		return ((AssociationResultPValuePair) o).p == this.p;
	}

	@Override
	public int hashCode() {
		return Objects.hash(p);
	}


	public double getP() {
		return p;
	}

	public AssociationResult getAssociationResult() {
		return r;
	}


}
