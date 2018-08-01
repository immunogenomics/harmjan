package nl.harmjanwestra.polypeak.containers;

import nl.harmjanwestra.utilities.enums.Chromosome;

/**
 * Created by hwestra on 1/7/15.
 */
public class Sequence {

	private final String name;
	private final Integer length;

	// this class provides a synomym list for the different HG19 builds..
	public Sequence(String name) {
		Chromosome chr = Chromosome.parseChr(name);
		if (chr.equals(Chromosome.NA)) {
			this.name = name;
		} else {
			this.name = chr.getName();
		}
		this.length = 0;
	}

	public Sequence(String name, Integer length) {
		Chromosome chr = Chromosome.parseChr(name);
		if (chr.equals(Chromosome.NA)) {
			this.name = name;
		} else {
			this.name = chr.getName();
		}
		this.length = length;
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (!(o instanceof Sequence)) return false;

		Sequence sequence = (Sequence) o;

		if (!name.equals(sequence.name)) return false;

		return true;
	}

	@Override
	public int hashCode() {
		return name.hashCode();
	}

	public String getName() {
		return name;
	}


	private boolean compareSequenceName(Sequence s) {
		String othername = s.getName();
		if (name == null && othername == null) {
			return true;
		}
		if (name != null && othername != null) {
			Chromosome chr = Chromosome.parseChr(name);
			Chromosome otherChr = Chromosome.parseChr(othername);
			if (chr.equals(Chromosome.NA) && otherChr.equals(Chromosome.NA)) {
				// do string compare
				return this.name.equals(othername);
			} else {
				return chr.equals(otherChr);
			}
		}
		return false;
	}

	public Integer getLength() {
		return length;
	}


	public boolean compareLength(Sequence other) {
		if (this.length == null && other.getLength() == null) {
			return true;
		} else if (this.length != null && other.getLength() != null) {
			return this.getLength().equals(other.getLength());
		} else {
			return false;
		}
	}
}
