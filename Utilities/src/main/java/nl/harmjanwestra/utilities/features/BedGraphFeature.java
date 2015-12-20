package nl.harmjanwestra.utilities.features;

/**
 * Created by hwestra on 12/3/15.
 */
public class BedGraphFeature extends Feature {

	private double value;

	public BedGraphFeature(Chromosome chr, Integer pos, int i) {
		super(chr, pos, i);
	}

	public void setValue(double v) {
		this.value = v;
	}

	public double getValue() {
		return value;
	}

}
