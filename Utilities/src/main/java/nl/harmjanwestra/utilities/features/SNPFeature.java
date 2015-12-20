package nl.harmjanwestra.utilities.features;

/**
 * Created by hwestra on 11/11/15.
 */
public class SNPFeature extends Feature {
	double p;

	public SNPFeature(){

	}

	public SNPFeature(SNPFeature f2) {
		super(f2);
		this.p = f2.getP();
	}

	public void setP(double p) {
		this.p = p;
	}

	public double getP() {
		return p;
	}
}
