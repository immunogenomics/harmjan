package hms.hwestra.harmonics.math;

import cern.jet.random.tdouble.engine.DRand;


/**
 * Created by hwestra on 6/28/15.
 */
public class NegativeBinomial {


	private final int[] data;

	public NegativeBinomial(int[] data) {
		this.data = data;
	}

	public double[] estimateParameters() {

		// determine R using Newton-Raphson
		NegativeBinomialLikelihoodFunction fx = new NegativeBinomialLikelihoodFunction(data);
		NegativeBinomialRDerivativeFunction fprime = new NegativeBinomialRDerivativeFunction(data);

		NewtonRaphson nr = new NewtonRaphson(fx, fprime, 0.5, 1000);
		double r = nr.findRoot();
		double p = fx.determineP(r);
		return new double[]{r, p};
	}

	public cern.jet.random.tdouble.NegativeBinomial getDist() {
		double[] params = estimateParameters();
		return new cern.jet.random.tdouble.NegativeBinomial(
				(int) Math.ceil(params[0]),
				params[1],
				new DRand());
	}


}
