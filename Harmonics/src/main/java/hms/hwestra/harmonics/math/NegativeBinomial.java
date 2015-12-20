package hms.hwestra.harmonics.math;

import cern.jet.random.tdouble.engine.DRand;


/**
 * Created by hwestra on 6/28/15.
 * WARNING: THIS CODE IS EXPERIMENTAL AND IS BOUND TO YIELD INCORRECT RESULTS
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

		// get the mean and variance over the data
		double mean = JSci.maths.ArrayMath.mean(data);
		double sd = JSci.maths.ArrayMath.standardDeviation(data);

		double start = mean + 3*sd;

		NewtonRaphson nr = new NewtonRaphson(fx, fprime, start, 200);
		double r = nr.findRoot();
		NegativeBinomialLikelihoodFunction nb = new NegativeBinomialLikelihoodFunction(data);
		double p = nb.determineP(r);
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
