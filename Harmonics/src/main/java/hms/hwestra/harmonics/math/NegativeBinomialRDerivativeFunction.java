package hms.hwestra.harmonics.math;

import JSci.maths.ArrayMath;
import org.apache.commons.math3.special.Gamma;

/**
 * Created by hwestra on 6/28/15.
 * Formulae graciously taken from https://en.wikipedia.org/wiki/Negative_binomial_distribution#Maximum_likelihood_estimation
 */
public class NegativeBinomialRDerivativeFunction implements Function {


	private final int N ;
	private final int[] data;

	public NegativeBinomialRDerivativeFunction(int[] data) {
		this.data = data;
		this.N = data.length;
	}

	@Override
	public double calc(double r) {

		double sum = 0;
		double meanK = ArrayMath.mean(data);
		for (int i = 0; i < N; i++) {
			sum += Gamma.digamma(data[i] + r)
					- (N * Gamma.digamma(r))
					+ (N * Math.log(r / r + meanK));
		}

		return sum;
	}


}
