package hms.hwestra.harmonics.math;

import org.apache.commons.math3.special.Gamma;

/**
 * Created by hwestra on 6/30/15.
 * * WARNING: THIS CODE IS EXPERIMENTAL AND IS BOUND TO YIELD INCORRECT RESULTS
 */
public class NegativeBinomialRSecondDerivativeFunction implements Function {
	private final int[] data;

	private final int N;

	public NegativeBinomialRSecondDerivativeFunction(int[] data) {
		this.data = data;
		this.N = data.length;
	}


	@Override
	public double calc(double x) {

		double suma = 0;
		double sumb = 0;

		for (int i = 0; i < data.length; i++) {
			suma += Gamma.trigamma(data[i] + x);
			sumb += ((double) data[i]);
		}


		double ll = suma - (N * Gamma.trigamma(x)) + (N * 1 / (x / (x + (sumb/N))));
		return ll;

	}
}
