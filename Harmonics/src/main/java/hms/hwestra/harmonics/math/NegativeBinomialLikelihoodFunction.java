package hms.hwestra.harmonics.math;


import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.CombinatoricsUtils;

/**
 * Created by hwestra on 6/28/15.
 * Formulae graciously taken from https://en.wikipedia.org/wiki/Negative_binomial_distribution#Maximum_likelihood_estimation
 * * WARNING: THIS CODE IS EXPERIMENTAL AND IS BOUND TO YIELD INCORRECT RESULTS
 */

public class NegativeBinomialLikelihoodFunction implements Function {


	private final int[] data;

	private final int N;

	public NegativeBinomialLikelihoodFunction(int[] data) {
		this.data = data;
		this.N = data.length;
	}

	@Override
	public double calc(double r) {

		// given r, calculate p first

		double p = determineP(r);
		double logp = Math.log(p);
		double log1minp = Math.log(1 - p);
		double suma = 0;
		double sumb = 0;
		double sumc = 0;
		for (int i = 0; i < N; i++) {
			suma += Gamma.logGamma(data[i] + r);
			sumb += CombinatoricsUtils.factorialLog(data[i]);
			sumc += (data[i] * logp);
		}

		double logLikelihood = suma - sumb - (N * Gamma.logGamma(r)) + sumc + (N * r * log1minp);
		System.out.println("ll NB: " + logLikelihood);
		return logLikelihood;
	}

	public double determineP(double r) {

		double sumK = 0;
		for (int i = 0; i < N; i++) {
			sumK += data[i];
		}
		double p = sumK / ((N * r) + sumK);
		System.out.println("r: " + r + "\tp: " + p);
		return p;
	}


}
