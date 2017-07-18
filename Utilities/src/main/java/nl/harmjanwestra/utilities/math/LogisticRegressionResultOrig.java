package nl.harmjanwestra.utilities.math;

/**
 * Created by hwestra on 11/10/15.
 */
public class LogisticRegressionResultOrig {


	private final double[] beta;
	private final double[] stderrs;
	private final double[] sigprms;
	private final double deviance;
	private final double[] loglike;
	private final double loglike0;

	public LogisticRegressionResultOrig(double[] beta, double[] stderrs, double[] sigprms, double deviance, double[] loglike, double loglike0) {
		this.beta = beta;
		this.stderrs = stderrs;
		this.sigprms = sigprms;
		this.deviance = deviance;
		this.loglike = loglike;
		this.loglike0 = loglike0;
	}

	public double[] getBeta() {
		return beta;
	}

	public double[] getStderrs() {
		return stderrs;
	}

	public double[] getSigprms() {
		return sigprms;
	}

	public double getDeviance() {
		return deviance;
	}

	public double[] getLoglike() {
		return loglike;
	}

	public double getLoglike0() {
		return loglike0;
	}
}
