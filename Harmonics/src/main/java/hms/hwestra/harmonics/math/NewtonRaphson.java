package hms.hwestra.harmonics.math;

/**
 * Created by hwestra on 6/28/15.
 * Freely implemented from: https://en.wikipedia.org/wiki/Newton%27s_method#Pseudocode
 */
public class NewtonRaphson {

	private double x0 = 1;
	private Function fx = null;
	private Function fprime = null;
	private double tolerance = 1E-17;
	private double epsilon = 10E-14;
	private int maxIterations = 100;

	public NewtonRaphson(Function fx, Function fxprime) {
		this.fx = fx;
		this.fprime = fxprime;
	}

	public NewtonRaphson(Function fx, Function fxprime, double x0) {
		this(fx, fxprime);
		this.x0 = x0;
	}

	public NewtonRaphson(Function fx, Function fxprime, double x0, int maxIterations) {
		this(fx, fxprime, x0);
		this.maxIterations = maxIterations;
	}

	double findRoot() {
		int niter = 0;
		while (niter < maxIterations) {


			double y = fx.calc(x0);
			double yprime = fprime.calc(x0);


			System.out.println("nr step1: " + niter + "\ty: " + y + "\typrime: " + yprime);

			if (Math.abs(yprime) < epsilon) {
				// already at minimum
				return x0;
			}

			double x1 = x0 - ((y / yprime));

			System.out.println("nr step2: " + niter + "\tx0: " + x0 + "\tx1: " + x1 + "\t" + (Math.abs(x1 - x0) / Math.abs(x1)));

			if ((Math.abs(x1 - x0) / Math.abs(x1)) < tolerance) {
				// found new minimum;
				return x0;
			}




			x0 = x1;

			niter++;
		}

		// we did not converge
		return Double.NaN;
	}
}
