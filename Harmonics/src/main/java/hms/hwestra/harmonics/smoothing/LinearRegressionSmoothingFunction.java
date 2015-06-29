package hms.hwestra.harmonics.smoothing;

import umcg.genetica.math.stats.Regression;

/**
 * Created by hwestra on 6/25/15.
 */
public class LinearRegressionSmoothingFunction implements SmoothingFunction {


	@Override
	public double kernel(int[] data) {
		double[] x = new double[data.length];
		double[] y = new double[data.length];
		for (int i = 0; i < x.length; i++) {
			x[i] = i;
			y[i] = data[i];
		}

		double[] coeff = Regression.getLinearRegressionCoefficients(x, y);
		double output = coeff[1] + (coeff[0] * data[data.length / 2]);
		return output;
	}
}
