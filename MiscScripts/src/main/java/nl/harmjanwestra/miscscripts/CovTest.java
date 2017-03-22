package nl.harmjanwestra.miscscripts;

/**
 * Created by hwestra on 3/20/17.
 */
public class CovTest {

	public static void main(String[] args) {

		double[] y = new double[10];

		for (int i = 0; i < y.length; i++) {
			y[i] = Math.random();
		}

		double[][] cov = new double[y.length][y.length];

		double avg = JSci.maths.ArrayMath.mean(y);

		for (int i = 0; i < y.length; i++) {
			for (int j = 0; j < y.length; j++) {

			}
		}

	}
}
