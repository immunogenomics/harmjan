package nl.harmjanwestra.harmonics.smoothing;

/**
 * Created by hwestra on 6/25/15.
 */
public class AverageSmoothingFunction implements SmoothingFunction {


	@Override
	public double kernel(int[] data) {
		double sum = 0;
		int vals = data.length;
		for (int i = 0; i < data.length; i++) {
			sum += data[i];
		}
		return sum / vals;
	}

	@Override
	public double kernel(double[] data) {
		double sum = 0;
		int vals = data.length;
		for (int i = 0; i < data.length; i++) {
			sum += data[i];
		}
		return sum / vals;
	}

}
