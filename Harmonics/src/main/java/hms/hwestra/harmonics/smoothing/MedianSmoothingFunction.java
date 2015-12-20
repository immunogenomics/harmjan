package hms.hwestra.harmonics.smoothing;

import java.util.Arrays;

/**
 * Created by hwestra on 6/25/15.
 */
public class MedianSmoothingFunction implements SmoothingFunction {
	@Override
	public double kernel(int[] data) {
		Arrays.sort(data);
		if (data.length < 2) {
			return data[0];
		}
		if (data.length % 2 == 0) {
			// even number
			int mid = data.length / 2;
			double s = ((double) data[mid] + data[mid - 1]) / 2;
			return s;
		} else {
			return data[(int) Math.floor((double) data.length / 2)];
		}

	}

	@Override
	public double kernel(double[] data) {
		Arrays.sort(data);
		if (data.length < 2) {
			return data[0];
		}
		if (data.length % 2 == 0) {
			// even number
			int mid = data.length / 2;
			double s = ((double) data[mid] + data[mid - 1]) / 2;
			return s;
		} else {
			return data[(int) Math.floor((double) data.length / 2)];
		}

	}
}
