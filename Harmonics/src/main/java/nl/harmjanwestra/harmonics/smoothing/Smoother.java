package nl.harmjanwestra.harmonics.smoothing;

/**
 * Created by hwestra on 6/25/15.
 */
public class Smoother {

	private final SmoothingFunction smoothingFunction;

	public Smoother(SmoothingFunction f) {
		this.smoothingFunction = f;
	}

	public int[][] smoothMatrix(int[][] data, int smoothingwindow) {
		int[][] output = new int[data.length][0];
		for (int i = 0; i < data.length; i++) {
			output[i] = this.smooth(output[i], smoothingwindow);
		}
		return output;
	}


	public int[] smooth(int[] data, int smoothingwindow) {
		int[] output = new int[data.length];
		int half = smoothingwindow / 2;
		for (int i = 0; i < output.length; i++) {

			// put the current i in the middle

			int start = i - half;
			int stop = i + half;
			if (start < 0) {
				start = 0;
				stop = i + smoothingwindow;
			}
			if (stop >= output.length) {
				stop = output.length - 1;
				start = stop - smoothingwindow;
				if (start < 0) {
					start = 0;
				}
			}


			int[] tmpData = new int[stop - start];
			System.arraycopy(data, start, tmpData, 0, tmpData.length);

			output[i] = (int) Math.ceil(smoothingFunction.kernel(tmpData));

		}
		return output;
	}

	public int[][] reduceMatrix(int[][] data, int newResolution) {
		int[][] output = new int[data.length][0];
		for (int i = 0; i < data.length; i++) {
			output[i] = this.reduce(output[i], newResolution);
		}
		return output;
	}


	public int[] reduce(int[] data, int newResolution) {
		int smoothingWindow = data.length / newResolution;
		int[] output = new int[newResolution];

		int windowCtr = 0;
		for (int i = 0; i < data.length; i += smoothingWindow) {
			int start = i;
			int stop = i + smoothingWindow;
			if (stop >= data.length) {
				stop = data.length - 1;
			}

			int[] tmpData = new int[stop - start];
			System.arraycopy(data, start, tmpData, 0, tmpData.length);
			output[windowCtr] = (int) Math.ceil(smoothingFunction.kernel(tmpData));
			windowCtr++;
		}
		return output;
	}


	public double[] reduce(double[] data, int newResolution) {
		double smoothingWindowSize = (double) data.length / newResolution;

		double[] output = new double[newResolution];

		System.out.println(data.length + "\t" + newResolution + "\t" + smoothingWindowSize);

		int windowCtr = 0;
		for (int window = 0; window < newResolution; window++) {


			int start = (int) Math.floor(window * smoothingWindowSize);
			int stop = (int) Math.floor(start + smoothingWindowSize);
			if (stop >= data.length) {
				stop = data.length - 1;
			}
			System.out.println("window " + windowCtr + "/" + newResolution + "\t" + window + "\t" + start + "\t" + stop);

			double[] tmpData = new double[stop - start];
			System.arraycopy(data, start, tmpData, 0, tmpData.length);
			output[windowCtr] = Math.ceil(smoothingFunction.kernel(tmpData));
			windowCtr++;
		}
		return output;
	}
}
