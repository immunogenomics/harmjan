package hms.hwestra.harmonics.smoothing;

/**
 * Created by hwestra on 6/28/15.
 */
public interface SmoothingFunction {

	double kernel(int[] data);

}
