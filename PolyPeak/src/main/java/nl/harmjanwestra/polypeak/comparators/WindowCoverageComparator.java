package nl.harmjanwestra.polypeak.comparators;

import nl.harmjanwestra.polypeak.containers.Window;

import java.util.Comparator;

/**
 * Created by hwestra on 1/16/15.
 */
public class WindowCoverageComparator implements Comparator<Window> {


	@Override
	public int compare(Window w1, Window w2) {
		if (w1.getMaxCoverage() > w2.getMaxCoverage()) {
			return -1;
		} else if (w1.getMaxCoverage() == w2.getMaxCoverage()) {
			return 0;
		}
		return 1;
	}
}
