package nl.harmjanwestra.utilities.bedfile;

import nl.harmjanwestra.utilities.features.Feature;

/**
 * Created by hwestra on 1/6/15.
 */
public class BedFileFeatureFilter {

	enum MODE {INCLUSIVE, EXCLUSIVE}

	;

	Feature[] features;
	MODE mode = MODE.INCLUSIVE;

	public BedFileFeatureFilter(Feature[] featuresToFilter, MODE mode) {
		this.mode = mode;
		this.features = featuresToFilter;
	}


	// does the feature f overlap with any feature in Feature[] if MODE==inclusive?
	public boolean passesFilter(Feature toCheck) {
		boolean overlap = false;
		for (Feature f : features) {

			if(f.overlaps(toCheck)){
				overlap = true;
				break;
			}

		}

		if(mode == MODE.INCLUSIVE){
			return overlap;
		} else {
			return !overlap;
		}
	}

}
