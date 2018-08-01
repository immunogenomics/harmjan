package nl.harmjanwestra.utilities.sets;

import nl.harmjanwestra.utilities.features.Feature;

import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 11/30/15.
 */
public class FeatureSets {


	public static HashSet<Feature> intersect(ArrayList<Feature> f1, ArrayList<Feature> f2) {
		HashSet<Feature> shared = new HashSet<Feature>();
		for (Feature feat1 : f1) {
			for (Feature feat2 : f2) {
				if (feat1.overlaps(feat2)) {
					shared.add(feat1);
				}
			}
		}
		return shared;

	}
}
