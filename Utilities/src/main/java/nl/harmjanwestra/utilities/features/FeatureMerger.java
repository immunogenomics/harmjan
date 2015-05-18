package nl.harmjanwestra.utilities.features;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by hwestra on 5/13/15.
 */
public class FeatureMerger {

	public static ArrayList<Feature> merge(ArrayList<Feature> f) {

		Feature[] featureArr = new Feature[f.size()];
		for (int i = 0; i < f.size(); i++) {
			featureArr[i] = f.get(i);
		}
		Arrays.sort(featureArr,new FeatureComparator(false));

		System.out.println(featureArr.length);
		int overlap = 0;
		while (overlap != 0) {
			overlap = 0;
			for (int i = 0; i < featureArr.length; i++) {
				Feature f1 = featureArr[i];
				if (f1 != null) {
					for (int j = i + 1; j < featureArr.length; j++) {
						Feature f2 = featureArr[j];
						if (f2 != null) {
							if(!f1.equals(f2.getChromosome())){
								break;// list is sorted: no need to continue merging
							}
							if (f1.overlaps(f2)) {
								if (f2.getStart() < f1.getStart()) {
									f1.setStart(f2.getStart());
								}
								if (f2.getStop() > f1.getStop()) {
									f1.setStop(f2.getStop());
								}

								featureArr[j] = null;
								overlap++;
							}
						}

					}
				}
			}

		}

		ArrayList<Feature> output = new ArrayList<Feature>();
		for (int i = 0; i < featureArr.length; i++) {
			if (featureArr[i] != null) {
				output.add(featureArr[i]);
			}
		}
		System.out.println(output.size());
		return output;
	}
}
