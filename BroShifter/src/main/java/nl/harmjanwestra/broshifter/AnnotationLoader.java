package nl.harmjanwestra.broshifter;


import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureMerger;
import nl.harmjanwestra.utilities.features.PeakFeature;
import nl.harmjanwestra.utilities.features.Track;
import nl.harmjanwestra.utilities.peakfiles.XLSFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 12/3/15.
 */
public class AnnotationLoader {

	public Track loadAnnotations(String annotation1,
								 boolean usePeakCenter,
								 int bpToExtendAnnotation,
								 boolean mergeoverlapping,
								 ArrayList<Feature> regions) throws IOException {

		ArrayList<Feature> allAnnotations = null;
		if (annotation1.endsWith(".xls")) {
			XLSFile xlsFile = new XLSFile();
			ArrayList<PeakFeature> pf = xlsFile.readAllPeaks(annotation1, false, 0.05);
			allAnnotations = new ArrayList<>(10000);
			for (Feature p : pf) {
				if (regions == null) {
					allAnnotations.add(p);
				} else {
					boolean overlap = false;
					for (Feature f : regions) {
						if (f.overlaps(p)) {
							overlap = true;
						}
					}
					if (overlap) {
						allAnnotations.add(p);
					}
				}
			}
		} else {
			BedFileReader bf = new BedFileReader();
			allAnnotations = new ArrayList<>(10000);
			if (annotation1.endsWith(".csv")) {
				bf = new BedFileReader(Strings.semicolon);
			}

			Track t = bf.readAsTrack(annotation1, annotation1);
			for (Feature f : regions) {
				allAnnotations.addAll(t.getFeatureSet(f));
			}
		}

		if (usePeakCenter) {
			updateAnnotationsExtendAroundPeakCenter(allAnnotations, bpToExtendAnnotation);
		}

		// merge any annotation that may overlap within this dataset
		if (mergeoverlapping) {
			allAnnotations = FeatureMerger.merge(allAnnotations, true);
		}
		Track t = new Track(annotation1);
		t.addFeatures(allAnnotations);
		return t;
	}

	private void updateAnnotationsExtendAroundPeakCenter(ArrayList<Feature> allAnnotations, int bpToExtend) {

		for (int i = 0; i < allAnnotations.size(); i++) {
			Feature p = allAnnotations.get(i);
			int start = p.getStart();
			int stop = p.getStop();
			int mid = start + ((stop - start) / 2);
			if (p instanceof PeakFeature) {
				// this thing has peak center info.. sweet!
				PeakFeature pf = (PeakFeature) p;
				int summit = pf.getSummit();
				if (summit > 0) {
					mid = summit;
				}
			}
			start = mid - (bpToExtend / 2);
			stop = mid + (bpToExtend / 2);
			p.setStart(start);
			p.setStop(stop);

		}
	}
}
