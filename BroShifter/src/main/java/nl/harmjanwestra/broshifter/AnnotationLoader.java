package nl.harmjanwestra.broshifter;


import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.*;
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
		ArrayList<Feature> tmpFeatures = null;
		if (annotation1.endsWith(".xls") || annotation1.endsWith(".xls.gz")) {
			XLSFile xlsFile = new XLSFile();
			ArrayList<PeakFeature> pf = xlsFile.readAllPeaks(annotation1, false, 0.05);
			tmpFeatures = new ArrayList<>();
			tmpFeatures.addAll(pf);
		} else {
			BedFileReader bf = new BedFileReader();
			allAnnotations = new ArrayList<>(10000);
			if (annotation1.endsWith(".csv")) {
				bf = new BedFileReader(Strings.semicolon);
			}
			tmpFeatures = bf.readAsList(annotation1);
		}

		// filter if regions set
		if (regions != null) {
			Track t = new Track("tmp");
			t.addFeatures(tmpFeatures);
			for (Feature f : regions) {
				allAnnotations.addAll(t.getFeatureSet(f));
			}
		} else {
			allAnnotations = tmpFeatures;
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
