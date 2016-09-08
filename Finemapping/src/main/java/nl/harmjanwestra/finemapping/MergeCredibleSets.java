package nl.harmjanwestra.finemapping;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.AssociationResultPValuePair;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Created by hwestra on 9/7/16.
 */
public class MergeCredibleSets {

	public void run(String bedregions, String[] assocFiles, String outfile) throws IOException {


		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedregions);

		// region variant1 pval1 posterior1 variant2 pval2 posterior2 variant3 pval3 posterior3
		AssociationFile f = new AssociationFile();


		// [dataset][region][variants]
		AssociationResult[][][] crediblesets = new AssociationResult[assocFiles.length][regions.size()][];
		AssociationResult[][][] data = new AssociationResult[assocFiles.length][regions.size()][];

		ArrayList<ArrayList<AssociationResult>> results = new ArrayList<>();
		ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
		ArrayList<Feature> regionsWithCredibleSets = new ArrayList<>();
		for (int d = 0; d < regions.size(); d++) {
			boolean hasSet = false;
			for (int i = 0; i < assocFiles.length; i++) {

				ArrayList<AssociationResult> allDatasetData = f.read(assocFiles[i], regions.get(d));
				data[i][d] = allDatasetData.toArray(new AssociationResult[0]);
				ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(results.get(0), 0.9);
				crediblesets[i][d] = credibleSet.toArray(new AssociationResult[0]);
				if (credibleSet.size() <= 5) {
					hasSet = true;
				}
			}
			if (hasSet) {
				regionsWithCredibleSets.add(regions.get(d));
			}
		}


		for (int d = 0; d < regions.size(); d++) {
			for (int i = 0; i < assocFiles.length; i++) {

			}
		}

		TextFile out = new TextFile(outfile, TextFile.R);
		for (int d = 0; d < regions.size(); d++) {
			int len = crediblesets[0][d].length;
			Feature region = regions.get(d);
			if (len <= 5) {
				// print

				String ln = "";
				for (int i = 0; i < data.length; i++) {
					ArrayList<AssociationResult> topResults = getTopVariants(data, i, d, len);

				}

			}
		}
		out.close();

	}

	private ArrayList<AssociationResult> getTopVariants(AssociationResult[][][] data, int i, int d, int len) {
		ArrayList<AssociationResultPValuePair> pairs = new ArrayList<AssociationResultPValuePair>();

		for (AssociationResult r : data[i][d]) {
			AssociationResultPValuePair p = new AssociationResultPValuePair(r, r.getPosterior(), false);
			if (!Double.isInfinite(p.getP()) && !Double.isNaN(p.getP())) {
				pairs.add(p);
			}
		}

		ArrayList<AssociationResult> credibleSet = new ArrayList<AssociationResult>();
		if (!pairs.isEmpty()) {
			Collections.sort(pairs);
			ArrayList<AssociationResult> output = new ArrayList<>();
			int ctr = 0;
			while (ctr < len) {
				output.add(pairs.get(ctr).getAssociationResult());
				ctr++;
			}
			return output;
		}

		return null;
	}
}
