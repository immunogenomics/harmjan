package nl.harmjanwestra.finemapping;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.AssociationResultPValuePair;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

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


		int len = 5;
		TextFile out = new TextFile(outfile, TextFile.R);
		String header2 = "";
		String header1 = "region";
		for (int i = 0; i < data.length; i++) {
			header1 += "\tNrVariantsInCredibleSet\tSumPosteriorTop5Variants\tVariants\tAlleles\tOR\tPval\tPosterior";
		}
		out.writeln(header1);
		for (int regionId = 0; regionId < regionsWithCredibleSets.size(); regionId++) {
			Feature region = regionsWithCredibleSets.get(regionId);

			// region nrCrediblesetVariants posteriorsumtop5 topvariants alleles or pval posterior
			ArrayList<ArrayList<AssociationResult>> resultsPerDs = new ArrayList<>();
			for (int i = 0; i < data.length; i++) {
				ArrayList<AssociationResult> topResults = getTopVariants(data, i, regionId, len);
				resultsPerDs.add(topResults);
			}

			double[] regionsums = new double[data.length];
			for (int snpId = 0; snpId < len; snpId++) {
				String ln = "";

				if (snpId == 0) {
					ln = region.toBedString();
					for (int datasetId = 0; datasetId < data.length; datasetId++) {
						double sum = 0;
						AssociationResult r = data[datasetId][regionId][snpId];
						for (int var = 0; var < 5; var++) {
							sum += data[datasetId][regionId][var].getPosterior();
						}
						ln += "\t" + crediblesets[datasetId][regionId].length
								+ "\t" + sum
								+ "\t" + r.getSnp().toString()
								+ "\t" + Strings.concat(r.getSnp().getAlleles(), Strings.comma)
								+ "\t" + Strings.concat(r.getORs(), Strings.semicolon)
								+ "\t" + r.getLog10Pval()
								+ "\t" + r.getPosterior();
						regionsums[datasetId] += r.getPosterior();
					}
				} else {
					ln = "";
					for (int datasetId = 0; datasetId < data.length; datasetId++) {

						AssociationResult r = data[datasetId][regionId][snpId];
						if (regionsums[datasetId] < 0.9) {
							ln += "\t"
									+ "\t"
									+ "\t" + r.getSnp().toString()
									+ "\t" + Strings.concat(r.getSnp().getAlleles(), Strings.comma)
									+ "\t" + Strings.concat(r.getORs(), Strings.semicolon)
									+ "\t" + r.getLog10Pval()
									+ "\t" + r.getPosterior();
						} else {
							ln += "\t"
									+ "\t"
									+ "\t"
									+ "\t"
									+ "\t"
									+ "\t"
									+ "\t";
						}
						regionsums[datasetId] += r.getPosterior();


					}

				}
				out.writeln(ln);
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
