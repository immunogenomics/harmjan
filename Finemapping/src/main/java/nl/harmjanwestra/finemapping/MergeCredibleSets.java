package nl.harmjanwestra.finemapping;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.AssociationResultPValuePair;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.CircularHeatmapPanel;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 9/7/16.
 */
public class MergeCredibleSets {


	public static void main(String[] args) {
		String bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		String genenames = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/AllLoci-GenesPerLocus.txt";


		try {
			MergeCredibleSets c = new MergeCredibleSets();
//			c.run(bedregions, assocfiles, datasetnames, genenames, outfile);


			String outfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/MergedCredibleSets/mergedCredibleSets.txt";
			String outplot = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/MergedCredibleSets/mergedCredibleSets-plot.pdf";
			String[] assocfiles = new String[]{
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/RA-assoc0.3-COSMO-merged-posterior.txt.gz",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/META-assoc0.3-COSMO-merged-posterior.txt.gz"

			};
			String[] datasetnames = new String[]{
					"RA",
					"T1D",
					"Combined"
			};
			double threshold = 7.5E-7;
			int nrVariantsInCredibleSet = 10;
			double maxPosteriorCredibleSet = 0.9;
			c.run(bedregions, assocfiles, datasetnames, genenames, outfile, maxPosteriorCredibleSet, nrVariantsInCredibleSet);
			c.makePlot(bedregions, assocfiles, datasetnames, genenames, outplot, threshold, maxPosteriorCredibleSet, nrVariantsInCredibleSet);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
	}

	public void run(String bedregions,
					String[] assocFiles,
					String[] datasetNames,
					String genenamefile,
					String outfile, double maxPosteriorCredibleSet, int maxNrVariantsInCredibleSet) throws IOException {


		HashMap<String, String> locusToGene = new HashMap<String, String>();
		TextFile genefiletf = new TextFile(genenamefile, TextFile.R);
		String[] genelems = genefiletf.readLineElems(TextFile.tab);
		while (genelems != null) {
			if (genelems.length > 1) {
				locusToGene.put(genelems[0], genelems[1]);
			} else {
				locusToGene.put(genelems[0], "");
			}
			genelems = genefiletf.readLineElems(TextFile.tab);
		}
		genefiletf.close();

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedregions);

		// region variant1 pval1 posterior1 variant2 pval2 posterior2 variant3 pval3 posterior3
		AssociationFile f = new AssociationFile();


		// [dataset][region][variants]
		AssociationResult[][][] crediblesets = new AssociationResult[assocFiles.length][regions.size()][];
		AssociationResult[][][] data = new AssociationResult[assocFiles.length][regions.size()][];


		ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
		ArrayList<Feature> regionsWithCredibleSets = new ArrayList<>();
		for (int d = 0; d < regions.size(); d++) {
			boolean hasSet = false;
			for (int i = 0; i < assocFiles.length; i++) {
				ArrayList<AssociationResult> allDatasetData = f.read(assocFiles[i], regions.get(d));
				data[i][d] = allDatasetData.toArray(new AssociationResult[0]);
				ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(allDatasetData, maxPosteriorCredibleSet);
				crediblesets[i][d] = credibleSet.toArray(new AssociationResult[0]);
				if (credibleSet.size() <= maxNrVariantsInCredibleSet) {
					hasSet = true;
				}
			}
			if (hasSet) {
				regionsWithCredibleSets.add(regions.get(d));
			}
		}

		HashSet<Feature> regionsWithCredibleSetsHash = new HashSet<Feature>();
		regionsWithCredibleSetsHash.addAll(regionsWithCredibleSets);

		int len = maxNrVariantsInCredibleSet;
		TextFile out = new TextFile(outfile, TextFile.W);
		TextFile outg = new TextFile(outfile + "-genes.txt", TextFile.W);
		String header2 = "\t\t";

		String header1 = "region\tgene";
		for (int i = 0; i < data.length; i++) {
			header2 += datasetNames[i]
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t";
			header1 += "\tNrVariantsInCredibleSet" +
					"\tSumPosteriorTop5Variants" +
					"\tVariants" +
					"\tAlleles" +
					"\tOR" +
					"\tPval" +
					"\tPosterior";
		}
		out.writeln(header2);
		out.writeln(header1);
		for (int regionId = 0; regionId < regions.size(); regionId++) {
			Feature region = regions.get(regionId);
			if (regionsWithCredibleSetsHash.contains(region)) {
				// region nrCrediblesetVariants posteriorsumtop5 topvariants alleles or pval posterior
				ArrayList<ArrayList<AssociationResult>> resultsPerDs = new ArrayList<>();
				for (int i = 0; i < data.length; i++) {
					ArrayList<AssociationResult> topResults = getTopVariants(data, i, regionId, len);
					resultsPerDs.add(topResults);
				}


				double[] sumsperregion = new double[datasetNames.length];
				for (int datasetId = 0; datasetId < datasetNames.length; datasetId++) {
					double sum = 0;
					for (int snpId = 0; snpId < len; snpId++) {
						AssociationResult r = resultsPerDs.get(datasetId).get(snpId); //data[datasetId][regionId][snpId];

						sum += r.getPosterior();

					}
					sumsperregion[datasetId] = sum;
				}


				outg.writeln(region.toString() + "\t" + locusToGene.get(region.toString()));
				double[] regionsums = new double[data.length];

				for (int snpId = 0; snpId < len; snpId++) {
					String ln = "";


					boolean allSNPsPrinted = true;
					for (int datasetId = 0; datasetId < data.length; datasetId++) {
						if (regionsums[datasetId] < maxPosteriorCredibleSet) {
							allSNPsPrinted = false;
						}
					}

					if (!allSNPsPrinted) {
						if (snpId == 0) {
							ln = region.toString() + "\t" + locusToGene.get(region.toString());
							for (int datasetId = 0; datasetId < data.length; datasetId++) {
								AssociationResult r = resultsPerDs.get(datasetId).get(snpId); //data[datasetId][regionId][snpId];
								ln += "\t" + crediblesets[datasetId][regionId].length
										+ "\t" + sumsperregion[datasetId]
										+ "\t" + r.getSnp().toString()
										+ "\t" + Strings.concat(r.getSnp().getAlleles(), Strings.comma)
										+ "\t" + Strings.concat(r.getORs(), Strings.semicolon)
										+ "\t" + r.getLog10Pval()
										+ "\t" + r.getPosterior();
								regionsums[datasetId] += r.getPosterior();
							}
						} else {

							for (int datasetId = 0; datasetId < data.length; datasetId++) {
								if (datasetId == 0) {
									ln = "\t\t\t\t";
								} else {
									ln += "\t\t\t";
								}

								AssociationResult r = resultsPerDs.get(datasetId).get(snpId); //data[datasetId][regionId][snpId];
								if (regionsums[datasetId] < maxPosteriorCredibleSet) {
									ln += r.getSnp().toString()
											+ "\t" + Strings.concat(r.getSnp().getAlleles(), Strings.comma)
											+ "\t" + Strings.concat(r.getORs(), Strings.semicolon)
											+ "\t" + r.getLog10Pval()
											+ "\t" + r.getPosterior();
								} else {
									ln += "\t"
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
			}


		}
		outg.close();
		out.close();

	}


	private void makePlot(String bedregions, String[] assocFiles, String[] datasetNames, String genenamefile, String outfile, double threshold, double maxPosteriorCredibleSet, int nrVariantsInCredibleSet) throws IOException, DocumentException {


		HashMap<String, String> locusToGene = new HashMap<String, String>();
		TextFile genefiletf = new TextFile(genenamefile, TextFile.R);
		String[] genelems = genefiletf.readLineElems(TextFile.tab);
		while (genelems != null) {
			if (genelems.length > 1) {
				locusToGene.put(genelems[0], genelems[1]);
			} else {
				locusToGene.put(genelems[0], "");
			}
			genelems = genefiletf.readLineElems(TextFile.tab);
		}
		genefiletf.close();

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedregions);

		// region variant1 pval1 posterior1 variant2 pval2 posterior2 variant3 pval3 posterior3
		AssociationFile f = new AssociationFile();


		// [dataset][region][variants]
		AssociationResult[][][] crediblesets = new AssociationResult[assocFiles.length][regions.size()][];
		AssociationResult[][][] data = new AssociationResult[assocFiles.length][regions.size()][];


		ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
		ArrayList<Feature> regionsWithCredibleSets = new ArrayList<>();
		for (int d = 0; d < regions.size(); d++) {
			boolean hasSet = false;
			for (int i = 0; i < assocFiles.length; i++) {
				ArrayList<AssociationResult> allDatasetData = f.read(assocFiles[i], regions.get(d));
				data[i][d] = allDatasetData.toArray(new AssociationResult[0]);
				ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(allDatasetData, maxPosteriorCredibleSet);

				boolean abovethresh = false;
				for (AssociationResult r : credibleSet) {
					if (r.getPval() < threshold) {
						abovethresh = true;
					}
				}

				crediblesets[i][d] = credibleSet.toArray(new AssociationResult[0]);
				if (credibleSet.size() <= nrVariantsInCredibleSet && abovethresh) {
					hasSet = true;
				}
			}
			if (hasSet) {
				regionsWithCredibleSets.add(regions.get(d));
			}
		}

		HashSet<Feature> regionsWithCredibleSetsHash = new HashSet<Feature>();
		regionsWithCredibleSetsHash.addAll(regionsWithCredibleSets);

		int len = nrVariantsInCredibleSet;

		double[][][] dataForPlotting = new double[data.length][regionsWithCredibleSets.size()][];
		double[][][] dataForPlotting2 = new double[data.length][regionsWithCredibleSets.size()][];

		int regionCtr = 0;
		ArrayList<Triple<Integer, Integer, String>> groups = new ArrayList<Triple<Integer, Integer, String>>();
		ArrayList<String> allColNames = new ArrayList<>();
		ArrayList<String> groupnames = new ArrayList<>();
		int prevCtr = 0;
		for (int regionId = 0; regionId < regions.size(); regionId++) {
			Feature region = regions.get(regionId);
			if (regionsWithCredibleSetsHash.contains(region)) {
				// region nrCrediblesetVariants posteriorsumtop5 topvariants alleles or pval posterior
				ArrayList<ArrayList<AssociationResult>> resultsPerDs = new ArrayList<>();
				HashMap<String, Integer> variantToInt = new HashMap<String, Integer>();
				ArrayList<String> variantsInRegion = new ArrayList<String>();
				for (int i = 0; i < data.length; i++) {
					ArrayList<AssociationResult> topResults = getTopVariants(data, i, regionId, len);
					resultsPerDs.add(topResults);

					double sum = 0;
					for (int s = 0; s < topResults.size(); s++) {
						double posterior = topResults.get(s).getPosterior();
						if (sum < maxPosteriorCredibleSet) {
							String variant = topResults.get(s).getSnp().getName();
							boolean isSignificant = (topResults.get(s).getPval() < threshold);
							if (!variantToInt.containsKey(variant) && isSignificant) {
								variantToInt.put(variant, variantToInt.size());
								variantsInRegion.add(variant);
							}
						}
						sum += posterior;
					}
				}

				allColNames.addAll(variantsInRegion);

				for (int i = 0; i < data.length; i++) {
					dataForPlotting[i][regionCtr] = new double[variantToInt.size()];
					dataForPlotting2[i][regionCtr] = new double[variantToInt.size()];

					ArrayList<AssociationResult> dsResults = resultsPerDs.get(i);
					for (int s = 0; s < dsResults.size(); s++) {
						String variant = dsResults.get(s).getSnp().getName();
						Integer variantId = variantToInt.get(variant);
						double posterior = dsResults.get(s).getPosterior();

						if (variantId != null) {
							if (variantId > dataForPlotting[i][regionCtr].length - 1) {
								System.out.println();
								System.out.println("WEIRDNESSSSSSSSSSSSSSSSS: " + variant);
								System.out.println();
							} else {
								if (dsResults.get(s).getPval() < threshold) {
									dataForPlotting[i][regionCtr][variantId] = 0.10 + (0.90 * posterior);
									dataForPlotting2[i][regionCtr][variantId] = 1;
								} else {
									dataForPlotting[i][regionCtr][variantId] = 0;
									dataForPlotting2[i][regionCtr][variantId] = 0;
								}
							}

						}
					}
				}

				DecimalFormat decimalFormat = new DecimalFormat("#");
				decimalFormat.setGroupingUsed(true);
				decimalFormat.setGroupingSize(3);

				String locusName = region.getChromosome().toString() + ":" + decimalFormat.format(region.getStart()) + "-" + decimalFormat.format(region.getStop());
				// locusName = region.toString();
				locusName += "; " + locusToGene.get(region.toString());
				groupnames.add(locusName);

				groups.add(new Triple<Integer, Integer, String>(regionCtr, regionCtr + 1, locusName));
				System.out.println("region\t" + region.toString() + "\t" + locusName);
				regionCtr++;
			}
		}


		Grid grid = new Grid(1000, 1000, 1, 1, 100, 0);
		CircularHeatmapPanel panel = new CircularHeatmapPanel(1, 1);
		panel.setData(datasetNames, groupnames.toArray(new String[0]), dataForPlotting);
		panel.setGroups(groups);
		grid.addPanel(panel);
		grid.draw(outfile);

		grid = new Grid(1000, 1000, 1, 1, 100, 0);
		panel = new CircularHeatmapPanel(1, 1);
		panel.setData(datasetNames, groupnames.toArray(new String[0]), dataForPlotting2);
		panel.setGroups(groups);
		grid.addPanel(panel);
		grid.draw(outfile + "-bin.pdf");

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
