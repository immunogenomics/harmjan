package nl.harmjanwestra.harmonics.posthoc;

import JSci.maths.ArrayMath;
import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.features.*;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.BoxPlotPanel;
import nl.harmjanwestra.utilities.graphics.panels.HeatmapPanel;
import nl.harmjanwestra.utilities.graphics.panels.HistogramPanel;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import nl.harmjanwestra.utilities.peakfiles.SafFile;
import nl.harmjanwestra.utilities.peakfiles.XLSFile;
import nl.harmjanwestra.utilities.peakfiles.ZinbaPeaksFile;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.*;

/**
 * Created by hwestra on 9/1/15.
 */
public class PeakMerge {


	public static void main(String[] args) {
		try {
			PeakMerge p = new PeakMerge();
			String sampleFileName = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/FileToSampleName.txt";
			String output = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/Coverage-Summary/summary-overAllChr.txt";

			//p.mergeAllChrFiles(sampleFileName, 0, 2, output);
			output = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/Coverage-Summary/summary-perChr.txt";
			//p.mergeCoverageFilesNew(true, sampleFileName, 0, 3, output);
//			output = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/Coverage-Summary/summary-perChr-WolMT.txt";
//			p.mergeCoverageFiles(false, sampleFileName, 0, 3, output);

//			// PEAK MERGER
			String outdir = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/Peaks-Summary/";
			Gpio.createDir(outdir);

			// p.summarizeAndMergePeaks(sampleFileName, 0, 1, outdir, false, 0.05, 10);

			outdir = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/Peaks-Summary/";
			Gpio.createDir(outdir);

			p.summarizeAndMergePeaks("/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/FileToSampleName.txt", 0, 1, outdir, false, 0.05, 2);
			p.basePairOverlapBetweenPeaksWithoutMergedPeaks("/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/FileToSampleName.txt", 0, 1, outdir + "mergedpeaks.saf", outdir);

			outdir = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/Peaks-Summary-ZinbaCompare/";
			Gpio.createDir(outdir);
			sampleFileName = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/FileToSampleName-zinba.txt";

		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	private void mergeAllChrFiles(String sampleFileName, int samplecol, int filecol, String output) throws IOException {
		Pair<String[], String[]> samplecombo = loadSamplePairs(sampleFileName, samplecol, filecol);
		String[] coverageFiles = samplecombo.getRight();
		String[] sampleNames = samplecombo.getLeft();

		TextFile tfout = new TextFile(output, TextFile.W);
		String header = "sampleName"
				+ "\tnrReads1"
				+ "\tnrReads2"
				+ "\tproperpair"
				+ "\tnrMappingToChr"
				+ "\tnrDups"
				+ "\tunmapped"
				+ "\tsecondary"
				+ "\tnrAlignments"
				+ "\ttotalReads"
				+ "\tpercProperPairs"
				+ "\tpercMapping"
				+ "\tpercDups"
				+ "\tpercUnmapped"
				+ "\tpercSecondaryAlignments";

		tfout.writeln(header);

		for (int d = 0; d < coverageFiles.length; d++) {
			String covFile = coverageFiles[d];
			String sampleName = sampleNames[d];
			/*
			 0 nrReads1
			 1 nrReads2
			 2 properpair
			 3 nrMappingToChr
			 4 nrDups
			 5 unmapped
			 6 secondary
			 7 nrAlignments
			  */

			TextFile tf = new TextFile(covFile, TextFile.R);
			tf.readLine();//skip header
			String[] elems = tf.readLineElems(Strings.tab);


			Integer nrReads1 = Integer.parseInt(elems[0]);
			Integer nrReads2 = Integer.parseInt(elems[1]);
			Integer properpair = Integer.parseInt(elems[2]);
			Integer nrMappingToChr = Integer.parseInt(elems[3]);
			Integer nrDups = Integer.parseInt(elems[4]);
			Integer unmapped = Integer.parseInt(elems[5]);
			Integer secondary = Integer.parseInt(elems[6]);
			Integer nrAlignments = Integer.parseInt(elems[7]);

			int totalReads = nrReads1 + nrReads2;
			double percProperPairs = (double) properpair / nrReads2;
			double percMapping = (double) nrMappingToChr / totalReads;
			double percDups = (double) nrDups / totalReads;
			double percUnmapped = (double) unmapped / totalReads;
			double percSecondaryAlignments = (double) secondary / nrAlignments;

			tfout.writeln(sampleName
							+ "\t" + nrReads1
							+ "\t" + nrReads2
							+ "\t" + properpair
							+ "\t" + nrMappingToChr
							+ "\t" + nrDups
							+ "\t" + unmapped

							+ "\t" + secondary
							+ "\t" + nrAlignments
							+ "\t" + totalReads
							+ "\t" + percProperPairs
							+ "\t" + percMapping
							+ "\t" + percDups
							+ "\t" + percUnmapped
							+ "\t" + percSecondaryAlignments
			);


			tf.close();

		}

		tfout.close();


	}

	/*

	-- overview of number of peaks (significant and not significant) per sample
	-- overview of coverage underneath peaks
	-- overview of peak size per sample
	-- number of peaks shared between samples
	   -- merge peaks that overlap

	-- IDR?

	Shared peaks:
	-- comparison of peak height between samples
	-- comparison of significance of shared peaks
	-- variation in peak size between samples

	 */


	public void mergeCoverageFilesNew(boolean includeMT,
									  String sampleFile,
									  int sampleNamecol,
									  int filecol,
									  String output) throws IOException {

		Pair<String[], String[]> samplecombo = loadSamplePairs(sampleFile, sampleNamecol, filecol);
		String[] datasetFiles = samplecombo.getRight();
		String[] datasetNames = samplecombo.getLeft();

		/*
		0 Seq
		1 NrBases
        2 nrBasesCovered
3 percBasesCovered
4 nrReads
5 nrReadsPassingFilter
6 nrDups
7 percDups
8 nrDupsPassingFilter
9 nrFragments
10 nrFragmentsDup
11 nrFragmentsPassingFilter
12 avgInsertSize
13 avgCoverageAllBases
14 stDevCoverageAllBases
15 cvCoverageAllBases
16 avgCoverageCoveredBases
17 stDevCoverageCoveredBases
18 cvCoverageCoveredBases

		 */

		TextFile out = new TextFile(output, TextFile.W);
		out.writeln("Sample" +
				"\tpercBasesCovered" +
				"\tpercReadsPassingFilter" +
				"\tdupsOnAut" +
				"\tpercReadsOnMT" +
				"\tdupsOnMT" +
				"\tavgInsert" +
				"\tavgCoverageAllBases");
		for (int d = 0; d < datasetFiles.length; d++) {
			String dfile = datasetFiles[d];
			TextFile file = new TextFile(dfile, TextFile.R);
			file.readLine(); // skip header
			String[] elems = file.readLineElems(TextFile.tab);

			long totalbases = 0;
			long basesCovered = 0;
			long nrReadsOnMT = 0;
			long nrReadsOnOtherChr = 0;
			long nrReadsOnOtherChrPassingFilter = 0;
			long nrDupsOnMT = 0;
			long nrDupsOnOtherChr = 0;
			long totalReads = 0;
			double avgInsertSize = 0;
			double avgCoverage = 0;
			double avgCoverageCoveredBases = 0;


			while (elems != null) {

				/*
				0 Seq
		1 NrBases
        2 nrBasesCovered
4 nrReads
5 nrReadsPassingFilter
6 nrDups
12 avgInsertSize
13 avgCoverageAllBases
				 */
				Chromosome chr = Chromosome.parseChr(elems[0]);

				totalbases += Long.parseLong(elems[1]);
				basesCovered += Long.parseLong(elems[2]);
				long nrReads = Long.parseLong(elems[4]);
				totalReads += nrReads;
				long nrReadsPassingFilter = Long.parseLong(elems[5]);
				long nrDups = Long.parseLong(elems[6]);
				if (!chr.equals(Chromosome.MT)) {
					avgInsertSize += Double.parseDouble(elems[12]);
					avgCoverage += Double.parseDouble(elems[13]);
					nrReadsOnOtherChr += nrReads;
					nrReadsOnOtherChrPassingFilter += nrReadsPassingFilter;
					nrDupsOnOtherChr += nrDups;
					avgCoverageCoveredBases += Double.parseDouble(elems[16]);
				} else {
					nrReadsOnMT += nrReads;
					nrDupsOnMT += nrDups;
				}
				elems = file.readLineElems(TextFile.tab);
			}
			file.close();

			avgCoverage /= 24;
			avgInsertSize /= 24;
			avgCoverageCoveredBases /= 24;
			double percReadsOnMT = (double) nrReadsOnMT / totalReads;
			double percBasesCovered = (double) basesCovered / totalbases;
			double percReadsPassingFilter = (double) nrReadsOnOtherChrPassingFilter / nrReadsOnOtherChr;
			double dupsOnMT = (double) nrDupsOnMT / nrReadsOnMT;
			double dupsOnAut = (double) nrDupsOnOtherChr / nrReadsOnOtherChr;


			out.writeln(datasetNames[d]
					+ "\t" + percBasesCovered
					+ "\t" + percReadsPassingFilter
					+ "\t" + dupsOnAut
					+ "\t" + percReadsOnMT
					+ "\t" + dupsOnMT
					+ "\t" + avgInsertSize
					+ "\t" + avgCoverage
					+ "\t" + avgCoverageCoveredBases);

		}

		out.close();
	}

	public void mergeCoverageFiles(boolean includeMT,
								   String sampleFile,
								   int sampleNamecol,
								   int filecol,
								   String output) throws IOException {


		/*
0   Seq
1   NrBases
2   nrBasesCovered
x 3   percBasesCovered
x 4   nrReads
x 5   nrReadsPassingFilter
6   nrDups
x 7   percDups
8   nrDupsPassingFilter
x 9   avgInsertSize
x 10  avgCoverageAllBases
x 11  stDevCoverageAllBases
12  cvCoverageAllBases
x 13  avgCoverageCoveredBases
x 14  stDevCoverageCoveredBases
15  cvCoverageCoveredBases
		 */

		Pair<String[], String[]> samplecombo = loadSamplePairs(sampleFile, sampleNamecol, filecol);
		String[] datasetFiles = samplecombo.getRight();
		String[] datasetNames = samplecombo.getLeft();

		double[] percBasesCovered = new double[datasetNames.length];
		int[] nrreads = new int[datasetNames.length];
		int[] nrreadspassingfilter = new int[datasetNames.length];
		double[] percdups = new double[datasetNames.length];
		double[] avginsert = new double[datasetNames.length];
		double[] avgCovAllBases = new double[datasetNames.length];
		double[] stdevCovAllBases = new double[datasetNames.length];
		double[] avgCovCovered = new double[datasetNames.length];
		double[] stdevCovCovered = new double[datasetNames.length];

		TextFile out = new TextFile(output, TextFile.W);
		String header = "Sample\t" +
				"percBasesCovered\t" +
				"nrReads\t" +
				"nrReadsPassingFilter\t" +
				"percDups\t" +
				"avgInsertSize\t" +
				"avgCoverageAllBases\t" +
				"stDevCoverageAllBases\t" +
				"avgCoverageCoveredBases\t" +
				"stDevCoverageCoveredBases";
		out.writeln(header);


		for (int d = 0; d < datasetNames.length; d++) {

			TextFile dffile = new TextFile(datasetFiles[d], TextFile.R);

			dffile.readLine();
			String[] elems = dffile.readLineElems(TextFile.tab);
			int chrCtr = 0;
			while (elems != null) {

				Chromosome chr = Chromosome.parseChr(elems[0]);
				if (chr.isAutosome() || (chr.equals(Chromosome.MT) && includeMT)) {
					double basesCovered = Double.parseDouble(elems[3]);
					int nrReads = Integer.parseInt(elems[4]);
					int nrReadsPassingFilter = Integer.parseInt(elems[5]);
					double percDups = Double.parseDouble(elems[7]);
					double avgInsertSize = Double.parseDouble(elems[9]);
					double avgCoverageAllBases = Double.parseDouble(elems[10]);
					double stDevCoverageAllBases = Double.parseDouble(elems[11]);
					double avgCoverageCoveredBases = Double.parseDouble(elems[13]);
					double stDevCoverageCoveredBases = Double.parseDouble(elems[14]);

					percBasesCovered[d] += basesCovered;
					nrreads[d] += nrReads;
					nrreadspassingfilter[d] += nrReadsPassingFilter;
					percdups[d] += percDups;
					avginsert[d] += avgInsertSize;
					avgCovAllBases[d] += avgCoverageAllBases;
					stdevCovAllBases[d] += stDevCoverageAllBases * stDevCoverageAllBases;
					avgCovCovered[d] += avgCoverageCoveredBases;
					stdevCovCovered[d] += stDevCoverageCoveredBases * stDevCoverageCoveredBases;
					chrCtr++;
				}


				elems = dffile.readLineElems(TextFile.tab);
			}

			dffile.close();

			// average
			percBasesCovered[d] /= chrCtr;


			percdups[d] /= chrCtr;
			avginsert[d] /= chrCtr;
			avgCovAllBases[d] /= chrCtr;
			stdevCovAllBases[d] /= chrCtr;
			avgCovCovered[d] /= chrCtr;
			stdevCovCovered[d] /= chrCtr;

			System.out.println(chrCtr + " chr in " + datasetFiles[d]);

			// sqrt
			stdevCovAllBases[d] = Math.sqrt(stdevCovAllBases[d]);
			stdevCovCovered[d] = Math.sqrt(stdevCovCovered[d]);

			// write
			String ln = datasetNames[d]
					+ "\t" + percBasesCovered[d]
					+ "\t" + nrreads[d]
					+ "\t" + nrreadspassingfilter[d]
					+ "\t" + percdups[d]
					+ "\t" + avginsert[d]
					+ "\t" + avgCovAllBases[d]
					+ "\t" + stdevCovAllBases[d]
					+ "\t" + avgCovCovered[d]
					+ "\t" + stdevCovCovered[d];
			out.writeln(ln);

		}

		out.close();


	}


	public Pair<String[], String[]> loadSamplePairs(String sampleFileName, int samplecol, int filecol) throws IOException {
		TextFile samplefile = new TextFile(sampleFileName, TextFile.R);
		samplefile.readLine(); // skip header
		String[] elems = samplefile.readLineElems(Strings.tab);
		ArrayList<String> files = new ArrayList<String>();
		ArrayList<String> samples = new ArrayList<String>();

		while (elems != null) {
			if (elems.length > filecol) {
				String sample = elems[samplecol];
				String file = elems[filecol];
				if (file.trim().length() > 0) {
					if (Gpio.exists(file)) {
						files.add(file);
						samples.add(sample);
					}

				}
			}

			elems = samplefile.readLineElems(Strings.tab);
		}
		samplefile.close();
		System.out.println(files.size() + " file locations loaded.");
		return new Pair<String[], String[]>(samples.toArray(new String[0]), files.toArray(new String[0]));
	}


	public ArrayList<ArrayList<PeakFeature>> readPeakFiles(String[] peakFiles) throws IOException {
		ArrayList<ArrayList<PeakFeature>> output = new ArrayList<ArrayList<PeakFeature>>();
		for (int i = 0; i < peakFiles.length; i++) {
			System.out.println("Reading peaks: " + peakFiles[i]);
			ArrayList<PeakFeature> samplePeaks1;
			if (peakFiles[i].endsWith(".peaks")) {
				ZinbaPeaksFile p = new ZinbaPeaksFile();
				samplePeaks1 = p.readAllPeaks(peakFiles[i], false, 0.05);
			} else {
				XLSFile x = new XLSFile();
				samplePeaks1 = x.readAllPeaks(peakFiles[i], false, 0.05);
			}
			ArrayList<PeakFeature> featuresForSample = samplePeaks1;
			output.add(featuresForSample);
		}
		return output;
	}

	public ArrayList<Feature> mergeFeatures(ArrayList<PeakFeature> set1, ArrayList<PeakFeature> set2) {

		// combine features
		ArrayList<Feature> allFeaturesArr = new ArrayList<Feature>();
		allFeaturesArr.addAll(set1);
		allFeaturesArr.addAll(set2);

		ArrayList<Feature> output = new ArrayList<Feature>();

		for (Chromosome chr : Chromosome.values()) {
			ArrayList<Feature> chrFeatures = new ArrayList<Feature>();
			for (Feature f : allFeaturesArr) {
				if (f.getChromosome().equals(chr)) {
					chrFeatures.add(f);
				}
			}

			HashMap<Feature, Feature> remappedFeatures = new HashMap<Feature, Feature>();

			if (!chrFeatures.isEmpty()) {
				Collections.sort(chrFeatures, new FeatureComparator(false));

				for (int z = 0; z < chrFeatures.size(); z++) {
					Feature f = chrFeatures.get(z);
				}

				int ctr = 0;

				Feature origcurrent = chrFeatures.get(0);
				Feature remapFeat = copyFeat(origcurrent);
				Feature current = remapFeat;
				remappedFeatures.put(origcurrent, remapFeat);

				ctr++;
				int nrRemapped = 1;
				while (current != null) {
					if (ctr == chrFeatures.size()) {
						current = null;
					} else {
						Feature next = chrFeatures.get(ctr);
						ctr++;
						if (next.overlaps(remapFeat)) {
							// update remapped feature
							if (next.getStop() > remapFeat.getStop()) {
								remapFeat.setStop(next.getStop());
							}
							if (next.getStart() < remapFeat.getStart()) {
								remapFeat.setStart(next.getStart());
							}
							remappedFeatures.put(next, remapFeat);
						} else {
							current = next;

							remapFeat = copyFeat(current);
							nrRemapped++;
							remappedFeatures.put(current, remapFeat);
						}
					}
				}
				remappedFeatures.put(current, remapFeat);
			}

			HashSet<Feature> featureHash = new HashSet<Feature>();
			Set<Map.Entry<Feature, Feature>> remappedFeatureSet = remappedFeatures.entrySet();
			for (Map.Entry<Feature, Feature> set : remappedFeatureSet) {
				Feature val = set.getValue();
				if (!featureHash.contains(val)) {
					output.add(val);
					featureHash.add(val);
				}
			}
		}

		return output;


	}

	public void basePairOverlapBetweenPeaksWithoutMergedPeaks(String sampleFileName,
															  int samplecol,
															  int filecol,
															  String mergedPeakFile,
															  String outdir) throws IOException {


		SafFile mergedPeaksFile = new SafFile(mergedPeakFile, SafFile.R);
		ArrayList<Feature> mergedPeaks = mergedPeaksFile.read();
		mergedPeaksFile.close();

		Pair<String[], String[]> samplecombo = loadSamplePairs(sampleFileName, samplecol, filecol);
		String[] peakFiles = samplecombo.getRight();
		String[] sampleNames = samplecombo.getLeft();

		ArrayList<ArrayList<PeakFeature>> featuresPerSample = readPeakFiles(peakFiles);

		ArrayList<Track> tracksPerSample = convertToTracks(featuresPerSample, sampleNames);
		Track mergedTrack = new Track("merged");
		mergedTrack.addFeatures(mergedPeaks);


		int nrCols = 3;
		int remain = peakFiles.length % nrCols;
		int nrRows = (int) Math.ceil(peakFiles.length / nrCols) + remain;

		System.out.println(nrRows + " rows");
		System.out.println(nrCols + " cols");

		Grid peakWidthVsPeakHeight = new Grid(300, 300, nrRows, nrCols, 100, 100);
		Grid peakFoldEnrichVsPeakHeight = new Grid(300, 300, nrRows, nrCols, 100, 100);

		for (int i = 0; i < peakFiles.length; i++) {

			// determine correlation between peak size and nr of reads underneath peaks
			ArrayList<PeakFeature> featuresForSample = featuresPerSample.get(i);
			double[] peakWidths = new double[featuresForSample.size()];
			double[] peakHeights = new double[featuresForSample.size()];
			double[] peakFoldEnrichments = new double[featuresForSample.size()];

			for (int j = 0; j < featuresForSample.size(); j++) {
				PeakFeature f = featuresForSample.get(j);
				double width = f.getStop() - f.getStart();
				double height = f.getPileup();

				peakWidths[j] = width;
				peakHeights[j] = height;
				peakFoldEnrichments[j] = f.getFoldenrich();
			}

			ScatterplotPanel panelWidthVsHeight = new ScatterplotPanel(1, 1);
			panelWidthVsHeight.setLabels("PeakWidth", "PeakHeight");
			panelWidthVsHeight.setTitle(sampleNames[i]);
			panelWidthVsHeight.setData(peakWidths, peakHeights);

			ScatterplotPanel panelFoldChangeVsHeight = new ScatterplotPanel(1, 1);
			panelFoldChangeVsHeight.setLabels("Fold Enrichment", "PeakHeight");
			panelFoldChangeVsHeight.setTitle(sampleNames[i]);
			panelFoldChangeVsHeight.setData(peakFoldEnrichments, peakHeights);


			int row = i / nrCols;
			int col = i % nrCols;

			System.out.println(i + "\t" + row + "\t" + col);
			peakWidthVsPeakHeight.addPanel(panelWidthVsHeight, row, col);
			peakFoldEnrichVsPeakHeight.addPanel(panelFoldChangeVsHeight, row, col);
		}


		int[][] mergedPeaksNrBasesTotal = new int[featuresPerSample.size()][featuresPerSample.size()];
		int[][] mergedPeaksNrBasesOverlap = new int[featuresPerSample.size()][featuresPerSample.size()];
		for (int chr = 0; chr < 23; chr++) {
			for (int i = 0; i < featuresPerSample.size(); i++) {
				ArrayList<PeakFeature> features1 = getFeaturesForChr(featuresPerSample.get(i), chr);
				for (int j = i; j < featuresPerSample.size(); j++) {
					// create a merged track

					System.out.println("Determining merged peaks: chr: " + chr + "\td1: " + i + " - d2: " + j);

					ArrayList<PeakFeature> features2 = getFeaturesForChr(featuresPerSample.get(j), chr);

					ArrayList<Feature> mergedFeatures = mergeFeatures(features1, features2);

					for (Feature f : mergedFeatures) {
						int total = f.getStop() - f.getStart();
						mergedPeaksNrBasesTotal[i][j] += total;
						mergedPeaksNrBasesTotal[j][i] += total;
						for (Feature q : features1) {
							int overlap = q.determineBaseOverlap(f);
							mergedPeaksNrBasesOverlap[i][j] += overlap;
						}

						for (Feature q : features2) {
							int overlap = q.determineBaseOverlap(f);
							mergedPeaksNrBasesOverlap[j][i] += overlap;
						}

					}

				}
			}
		}

		double[][] mergedPeaksPercOverlap = new double[featuresPerSample.size()][featuresPerSample.size()];
		for (int i = 0; i < featuresPerSample.size(); i++) {
			for (int j = 0; j < featuresPerSample.size(); j++) {
				mergedPeaksPercOverlap[i][j] = (double) mergedPeaksNrBasesOverlap[i][j] / mergedPeaksNrBasesTotal[j][i];
			}
		}

		try {
			heatmapSharedPeaks(mergedPeaksPercOverlap, outdir + "pairwisePeakMerge.pdf", sampleNames);
		} catch (DocumentException e) {
			e.printStackTrace();
		}

		double[][] comparisonJaccardSums = new double[featuresPerSample.size()][featuresPerSample.size()];
		int[][] comparisonBasesOverlap = new int[featuresPerSample.size()][featuresPerSample.size()];
		int[][] nrComparisons = new int[featuresPerSample.size()][featuresPerSample.size()];
		Grid peakDistanceDistributions = new Grid(300, 300, featuresPerSample.size(), featuresPerSample.size(), 100, 100);
		for (int i = 0; i < featuresPerSample.size(); i++) {

			int[][] peakDistanceHistogram = new int[featuresPerSample.size()][500];

			for (int chr = 0; chr < 23; chr++) {
				// would be fastest with treeSet, but meh..
				ArrayList<PeakFeature> features1 = getFeaturesForChr(featuresPerSample.get(i), chr);
				for (int j = 0; j < featuresPerSample.size(); j++) {

					System.out.println("chr" + chr + "\tds1 " + i + "\tds2 " + j);
					ArrayList<PeakFeature> features2 = getFeaturesForChr(featuresPerSample.get(j), chr);

					for (int f1 = 0; f1 < features1.size(); f1++) {
						PeakFeature feat1 = features1.get(f1);

						for (int f2 = 0; f2 < features2.size(); f2++) {
							PeakFeature feat2 = features2.get(f2);

							if (feat1.overlaps(feat2)) {

								int summit1 = feat1.getSummit();
								int summit2 = feat2.getSummit();

								int distance = Math.abs(summit1 - summit2);
								if (distance >= 500) {
									distance = 499;
								}
								peakDistanceHistogram[j][distance]++;

								int sta = feat1.getStart();
								int sto = feat1.getStop();
								if (feat2.getStart() < sta) {
									sta = feat2.getStart();
								}
								if (feat2.getStop() > sto) {
									sto = feat2.getStop();
								}

								int total = sto - sta;

								// count the number of overlapping bases :)
								int overlap = feat1.determineBaseOverlap(feat2);

								double jaccard = (double) overlap / total;

								comparisonBasesOverlap[i][j] += overlap;

								comparisonJaccardSums[i][j] += jaccard;
								nrComparisons[i][j]++;
							}
						}
					}
				}
			}

			for (int j = 0; j < featuresPerSample.size(); j++) {
				HistogramPanel panel = new HistogramPanel(1, 1);
				panel.setTitle(sampleNames[i] + " - " + sampleNames[j]);
				panel.setLabels("Distance", "Log10(Frequency)");

				panel.setData(logTransform(peakDistanceHistogram[j]));
				peakDistanceDistributions.addPanel(panel, i, j);
			}
		}

		// compute overlap with merged track
		TextFile jaccardbaseoverlap = new TextFile(outdir + "", TextFile.W);
		for (int i = 0; i < sampleNames.length; i++) {
			Track sampleTrack = tracksPerSample.get(i);
			int overlap = 0;
			int nr = 0;
			int total = 0;
			for (int p = 0; p < mergedPeaks.size(); p++) {
				Feature feat1 = mergedPeaks.get(p);

				NavigableSet<Feature> set = sampleTrack.getFeatureSet(feat1);

				// determine overlap with all peaks overlapping feat 1
				if (set.size() > 0) {
					for (Feature f : set) {
						overlap += f.determineBaseOverlap(feat1);
					}
					nr += set.size();
					total += feat1.getStop() - feat1.getStart();
				}
			}

			double jaccard = (double) overlap / total;
			jaccardbaseoverlap.writeln(sampleNames[i] + "\t" + jaccard);
		}

		jaccardbaseoverlap.close();

		try {
			peakDistanceDistributions.draw(outdir + "PeakDistanceDistributions.png");
		} catch (DocumentException e) {
			e.printStackTrace();
		}


		for (int d = 0; d < comparisonJaccardSums.length; d++) {
			for (int d1 = 0; d1 < comparisonJaccardSums.length; d1++) {
				comparisonJaccardSums[d][d1] /= nrComparisons[d][d1];
			}
		}

		TextFile overlapout = new TextFile(outdir + "PeakOverlapBasesAverageJaccard.txt", TextFile.W);
		String overlapheader = "-\t" + Strings.concat(sampleNames, Strings.tab);
		overlapout.writeln(overlapheader);
		for (int d = 0; d < sampleNames.length; d++) {
			String ln = sampleNames[d];
			for (int d1 = 0; d1 < sampleNames.length; d1++) {
				ln += "\t" + comparisonJaccardSums[d][d1];
			}
			overlapout.writeln(ln);
		}
		overlapout.close();

		try {
			heatmapSharedPeaks(comparisonJaccardSums, outdir + "PeakOverlapBasesAverageJaccard.pdf", sampleNames);
		} catch (DocumentException e) {
			e.printStackTrace();
		}


		try {
			peakWidthVsPeakHeight.draw(outdir + "peakWidthVsPeakHeight.png");
		} catch (DocumentException e) {
			e.printStackTrace();
		}

		try {
			peakFoldEnrichVsPeakHeight.draw(outdir + "peakFoldEnrichVsPeakHeight.png");
		} catch (DocumentException e) {
			e.printStackTrace();
		}
	}


	private ArrayList<Track> convertToTracks(ArrayList<ArrayList<PeakFeature>> featuresPerSample, String[] sampleNames) {
		ArrayList<Track> output = new ArrayList<Track>();
		for (int i = 0; i < featuresPerSample.size(); i++) {
			Track sampleTrack = new Track(sampleNames[i]);
			ArrayList<PeakFeature> features = featuresPerSample.get(i);
			for (Feature f : features) {
				sampleTrack.addFeature(f);
			}
			output.add(sampleTrack);
		}
		return output;
	}

	private double[] logTransform(int[] peakDistanceHistogram) {

		double[] output = new double[peakDistanceHistogram.length];
		for (int i = 0; i < peakDistanceHistogram.length; i++) {
			if (peakDistanceHistogram[i] == 0) {
				output[i] = 0;
			} else {
				output[i] = Math.log10(peakDistanceHistogram[i]);
			}
		}
		return output;

	}

	private ArrayList<PeakFeature> getFeaturesForChr(ArrayList<PeakFeature> peakFeatures, int chr) {
		ArrayList<PeakFeature> output = new ArrayList<PeakFeature>();
		for (PeakFeature f : peakFeatures) {
			if (f.getChromosome().getNumber() == chr) {
				output.add(f);
			}
		}
		return output;
	}


	public void summarizeAndMergePeaks(String sampleFileName,
									   int samplecol,
									   int filecol,
									   String outdir,
									   boolean filterfdr,
									   double fdrthreshold,
									   int peakPresenceThreshold) throws IOException {

		Pair<String[], String[]> samplecombo = loadSamplePairs(sampleFileName, samplecol, filecol);

		String[] peakFiles = samplecombo.getRight();
		String[] sampleNames = samplecombo.getLeft();

		ArrayList<ArrayList<PeakFeature>> features = new ArrayList<ArrayList<PeakFeature>>();


		TextFile tf = new TextFile(outdir + "nrOfPeaks.txt", TextFile.W);
		String header = "Sample" +
				"\tNrPeaks" +
				"\tNrPeaksGt10" +
				"\tNrPeaks-bonferroni" +
				"\tNrPeaks-qval<0.05" +
				"\tAvgSize" +
				"\tstDevSize" +
				"\tAvgSize-qval<0.05" +
				"\tstDevSize-qval<0.05";
		tf.writeln(header);
		int covDistSize = 1000;
		int peakDistSize = 5000;
		int[][] coverageDist = new int[peakFiles.length][covDistSize];
		int[][] coverageDistSign = new int[peakFiles.length][covDistSize];
		int[][] peaksizeDist = new int[peakFiles.length][peakDistSize];
		int[][] peaksizeDistSign = new int[peakFiles.length][peakDistSize];


		double tenlogQval = -Math.log10(fdrthreshold);

		String safOutdir = outdir + "/safFiles/";
		Gpio.createDir(safOutdir);

		for (int i = 0; i < peakFiles.length; i++) {
			System.out.println("Reading peaks: " + peakFiles[i]);
			ArrayList<PeakFeature> samplePeaks;
			if (peakFiles[i].endsWith(".peaks")) {
				ZinbaPeaksFile p = new ZinbaPeaksFile();
				samplePeaks = p.readAllPeaks(peakFiles[i], filterfdr, fdrthreshold);
			} else {
				XLSFile x = new XLSFile();
				samplePeaks = x.readAllPeaks(peakFiles[i], filterfdr, fdrthreshold);
			}
			ArrayList<PeakFeature> featuresForSample = samplePeaks;
			saveAsSaf(safOutdir, peakFiles[i], featuresForSample);

			System.out.println(featuresForSample.size() + " peaks in file");
			features.add(featuresForSample);


			double tenlogbonferroni = -Math.log10(0.05 / featuresForSample.size());

			int nrPeaksBonferroni = 0;
			int nrPeaksQVal = 0;
			int nrPeaksGt10 = 0;


			double[] sizeArr = new double[featuresForSample.size()];
			ArrayList<Double> sizeArrQVal = new ArrayList<Double>();


			for (int q = 0; q < featuresForSample.size(); q++) {
				PeakFeature f = featuresForSample.get(q);
				int start = f.getStart();
				int end = f.getStop();
				int size = end - start;
				sizeArr[q] = size;
				double pileup = f.getPileup();
				double pval = f.getPval();
				double qval = f.getQval();

				int adjSize = size;
				if (adjSize >= peakDistSize) {
					adjSize = peakDistSize - 1;
				}
				peaksizeDist[i][adjSize]++;

				int pileupadj = (int) Math.ceil(pileup);
				if (pileupadj >= covDistSize) {
					pileupadj = covDistSize - 1;
				}
				coverageDist[i][pileupadj]++;

				if (qval > tenlogQval) {
					nrPeaksQVal++;
					sizeArrQVal.add((double) size);
					peaksizeDistSign[i][adjSize]++;
					coverageDistSign[i][pileupadj]++;
				}
				if (pval > tenlogbonferroni) {
					nrPeaksBonferroni++;

				}

				if (pileup > 10) {
					nrPeaksGt10++;
				}

			}

			double[] sizeArrQValPrim = Primitives.toPrimitiveArr(sizeArrQVal.toArray(new Double[0]));

			double meanSize = ArrayMath.mean(sizeArr);
			double sdSize = ArrayMath.standardDeviation(sizeArr);


			double meanSizeQVal = ArrayMath.mean(sizeArrQValPrim);
			double sdSizeQVal = ArrayMath.standardDeviation(sizeArrQValPrim);

			String out = sampleNames[i]
					+ "\t" + featuresForSample.size()
					+ "\t" + nrPeaksGt10
					+ "\t" + nrPeaksBonferroni
					+ "\t" + nrPeaksQVal
					+ "\t" + meanSize
					+ "\t" + sdSize
					+ "\t" + meanSizeQVal
					+ "\t" + sdSizeQVal;

			tf.writeln(out);

		}
		tf.close();

		try {
			boxPlotPeakSize(features, outdir + "peakSizeBoxplot.pdf", sampleNames);
		} catch (DocumentException e) {
			e.printStackTrace();
		}

		writedist(outdir + "peaksizedists.txt", peaksizeDist, peakDistSize, sampleNames);
		writedist(outdir + "peaksizedistsSign.txt", peaksizeDistSign, peakDistSize, sampleNames);
		writedist(outdir + "coveragedists.txt", coverageDist, covDistSize, sampleNames);

		try {
			boxPlotPeakCoverage(coverageDist, outdir + "peakCoverageBoxplot.pdf", sampleNames);
		} catch (DocumentException e) {
			e.printStackTrace();
		}

		writedist(outdir + "coveragedistsSign.txt", coverageDistSign, covDistSize, sampleNames);

		// get unique features
		HashSet<Feature> allFeaturesHash = new HashSet<Feature>();

		for (ArrayList<PeakFeature> sampleFeatures : features) {
			if (filterfdr) {

				for (PeakFeature f : sampleFeatures) {
					if (f.getQval() > tenlogQval) {
						allFeaturesHash.add(f);
					}
				}

			} else {
				allFeaturesHash.addAll(sampleFeatures);
			}
		}

		ArrayList<Feature> allFeaturesArr = new ArrayList<Feature>();
		allFeaturesArr.addAll(allFeaturesHash);

		// merge overlapping features
		int[][] sharedPeaks = new int[sampleNames.length][sampleNames.length];
		double[][] sharedPeaksJaccardUnion = new double[sampleNames.length][sampleNames.length];


		int[][] peaksPerChr = new int[sampleNames.length][Chromosome.values().length];

		int chrCtr = 0;

		TextFile remappedFeaturesOut = new TextFile(outdir + "mergedpeaks.saf", TextFile.W);
		TextFile allFeaturesOut = new TextFile(outdir + "allpeaks.bed", TextFile.W);
		TextFile mergedFeatsPerDsOut = new TextFile(outdir + "mergedpeaksperDs.txt", TextFile.W);
		TextFile peaksPresentInNDsOut = new TextFile(outdir + "mergedPeaksPresentIn" + peakPresenceThreshold + "Datasets.txt", TextFile.W);


		String sampleNamesCombined = Strings.concat(sampleNames, Strings.tab);
		header = "Peak\tTotalShared\t" + sampleNamesCombined;
		mergedFeatsPerDsOut.writeln(header);
		peaksPresentInNDsOut.writeln(header);

		ArrayList<Feature> remappedFeaturesPresenceAboveThreshold = new ArrayList<Feature>();
		int totalremapctr = 0;
		for (Chromosome chr : Chromosome.values()) {

			ArrayList<Feature> chrFeatures = new ArrayList<Feature>();
			for (Feature f : allFeaturesArr) {
				if (f.getChromosome().equals(chr)) {
					chrFeatures.add(f);
				}
			}

			System.out.println("chr: " + chr.getName() + " - peaks: " + chrFeatures.size());

			HashMap<Feature, Feature> remappedFeatures = new HashMap<Feature, Feature>();

			if (!chrFeatures.isEmpty()) {
				Collections.sort(chrFeatures, new FeatureComparator(false));

				for (int z = 0; z < chrFeatures.size(); z++) {
					Feature f = chrFeatures.get(z);
					allFeaturesOut.writeln(f.getChromosome().getName() + "\t" + f.getStart() + "\t" + f.getStop());
				}

				int ctr = 0;

				Feature origcurrent = chrFeatures.get(0);
				Feature remapFeat = copyFeat(origcurrent);
				Feature current = remapFeat;
				remappedFeatures.put(origcurrent, remapFeat);

				ctr++;
				int nrRemapped = 1;
				totalremapctr++;
				while (current != null) {
					if (ctr == chrFeatures.size()) {
						current = null;
					} else {
						Feature next = chrFeatures.get(ctr);
						ctr++;
						if (next.overlaps(remapFeat)) {
							// update remapped feature
							if (next.getStop() > remapFeat.getStop()) {
								remapFeat.setStop(next.getStop());
							}
							if (next.getStart() < remapFeat.getStart()) {
								remapFeat.setStart(next.getStart());
							}
							remappedFeatures.put(next, remapFeat);
						} else {
							current = next;
							remappedFeaturesOut.writeln("remap_" + totalremapctr + "\t" + remapFeat.getChromosome().getName() + "\t" + remapFeat.getStart() + "\t" + remapFeat.getStop() + "\t" + Strand.POS.toString());
							remapFeat = copyFeat(current);
							nrRemapped++;
							totalremapctr++;
							remappedFeatures.put(current, remapFeat);
						}
					}
				}

				// write the last remapped feature to disk
				remappedFeaturesOut.writeln("remap_" + totalremapctr + "\t" + remapFeat.getChromosome().getName() + "\t" + remapFeat.getStart() + "\t" + remapFeat.getStop() + "\t" + Strand.POS.toString());

				// index
				HashSet<Feature> featureHash = new HashSet<Feature>();
				Set<Map.Entry<Feature, Feature>> remappedFeatureSet = remappedFeatures.entrySet();
				ArrayList<Feature> remapFeatArr = new ArrayList<Feature>();
				for (Map.Entry<Feature, Feature> set : remappedFeatureSet) {
					Feature val = set.getValue();
					if (!featureHash.contains(val)) {
						remapFeatArr.add(val);
						featureHash.add(val);
					}
				}

				HashMap<Feature, Integer> remappedFeatureIndex = new HashMap<Feature, Integer>();
				for (int i = 0; i < remapFeatArr.size(); i++) {
					Feature f = remapFeatArr.get(i);
					remappedFeatureIndex.put(f, i);
				}


				// now count the number of remapped peaks per dataset
				byte[][] peakPresence = new byte[remappedFeatureIndex.size()][sampleNames.length];
				System.out.println("Nr peaks remapped: " + nrRemapped + " merged peaks\t" + remappedFeatures.size() + " unique original peaks");
				System.out.println(remappedFeatureIndex.size() + " merged features indexed");
				ArrayList<HashSet<Feature>> remappedFeatsPerDs = new ArrayList<HashSet<Feature>>();
				for (int d = 0; d < sampleNames.length; d++) {
					ArrayList<PeakFeature> f = features.get(d);
					HashSet<Feature> remaps = new HashSet<Feature>();

					if (filterfdr) {
						for (int q = 0; q < f.size(); q++) {
							PeakFeature sampleFeat = f.get(q);
							if (sampleFeat.getQval() > tenlogQval) {
								Feature remap = remappedFeatures.get(sampleFeat);
								if (remap != null) {
									Integer remapIndex = remappedFeatureIndex.get(remap);
									if (remapIndex == null) {
										System.out.println(remap.toString() + " does not have an index");
									} else {
										peakPresence[remapIndex][d] = 1;
									}
									remaps.add(remap);
								}
							}
						}
					} else {
						for (int q = 0; q < f.size(); q++) {
							Feature sampleFeat = f.get(q);
							Feature remap = remappedFeatures.get(sampleFeat);
							if (remap != null) {
								Integer remapIndex = remappedFeatureIndex.get(remap);
								if (remapIndex == null) {
									System.out.println(remap.toString() + " does not have an index");
								} else {
									peakPresence[remapIndex][d] = 1;
								}

								remaps.add(remap);
							}
						}
					}
					remappedFeatsPerDs.add(remaps);
					peaksPerChr[d][chrCtr] = remaps.size();
				}


				// write presence of peaks to file...
				for (int i = 0; i < remapFeatArr.size(); i++) {
					int sum = 0;
					Feature remap = remapFeatArr.get(i);
					String featLn = remap.toString();
					String sampleLn = "";
					for (int j = 0; j < sampleNames.length; j++) {
						sampleLn += "\t" + peakPresence[i][j];
						sum += peakPresence[i][j];
					}

					if (sum >= peakPresenceThreshold) {
						peaksPresentInNDsOut.writeln(featLn + "\t" + sum + sampleLn);
						remappedFeaturesPresenceAboveThreshold.add(remap);
					}
					mergedFeatsPerDsOut.writeln(featLn + "\t" + sum + sampleLn);

				}


				// and perform pairwise comparisons
				for (int d = 0; d < sampleNames.length; d++) {
					HashSet<Feature> d1Remaps = remappedFeatsPerDs.get(d);


					for (int d2 = 0; d2 < sampleNames.length; d2++) {
						HashSet<Feature> merger = new HashSet<Feature>();
						for (Feature f : d1Remaps) {
							merger.add(f);
						}

						HashSet<Feature> d2Remaps = remappedFeatsPerDs.get(d2);
						for (Feature f : d2Remaps) {
							merger.add(f);
						}
						int shared = sharedPeaks[d][d2];
						merger.addAll(d2Remaps);
						for (Feature f : d1Remaps) {
							if (d2Remaps.contains(f)) {
								shared++;
							}
						}
						sharedPeaksJaccardUnion[d][d2] += merger.size();
						sharedPeaks[d][d2] = shared;
					}
				}
			}
			chrCtr++;
		}
		peaksPresentInNDsOut.close();
		mergedFeatsPerDsOut.close();

		double[][] sharedPeaksJaccard = new double[sampleNames.length][sampleNames.length];
		for (int d = 0; d < sharedPeaksJaccard.length; d++) {
			for (int d2 = 0; d2 < sharedPeaksJaccard[d].length; d2++) {
				sharedPeaksJaccard[d][d2] = (double) sharedPeaks[d][d2] / sharedPeaksJaccardUnion[d][d2];
			}
		}

		allFeaturesOut.close();
		remappedFeaturesOut.close();

		TextFile overlapout = new TextFile(outdir + "overlappingpeaks.txt", TextFile.W);
		String overlapheader = "-\t" + Strings.concat(sampleNames, Strings.tab);
		overlapout.writeln(overlapheader);
		for (int d = 0; d < sampleNames.length; d++) {
			String ln = sampleNames[d];
			for (int d1 = 0; d1 < sampleNames.length; d1++) {
				ln += "\t" + sharedPeaks[d][d1];
			}
			overlapout.writeln(ln);
		}
		overlapout.close();

		TextFile overlapoutjaccard = new TextFile(outdir + "overlappingpeaksJaccard.txt", TextFile.W);

		overlapoutjaccard.writeln(overlapheader);
		for (int d = 0; d < sampleNames.length; d++) {
			String ln = sampleNames[d];
			for (int d1 = 0; d1 < sampleNames.length; d1++) {
				ln += "\t" + sharedPeaksJaccard[d][d1];
			}
			overlapoutjaccard.writeln(ln);
		}
		overlapoutjaccard.close();

		TextFile overlapoutPercOfTotal = new TextFile(outdir + "overlappingpeaksPercOfTotal.txt", TextFile.W);

		overlapoutPercOfTotal.writeln(overlapheader);
		double[][] percShared = new double[sampleNames.length][sampleNames.length];
		for (int d = 0; d < sampleNames.length; d++) {
			String ln = sampleNames[d];
			for (int d1 = 0; d1 < sampleNames.length; d1++) {
				ln += "\t" + ((double) sharedPeaks[d][d1] / sharedPeaks[d][d]);
				percShared[d][d1] = (double) sharedPeaks[d][d1] / sharedPeaks[d][d];
			}
			overlapoutPercOfTotal.writeln(ln);
		}
		overlapoutPercOfTotal.close();

		try {
			heatmapSharedPeaks(percShared, outdir + "percSharedPeaks.pdf", sampleNames);
		} catch (DocumentException e) {
			e.printStackTrace();
		}

		try {
			heatmapSharedPeaks(sharedPeaksJaccard, outdir + "jaccardSharedPeaks.pdf", sampleNames);
		} catch (DocumentException e) {
			e.printStackTrace();
		}


		TextFile nrpeaksperds = new TextFile(outdir + "nrmergedpeaksperds.txt", TextFile.W);
		nrpeaksperds.writeln(overlapheader);
		chrCtr = 0;
		for (Chromosome chr : Chromosome.values()) {
			String ln = chr.getName();
			for (int d1 = 0; d1 < sampleNames.length; d1++) {
				ln += "\t" + peaksPerChr[d1][chrCtr];
			}
			nrpeaksperds.writeln(ln);
			chrCtr++;
		}

		String totalLn = "Total";
		for (int d1 = 0; d1 < sampleNames.length; d1++) {
			chrCtr = 0;
			int sum = 0;
			for (Chromosome chr : Chromosome.values()) {
				sum += peaksPerChr[d1][chrCtr];
				chrCtr++;
			}

			totalLn += "\t" + sum;
		}
		nrpeaksperds.writeln(totalLn);
		nrpeaksperds.close();


		TextFile importantpeaksout = new TextFile(outdir + "peaksPresentInAtLeast" + peakPresenceThreshold + "samples.bed", TextFile.W);
		TextFile importantpeaksoutSAF = new TextFile(outdir + "peaksPresentInAtLeast" + peakPresenceThreshold + "samples.saf", TextFile.W);
		importantpeaksoutSAF.writeln("GeneID\tChr\tStart\tEnd\tStrand");

		Collections.sort(remappedFeaturesPresenceAboveThreshold, new FeatureComparator(false));
		for (int i = 0; i < remappedFeaturesPresenceAboveThreshold.size(); i++) {
			Feature f = remappedFeaturesPresenceAboveThreshold.get(i);
			importantpeaksoutSAF.writeln("Peak-" + i + "\t" + f.getChromosome().getName().toLowerCase() + "\t" + f.getStart() + "\t" + f.getStop() + "\t" + Strand.POS.toString());
			importantpeaksout.writeln(f.getChromosome().getName() + "\t" + f.getStart() + "\t" + f.getStop());
		}
		importantpeaksout.close();
		importantpeaksoutSAF.close();

	}

	private void boxPlotPeakCoverage(int[][] coverageDist, String outfilename, String[] sampleNames) throws IOException, DocumentException {
		Grid grid = new Grid(1200, 600, 1, 1, 100, 100);
		BoxPlotPanel boxplotpanel = new BoxPlotPanel(1, 1);
		Range r = new Range(0, 0, 0, 150);
		boxplotpanel.setDrawDataPoints(false);
		boxplotpanel.setLabels(sampleNames);
		boxplotpanel.setOutputIQRS(outfilename + "-iqrs.txt");
		boxplotpanel.setRange(r);
		boxplotpanel.setData(convertIntMatToDouble(coverageDist));
		grid.addPanel(boxplotpanel);
		grid.draw(outfilename);
	}

	private double[][] convertIntMatToDouble(int[][] coverageDist) {
		double[][] data = new double[coverageDist.length][0];
		for (int i = 0; i < data.length; i++) {
			int nrVals = 0;
			for (int j = 0; j < coverageDist[i].length; j++) {

				nrVals += coverageDist[i][j];
			}
			double[] datadata = new double[nrVals];
			int val = 0;
			for (int j = 0; j < coverageDist[i].length; j++) {
				int ct = coverageDist[i][j];
				for (int q = 0; q < ct; q++) {
					datadata[val] = j;
					val++;
				}
			}
			data[i] = datadata;
		}
		return data;
	}

	private void heatmapSharedPeaks(double[][] percShared, String outfilename, String[] sampleNames) throws IOException, DocumentException {
		Grid grid = new Grid(1200, 1200, 1, 1, 100, 100);
		HeatmapPanel panel = new HeatmapPanel(1, 1);
		panel.setRange(new Range(0, 0, 0, 1));
		panel.setData(percShared, sampleNames, sampleNames);

		grid.addPanel(panel);
		grid.draw(outfilename);
	}

	private void boxPlotPeakSize(ArrayList<ArrayList<PeakFeature>> features, String outfilename, String[] sampleNames) throws IOException, DocumentException {

		Grid grid = new Grid(1200, 600, 1, 1, 100, 100);
		BoxPlotPanel boxplotpanel = new BoxPlotPanel(1, 1);
		Range r = new Range(0, 0, 0, 1250);
		boxplotpanel.setDrawDataPoints(false);
		boxplotpanel.setUseMeanAndSd(true);
		boxplotpanel.setLabels(sampleNames);
		boxplotpanel.setRange(r);
		boxplotpanel.setOutputIQRS(outfilename + "-iqrs.txt");
		boxplotpanel.setData(convertFeaturesToDouble(features));
		grid.addPanel(boxplotpanel);
		grid.draw(outfilename);

	}

	private double[][] convertFeaturesToDouble(ArrayList<ArrayList<PeakFeature>> features) {
		double[][] data = new double[features.size()][0];
		for (int i = 0; i < features.size(); i++) {
			ArrayList<PeakFeature> f = features.get(i);

			double[] datai = new double[f.size()];
			for (int j = 0; j < f.size(); j++) {
				PeakFeature p = f.get(j);
				datai[j] = p.getStop() - p.getStart();
			}
			data[i] = datai;
			System.out.println("dataset " + i + " " + datai.length + " peaks");
		}

		return data;

	}

	private void saveAsSaf(String safOutdir, String sampleFileName, ArrayList<PeakFeature> featuresForSample) throws IOException {
		Collections.sort(featuresForSample, new FeatureComparator(false));

		String[] fname = sampleFileName.split("/");
		String sampleName = fname[fname.length - 1];

		TextFile out = new TextFile(safOutdir + sampleName + ".saf", TextFile.W);

		for (int i = 0; i < featuresForSample.size(); i++) {
			Feature f = featuresForSample.get(i);
			out.writeln("Peak_" + i
					+ "\t" + f.getChromosome().getName().toLowerCase()
					+ "\t" + f.getStart()
					+ "\t" + f.getStop()
					+ "\t" + Strand.POS.toString());
		}

		out.close();


	}

	private Feature copyFeat(Feature current) {
		Feature output = new Feature();

		output.setChromosome(current.getChromosome());
		output.setName(current.getName());
		output.setStart(current.getStart());
		output.setStop(current.getStop());
		output.setStrand(current.getStrand());
		return output;
	}

	private void writedist(String name, int[][] dist, int size, String[] sampleNames) throws IOException {
		TextFile peakdistoutsign = new TextFile(name, TextFile.W);
		String headerpeak = "Bin";
		String samples = Strings.concat(sampleNames, Strings.tab);
		peakdistoutsign.writeln(headerpeak + "\t" + samples);
		for (int i = 0; i < size; i++) {
			String ln = "" + i;
			for (int j = 0; j < sampleNames.length; j++) {
				ln += "\t" + dist[j][i];
			}
			peakdistoutsign.writeln(ln);
		}
		peakdistoutsign.close();
	}
}
