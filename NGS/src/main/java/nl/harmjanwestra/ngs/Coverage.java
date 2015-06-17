package nl.harmjanwestra.ngs;

import com.xeiam.xchart.*;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import nl.harmjanwestra.utilities.bamfile.BamFileReader;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.FeatureMerger;
import nl.harmjanwestra.utilities.graphics.ColorGenerator;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.broadinstitute.gatk.engine.filters.*;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.HCMappingQualityFilter;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Created by hwestra on 1/29/15.
 */
public class Coverage {

	public static void main(String[] args) {

//		Coverage c = new Coverage();
//		try {
////			// c.bamToBed(bamfile, outdir);
////
//			String coveragefile = "/Data/Projects/2014-FR-Reseq/2015-01-31-HaplotypeCallerVariants-762Samples/LowQualBamFiles/region-coverage.txt";
//			String sampefile = "/Data/Projects/2014-FR-Reseq/2015-02-CoverageMetrics/runsAndSampleNames.txt";
//			String outdir = "/Data/Projects/2014-FR-Reseq/2015-01-31-HaplotypeCallerVariants-762Samples/LowQualBamFiles/plots/";
//			c.createCoveragePlots(coveragefile, sampefile, outdir);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//
		Coverage c = new Coverage();

//		try {
//
//			String regions = args[0];
//			String inputdir = args[1];
//			String outputdir = args[2];
//			c.correlateCoverageWithinRegionFiles(regions, inputdir, outputdir);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

		try {
			if (args.length < 5) {
				System.out.println("Usage: bamfilelist outdir regionfile makebedgraphfiles threads");
			} else {
				c.bamToBedWithinRegionsForList(args[0], args[1], args[2], Boolean.parseBoolean(args[3]), Integer.parseInt(args[4]));
////					c.bamToBedWithinRegions(file, sampleOutDir, targetregions);
				//String bamfile, String outdir, String regionFile, boolean outputcoverageperregion
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

	}




	public void bamToBedWithinRegionsForList(String listFile, String outdir, String targetregions, boolean outputcoverageperregion, int threads) throws IOException {
		TextFile lf = new TextFile(listFile, TextFile.R);
		ArrayList<String> samples = new ArrayList<String>();
		ArrayList<String> files = new ArrayList<String>();

		String[] lfelems = lf.readLineElems(TextFile.tab);
		while (lfelems != null) {

			samples.add(lfelems[0]);
			files.add(lfelems[1]);

			lfelems = lf.readLineElems(TextFile.tab);
		}

		lf.close();


		int filectr = 0;

		if (threads <= 0) {
			threads = 1;
		}

		System.out.println("Loading regions from: " + targetregions);
		ArrayList<Feature> features = loadRegions(targetregions);
		System.out.println(features.size() + " target regions loaded");
		if (outputcoverageperregion) {
			System.out.println("Merging overlapping regions for .bedGraph output");
			features = FeatureMerger.merge(features);
			System.out.println(features.size() + " regions remain after mergin overlapping regions");
		}

		System.out.println("Opening threadpool for " + threads + " threads.");
		ExecutorService threadPool = Executors.newFixedThreadPool(threads);
		CompletionService<Boolean> pool = new ExecutorCompletionService<Boolean>(threadPool);

		for (int i = 0; i < samples.size(); i++) {

			String sample = samples.get(i);
			String file = files.get(i);


			String sampleOutDir = outdir + "samples/" + sample + "/";
			Gpio.createDir(outdir);
			Gpio.createDir(sampleOutDir);


			CoverageTask t = new CoverageTask(file, sampleOutDir, features, outputcoverageperregion);
			pool.submit(t);

			filectr++;
		}


		int returned = 0;
		while (returned < files.size()) {

			try {
				Boolean result = pool.take().get();
				if (result) {
					returned++;
				}

				System.out.println(returned + " / " + files.size() + " returned");
			} catch (Exception e) {
				e.printStackTrace();
			}

		}

		threadPool.shutdown();

		// concatenate everything
		HashMap<String, Integer> regionOrder = null;
		TextFile out = new TextFile(outdir + "region-coverage.txt", TextFile.W);
		TextFile out2 = new TextFile(outdir + "region-summary.txt", TextFile.W);

		boolean headerWritten = false;


		for (int i = 0; i < samples.size(); i++) {

			String sample = samples.get(i);
			String file = files.get(i);


			String sampleOutDir = outdir + "samples/" + sample + "/";
			String regionCoverage = sampleOutDir + "/region-coverage.txt.gz";
			TextFile tf = new TextFile(regionCoverage, TextFile.R);
			String[] header = tf.readLineElems(TextFile.tab);
			if (regionOrder == null) {
				regionOrder = new HashMap<String, Integer>();
				out.writeln(Strings.concat(header, Strings.tab));
				for (int j = 0; j < header.length; j++) {
					regionOrder.put(header[j], j);
				}
			}
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				String[] lineOut = new String[header.length];
				lineOut[0] = elems[0];
				for (int j = 1; j < header.length; j++) {
					Integer newPos = regionOrder.get(header[j]);
					lineOut[newPos] = elems[j];
				}
				out.writeln(Strings.concat(lineOut, Strings.tab));
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

			String regionSummary = sampleOutDir + "/summary.txt.gz";


			TextFile q = new TextFile(regionSummary, TextFile.R);
			elems = q.readLineElems(TextFile.tab);
			int lnctr = 0;
			while (elems != null) {

				if (!headerWritten) {
					out2.writeln(Strings.concat(elems, Strings.tab));
					headerWritten = true;
				} else if (lnctr > 0) {
					out2.writeln(Strings.concat(elems, Strings.tab));
				}

				lnctr++;
				elems = q.readLineElems(TextFile.tab);
			}

			q.close();
		}
		out.close();
		out2.close();

		if (outputcoverageperregion) {
			// concatenate region files...

			//

			HashMap<String, Integer> featureMap = new HashMap<String, Integer>();
			ArrayList<String> featureStr = new ArrayList<String>();
			for (int i = 0; i < features.size(); i++) {

				featureMap.put(features.get(i).toString(), i);
				featureStr.add(features.get(i).toString());
			}

			double[][] output1 = new double[features.size()][files.size()];
			double[][] output2 = new double[500][files.size()];

			ArrayList<String> rowNames2 = new ArrayList<String>();
			for (int i = 0; i < output2.length; i++) {
				rowNames2.add("" + i);
			}

			int fctr = 0;
			ArrayList<String> sampleNames = new ArrayList<String>();
			for (int i = 0; i < samples.size(); i++) {

				String sample = samples.get(i);
				String file = files.get(i);

				String sampleOutDir = outdir + "samples/" + sample + "/";

				DoubleMatrixDataset<String, String> ds1 = new DoubleMatrixDataset<String, String>(sampleOutDir + "region-percreadsgt20.txt.gz");
				DoubleMatrixDataset<String, String> ds2 = new DoubleMatrixDataset<String, String>(sampleOutDir + "histogramCumulative.txt.gz");


				for (int col = 0; col < ds1.nrCols; col++) {
					String feat = ds1.colObjects.get(col);
					Integer id = featureMap.get(feat);
					if (id != null) {
						output1[id][fctr] = ds1.rawData[0][col];
					}
				}

				sampleNames.add(ds2.colObjects.get(0));
				for (int j = 0; j < ds2.nrRows; j++) {
					output2[j][fctr] = ds2.rawData[j][0];
				}
				fctr++;
			}

			DoubleMatrixDataset<String, String> dsout = new DoubleMatrixDataset<String, String>();
			dsout.rawData = output1;
			dsout.colObjects = sampleNames;
			dsout.rowObjects = featureStr;
			dsout.recalculateHashMaps();

			dsout.save(outdir + "regionBasesCoveredGt20xPerRegion.txt");

			dsout = new DoubleMatrixDataset<String, String>();
			dsout.rawData = output2;
			dsout.colObjects = sampleNames;
			dsout.rowObjects = rowNames2;
			dsout.recalculateHashMaps();
			dsout.save(outdir + "sampleCoverageHistogramsMerged.txt");


//
//			String regionSummaryOutput = outdir + "regionOutput/";
//
//			Gpio.createDir(regionSummaryOutput);
//
//			for (Feature f : features) {
//				String filename = returnFileNameForRegion(f);
//				TextFile regionout = new TextFile(regionSummaryOutput + filename, TextFile.W);
//
//					String regionOutdir = sampleOutDir + "regions/";
//					String specificRegionFile = regionOutdir + filename;
//					// does this thing have a header
////					TextFile tf = new TextFile(specificRegionFile, TextFile.R);
////					String ln = tf.readLine();
////					while (ln != null) {
////						regionout.writeln(ln);
////						ln = tf.readLine();
////					}
////					tf.close();
//				}
//				regionout.close();
//			}


		}
	}

	private String returnFileNameForRegion(Feature f) {
		Chromosome c = f.getChromosome();
		int start = f.getStart();
		int stop = f.getStop();
		return c.getName() + "-" + start + "-" + stop + ".txt";
	}


	public void bamToBed(String bamfile, String outdir) throws IOException {

		BamFileReader reader = new BamFileReader(new File(bamfile));
		List<SAMReadGroupRecord> readGroups = reader.getReadGroups();
		HashMap<SAMReadGroupRecord, Integer> readgroupMap = new HashMap<SAMReadGroupRecord, Integer>();
		String[] samples = null;
		TextFile[] outfiles = null;
		if (readGroups.size() > 0) {
			samples = new String[readGroups.size()];
			outfiles = new TextFile[readGroups.size()];
			int rgctr = 0;
			for (SAMReadGroupRecord r : readGroups) {
				readgroupMap.put(r, rgctr);
				String sample = r.getSample();
				sample = sample.replaceAll("/", "-");
				samples[rgctr] = sample;
				outfiles[rgctr] = new TextFile(outdir + sample + ".txt", TextFile.W);
				rgctr++;
			}
		} else {
			System.err.println("WARNING: no readgroups found");
			samples = new String[1];
			samples[0] = bamfile;

		}

		ArrayList<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();
		AggregateFilter filter = new AggregateFilter(filters);

		filters.add(new NotPrimaryAlignmentFilter());
		filters.add(new FailsVendorQualityCheckFilter());
		filters.add(new DuplicateReadFilter());
		filters.add(new UnmappedReadFilter());
		filters.add(new MappingQualityUnavailableFilter());
		filters.add(new HCMappingQualityFilter());
		// ilters.add(new MalformedReadFilter());

		SAMRecordIterator it = reader.iterator();
		SAMRecord record = it.next();

		int[][] data = new int[samples.length][15];
		int[][] mapqDist = new int[samples.length][50];

		int ln = 0;
		int ln2 = 0;


		while (it.hasNext()) {

			SAMReadGroupRecord rg = record.getReadGroup();
			Integer sampleId = readgroupMap.get(rg);

			if (record.getDuplicateReadFlag()) {
// 1
				data[sampleId][0]++;
			}
			if (record.getFirstOfPairFlag()) {
// 2
				data[sampleId][1]++;
			}

			if (record.getMateNegativeStrandFlag()) {
// 3
				data[sampleId][2]++;
			}

			if (record.getMateUnmappedFlag()) {
// 4
				data[sampleId][3]++;
			}

			if (record.getNotPrimaryAlignmentFlag()) {
// 5
				data[sampleId][4]++;
			}

			if (record.getProperPairFlag()) {
// 6
				data[sampleId][5]++;
			}

			if (record.getReadFailsVendorQualityCheckFlag()) {
// 7
				data[sampleId][6]++;
			}

			if (record.getReadNegativeStrandFlag()) {
// 8
				data[sampleId][7]++;
			}

			if (record.getReadPairedFlag()) {
// 9
				data[sampleId][8]++;
			}

			if (record.getReadUnmappedFlag()) {
// 10
				data[sampleId][9]++;
			}

			if (record.getSecondOfPairFlag()) {
// 11
				data[sampleId][11]++;
			}

			if (filter.filterOut(record)) {
// 12
				data[sampleId][12]++; // nr of reads filtered out by HC filter
			}

			data[sampleId][13]++; // total nr of reads

			if (record.getFirstOfPairFlag() && record.getReferenceIndex() != null) {
				if (!record.getDuplicateReadFlag()) { // no dups
					if (!record.getNotPrimaryAlignmentFlag()) { // primary only
						if (record.getProperPairFlag()) { // proper pair
							if (record.getReadPairedFlag()) { // paired
								if (!record.getReadUnmappedFlag() && !record.getMateUnmappedFlag()) { // both mates mapped
									data[sampleId][14]++;
								}
							}
						}
					}
				}
			}

			String refname = record.getReferenceName();
			int len = record.getReadLength();
			int sta = record.getAlignmentStart();
			int sto = record.getAlignmentEnd();
			int insert = record.getInferredInsertSize();
			int mapq = record.getMappingQuality();

			if (!record.getReadUnmappedFlag()) {
				String outputStr = refname + "\t" + sta + "\t" + sto + "\t" + record.getFlags() + "\t" + mapq + "\t" + insert + "\t" + len;
				outfiles[sampleId].writeln(outputStr);
			}

			if (!filter.filterOut(record)) {
				if (mapq > 49) {
					mapq = 49;
				}
				mapqDist[sampleId][mapq]++;
			}

			record = it.next();
			ln++;
			if (ln % 1000000 == 0) {
				System.out.println(ln2 + "\t" + ln + " processed");
				ln2++;
			}
		}


		// chr sta sto dup passf

		for (int i = 0; i < outfiles.length; i++) {
			outfiles[i].close();
		}

		TextFile summary = new TextFile(outdir + "summary.txt", TextFile.W);

		String header = "-\tgetDuplicateReadFlag" +
				"\tgetFirstOfPairFlag" +
				"\tgetMateNegativeStrandFlag" +
				"\tgetMateUnmappedFlag" +
				"\tgetNotPrimaryAlignmentFlag" +
				"\tgetProperPairFlag" +
				"\tgetReadFailsVendorQualityCheckFlag" +
				"\tgetReadPairedFlag" +
				"\tgetReadUnmappedFlag" +
				"\tgetSecondOfPairFlag" +
				"\tHCfilterOut" +
				"\tTotalReads" +
				"\tpairsNoDupPrimaryProperPairPairedBothMapped";
		summary.writeln(header);
		for (int s = 0; s < data.length; s++) {
			String outputStr = samples[s];
			for (int d = 0; d < data[s].length; d++) {
				outputStr += "\t" + data[s][d];
			}
			summary.writeln(outputStr);
		}
		summary.close();

		TextFile summarymapq = new TextFile(outdir + "mapq.txt", TextFile.W);
		header = "mapq";
		for (int s = 0; s < data.length; s++) {
			header += "\t" + samples[s];
		}
		summarymapq.writeln(header);
		for (int d = 0; d < 50; d++) {
			String outln = d + "";
			for (int s = 0; s < data.length; s++) {
				outln += "\t" + mapqDist[s][d];
			}
			summarymapq.writeln(outln);
		}

		summarymapq.close();

	}

	public void correlateCoverageWithinRegionFiles(String regionFile, String inputdir, String output) throws IOException {
		ArrayList<Feature> regions = loadRegions(regionFile);

		TextFile outputF = new TextFile(output + "averageCorrelationPerRegion.txt", TextFile.W);

		for (Feature f : regions) {
			String filename = returnFileNameForRegion(f);
			TextFile tf = new TextFile(inputdir + filename, TextFile.R);

			ArrayList<double[]> cts = new ArrayList<double[]>();
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {

				double[] counts = new double[elems.length - 1];
				for (int i = 1; i < counts.length; i++) {
					counts[i - 1] = Double.parseDouble(elems[i]);
				}
				cts.add(counts);
				elems = tf.readLineElems(TextFile.tab);
			}

			tf.close();

			ArrayList<Double> correlations = new ArrayList<Double>();
			PearsonsCorrelation corr = new PearsonsCorrelation();
			for (int i = 0; i < cts.size(); i++) {
				double[] vals1 = cts.get(i);
				for (int j = i + 1; j < cts.size(); j++) {
					double[] vals2 = cts.get(j);
					double pearson = corr.correlation(vals1, vals2);
					correlations.add(pearson);
				}
			}

			double[] correlationsArr = Primitives.toPrimitiveArr(correlations.toArray(new Double[0]));
			double mean = Descriptives.mean(correlationsArr);
			double sd = Math.sqrt(Descriptives.variance(correlationsArr));
			outputF.writeln(f.getChromosome().getName() + "-" + f.getStart() + "-" + f.getStop() + "\t" + mean + "\t" + sd);
		}
		outputF.close();
	}

	public ArrayList<Feature> loadRegions(String regionfile) throws IOException {
		ArrayList<Feature> regions = new ArrayList<Feature>();

		TextFile featurefile = new TextFile(regionfile, TextFile.R);
		String[] felems = featurefile.readLineElems(TextFile.tab);
		while (felems != null) {

			if (felems.length >= 3) {
				Feature f = new Feature();
				Chromosome c = Chromosome.parseChr(felems[0]);
				f.setChromosome(c);
				f.setStart(Integer.parseInt(felems[1]));
				f.setStop(Integer.parseInt(felems[2]));
				regions.add(f);
			}

			felems = featurefile.readLineElems(TextFile.tab);
		}
		featurefile.close();


		// now sort them
		Collections.sort(regions, new FeatureComparator(false));
		return regions;
	}


	public void createCoveragePlots(String coverageFile, String sampleFile, String outdir) throws IOException {

		Gpio.createDir(outdir);


		HashMap<String, String> sampleToRun = null;
		HashMap<String, Color> runToColor = null;
		HashSet<String> runIds = null;
		if (sampleFile != null && Gpio.exists(sampleFile)) {

			TextFile sf = new TextFile(sampleFile, TextFile.R);
			sampleToRun = (HashMap<String, String>) sf.readAsHashMap(1, 0);
			sf.close();
			sf.open();
			runIds = (HashSet<String>) sf.readAsSet(0, TextFile.tab);
			sf.close();

			runToColor = new HashMap<String, Color>();

			for (String s : runIds) {
				runToColor.put(s, ColorGenerator.generate());
			}


		}

		// HashMap<String, String> sampleToRun = loadSampleRuns(sampleFile);

		Triple<String[], String[], double[][]> coverageObj = loadMatrix(coverageFile, true);
		String[] samples = coverageObj.getLeft();
		String[] regions = coverageObj.getMiddle();
		double[][] coverage = coverageObj.getRight();

		HashMap<String, Integer> rowIndex = index(coverageObj.getLeft());
		HashMap<String, Integer> colIndex = index(coverageObj.getMiddle());


		int maxBinCoverage = 300;
		int increment = 5;
		double[][] bins = new double[samples.length][maxBinCoverage / increment];
		int binNo = 0;

		double maxCoverage = 0;
		// for each sample,
		// for 5x, 10x, 15x etc, determine the number of regions
		// that pass this threshold
		for (int i = 0; i < maxBinCoverage; i += increment) {
			for (int s = 0; s < samples.length; s++) {
				for (int r = 0; r < regions.length; r++) {
					double d = coverage[s][r];
					if (d >= i) {
						bins[s][binNo]++;
					}
					if (d > maxCoverage) {
						maxCoverage = d;
					}
				}
			}
			binNo++;
		}

		// convert to frequencies
		double[] x = new double[bins[0].length];
		for (int s = 0; s < samples.length; s++) {
			for (int d = 0; d < bins[s].length; d++) {
				bins[s][d] /= regions.length;
				if (s == 0) {
					System.out.println(d + "\t" + bins[s][d]);
				}
				x[d] = (d) * increment;
			}
		}


		// plotVariantsUniqueIneachDataset
		makeCoveragePlot(samples, null, x, bins, runIds, outdir + "coveragePerSampleCumulative-", 100);

		int nrBins = (int) Math.ceil(maxBinCoverage / increment);
		double[][] bins2 = new double[samples.length][nrBins];

		// for each sample,
		// make a distribution of coverage over all regions
		for (int s = 0; s < samples.length; s++) {
			for (int r = 0; r < regions.length; r++) {
				double d = coverage[s][r];
				int bin = (int) Math.ceil(d / increment);
				if (bin >= nrBins) {
					bin = nrBins - 1;
				}
				bins2[s][bin]++;
			}
		}

		makeCoveragePlot(samples, null, x, bins2, runIds, outdir + "coveragePerSamplePerRegion-", 100);


		int threshold = 20;

		double[] bins3 = new double[100];

		// percentage of
		// y : proportion of samples
		// x : percent of regions covered

		for (int s = 0; s < samples.length; s++) {
			int nrHigherThenThreshold = 0;
			for (int r = 0; r < regions.length; r++) {
				if (coverage[s][r] > threshold) {
					nrHigherThenThreshold++;
				}
			}
			double perc = (double) nrHigherThenThreshold / regions.length;

			perc -= .8; // move > 80 to 0
			if (perc < 0) {
				perc = 0;
			}
			perc /= .20; // move back to percentage space


			int bin = (int) Math.floor(perc * 100);
			bins3[bin]++;
		}

		double[] xvals = new double[bins3.length];
		for (int bin = 0; bin < bins3.length; bin++) {
			bins3[bin] /= samples.length;
			xvals[bin] = 80 + (bin * 20 / 100); // bin 80 + (bin * 20/100)
		}

		Chart barchart = initBarChart("Coverage > " + threshold + "x", "Percent of targeted regions covered", "Proportion of samples", 1200, 800, true, 50);
		barchart.addSeries("Series1", xvals, bins3);
		barchart.getStyleManager().setLegendVisible(false);
		BitmapEncoder.saveBitmapWithDPI(barchart, outdir + "regionsWithCoverageLt" + threshold + "x", BitmapEncoder.BitmapFormat.PNG, 300);


		// determine average coverage per region
		TextFile out = new TextFile(outdir + "averageCoveragePerRegion.txt", TextFile.W);
		double[] plotVals = new double[regions.length];
		double[] plotValsX = new double[regions.length];
		for (int r = 0; r < regions.length; r++) {
			plotValsX[r] = r;
			double sum = 0;
			double[] vals = new double[samples.length];
			for (int s = 0; s < samples.length; s++) {
				sum += coverage[s][r];
				vals[s] = coverage[s][r];
			}
			sum /= samples.length;
			plotVals[r] = sum;
			double sd = JSci.maths.ArrayMath.standardDeviation(vals);
			out.writeln(regions[r] + "\t" + sum + "\t" + sd);
		}
		out.close();

		Arrays.sort(plotVals);

		barchart = initBarChart("Average coverage per targeted region", "Targeted region", "Average coverage (bp)", 1200, 800, true, 50);
		barchart.getStyleManager().setXAxisTicksVisible(false);
		barchart.addSeries("Series1", plotValsX, plotVals);
		barchart.getStyleManager().setLegendVisible(false);
		BitmapEncoder.saveBitmapWithDPI(barchart, outdir + "averageCoveragePerRegion", BitmapEncoder.BitmapFormat.PNG, 300);

		// determine average coverage per sample
		out = new TextFile(outdir + "averageCoveragePerSample.txt", TextFile.W);
		plotVals = new double[samples.length];
		plotValsX = new double[samples.length];
		for (int s = 0; s < samples.length; s++) {
			double sum = 0;
			plotValsX[s] = s;
			double[] vals = new double[regions.length];
			for (int r = 0; r < regions.length; r++) {
				sum += coverage[s][r];
				vals[s] = coverage[s][r];
			}
			sum /= regions.length;
			plotVals[s] = sum;
			double sd = JSci.maths.ArrayMath.standardDeviation(vals);
			out.writeln(samples[s] + "\t" + sum + "\t" + sd);
		}
		out.close();

		Arrays.sort(plotVals);

		barchart = initBarChart("Average coverage per sample", "Sample", "Average coverage (bp)", 1200, 800, true, 50);
		barchart.getStyleManager().setXAxisTicksVisible(false);
		barchart.addSeries("Series1", plotValsX, plotVals);
		barchart.getStyleManager().setLegendVisible(false);
		BitmapEncoder.saveBitmapWithDPI(barchart, outdir + "averageCoveragePerSample", BitmapEncoder.BitmapFormat.PNG, 300);

		System.out.println();

		for (int s = 0; s < samples.length; s++) {
			int count = 0;

			for (int r = 0; r < regions.length; r++) {
				if (coverage[s][r] > threshold) {
					count++;
				}
			}

			System.out.println(samples[s] + "\t" + count + "\t" + ((double) count / regions.length));

		}
	}

	private void makeCoveragePlot(String[] samples, HashMap<String, String> sampleToRun,
								  double[] x, double[][] bins, HashSet<String> runIds, String outdir, int ticksEvery) throws IOException {
		if (sampleToRun != null) {

			HashMap<String, Chart> chartsForRun = new HashMap<String, Chart>();
			for (int s = 0; s < samples.length; s++) {
				String sample = samples[s];
				Color seriesColor = ColorGenerator.generate();
				String seriesName = sample;

				String run = sampleToRun.get(sample);
				Chart chart = chartsForRun.get(run);
				if (chart == null) {
					chart = initLineChart("Coverage per sample", "Coverage (bp)", "Proportion of targeted regions", 1200, 800, true, ticksEvery);
					chart.getStyleManager().setXAxisMin(0);
				}

				double[] y = bins[s];
				Series series = chart.addSeries(seriesName, x, y);
				series.setLineColor(seriesColor);
				series.setMarker(SeriesMarker.NONE);
				series.setLineStyle(SeriesLineStyle.SOLID);

				chartsForRun.put(run, chart);
			}

			for (String run : runIds) {
				Chart chart = chartsForRun.get(run);
				BitmapEncoder.saveBitmapWithDPI(chart, outdir + run, BitmapEncoder.BitmapFormat.PNG, 300);
			}

		} else {
			Chart chart = initLineChart("Coverage per sample", "Coverage (bp)", "Proportion of targeted regions", 1200, 800, false, ticksEvery);
			chart.getStyleManager().setXAxisMin(0);
			for (int s = 0; s < samples.length; s++) {
				String sample = samples[s];
				Color seriesColor = Color.black;
				String seriesName = sample;


				double[] y = bins[s];
				Series series = chart.addSeries(seriesName, x, y);
				series.setLineColor(seriesColor);
				series.setMarker(SeriesMarker.NONE);
				series.setLineStyle(SeriesLineStyle.SOLID);
			}

			BitmapEncoder.saveBitmapWithDPI(chart, outdir + "overall", BitmapEncoder.BitmapFormat.PNG, 300);
		}
	}


	private Chart initLineChart(String title, String xaxis, String yaxis, int width, int height, boolean showLegend, int ticksEvery) {
		Chart chart = new ChartBuilder().width(width).height(height).build();

		chart.setChartTitle(title);
		chart.setXAxisTitle(xaxis);
		chart.setYAxisTitle(yaxis);

		chart.getStyleManager().setPlotBorderVisible(false);
		chart.getStyleManager().setPlotBackgroundColor(Color.WHITE);
		chart.getStyleManager().setPlotGridLinesColor(Color.WHITE);
		chart.getStyleManager().setPlotGridLinesVisible(false);
		chart.getStyleManager().setPlotPadding(20);

		chart.getStyleManager().setChartFontColor(Color.BLACK);
		chart.getStyleManager().setChartBackgroundColor(Color.WHITE);

		chart.getStyleManager().setChartTitleFont(new Font(Font.SANS_SERIF, Font.BOLD, 24));

		chart.getStyleManager().setLegendVisible(showLegend);
		chart.getStyleManager().setLegendFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));

		chart.getStyleManager().setLegendBackgroundColor(Color.white);
		chart.getStyleManager().setLegendBorderColor(Color.white);
		chart.getStyleManager().setLegendPosition(StyleManager.LegendPosition.OutsideE);

		chart.getStyleManager().setAxisTickPadding(10);
		chart.getStyleManager().setAxisTickMarkLength(5);
		chart.getStyleManager().setAxisTitleFont(new Font(Font.SANS_SERIF, Font.PLAIN, 18));
		chart.getStyleManager().setAxisTickLabelsFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));

		chart.getStyleManager().setXAxisTickMarkSpacingHint(ticksEvery);
		return chart;
	}

	public Chart initBarChart(String title, String xaxis, String yaxis, int width, int height, boolean showLegend, int ticksEvery) {
		Chart chart = new ChartBuilder().chartType(StyleManager.ChartType.Bar).width(width).height(height).build();

		chart.setChartTitle(title);
		chart.setXAxisTitle(xaxis);
		chart.setYAxisTitle(yaxis);

		chart.getStyleManager().setPlotBorderVisible(false);
		chart.getStyleManager().setPlotBackgroundColor(Color.WHITE);
		chart.getStyleManager().setPlotGridLinesColor(Color.WHITE);
		chart.getStyleManager().setPlotGridLinesVisible(false);
		chart.getStyleManager().setPlotPadding(20);

		chart.getStyleManager().setChartFontColor(Color.BLACK);
		chart.getStyleManager().setChartBackgroundColor(Color.WHITE);

		chart.getStyleManager().setChartTitleFont(new Font(Font.SANS_SERIF, Font.BOLD, 24));

		chart.getStyleManager().setLegendVisible(showLegend);
		chart.getStyleManager().setLegendFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));

		chart.getStyleManager().setLegendBackgroundColor(Color.white);
		chart.getStyleManager().setLegendBorderColor(Color.white);
		chart.getStyleManager().setLegendPosition(StyleManager.LegendPosition.OutsideE);

		chart.getStyleManager().setAxisTickPadding(10);
		chart.getStyleManager().setAxisTickMarkLength(5);
		chart.getStyleManager().setAxisTitleFont(new Font(Font.SANS_SERIF, Font.PLAIN, 18));
		chart.getStyleManager().setAxisTickLabelsFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));

		chart.getStyleManager().setXAxisTickMarkSpacingHint(ticksEvery);
		return chart;
	}

	private Triple<String[], String[], double[][]> loadMatrix(String f, boolean hasHeader) throws IOException {

		ArrayList<String> rows = new ArrayList<String>(); // samples
		ArrayList<String> cols = new ArrayList<String>(); // regions
		ArrayList<double[]> data = new ArrayList<double[]>();

		TextFile tf = new TextFile(f, TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);

		if (hasHeader) {
			for (int i = 1; i < elems.length; i++) {
				cols.add(elems[i - 1].intern());
			}
			elems = tf.readLineElems(TextFile.tab);
		} else {
			for (int i = 0; i < elems.length - 1; i++) {
				cols.add(("Col-" + (i + 1)).intern());
			}
		}

		HashSet<String> visitedSamples = new HashSet<String>();
		while (elems != null) {
			int nrCols = elems.length - 1;
			double[] datarow = new double[nrCols];
			String sample = elems[0].intern();
			String finalSample = sample;

			// prevent duplicates

			if (visitedSamples.contains(sample)) {
				int q = 1;
				while (visitedSamples.contains(finalSample)) {
					finalSample = sample + "_" + q;
					q++;
				}
				System.out.println("duplicate sample:\t" + sample + "\treplaced with:\t" + finalSample);
				sample = finalSample;
			}
			visitedSamples.add(sample);

			rows.add(sample);
			for (int i = 1; i < elems.length; i++) {
				datarow[i - 1] = Double.parseDouble(elems[i]);
			}
			data.add(datarow);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		double[][] datad = new double[rows.size()][cols.size()];
		int rc = 0;
		for (double[] row : data) {
			datad[rc] = row;
			rc++;
		}

		return new Triple<String[], String[], double[][]>(
				rows.toArray(new String[0]),
				cols.toArray(new String[0]),
				datad
		);

	}

	private HashMap<String, Integer> index(String[] s) {
		int c = 0;
		HashMap<String, Integer> index = new HashMap<String, Integer>();
		for (String s1 : s) {
			index.put(s1, c);
			c++;
		}
		return index;
	}

}


