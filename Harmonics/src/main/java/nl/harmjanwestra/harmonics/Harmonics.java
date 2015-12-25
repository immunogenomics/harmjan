/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.harmonics;

import JSci.maths.ArrayMath;
import com.itextpdf.text.DocumentException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import nl.harmjanwestra.harmonics.smoothing.AverageSmoothingFunction;
import nl.harmjanwestra.harmonics.smoothing.Smoother;
import nl.harmjanwestra.utilities.bamfile.BamFileReader;
import nl.harmjanwestra.utilities.bamfile.QueryableMergingSamRecordIterator;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.coverage.Coverage;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureMerger;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.features.Strand;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.CoveragePanel;
import nl.harmjanwestra.utilities.graphics.panels.GenePanel;
import nl.harmjanwestra.utilities.graphics.panels.HistogramPanel;
import nl.harmjanwestra.utilities.gtf.GTFAnnotation;
import org.apache.commons.math3.distribution.PoissonDistribution;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umontreal.iro.lecuyer.probdist.NegativeBinomialDist;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * @author hwestra
 */
public class Harmonics {

	public void run(String bedfile,
	                String[] bamFileLocs,
	                String[] sampleNames,
	                int posShift,
	                int negShift,
	                int maxDistSize,
	                int smoothingWindow,
	                String outputdir) throws IOException {

		// readAsTrack bedfile
		BedFileReader bedFileReader = new BedFileReader();
		ArrayList<Feature> regions = bedFileReader.readAsList(bedfile);
		FeatureMerger featureMerger = new FeatureMerger();
		regions = featureMerger.merge(regions, false);

		// initialize Bamfile readers
		// TODO: replace this with a more useful iterator, that is also able to run queries
		// TODO: split regions up in 1mb windows.
		ArrayList<BamFileReader> readers = new ArrayList<BamFileReader>();

		HashMap<String, Integer> sampleToPosition = new HashMap<String, Integer>();
		int sampleCounter = 0;

		ArrayList<String> allSampleNames = new ArrayList<String>();

		for (int bamfileNr = 0; bamfileNr < bamFileLocs.length; bamfileNr++) {
			String f = bamFileLocs[bamfileNr];

			String[] felems = f.split("/");
			String filename = felems[felems.length - 1];

			BamFileReader reader = new BamFileReader(new File(f), true);
			SAMFileHeader header = reader.getHeader();
			List<SAMReadGroupRecord> readgroups = header.getReadGroups();
			if (readgroups.isEmpty()) {
				System.out.println("Adding: " + filename);
				header.setAttribute("filename", filename);
				sampleToPosition.put(filename, sampleCounter);
				allSampleNames.add(filename);
				sampleCounter++;
			} else {
				// get readgroups
				for (SAMReadGroupRecord rg : readgroups) {
					sampleToPosition.put(rg.getId(), sampleCounter);
					allSampleNames.add(rg.getId());
					sampleCounter++;
				}
			}
			readers.add(reader);
		}

		QueryableMergingSamRecordIterator it = new QueryableMergingSamRecordIterator(readers, regions.get(0));
//		while (it.hasNext()) {
//			it.next().getFileSource();
//		}

		Coverage c = new Coverage();
		c.setReadGroupMap(sampleToPosition);

		String gtf = "/Data/Annotation/UCSC/genes.gtf";
		GTFAnnotation annot = new GTFAnnotation(gtf);
		String plotOut = outputdir + "plots";

		TextFile tmpOut = new TextFile(outputdir + "Coverage.txt", TextFile.W);

		for (Feature region : regions) {

			// TODO: query-able merging iterator

			c.calculate(it, region, posShift, negShift); // calculate coverage for region
//			int[][][] coverage = c.getCoverageStranded();
			int[][] unstranded = c.getCoverageUnStranded();

			// smooth the data
			Smoother s = new Smoother(new AverageSmoothingFunction());

			//	int[][] smoothedCoverage = s.smoothMatrix(unstranded, smoothingWindow);
			int[] overallCoverage = sumRows(unstranded);

			// make readAsTrack coverage distributions

//			// TODO: plot the distributions
//			// fit a distribution model (guassian?)
//			for (int i = 0; i < distributions.length; i++) {
//				fitPoisson(distributions[i], smoothedCoverage[i]);
//			}

			outputdir = Gpio.formatAsDirectory(outputdir);

			Gpio.createDir(outputdir);
			Grid grid = new Grid(640, 480, unstranded.length + 2, 3, 100, 100);

			Grid gridHist = new Grid(640 * 3, 480, unstranded.length + 1, 1, 100, 100);

			GenePanel genePanel = new GenePanel(1, 2);
			genePanel.setRegion(region);

			TreeSet<Gene> genes = annot.getGeneTree();
			Gene geneStart = new Gene("", region.getChromosome(), Strand.POS, region.getStart(), region.getStart());
			Gene geneStop = new Gene("", region.getChromosome(), Strand.NEG, region.getStop(), region.getStop());
			SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);
			genePanel.setGenes(overlappingGenes);
			grid.addPanel(genePanel, 0, 0);

			for (int i = 0; i < unstranded.length; i++) {
//			for (int i = 0; i < 2; i++) {
				System.out.println();
				double mean = ArrayMath.mean(unstranded[i]);

				System.out.println("Fitting NB: " + i + "\tMean: " + mean);


				Triple<Pair<Double, Double>, Pair<double[], double[]>, Double> fitted = fitNegativeBinomial(unstranded[i], maxDistSize);
				Pair<Double, Double> params = fitted.getLeft();
				Pair<double[], double[]> dists = fitted.getMiddle();
				double kl = fitted.getRight();

				System.out.println("params: n" + params.getLeft() + "\tp" + params.getRight());
				System.out.println("kl: " + kl);

//				int plotwidth = 1000;
//				int plotheight = 1000;
//				int margin = 100;


				double[][] distArr = new double[3][0];
				distArr[0] = dists.getLeft();
				distArr[1] = dists.getLeft();
				distArr[2] = dists.getRight();

				HistogramPanel.PLOTTYPE[] types = new HistogramPanel.PLOTTYPE[3];
				types[0] = HistogramPanel.PLOTTYPE.BAR;
				types[1] = HistogramPanel.PLOTTYPE.POLY;
				types[2] = HistogramPanel.PLOTTYPE.POLY;

//				String outfilename = plotOut + allSampleNames.get(i) + ".pdf";
//				try {
//					HistogramPlot plot = new HistogramPlot(outfilename, plotwidth, plotheight);
//					plot.setMargin(margin);
//					plot.setData(distArr, types);
//					plot.draw();
//				} catch (DocumentException e) {
//					e.printStackTrace();
//				}

				// add panel for coverage.
				CoveragePanel coveragePanel = new CoveragePanel(1, 2);
				coveragePanel.setRegion(region);
				coveragePanel.setCoverageData(unstranded[i]);
//				coveragePanel.setReads(c.getReadPositionsPerSample().get(i));


				grid.addPanel(coveragePanel, i + 1, 0);

				HistogramPanel histogram = new HistogramPanel(1, 1);
				histogram.setData(distArr);
				histogram.setPlotTypes(types);

				grid.addPanel(histogram, i + 1, 2);

				// TODO: plot the data for all samples

				StringBuilder lnout = new StringBuilder();
				lnout.append("Sample_" + i);
				for (int q = 0; q < unstranded[i].length; q++) {
					lnout.append("\t").append(unstranded[i][q]);
				}
				tmpOut.writeln(lnout.toString());

				// split dist and fit two negative binomials
				Pair<Pair<double[], double[]>, Pair<double[], double[]>> fitted2 = fitTwoNegativeBinomial(unstranded[i], maxDistSize * 4, 0.05);

				HistogramPanel histogram2 = new HistogramPanel(1, 1);
				double[][] dataForHist = new double[2][0];
				dataForHist[0] = fitted2.getLeft().getLeft();
				dataForHist[1] = fitted2.getRight().getLeft();

//				dataForHist[2] = fitted2.getLeft().getRight();
//				dataForHist[3] = fitted2.getRight().getRight();

				histogram2.setData(dataForHist);

				HistogramPanel.PLOTTYPE[] types2 = new HistogramPanel.PLOTTYPE[4];
				types2[0] = HistogramPanel.PLOTTYPE.BAR;
				types2[1] = HistogramPanel.PLOTTYPE.BAR;

//				types2[2] = HistogramPanel.PLOTTYPE.POLY;
//				types2[3] = HistogramPanel.PLOTTYPE.POLY;

				gridHist.addPanel(histogram2, i, 0);
			}


			System.out.println("");
			System.out.println("Fitting poisson");
			System.out.println(overallCoverage.length + " length of array");
			CoveragePanel coveragePanel = new CoveragePanel(1, 2);
			coveragePanel.setRegion(region);
			coveragePanel.setCoverageData(overallCoverage);
			grid.addPanel(coveragePanel, unstranded.length + 1, 0);

			Triple<Pair<Double, Double>, Pair<double[], double[]>, Double> fitted = fitNegativeBinomial(overallCoverage, maxDistSize);

			//double mean = fitted.getLeft();
			Pair<double[], double[]> dists = fitted.getMiddle();
			double kl = fitted.getRight();

			// System.out.println("params: mean: " + mean);
			System.out.println("kl: " + kl);

//				int plotwidth = 1000;
//				int plotheight = 1000;
//				int margin = 100;

			double[][] distArr = new double[3][0];
			distArr[0] = dists.getLeft();
			distArr[1] = dists.getLeft();
			distArr[2] = dists.getRight();

			HistogramPanel.PLOTTYPE[] types = new HistogramPanel.PLOTTYPE[3];
			types[0] = HistogramPanel.PLOTTYPE.BAR;
			types[1] = HistogramPanel.PLOTTYPE.POLY;
			types[2] = HistogramPanel.PLOTTYPE.POLY;

			HistogramPanel histogram = new HistogramPanel(1, 1);
			histogram.setData(distArr);
			histogram.setPlotTypes(types);

			grid.addPanel(histogram, unstranded.length + 1, 2);

			StringBuilder lnout = new StringBuilder();
			lnout.append("Sample_Sum");
			for (int q = 0; q < overallCoverage.length; q++) {
				lnout.append("\t").append(overallCoverage[q]);
			}
			tmpOut.writeln(lnout.toString());


			// split dist and fit two negative binomials
			Pair<Pair<double[], double[]>, Pair<double[], double[]>> fitted2 = fitTwoNegativeBinomial(overallCoverage, maxDistSize * 4, 0.05);

			HistogramPanel histogram2 = new HistogramPanel(1, 1);
			double[][] dataForHist = new double[2][0];
			dataForHist[0] = fitted2.getLeft().getLeft();
			dataForHist[1] = fitted2.getRight().getLeft();

//			dataForHist[2] = fitted2.getLeft().getRight();
//			dataForHist[3] = fitted2.getRight().getRight();

			histogram2.setData(dataForHist);

			HistogramPanel.PLOTTYPE[] types2 = new HistogramPanel.PLOTTYPE[4];
			types2[0] = HistogramPanel.PLOTTYPE.BAR;
			types2[1] = HistogramPanel.PLOTTYPE.BAR;

//			types2[2] = HistogramPanel.PLOTTYPE.POLY;
//			types2[3] = HistogramPanel.PLOTTYPE.POLY;

			gridHist.addPanel(histogram2, unstranded.length, 0);


			try {
				grid.draw(plotOut + "test.pdf", DefaultGraphics.Output.PDF);
				gridHist.draw(plotOut + "hist.pdf", DefaultGraphics.Output.PDF);
			} catch (DocumentException e) {
				e.printStackTrace();
			}


			// assign p-values to each position (also for overall coverage)
			// meta-analyze the p-values?
			// assign peak regions at certain p-value cutoff


			// plot both smoothed and unsmoothed data
			// also plot the peak regions per sample
		}
		tmpOut.close();
	}

	private Pair<Pair<double[], double[]>, Pair<double[], double[]>> fitTwoNegativeBinomial(int[] observations, int maxDistSize, double threshold) {
		Pair<int[], int[]> splitObservations = splitObservations(observations, threshold);

		int[] above = splitObservations.getRight();
		int[] below = splitObservations.getLeft();

//		NegativeBinomialDist nbAbove = NegativeBinomialDist.getInstanceFromMLE(above, above.length);
//		NegativeBinomialDist nbBelow = NegativeBinomialDist.getInstanceFromMLE(below, below.length);

		double[] probDistAbove = new double[maxDistSize];
		double[] probDistBelow = new double[maxDistSize];
		int[] obsHistAbove = makeDistribution(above, maxDistSize);
		long sumAbove = sum(obsHistAbove);

		int[] obsHistBelow = makeDistribution(below, maxDistSize);
		long sumBelow = sum(obsHistBelow);

		double[] obsDistAbove = new double[maxDistSize];
		double[] obsDistBelow = new double[maxDistSize];

		for (int i = 0; i < obsHistAbove.length; i++) {
			obsDistAbove[i] = (double) obsHistAbove[i] / sumAbove;
//			probDistAbove[i] = nbAbove.prob(i);
		}

		for (int i = 0; i < obsDistBelow.length; i++) {
			obsDistBelow[i] = (double) obsHistBelow[i] / sumBelow;
//			probDistBelow[i] = nbBelow.prob(i);
		}


		return new Pair<Pair<double[], double[]>, Pair<double[], double[]>>(
				new Pair<double[], double[]>(obsDistAbove, probDistAbove),
				new Pair<double[], double[]>(obsDistBelow, probDistBelow)
		);
	}

	// split observations in upper and lower $threshold%.
	private Pair<int[], int[]> splitObservations(int[] observations, double threshold) {

		int[] observationsSorted = new int[observations.length];
		System.arraycopy(observations, 0, observationsSorted, 0, observations.length);
		Arrays.sort(observationsSorted);


		double max = observationsSorted[observationsSorted.length - 1];

		int[] obsHist = makeDistribution(observations, (int) max);
		double sumObs = sum(obsHist);

		double probSum = 0;
		int cutoff = 0;

		int nrObs = 0;
		double upperthres = 1 - threshold;

		while (probSum < upperthres) {
			probSum += (double) obsHist[cutoff] / sumObs;
			nrObs += obsHist[cutoff];
			cutoff++;
		}


		int[] aboveThreshold = new int[observations.length - nrObs];
		int[] belowThreshold = new int[nrObs];

		System.out.println(nrObs + "\t" + observations.length);

		System.arraycopy(observationsSorted, 0, belowThreshold, 0, nrObs);
		System.arraycopy(observationsSorted, nrObs, aboveThreshold, 0, aboveThreshold.length);


		return new Pair<int[], int[]>(belowThreshold, aboveThreshold);

//		return new Pair<int[], int[]>(Primitives.toPrimitiveArr(belowThreshold.toArray(new Integer[0])), Primitives.toPrimitiveArr(aboveThreshold.toArray(new Integer[0])));
	}


	private Triple<Double, Pair<double[], double[]>, Double> fitPoisson(int[] observations, int maxDistSize) {
		double mean = JSci.maths.ArrayMath.mean(observations);
		PoissonDistribution poissonDistribution = new PoissonDistribution(mean);

		// convert to frequency dist
		int[] obsHist = makeDistribution(observations, maxDistSize);
		long sum = sum(obsHist);

		// for each bin, get expected from distribution
		double[] probDist = new double[maxDistSize];
		double[] obsDist = new double[maxDistSize];
		for (int i = 0; i < obsDist.length; i++) {
			obsDist[i] = (double) obsHist[i] / sum;
			probDist[i] = poissonDistribution.probability(i);
		}

		double kl = kldivergence(probDist, obsDist);

		double cumulative = 0;


		return new Triple<Double, Pair<double[], double[]>, Double>(mean, new Pair<double[], double[]>(obsDist, probDist), kl);


	}

	private double kldivergence(double[] p, double[] q) {
		double kl = 0;
		for (int i = 0; i < p.length; i++) {
			if (p[i] != 0d) {
				kl += (p[i] * Math.log(p[i] / q[i]));
			}
		}
		return kl;
	}

	private long sum(int[] data) {
		long sum = 0;
		for (int i = 0; i < data.length; i++) {
			sum += data[i];
		}
		return sum;
	}

	private double[] toFreqDist(int[] data) {
		double[] distout = new double[data.length];

		long sum = sum(data);
		for (int i = 0; i < data.length; i++) {
			distout[i] = (double) data[i] / sum;
		}
		return distout;
	}

	private Triple<Pair<Double, Double>, Pair<double[], double[]>, Double> fitNegativeBinomial(int[] observations, int maxDistSize) {


		double mean = ArrayMath.mean(observations);
		double sd = ArrayMath.standardDeviation(observations);
		double nest = (mean * mean) / ((sd * sd) - mean);

		System.out.println(nest);

		NegativeBinomialDist umontreal = NegativeBinomialDist.getInstanceFromMLE(observations, observations.length);


		double p = umontreal.getP();
		double n = umontreal.getN();

		System.out.println("n " + n + "\tp " + p + "\tnest " + nest);


		Pair<Double, Double> params = new Pair<Double, Double>(n, p);

		double[] probDist = new double[maxDistSize];
		int[] obsHist = makeDistribution(observations, maxDistSize);
		long sum = sum(obsHist);


		double[] obsDist = new double[maxDistSize];
		for (int i = 0; i < obsDist.length; i++) {
			obsDist[i] = (double) obsHist[i] / sum;
			probDist[i] = umontreal.prob(i);
		}


		Pair<double[], double[]> dists = new Pair<double[], double[]>(obsDist, probDist);

		double kl = kldivergence(probDist, obsDist);


		Triple<Pair<Double, Double>, Pair<double[], double[]>, Double> output = new Triple<Pair<Double, Double>, Pair<double[], double[]>, Double>(params, dists, kl);

		return output;
	}

	public int[] makeDistribution(int[] data, int maxDistSize) {
		int[] dist = new int[maxDistSize];
		for (int i = 0; i < data.length; i++) {
			int bin = data[i];
			if (bin >= maxDistSize) {
				bin = maxDistSize - 1;
			}
			dist[bin]++;
		}
		return dist;
	}

	public int[][] makeDistributions(int[][] data, int maxDistSize) {
		int[][] distributions = new int[data.length][0];
		for (int i = 0; i < data.length; i++) {
			distributions[i] = makeDistribution(data[i], maxDistSize);
		}
		return distributions;
	}

	public int[] sumRows(int[][] data) {
		int[] output = new int[data[0].length];

		for (int j = 0; j < data[0].length; j++) {
			for (int i = 0; i < data.length; i++) {
				output[j] += data[i][j];
			}
		}

		return output;
	}
}
