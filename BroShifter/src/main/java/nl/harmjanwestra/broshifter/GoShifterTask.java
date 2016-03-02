package nl.harmjanwestra.broshifter;


import nl.harmjanwestra.broshifter.CLI.BroShifterOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.*;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.math.stats.ZScores;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.Callable;

/**
 * Created by hwestra on 11/17/15.
 * Parallelize broshifting by running multiple annotations in different threads
 * Returns: all output for a certain (combination of) annotation(s)
 */
public class GoShifterTask implements Callable<Pair<String, ArrayList<String>>> {

	boolean DEBUG = false;
	private int threadNum;
	private String annotation1;
	private String annotation2;
	private BroShifterOptions options;


	public GoShifterTask(int threadNum,
						 String annotation1,
						 String annotation2,
						 BroShifterOptions options) {
		this.threadNum = threadNum;
		this.annotation1 = annotation1;
		this.annotation2 = annotation2;
		this.options = options;
	}

	public GoShifterTask() {

	}

	public Pair<String, ArrayList<String>> call() throws IOException {
		// load annotations
		AnnotationLoader loader = new AnnotationLoader();
		System.out.println("Thread " + threadNum + " | Testing: " + annotation1);
		Track annotation1Track = loader.loadAnnotations(annotation1, options.usePeakCenter, options.bpToExtendAnnotation, true, null);
		String annotation1name = new File(annotation1).getName();
		Track annotation2Track = null;
		String annotation2name = null;
		if (annotation2 != null) {
			annotation2Track = loader.loadAnnotations(annotation2, options.usePeakCenter, options.bpToExtendAnnotation, true, null);
			annotation2name = new File(annotation2).getName();
		}


		int medianAnnotationLength = determineMedianAnnotationLength(annotation1Track);

		BedFileReader bf = new BedFileReader();

		// load snp data
		ArrayList<Pair<SNPFeature, ArrayList<SNPFeature>>> allSNPs = loadSNPs(options.snpfile);

		// define regions using snp data
		ArrayList<Feature> regions = new ArrayList<>();

		for (int fctr = 0; fctr < allSNPs.size(); fctr++) {
			// load proxies
			Pair<SNPFeature, ArrayList<SNPFeature>> proxies = allSNPs.get(fctr);
			Feature region = defineRegion(proxies, medianAnnotationLength);
			regions.add(region);
		}


		System.out.println("Thread " + threadNum + " | " + allSNPs.size() + " SNPs in: " + options.snpfile);

		// total nr of overlapping variants

		ArrayList<String> outputLines = new ArrayList<>();

		// iterate all regions
		int nrOfLociWithOverlap = 0;
		boolean[][] overlapMatrix = new boolean[allSNPs.size()][options.nrIterations];

		int totalNrOfVariants = 0;
		int totalNrOfOverlappingVariants = 0;

		BroShifterTask broshifter = new BroShifterTask();

		for (int fctr = 0; fctr < allSNPs.size(); fctr++) {
			Feature region = regions.get(fctr);

			String origRegion = region.toString();

			// read posteriors
			double locusScore = 0;
			int annotation2size = 0;

			Pair<SNPFeature, ArrayList<SNPFeature>> proxies = allSNPs.get(fctr);
			ArrayList<SNPFeature> snps = proxies.getRight();


			if (snps.size() > 0) {

				if (options.trimRegions) {
					region = broshifter.trimRegion(region, snps);
				}

				// get subset of annotations
				Track subsetOfAnnotation1 = annotation1Track.getSubset(region.getChromosome(), region.getStart(), region.getStop());

				int nrOverlapping = 0;

				if (annotation2 != null) {
					// perform conditional analysis
					Track subsetOfAnnotation2 = annotation2Track.getSubset(region.getChromosome(), region.getStart(), region.getStop());
					annotation2size = subsetOfAnnotation2.getAllFeatures().size();
					// split up the region in Y and Y-hat
					Pair<Triple<Feature, ArrayList<Feature>, ArrayList<SNPFeature>>, Triple<Feature, ArrayList<Feature>, ArrayList<SNPFeature>>>
							conditionalRegions = broshifter.split(region, snps, subsetOfAnnotation1, subsetOfAnnotation2);


					Triple<Feature, ArrayList<Feature>, ArrayList<SNPFeature>> y = conditionalRegions.getLeft();
					Triple<Feature, ArrayList<Feature>, ArrayList<SNPFeature>> ybar = conditionalRegions.getRight();

					Pair<Double, Integer> signalUncond = broshifter.getOverlap(subsetOfAnnotation1, snps);
					int nrOverlapUncond1 = signalUncond.getRight();

					Pair<Double, Integer> signalUncond2 = broshifter.getOverlap(subsetOfAnnotation2, snps);
					int nrOverlapUncond2 = signalUncond2.getRight();


					Track ty = new Track("Y");
					ty.addFeatures(y.getMiddle());

					Track tybar = new Track("Ybar");
					tybar.addFeatures(ybar.getMiddle());

					Pair<Double, Integer> signal1 = broshifter.getOverlap(ty, y.getRight());
					Pair<Double, Integer> signal2 = broshifter.getOverlap(tybar, ybar.getRight());

					if (DEBUG) {
						System.out.println(y.getRight().size() + "\t" + ybar.getRight().size() + "\t" + nrOverlapUncond1 + "\t" + nrOverlapUncond2);
						System.out.println(signal1.getRight() + "\t" + signal2.getRight());
					}

					nrOverlapping = signal1.getRight() + signal2.getRight();
					if (nrOverlapping > 0) {
						nrOfLociWithOverlap++;
					}

					// locus score == % of iterations with >= 1 overlap
					for (int i = 0; i < options.nrIterations; i++) {
						broshifter.shift(y.getRight(), y.getLeft());
						broshifter.shift(ybar.getRight(), ybar.getLeft());
						Pair<Double, Integer> null1 = broshifter.getOverlap(ty, y.getRight());
						Pair<Double, Integer> null2 = broshifter.getOverlap(tybar, ybar.getRight());

						int nrOverlappingNull = null1.getRight() + null2.getRight();
						if (nrOverlappingNull > 0) {
							overlapMatrix[fctr][i] = true;
							locusScore++;
						}
					}
				} else {
					// run unconditional analysis
					// determine sigmaPosterior
					Pair<Double, Integer> signal = broshifter.getOverlap(subsetOfAnnotation1, snps);
					nrOverlapping = signal.getRight();
					if (nrOverlapping > 0) {
						nrOfLociWithOverlap++;
					}
					for (int i = 0; i < options.nrIterations; i++) {
						broshifter.shift(snps, region);
						Pair<Double, Integer> output = broshifter.getOverlap(subsetOfAnnotation1, snps);
						int nrNullOverlapping = output.getRight();
						if (nrNullOverlapping > 0) {
							overlapMatrix[fctr][i] = true;
							locusScore++;
						}
					}

				}

				locusScore /= options.nrIterations;

				String regionName = region.toString();
				totalNrOfVariants += snps.size();
				totalNrOfOverlappingVariants += nrOverlapping;

				StringBuilder builder = new StringBuilder(500);

				if (nrOverlapping == 0) {
					builder.append(locusScore)
							.append("\t").append(regionName)
							.append("\t").append(origRegion)
							.append("\t").append(nrOverlapping)
							.append("\t").append(snps.size())
							.append("\t").append(0)
							.append("\t").append(annotation1name)
							.append("\t").append(subsetOfAnnotation1.getSize());
				} else {
					builder.append(locusScore)
							.append("\t").append(regionName)
							.append("\t").append(origRegion)
							.append("\t").append(nrOverlapping)
							.append("\t").append(snps.size())
							.append("\t").append(((double) nrOverlapping / snps.size()))
							.append("\t").append(annotation1name)
							.append("\t").append(subsetOfAnnotation1.getAllFeatures().size());
				}
				if (annotation2 != null) {
					builder.append("\t").append(annotation2name)
							.append("\t").append(annotation2size);
				}
				outputLines.add(builder.toString());
			}

			if (fctr % 10 == 0) {
				System.out.println("Thread " + threadNum + " | " + fctr + " out of " + regions.size() + " regions processed.");
			}
		}

		double pval = 0;
		int[] nrLociWithOverlap = new int[options.nrIterations];
		for (int i = 0; i < options.nrIterations; i++) {
			int nrOverlappingnull = 0;
			for (int fctr = 0; fctr < allSNPs.size(); fctr++) {
				if (overlapMatrix[fctr][i]) {
					nrOverlappingnull++;
				}
			}
			nrLociWithOverlap[i] = nrOverlappingnull;
		}
		Arrays.sort(nrLociWithOverlap);

		// count nr of iterations with fewer nr loci of overlap


		System.out.println("Thread " + threadNum + " | " + regions.size() + " out of " + regions.size() + " regions processed.");

		StringBuilder builder = new StringBuilder(500);
		builder.append(pval).
				append("\t").append(nrOfLociWithOverlap).
				append("\t").append(meanNrOfLociWithOverlapNull).
				append("\t").append(totalNrOfOverlappingVariants).
				append("\t").append(totalNrOfVariants).
				append("\t").append(((double) totalNrOfOverlappingVariants / totalNrOfVariants)).
				append("\t").append(annotation1name);
		if (annotation2 != null) {
			builder.append("\t").append(annotation2name);
		}


		System.out.println("Thread " + threadNum + " | done testing.");
		return new Pair<>(builder.toString(), outputLines);
	}

	private Feature defineRegion(Pair<SNPFeature, ArrayList<SNPFeature>> snpFeature, int medianAnnotationLength) {
		int maxPos = 0;
		int minPos = Integer.MAX_VALUE;

		for (SNPFeature feature : snpFeature.getRight()) {

			int pos = feature.getStart();
			if (pos > maxPos) {
				maxPos = pos;
			}
			if (pos < minPos) {
				minPos = pos;
			}
		}

		if (minPos - medianAnnotationLength > 0) {
			minPos -= medianAnnotationLength;
		} else {
			minPos = 0;
		}

		if (maxPos + medianAnnotationLength > Integer.MAX_VALUE) {
			maxPos = Integer.MAX_VALUE;
		} else {
			maxPos += medianAnnotationLength;
		}

		Feature output = new Feature();
		output.setName(snpFeature.getLeft().getName());
		output.setChromosome(snpFeature.getLeft().getChromosome());
		output.setStart(minPos);
		output.setStop(maxPos);


		return output;
	}

	private ArrayList<SNPFeature> loadSNPs(String s) throws IOException {
		return null;
	}

	private int determineMedianAnnotationLength(Track annotation1Track) {
		return 0;
	}

}
