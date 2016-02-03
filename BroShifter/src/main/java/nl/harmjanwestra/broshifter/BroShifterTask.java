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
public class BroShifterTask implements Callable<Pair<String, ArrayList<String>>> {

	private final int threadNum;

	private final String annotation1;
	private final String annotation2;
	private final BroShifterOptions options;
	boolean DEBUG = false;


	public BroShifterTask(int threadNum,
						  String annotation1,
						  String annotation2,
						  BroShifterOptions options) {
		this.threadNum = threadNum;
		this.annotation1 = annotation1;
		this.annotation2 = annotation2;
		this.options = options;
	}

	public Pair<String, ArrayList<String>> call() throws IOException {
		// load bed regions to test
		BedFileReader bf = new BedFileReader();
		ArrayList<Feature> regions = bf.readAsList(options.regionFile);
		System.out.println("Thread " + threadNum + " | " + regions.size() + " regions in: " + options.regionFile);

		// load annotations
		AnnotationLoader loader = new AnnotationLoader();
		System.out.println("Thread " + threadNum + " | Testing: " + annotation1);
		Track annotation1Track = loader.loadAnnotations(annotation1, options.usePeakCenter, options.bpToExtendAnnotation, true);
		String annotation1name = new File(annotation1).getName();
		Track annotation2Track = null;
		String annotation2name = null;
		if (annotation2 != null) {
			annotation2Track = loader.loadAnnotations(annotation2, options.usePeakCenter, options.bpToExtendAnnotation, true);
			annotation2name = new File(annotation2).getName();
		}

		// keep track of some scores
		// sum of overlapping posteriors overall loci
		double sigmaSigmaPosterior = 0;
		double[] sigmaSigmaPosteriorNull = new double[options.nrIterations];

		// number of loci with >= 1 overlapping variant
		double sigmaLociWithAtLeastOneOverlap = 0;
		double[] sigmaLociWithAtLeastOneOverlapNull = new double[options.nrIterations];

		// total nr of overlapping variants
		int totalNrOfOverlappingVariants = 0;
		int totalNrOfVariants = 0;

		// TODO: replace this with a linkedblocking (FIFO) queue to save some memory
		ArrayList<String> outputLines = new ArrayList<>();

		// iterate all regions
		for (int fctr = 0; fctr < regions.size(); fctr++) {
			Feature region = regions.get(fctr);

			String origRegion = region.toString();

			// read posteriors
			String regionStr = region.toString();
			double locusScore = 0;

			int annotation2size = 0;
			if (Gpio.exists(options.posteriorFile) && !region.getChromosome().equals(Chromosome.X)) {

				ArrayList<SNPFeature> snps = readPosteriors(options.posteriorFile, region);

				if (snps.size() > 0) {

					if (options.trimRegions) {
						region = trimRegion(region, snps);
					}

					// get subset of annotations
					Track subsetOfAnnotation1 = annotation1Track.getSubset(region.getChromosome(), region.getStart(), region.getStop());

//					String ln = "";
					double sigmaPosterior = 0;
					double[] sigmaPosteriorNull = new double[options.nrIterations];
					int nrOverlapping = 0;

					if (annotation2 != null) {
						Track subsetOfAnnotation2 = annotation2Track.getSubset(region.getChromosome(), region.getStart(), region.getStop());
						annotation2size = subsetOfAnnotation2.getAllFeatures().size();
						// split up the region in Y and Y-hat
						Pair<Triple<Feature, ArrayList<Feature>, ArrayList<SNPFeature>>, Triple<Feature, ArrayList<Feature>, ArrayList<SNPFeature>>>
								conditionalRegions = split(region, snps, subsetOfAnnotation1, subsetOfAnnotation2);


						Triple<Feature, ArrayList<Feature>, ArrayList<SNPFeature>> y = conditionalRegions.getLeft();
						Triple<Feature, ArrayList<Feature>, ArrayList<SNPFeature>> ybar = conditionalRegions.getRight();

						Pair<Double, Integer> signalUncond = getOverlap(subsetOfAnnotation1, snps);
						int nrOverlapUncond1 = signalUncond.getRight();

						Pair<Double, Integer> signalUncond2 = getOverlap(subsetOfAnnotation2, snps);
						int nrOverlapUncond2 = signalUncond2.getRight();


						Track ty = new Track("");
						ty.addFeatures(y.getMiddle());

						Track tybar = new Track("");
						tybar.addFeatures(ybar.getMiddle());

						Pair<Double, Integer> signal1 = getOverlap(ty, y.getRight());
						Pair<Double, Integer> signal2 = getOverlap(tybar, ybar.getRight());

						if (DEBUG) {
							System.out.println(y.getRight().size() + "\t" + ybar.getRight().size() + "\t" + nrOverlapUncond1 + "\t" + nrOverlapUncond2);
							System.out.println(signal1.getRight() + "\t" + signal2.getRight());
						}

						sigmaPosterior = signal1.getLeft() + signal2.getLeft();
						sigmaSigmaPosterior += sigmaPosterior;

						nrOverlapping = signal1.getRight() + signal2.getRight();

						// locus score == % of iterations with >= 1 overlap
						for (int i = 0; i < options.nrIterations; i++) {
							shift(y.getRight(), y.getLeft());
							shift(ybar.getRight(), ybar.getLeft());
							Pair<Double, Integer> null1 = getOverlap(ty, y.getRight());
							Pair<Double, Integer> null2 = getOverlap(tybar, ybar.getRight());

							double posteriorsum = null1.getLeft() + null2.getLeft();
							sigmaPosteriorNull[i] = posteriorsum;
							sigmaSigmaPosteriorNull[i] += posteriorsum;

							int nrOverlappingNull = null1.getRight() + null2.getRight();
							if (nrOverlappingNull > 0) {
								sigmaLociWithAtLeastOneOverlapNull[i]++;
								locusScore++;
							}
						}
					} else {
						// run unconditional analysis
						// determine sigmaPosterior
						Pair<Double, Integer> signal = getOverlap(subsetOfAnnotation1, snps);
						sigmaPosterior = signal.getLeft();
						sigmaSigmaPosterior += sigmaPosterior;
						nrOverlapping = signal.getRight();

						for (int i = 0; i < options.nrIterations; i++) {
							shift(snps, region);
							Pair<Double, Integer> output = getOverlap(subsetOfAnnotation1, snps);
							sigmaPosteriorNull[i] = output.getLeft();
							sigmaSigmaPosteriorNull[i] += output.getLeft();
							int nrNullOverlapping = output.getRight();
							if (nrNullOverlapping > 0) {
								sigmaLociWithAtLeastOneOverlapNull[i]++;
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
						builder.append(1)
								.append("\t").append(1)
								.append("\t").append(0)
								.append("\t").append(sigmaPosterior)
								.append("\t").append(JSci.maths.ArrayMath.mean(sigmaPosteriorNull))
								.append("\t").append(JSci.maths.ArrayMath.standardDeviation(sigmaPosteriorNull))
								.append("\t").append(locusScore)
								.append("\t").append(regionName)
								.append("\t").append(origRegion)
								.append("\t").append(nrOverlapping)
								.append("\t").append(snps.size())
								.append("\t").append(0)
								.append("\t").append(annotation1name)
								.append("\t").append(subsetOfAnnotation1.getSize());
					} else {
						sigmaLociWithAtLeastOneOverlap++;

						double[] stats = getP(sigmaPosteriorNull, sigmaPosterior);

						builder.append(stats[0])
								.append("\t").append(stats[1])
								.append("\t").append(stats[2])
								.append("\t").append(sigmaPosterior)
								.append("\t").append(stats[3])
								.append("\t").append(stats[4])
								.append("\t").append(locusScore)
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
			} else {
				System.err.println("Thread " + threadNum + " | Cannot find region file: " + options.posteriorFile);
			}
			if (fctr % 10 == 0) {
				System.out.println("Thread " + threadNum + " | " + fctr + " out of " + regions.size() + " regions processed.");
			}
		}
		System.out.println("Thread " + threadNum + " | " + regions.size() + " out of " + regions.size() + " regions processed.");

		// calculate the overall score...
		double[] stats = getP(sigmaSigmaPosteriorNull, sigmaSigmaPosterior);
		double[] goshifterStats = getP(sigmaLociWithAtLeastOneOverlapNull, sigmaLociWithAtLeastOneOverlap);

		StringBuilder builder = new StringBuilder(500);
		builder.append(stats[0])
				.append("\t").append(stats[1])
				.append("\t").append(stats[2])
				.append("\t").append(sigmaSigmaPosterior)
				.append("\t").append(stats[3])
				.append("\t").append(stats[4])
				.append("\t").append(goshifterStats[0])
				.append("\t").append(goshifterStats[1])
				.append("\t").append(goshifterStats[2])
				.append("\t").append(sigmaLociWithAtLeastOneOverlap)
				.append("\t").append(goshifterStats[3])
				.append("\t").append(goshifterStats[4])
				.append("\t").append(totalNrOfOverlappingVariants)
				.append("\t").append(totalNrOfVariants)
				.append("\t").append(((double) totalNrOfOverlappingVariants / totalNrOfVariants))
				.append("\t").append(annotation1name);
		if (annotation2 != null) {
			builder.append("\t").append(annotation2name);
		}


		System.out.println("Thread " + threadNum + " | done testing.");
		return new Pair<>(builder.toString(), outputLines);
	}

	// determine p-values from permutation results
	private double[] getP(double[] overlapnull, double overlap) {

		Double[] nullcopy = new Double[overlapnull.length];
		for (int i = 0; i < overlapnull.length; i++) {
			nullcopy[i] = overlapnull[i];
		}

		Arrays.sort(nullcopy, Collections.reverseOrder());
		int nrAboveSignal = 0;
		for (int i = 0; i < nullcopy.length; i++) {
			double nulloverlap = nullcopy[i];
			if (nulloverlap >= overlap) {
				nrAboveSignal++;
			} else {
				break;
			}
		}

		// determine significance
		double mean = JSci.maths.ArrayMath.mean(overlapnull);
		double sd = JSci.maths.ArrayMath.standardDeviation(overlapnull);
		double z = (overlap - mean) / sd;
		double pz = ZScores.zToP(z);
		double p = (double) nrAboveSignal / options.nrIterations;

		if (mean == 0) {
			sd = 0;
			z = 0;
			pz = 1;
			p = 1;
		}

		return new double[]{p, pz, z, mean, sd};

	}


	// split a region into multiple smaller regions dependent upon annotation overlap
	// creates one region of annotation 1 that overlaps annotation2, and a region of annotation 1 without overlap with annotation 2
	private Pair<Triple<Feature, ArrayList<Feature>, ArrayList<SNPFeature>>, Triple<Feature, ArrayList<Feature>, ArrayList<SNPFeature>>> split(Feature queryRegion,
																																			   ArrayList<SNPFeature> snps,
																																			   Track subsetOfAnnotation1,
																																			   Track subsetOfAnnotation2) {
		// make features for those regions that are not overlapping annotation 2
		ArrayList<Feature> annotation2 = subsetOfAnnotation2.getAllFeatures();
		Collections.sort(annotation2, new FeatureComparator(false));

		ArrayList<Feature> ybarRegions = new ArrayList<>();
		int currentStart = queryRegion.getStart();

		for (int i = 0; i < annotation2.size(); i++) {
			Feature f1 = annotation2.get(i);
			if (queryRegion.overlaps(f1)) {
				int start = f1.getStart();
				int stop = f1.getStop();

				int startdiff = start - currentStart;
				if (startdiff > 0) {
					Feature ybarregion = new Feature(queryRegion.getChromosome(), currentStart, start);
					ybarRegions.add(ybarregion);
				}
				currentStart = stop;
			}
		}

		// last region not closed
		int diffstop = queryRegion.getStop() - currentStart;
		if (diffstop > 0) {
			Feature ybarregion = new Feature(queryRegion.getChromosome(), currentStart, queryRegion.getStop());
			ybarRegions.add(ybarregion);
		}

		// determine which parts of the annotations overlap and which ones don't
		ArrayList<Pair<Feature, ArrayList<Feature>>> partsOfAnnotation1OverlappingAnnotation2 =
				FeatureMerger.makeSubregionsUsingOverlap(subsetOfAnnotation1.getAllFeatures(), subsetOfAnnotation2.getAllFeatures());

		ArrayList<Pair<Feature, ArrayList<Feature>>> partsOfAnnotation1NotOverlappingAnnotation2 =
				FeatureMerger.makeSubregionsUsingOverlap(subsetOfAnnotation1.getAllFeatures(), ybarRegions);

		if (DEBUG) {
			// check split regions.
			int querysize = queryRegion.getSize();
			int size1 = 0;
			int size2 = 0;

			int sizex1 = 0;
			int sizex2 = 0;

			for (Pair<Feature, ArrayList<Feature>> p : partsOfAnnotation1OverlappingAnnotation2) {
				size1 += p.getLeft().getSize();
				ArrayList<Feature> q = p.getRight();
				for (Feature f : q) {
					int size = f.getSize();
					sizex1 += size;
				}
			}

			for (Pair<Feature, ArrayList<Feature>> p : partsOfAnnotation1NotOverlappingAnnotation2) {
				ArrayList<Feature> q = p.getRight();
				size1 += p.getLeft().getSize();
				for (Feature f : q) {
					int size = f.getSize();
					sizex2 += size;
				}
			}

			int totalsize = size1 + size2;
			int totalsizex = sizex1 + sizex2;

			ArrayList<Feature> x = subsetOfAnnotation1.getAllFeatures();
			int xsum = 0;
			for (int i = 0; i < x.size(); i++) {
				xsum += x.get(i).getSize();
			}

			System.out.println(queryRegion.toString()
					+ "\t" + querysize
					+ "\t" + totalsize
					+ "\t" + (totalsize == querysize)
					+ "\t" + xsum
					+ "\t" + totalsizex
					+ "\t" + (xsum == totalsizex));
		}

		// stitch regions together
		Triple<Feature, ArrayList<Feature>, ArrayList<SNPFeature>> yTriple = stitchRegion(queryRegion, partsOfAnnotation1OverlappingAnnotation2, snps);
		Triple<Feature, ArrayList<Feature>, ArrayList<SNPFeature>> yBarTriple = stitchRegion(queryRegion, partsOfAnnotation1NotOverlappingAnnotation2, snps);

		if (DEBUG) {
			// check if all SNPs are still there
			int nrSNPs = snps.size();
			int nrSNPs1 = yBarTriple.getRight().size();
			int nrSNPs2 = yTriple.getRight().size();
			System.out.println("nrSNPs: " + queryRegion.toString() + "\t" + nrSNPs + "\t" + nrSNPs1 + "\t" + nrSNPs2 + "\t" + (nrSNPs1 + nrSNPs2));
		}

		return new Pair<>(yTriple, yBarTriple);

	}


	// merge multiple smaller regions into one bigger one. rearrange SNPs and other annotations accordingly
	private Triple<Feature, ArrayList<Feature>, ArrayList<SNPFeature>> stitchRegion(Feature queryRegion,
																					ArrayList<Pair<Feature, ArrayList<Feature>>> partsOfAnnotation1OverlappingAnnotation2,
																					ArrayList<SNPFeature> snps) {
		int ywidth = 0;
		ArrayList<Feature> annotationXWithinAnnotationY = new ArrayList<>();
		ArrayList<SNPFeature> snpsWithinY = new ArrayList<SNPFeature>();

		for (int r = 0; r < partsOfAnnotation1OverlappingAnnotation2.size(); r++) {
			Pair<Feature, ArrayList<Feature>> p = partsOfAnnotation1OverlappingAnnotation2.get(r);
			Feature yRegion = p.getLeft();
			ArrayList<Feature> overlappingParts = p.getRight();

			if (yRegion.overlaps(queryRegion)) {
				int sta = yRegion.getStart();
				int sto = yRegion.getStop();

				if (sta < queryRegion.getStart()) {
					sta = queryRegion.getStart();
				}
				if (sto > queryRegion.getStop()) {
					sto = queryRegion.getStop();
				}

				int newYStart = ywidth;
				// make X relative to stitched together Y
				for (int r2 = 0; r2 < overlappingParts.size(); r2++) {
					Feature f2 = overlappingParts.get(r2);
					int f2start = f2.getStart();
					if (f2start < yRegion.getStart()) {
						f2start = yRegion.getStart();
					}
					int f2stop = f2.getStop();
					if (f2stop > yRegion.getStop()) {
						f2stop = yRegion.getStop();
					}
					int width = f2stop - f2start;

					int relativestart = f2start - yRegion.getStart();
					int newStart = newYStart + relativestart;
					int newStop = newStart + width;
					Feature f3 = new Feature(f2);
					f3.setStart(newStart);
					f3.setStop(newStop);
					annotationXWithinAnnotationY.add(f3);
				}

				// now we need to accomodate the SNPs as well ....
				for (SNPFeature f2 : snps) {
					if (f2.overlaps(yRegion) && f2.overlaps(queryRegion)) {
						// remap
						int pos = f2.getStart();
						int relativestart = pos - yRegion.getStart();
						int newStart = newYStart + relativestart;
						SNPFeature snpCopy = new SNPFeature(f2);
						snpCopy.setStart(newStart);
						snpCopy.setStop(newStart);
						snpsWithinY.add(snpCopy);
					}
				}
				int size = sto - sta;
				ywidth += size;
			}
		}

		Feature yfeature = new Feature(queryRegion.getChromosome(), 0, ywidth);

		return new Triple<Feature, ArrayList<Feature>, ArrayList<SNPFeature>>(yfeature, annotationXWithinAnnotationY, snpsWithinY);
	}

	private Feature trimRegion(Feature region, ArrayList<SNPFeature> snps) {
		// determine region size using variants
		int minStart = Integer.MAX_VALUE;
		int maxStop = 0;
		for (SNPFeature f : snps) {
			int pos = f.getStart();
			if (pos < minStart) {
				minStart = pos;
			}
			if (pos > maxStop) {
				maxStop = pos;
			}
		}

		minStart -= options.defaultRegionExtend;
		maxStop += options.defaultRegionExtend;
		region.setStart(minStart);
		region.setStop(maxStop);
		return region;
	}


	private void shift(ArrayList<SNPFeature> snps, Feature region) {
		int start = region.getStart();
		int stop = region.getStop();
		int regionSize = stop - start;

		int direction = 1;
		if (Math.random() > 0.5) {
			direction = -1;
		}
		int shift = direction * (int) Math.ceil(Math.random() * regionSize);

		for (SNPFeature snp : snps) {
			int pos = snp.getStart();
			int newpos = pos - shift;
			if (newpos > stop) {
				int diff = newpos - stop;
				newpos = start + diff;
			} else if (newpos < start) {
				int diff = start - newpos;
				newpos = stop - diff;
			}
			snp.setStart(newpos);
			snp.setStop(newpos);
		}
	}

	private Pair<Double, Integer> getOverlap(Track subsetOfAnnotations, ArrayList<SNPFeature> snps) {
		double sum = 0;
		int nrOverlap = 0;
		int maxalloweddist = options.getMaxAllowedDistance();
		for (SNPFeature snp : snps) {

			Iterable<Feature> subset = subsetOfAnnotations.getFeatures();
			ArrayList<Feature> overlappingFeatures = new ArrayList<Feature>();
			for (Feature f : subset) {
				if (f.overlaps(snp)) {
					overlappingFeatures.add(f);
				}
			}

			if (!overlappingFeatures.isEmpty()) {
				if (options.distanceweight.equals(BroShifterOptions.DISTANCEWEIGHT.NONE)) {
					sum += snp.getP();
					nrOverlap++;
				} else {
					// reweight distance
					// get closest annotation
					int minDist = Integer.MAX_VALUE;
					int snpPos = snp.getStart();

					PeakFeature overlapPeak = null;
					for (Feature f : overlappingFeatures) {
						int d = Integer.MAX_VALUE;
						if (f instanceof PeakFeature) {
							PeakFeature peak = (PeakFeature) f;
							overlapPeak = peak;
							int summit = peak.getSummit();
							d = Math.abs(summit - snpPos);
						} else {
							int midpoint = (f.getStop() + f.getStart()) / 2;
							d = Math.abs(midpoint - snpPos);
						}
						if (d < minDist) {
							minDist = d;
						}
					}
					double weight = 1;

					if (minDist == 0) {
						sum += snp.getP();
						nrOverlap++;
					} else if (minDist < maxalloweddist) {
						if (options.distanceweight.equals(BroShifterOptions.DISTANCEWEIGHT.SQUAREROOT)) {
							weight = 1d / (Math.sqrt(minDist));
						} else if (options.distanceweight.equals(BroShifterOptions.DISTANCEWEIGHT.INVERSE)) {
							weight = 1d / minDist;
						} else if (options.distanceweight.equals(BroShifterOptions.DISTANCEWEIGHT.EXPONENT)) {
							weight = 1d / minDist;
						} else if (options.distanceweight.equals(BroShifterOptions.DISTANCEWEIGHT.HEIGHTOVERDISTANCE)) {
							weight = minDist / overlapPeak.getSize();
						} else if (options.distanceweight.equals(BroShifterOptions.DISTANCEWEIGHT.LINEAR)) {
							weight = 1 - (minDist / maxalloweddist);

						}
						sum += (snp.getP() * weight);
						nrOverlap++;
					}
				}
			}
		}


		return new Pair<Double, Integer>(sum, nrOverlap);
	}

	private ArrayList<SNPFeature> readPosteriors(String file, Feature region) throws IOException {

		AssociationFile assocFile = new AssociationFile();
		ArrayList<AssociationResult> results = assocFile.read(file, region);

		ArrayList<SNPFeature> features = new ArrayList<>(results.size());
		for (AssociationResult r : results) {
			Feature f = r.getSnp();
			if (region.overlaps(f)) {
				SNPFeature f2 = new SNPFeature();
				f2.setChromosome(f.getChromosome());
				f2.setStart(f.getStart());
				f2.setStop(f.getStop());
				f2.setName(f.getName());

				f2.setP(r.getPosterior());
				features.add(f2);
			}
		}
		return features;
	}
}
