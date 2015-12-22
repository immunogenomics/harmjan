package nl.harmjanwestra.broshifter;


import nl.harmjanwestra.broshifter.CLI.BroShifterOptions;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.*;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.Callable;

/**
 * Created by hwestra on 11/17/15.
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

		// for each region
		ArrayList<String> outputLines = new ArrayList<>();

		for (int fctr = 0; fctr < regions.size(); fctr++) {
			Feature region = regions.get(fctr);

			String origRegion = region.toString();

			// read posteriors
			String regionStr = region.toString();
			double locusScore = 0;

			// TODO: make this some kind of default format
			// ALSO: make it independent of the region format
			String file = options.directoryWithPosteriors + regionStr + ".txt";
			int annotation2size = 0;
			if (Gpio.exists(file) && !region.getChromosome().equals(Chromosome.X)) {

				ArrayList<SNPFeature> snps = readPosteriors(file, region);

				if (snps.size() > 0) {

					if (options.trimRegions) {
						region = trimRegion(region, snps);
					}

					// get subset of annotations
					Track subsetOfAnnotation1 = annotation1Track.getSubset(region.getChromosome(), region.getStart(), region.getStop());

					String ln = "";
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

					if (nrOverlapping == 0) {
						ln = 1
								+ "\t" + 1
								+ "\t" + 0
								+ "\t" + sigmaPosterior
								+ "\t" + JSci.maths.ArrayMath.mean(sigmaPosteriorNull)
								+ "\t" + JSci.maths.ArrayMath.standardDeviation(sigmaPosteriorNull)
								+ "\t" + locusScore
								+ "\t" + regionName
								+ "\t" + origRegion
								+ "\t" + nrOverlapping
								+ "\t" + snps.size()
								+ "\t0"
								+ "\t" + annotation1name
								+ "\t" + subsetOfAnnotation1.getSize();
						if (annotation2 != null) {
							ln += "\t" + annotation2name;
							ln += "\t" + annotation2size;
						}
					} else {
						sigmaLociWithAtLeastOneOverlap++;

						double[] stats = getP(sigmaPosteriorNull, sigmaPosterior);

						ln = stats[0]
								+ "\t" + stats[1]
								+ "\t" + stats[2]
								+ "\t" + sigmaPosterior
								+ "\t" + stats[3]
								+ "\t" + stats[4]
								+ "\t" + locusScore
								+ "\t" + regionName
								+ "\t" + origRegion
								+ "\t" + nrOverlapping
								+ "\t" + snps.size()
								+ "\t" + ((double) nrOverlapping / snps.size())
								+ "\t" + annotation1name
								+ "\t" + subsetOfAnnotation1.getAllFeatures().size();

						if (annotation2 != null) {
							ln += "\t" + annotation2name;
							ln += "\t" + annotation2size;
						}
					}
					outputLines.add(ln);
				}
			} else {
				System.err.println("Thread " + threadNum + " | Cannot find region file: " + file);
			}
			if (fctr % 10 == 0) {
				System.out.println("Thread " + threadNum + " | " + fctr + " out of " + regions.size() + " regions processed.");
			}
		}
		System.out.println("Thread " + threadNum + " | " + regions.size() + " out of " + regions.size() + " regions processed.");

		// calculate the overall score...
		double[] stats = getP(sigmaSigmaPosteriorNull, sigmaSigmaPosterior);
		double[] goshifterStats = getP(sigmaLociWithAtLeastOneOverlapNull, sigmaLociWithAtLeastOneOverlap);

		String ln = stats[0]
				+ "\t" + stats[1]
				+ "\t" + stats[2]
				+ "\t" + sigmaSigmaPosterior
				+ "\t" + stats[3]
				+ "\t" + stats[4]
				+ "\t" + goshifterStats[0]
				+ "\t" + goshifterStats[1]
				+ "\t" + goshifterStats[2]
				+ "\t" + sigmaLociWithAtLeastOneOverlap
				+ "\t" + goshifterStats[3]
				+ "\t" + goshifterStats[4]
				+ "\t" + totalNrOfOverlappingVariants
				+ "\t" + totalNrOfVariants
				+ "\t" + ((double) totalNrOfOverlappingVariants / totalNrOfVariants)
				+ "\t" + annotation1name;
		if (annotation2 != null) {
			ln += "\t" + annotation2name;
		}


		System.out.println("Thread " + threadNum + " | done testing.");
		return new Pair<>(ln, outputLines);
	}

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
					// get closest annotation
					int maxDist = Integer.MAX_VALUE;
					int snpPos = snp.getStart();
					;
					for (Feature f : overlappingFeatures) {
						int d = Integer.MAX_VALUE;
						if (f instanceof PeakFeature) {
							PeakFeature peak = (PeakFeature) f;
							int summit = peak.getSummit();
							d = Math.abs(summit - snpPos);
						} else {
							int midpoint = (f.getStop() + f.getStart()) / 2;
							d = Math.abs(midpoint - snpPos);
						}
						if (d < maxDist) {
							d = maxDist;
						}
					}
					double weight = 1;
					if (options.distanceweight.equals(BroShifterOptions.DISTANCEWEIGHT.SQUAREROOT)) {
						weight = 1d / (Math.sqrt(maxDist + 1));
					} else if (options.distanceweight.equals(BroShifterOptions.DISTANCEWEIGHT.INVERSE)) {
						weight = 1d / maxDist;
					} else if (options.distanceweight.equals(BroShifterOptions.DISTANCEWEIGHT.LINEAR)) {
						weight = (1 - (1d / 500)) * maxDist;

					}
					sum += (snp.getP() * weight);
					nrOverlap++;
				}
			}
		}


		return new Pair<Double, Integer>(sum, nrOverlap);
	}


	// TODO: replace with posteriorfile in utils
	private ArrayList<SNPFeature> readPosteriors(String file, Feature region) throws IOException {
		ArrayList<SNPFeature> output = new ArrayList<SNPFeature>();

		TextFile tf = new TextFile(file, TextFile.R);

		// Chr	Pos	Id	CombinedId	Beta	Se	PVal	Posterior
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

			SNPFeature f = new SNPFeature();
			f.setChromosome(region.getChromosome());

			Integer pos = Integer.parseInt(elems[1]);
			f.setStart(pos);
			f.setStop(pos);

			if (region.overlaps(f)) {
				Double posterior = Double.parseDouble(elems[7]);
				f.setP(posterior);
				output.add(f);
			}

			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();
		return output;
	}
}
