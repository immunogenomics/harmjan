package nl.harmjanwestra.utilities.graphics.panels;

import JSci.maths.ArrayMath;
import nl.harmjanwestra.utilities.coverage.Coverage;
import nl.harmjanwestra.utilities.bamfile.MinimalRead;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.graphics.Range;
import umcg.genetica.containers.Pair;

import java.awt.*;
import java.util.*;

/**
 * Created by hwestra on 7/19/15.
 */
public class CoveragePanel extends Panel {

	double[] coverageData;
	ArrayList<MinimalRead> reads;
	ArrayList<Feature> annotations;
	Feature region;

	Range dataRange;

	public CoveragePanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	public double[] getCoverageData() {
		return coverageData;
	}

	public void setCoverageData(int[] coverageDataInt) {
		this.coverageData = new double[coverageDataInt.length];
		for (int i = 0; i < coverageDataInt.length; i++) {
			coverageData[i] = coverageDataInt[i];
		}
	}

	public void setCoverageData(double[] coverageData) {
		this.coverageData = coverageData;
	}

	public ArrayList<MinimalRead> getReads() {
		return reads;
	}

	public void setReads(ArrayList<MinimalRead> reads) {
		this.reads = reads;
	}

	public void setReads(HashMap<String, MinimalRead> f) {
		Set<Map.Entry<String, MinimalRead>> readset = f.entrySet();

		reads = new ArrayList<MinimalRead>();

		for (Map.Entry<String, MinimalRead> e : readset) {
			reads.add(e.getValue());
		}

	}

	public ArrayList<Feature> getAnnotations() {
		return annotations;
	}

	public void setAnnotations(ArrayList<Feature> annotations) {
		this.annotations = annotations;
	}

	public Feature getRegion() {
		return region;
	}

	public void setRegion(Feature region) {
		this.region = region;
	}

	public Range getDataRange() {
		return dataRange;
	}

	public void setDataRange(Range dataRange) {
		this.dataRange = dataRange;
	}

	@Override
	public void draw(DefaultGraphics g) {

		if (dataRange == null) {
			dataRange = new Range(region.getStart(), 0, region.getStop(), 1);
		}

		Graphics2D g2d = g.getG2d();

		int annotationheight = 10;
		int regionsize = region.getStop() - region.getStart();

		int plotWidth = width - (2 * marginX);
		int plotHeight = height - (2 * marginY);

		if (annotations != null) {

			int y = y0 + marginY - annotationheight - 5;
			int ymax = plotHeight + annotationheight + 5;
			int x = x0 + marginX;
			int xmax = x0 + marginX + plotWidth;

			for (Feature f : annotations) {
				if (region.overlaps(f)) {
					// draw box indicating feature
					Pair<Integer, Integer> positions = getPlotPosition(f, region, regionsize, plotWidth, 0);
					int xpos1 = positions.getLeft();
					int xpos2 = positions.getRight();
					if (xpos1 < x) {
						xpos1 = x;
					}
					if (xpos2 > xmax) {
						xpos2 = xmax;
					}

					g2d.setColor(theme.getColor(1));
					g2d.fillRect(xpos1, y, (xpos2 - xpos1), annotationheight);
					// fill rect to bottom of plot
					g2d.setColor(theme.getColorSetOpacity(1, 0.75f));
					g2d.fillRect(xpos1, y, (xpos2 - xpos1), ymax);
				}
			}

		}

		if (reads != null) {
			// convert reads to features
			ArrayList<Feature> readFeatures = new ArrayList<Feature>();
			Coverage c = new Coverage();
			int ctr = 0;
			ArrayList<MinimalRead> allReads = new ArrayList<MinimalRead>();
			for (MinimalRead r : reads) {
				Pair<Integer, Integer> pos1 = c.getReadStartAndEnd(r, region);
				Feature f1 = new Feature(Chromosome.ONE, pos1.getLeft(), pos1.getRight());
				f1.setName("" + ctr);
				f1.setStrand(r.getStrand());
				readFeatures.add(f1);
				MinimalRead mate = r.getMate();
				ctr++;
				allReads.add(r);
				if (mate != null) {
					Pair<Integer, Integer> pos2 = c.getReadStartAndEnd(r, region);
					Feature f2 = new Feature(Chromosome.ONE, pos2.getLeft(), pos2.getRight());
					f2.setName("" + ctr);
					f2.setStrand(mate.getStrand());
					readFeatures.add(f2);
					allReads.add(mate);
					ctr++;
				}
			}


			Collections.sort(readFeatures, new FeatureComparator(false));

			// get y-lebels
			Pair<Integer, HashMap<Feature, Integer>> ylevels = getYLevels(readFeatures, region);

			int nrYLevels = ylevels.getLeft();
			HashMap<Feature, Integer> readYLevels = ylevels.getRight();
			g2d.setColor(theme.getLightGrey());

			int x = x0 + marginX;
			int xmax = x + plotWidth;
			int distanceYbetweenReads = 2;
			int y = y0 + marginY;
			int heightPerRead = (int) Math.floor(((double) plotHeight / (nrYLevels + (nrYLevels - 1 * distanceYbetweenReads))));
			if (heightPerRead > 10) {
				heightPerRead = 10;
			}
			if (heightPerRead < 1) {
				heightPerRead = 1;
			}

			int ymax = y0 + marginY + plotHeight;
			// TODO: compare minimalreads as features, and make sure to take CIGAR string into account when plotting.
			for (int i = 0; i < readFeatures.size(); i++) {


				Feature read = readFeatures.get(i);
				int ylevel = readYLevels.get(read);

				Pair<Integer, Integer> plotPositions = getPlotPosition(read, region, regionsize, plotWidth, 0);
				int xpos1 = plotPositions.getLeft();
				int xpos2 = plotPositions.getRight();

				if (xpos1 < x) {
					xpos1 = x;
				}
				if (xpos2 > xmax) {
					xpos2 = xmax;
				}

				int ypos = y + (ylevel * heightPerRead) + (ylevel * distanceYbetweenReads);

				if (ypos < ymax) {
					g2d.fillRect(xpos1, ypos, xpos2 - xpos1, heightPerRead);
				}


			}

			// draw axis..
		}


		if (coverageData != null) {

			// draw histogram
			int widthPerBar = (int) Math.ceil(plotWidth / dataRange.getRangeX());
			if (widthPerBar < 1) {
				widthPerBar = 1;
			}

			g2d.setColor(theme.getDarkGrey());
			double maxCoverage = ArrayMath.max(coverageData);
			for (int i = 0; i < coverageData.length; i++) {
				double c = coverageData[i];
				double perc = c / maxCoverage;
				int barheight = (int) Math.ceil(perc * plotHeight);


				double xperc = (double) i / coverageData.length;
				int xpos = x0 + marginX + (int) Math.ceil((xperc * plotWidth));

				int ypos = y0 + marginY + plotHeight - barheight;
				g2d.fillRect(xpos, ypos, widthPerBar, barheight);
			}

		}


	}

	private Pair<Integer, Integer> getPlotPosition(Feature g, Feature region, int regionSize, int nrPixels, int nrPixelsBetweenGene, FontMetrics metrics) {
		int start = g.getStart();
		int featurewidth = g.getStop() - start;
		int relativeStart = start - region.getStart();
		if (relativeStart < 0) {
			featurewidth -= Math.abs(relativeStart);
			relativeStart = 0;
		}


		int relativeStop = relativeStart + featurewidth;
		if (relativeStop > region.getStop()) {
			relativeStop = regionSize;
		}

		double percStart = (double) relativeStart / regionSize;
		double percStop = (double) relativeStop / regionSize;

		int genenamewidth = 0;
		if (metrics != null) {
			genenamewidth = metrics.stringWidth(g.getName());
		}
		int pixelStart = (int) Math.ceil(percStart * nrPixels);
		int pixelStop = (int) Math.ceil(percStop * nrPixels) + nrPixelsBetweenGene + genenamewidth + nrPixelsBetweenGene;

		return new Pair<Integer, Integer>(pixelStart, pixelStop);
	}

	private Pair<Integer, Integer> getPlotPosition(Feature g, Feature region, int regionSize, int nrPixels, int marginbetweenelements) {
		return getPlotPosition(g, region, regionSize, nrPixels, marginbetweenelements, null);
	}

	private Pair<Integer, HashMap<Feature, Integer>> getYLevels(ArrayList<Feature> features, Feature region) {

		int regionWidth = region.getStop() - region.getStart();
		int plotWidth = width - (2 * marginX);

		// determine y-plot levels
		HashMap<Feature, Integer> ylevels = new HashMap<Feature, Integer>();
		ArrayList<Feature> genesAtCurrentLevel = new ArrayList<Feature>();
		ArrayList<Feature> genesAtNextLevel = new ArrayList<Feature>();
		genesAtCurrentLevel.addAll(features);

		// get gene
		int ylevel = 0;

		int marginXBetweenGenes = 10;
		int marginYBetweenGenes = 10;// determine y-plot levels

		// get gene
		while (ylevels.size() != features.size()) {

			int currentGeneIndex = 0;
			int nextGeneIndex = 1;
			// System.out.println(ylevel + " ylevel: " + genesAtCurrentLevel.size() + " genes left. " + currentGeneIndex + " current gene");

			if (genesAtCurrentLevel.size() == 1) {
				Feature currentGene = genesAtCurrentLevel.get(currentGeneIndex);
				ylevels.put(currentGene, ylevel);
			} else {
				while (currentGeneIndex < genesAtCurrentLevel.size()) {
					Feature currentGene = genesAtCurrentLevel.get(currentGeneIndex);
					ylevels.put(currentGene, ylevel);

					Pair<Integer, Integer> pixels1 = getPlotPosition(currentGene, region, regionWidth, plotWidth, marginXBetweenGenes);
					boolean overlap = true;
					if (currentGeneIndex == genesAtCurrentLevel.size() - 1) {
						break;
					}
					if (nextGeneIndex == genesAtCurrentLevel.size()) {
						// reached the end
						break;
					}
					while (overlap && nextGeneIndex < genesAtCurrentLevel.size()) {
						Feature nextGene = genesAtCurrentLevel.get(nextGeneIndex);
						Pair<Integer, Integer> pixels2 = getPlotPosition(nextGene, region, regionWidth, plotWidth, marginXBetweenGenes);

						// check overlap
						overlap = overlap(pixels1, pixels2);


						// if overlap, place to next y level
						if (overlap) {
							genesAtNextLevel.add(nextGene);
						} else {
							currentGeneIndex = nextGeneIndex;
						}

						// if no overlap, continue with next gene
						//System.out.println(currentGeneIndex + " - " + nextGeneIndex + "\t" + overlap + "\t" + pixels1 + "\t" + pixels2);
						nextGeneIndex++;
					}

				}
			}


			if (!genesAtNextLevel.isEmpty()) {
				//System.out.println(genesAtNextLevel.size() + " genes at next level.");
				genesAtCurrentLevel = genesAtNextLevel;
				genesAtNextLevel = new ArrayList<Feature>();
				ylevel++;
			}
			//System.out.println(ylevel + " next y level. " + ylevels.size() + " vs " + features.size() + "\t" + genesAtNextLevel.size());


		}
		return new Pair<Integer, HashMap<Feature, Integer>>(ylevel, ylevels);
	}

	private boolean overlap(Pair<Integer, Integer> pixel1, Pair<Integer, Integer> pixels2) {

		Feature f1 = new Feature();
		f1.setChromosome(Chromosome.ONE);
		f1.setStart(pixel1.getLeft());
		f1.setStop(pixel1.getRight());

		Feature f2 = new Feature();
		f2.setChromosome(Chromosome.ONE);
		f2.setStart(pixels2.getLeft());
		f2.setStop(pixels2.getRight());

		return f1.overlaps(f2);
	}
}
