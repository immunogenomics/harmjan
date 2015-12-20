package nl.harmjanwestra.utilities.graphics.panels;

import nl.harmjanwestra.utilities.features.*;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import umcg.genetica.containers.Pair;

import java.awt.*;
import java.util.*;

/**
 * Created by hwestra on 7/17/15.
 */
public class GenePanel extends Panel {

	public GenePanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	private ArrayList<Gene> genes;
	private Feature region;

	public void setData(Feature region, ArrayList<Gene> genes) {
		this.region = region;
		this.genes = genes;
	}

	public ArrayList<Gene> getGenes() {
		return genes;
	}

	public void setGenes(Set<Gene> geneSet) {
		ArrayList<Gene> tmpGenes = new ArrayList<Gene>();
		tmpGenes.addAll(geneSet);
		this.genes = tmpGenes;
	}

	public void setGenes(ArrayList<Gene> genes) {
		this.genes = genes;
	}

	public Feature getRegion() {
		return region;
	}

	public void setRegion(Feature region) {
		this.region = region;
	}

	@Override
	public void draw(DefaultGraphics g) {
		Graphics2D g2d = g.getG2d();
		g2d.setFont(theme.getSmallFont());


		// sort
		Collections.sort(genes, new FeatureComparator(false));

		int plotWidth = width - (2 * marginX);
		int plotHeight = height - (2 * marginY);
		int regionWidth = region.getStop() - region.getStart();

		// get the gene plot X start and end positions
		// remove duplicate genenames.
		HashSet<String> geneNames = new HashSet<String>();
		for (Gene gene : genes) {
			int ctr = 2;
			if (geneNames.contains(gene.getName())) {
				String name = gene.getName() + "_" + ctr;
				while (geneNames.contains(name)) {
					name = gene.getName() + "_" + ctr;
					ctr++;
				}
				gene.setName(name);
			}
			geneNames.add(gene.getName());
		}

		System.out.println(geneNames.size() + " genes to plot");
		// determine y-plot levels
		HashMap<Gene, Integer> ylevels = new HashMap<Gene, Integer>();
		ArrayList<Gene> genesAtCurrentLevel = new ArrayList<Gene>();
		ArrayList<Gene> genesAtNextLevel = new ArrayList<Gene>();
		genesAtCurrentLevel.addAll(genes);

		// get gene
		int ylevel = 0;

		FontMetrics metrics = g2d.getFontMetrics();
		int marginXBetweenGenes = 10;
		int marginYBetweenGenes = 10;

		while (ylevels.size() != genes.size()) {

			int currentGeneIndex = 0;
			int nextGeneIndex = 1;
//			System.out.println(ylevel + " ylevel: " + genesAtCurrentLevel.size() + " genes left. " + currentGeneIndex + " current gene");

			if (genesAtCurrentLevel.size() == 1) {
				Gene currentGene = genesAtCurrentLevel.get(currentGeneIndex);
				ylevels.put(currentGene, ylevel);
			} else {
				while (currentGeneIndex < genesAtCurrentLevel.size()) {
					Gene currentGene = genesAtCurrentLevel.get(currentGeneIndex);
					ylevels.put(currentGene, ylevel);
//					System.out.println(currentGene.getGeneId() + "\ty: " + ylevel);
					Pair<Integer, Integer> pixels1 = getGenePlotPositions(currentGene, region, regionWidth, plotWidth, marginXBetweenGenes, metrics);
					boolean overlap = true;
					if (currentGeneIndex == genesAtCurrentLevel.size() - 1) {
						break;
					}
					if (nextGeneIndex == genesAtCurrentLevel.size()) {
						// reached the end
						break;
					}
					while (overlap && nextGeneIndex < genesAtCurrentLevel.size()) {
						Gene nextGene = genesAtCurrentLevel.get(nextGeneIndex);
						Pair<Integer, Integer> pixels2 = getGenePlotPositions(nextGene, region, regionWidth, plotWidth, marginXBetweenGenes, metrics);

						// check overlap
						overlap = overlap(pixels1, pixels2);


						// if overlap, place to next y level
						if (overlap) {
							genesAtNextLevel.add(nextGene);
						} else {
							currentGeneIndex = nextGeneIndex;
						}

						// if no overlap, continue with next gene
//						System.out.println(currentGeneIndex + " - " + nextGeneIndex + "\t" + overlap + "\t" + pixels1 + "\t" + pixels2);
						nextGeneIndex++;
					}

				}
			}


			if (!genesAtNextLevel.isEmpty()) {
//				System.out.println(genesAtNextLevel.size() + " genes at next level.");
				genesAtCurrentLevel = genesAtNextLevel;
				genesAtNextLevel = new ArrayList<Gene>();
				ylevel++;
			}
//			System.out.println(ylevel + " next y level. " + ylevels.size() + " vs " + genes.size() + "\t" + genesAtNextLevel.size());


		}

//		System.out.println("Done placing genes.");

		g2d.setColor(theme.getDarkGrey());

		int heightPerGene = 10; // (int) Math.floor(plotHeight / (ylevel + 1)) - (ylevel * marginYBetweenGenes);

		g2d.drawString("Genes", x0, y0 - 10);

		for (Gene gene : genes) {
			Integer geneYlevel = ylevels.get(gene);
			if (geneYlevel == null) {
				System.out.println(gene.getGeneId() + " not mapped to y-level?");
			}
			int y1 = y0 + marginY + plotHeight - (geneYlevel * marginYBetweenGenes) - (geneYlevel * heightPerGene);

			ArrayList<Transcript> transcripts = gene.getTranscripts();
			HashSet<Exon> exons = new HashSet<Exon>();
			for (Transcript t : transcripts) {
				exons.addAll(t.getExons());
			}

			for (Exon f : exons) {
				if (region.overlaps(f)) {
					int start = f.getStart();
					int featurewidth = f.getStop() - start;
					int relativeStart = start - region.getStart();
					if (relativeStart < 0) {
						featurewidth -= Math.abs(relativeStart);
						relativeStart = 0;
					}

					int relativeStop = relativeStart + featurewidth;
					if (relativeStop > regionWidth) {
						relativeStop = regionWidth;
					}

					double percStart = (double) relativeStart / regionWidth;

					double percStop = (double) relativeStop / regionWidth;

					int pixelStart = (int) Math.ceil(percStart * plotWidth);
					int pixelStop = (int) Math.ceil(percStop * plotWidth);

					int exonwidth = pixelStop - pixelStart;
					if (pixelStop - pixelStart <= 0) {
						exonwidth = 1;
					}

					g2d.fillRect(x0 + marginX + pixelStart, y1, exonwidth, heightPerGene);
				}
			}

			Pair<Integer, Integer> pixelPair = getGenePlotPositions(gene, region, regionWidth, plotWidth);
			int pixelStart = pixelPair.getLeft();
			int pixelStop = pixelPair.getRight();

			
			int stop = x0 + marginX + pixelStop;
			if (x0 + marginX + pixelStop > x0 + marginX + plotWidth) {
				stop = x0 + marginX + plotWidth;
			}


			g2d.drawLine(x0 + marginX + pixelStart, y1 + 5, stop, y1 + 5);
			g2d.drawString(gene.getGeneId(), stop + marginXBetweenGenes, y1 + 10);
		}
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

	private Pair<Integer, Integer> getGenePlotPositions(Gene g, Feature region, int regionSize, int nrPixels, int nrPixelsBetweenGene, FontMetrics metrics) {
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
			genenamewidth = metrics.stringWidth(g.getGeneId());
		}
		int pixelStart = (int) Math.ceil(percStart * nrPixels);
		int pixelStop = (int) Math.ceil(percStop * nrPixels) + nrPixelsBetweenGene + genenamewidth + nrPixelsBetweenGene;

		return new Pair<Integer, Integer>(pixelStart, pixelStop);
	}

	private Pair<Integer, Integer> getGenePlotPositions(Gene g, Feature region, int regionSize, int nrPixels) {
		return getGenePlotPositions(g, region, regionSize, nrPixels, 0, null);
	}

}
