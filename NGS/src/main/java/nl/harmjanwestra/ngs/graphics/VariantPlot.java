package nl.harmjanwestra.ngs.graphics;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.features.*;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.vcf.VCFFunctions;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.vcf.filter.GenotypeQualityFilter;
import nl.harmjanwestra.utilities.vcf.filter.ReadDepthFilter;
import nl.harmjanwestra.utilities.vcf.filter.VCFGenotypeFilter;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.awt.*;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Created by hwestra on 4/23/15.
 */
public class VariantPlot extends DefaultGraphics {

	private int margin;

	public VariantPlot() {

	}

	public VariantPlot(String name, int pageWidth, int pageHeight) throws FileNotFoundException, DocumentException {
		super(name, pageWidth, pageHeight);
		figureWidth = pageWidth;
	}

	public void setMargin(int margin) {
		this.margin = margin;
	}

	public void plotVariantsUniqueIneachDataset(String[] variantFiles, String[] variantFileNames, String sequencedRegionFile, Feature region, String gtf) throws IOException {
		TextFile tf1 = new TextFile(sequencedRegionFile, TextFile.R);

		HashSet<Feature> sequencedRegions = new HashSet<Feature>();

		String[] elems = tf1.readLineElems(TextFile.tab);
		while (elems != null) {
			Feature f = new Feature();
			if (elems.length > 2) {
				f.setChromosome(Chromosome.parseChr(elems[0]));
				f.setStart(Integer.parseInt(elems[1]));
				f.setStop(Integer.parseInt(elems[2]));
				sequencedRegions.add(f);
			}
			elems = tf1.readLineElems(TextFile.tab);
		}
		tf1.close();

		int regionSize = region.getStop() - region.getStart();
		int nrPixels = figureWidth - (2 * margin);

		System.out.println(nrPixels);

		drawOverlappingGenes(gtf, region, regionSize, nrPixels);

		g2d.setColor(new Color(255, 1, 0, 202));
		for (Feature f : sequencedRegions) {
			if (f.overlaps(region)) {
				int start = f.getStart();
				int featurewidth = f.getStop() - start;
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

				int pixelStart = (int) Math.ceil(percStart * nrPixels);
				int pixelStop = (int) Math.ceil(percStop * nrPixels);

				System.out.println(f.toString());
				System.out.println(pixelStart + "\t" + pixelStop);
				int y1 = margin;
				int y2 = margin + 10;
				g2d.fillRect(margin + pixelStart, y1, pixelStop - pixelStart, 10);
			}
		}

		int betweenmargin = 25;
		int starty = margin + betweenmargin + 10;
		int dotsize = 50;
		int minDotSize = 2;

		VCFFunctions func = new VCFFunctions();
		HashMap<Feature, Integer> ctr = new HashMap<Feature, Integer>();
		for (int i = 0; i < variantFiles.length; i++) {
			String variantFile = variantFiles[i];
			TextFile tf = new TextFile(variantFile, TextFile.R);
			String ln = tf.readLine();
			while (ln != null) {

				if (!ln.startsWith("#")) {

					ArrayList<VCFGenotypeFilter> filters = new ArrayList<>();
					filters.add(new GenotypeQualityFilter(30));
					filters.add(new ReadDepthFilter(10));
					VCFVariant v = new VCFVariant(ln, filters, false);

					Feature f = new Feature();

					f.setChromosome(Chromosome.parseChr(v.getChr()));
					f.setStart(v.getPos());
					f.setStop(v.getPos());
					Integer ct = ctr.get(f);
					if (ct == null) {
						ct = 0;
					}
					ct++;
					ctr.put(f, ct);
				}
				ln = tf.readLine();
			}
			tf.close();

		}
		System.out.println(ctr.size() + " variants loaded");

		for (int i = 0; i < variantFiles.length; i++) {
			starty += betweenmargin + dotsize;
			String variantFile = variantFiles[i];
			String variantFileName = variantFileNames[i];

			TextFile tf = new TextFile(variantFile, TextFile.R);

			String ln = tf.readLine();

			Color grey = new Color(0, 0, 0, 28);
			g2d.setColor(grey);

			while (ln != null) {

				if (!ln.startsWith("#")) {
					ArrayList<VCFGenotypeFilter> filters = new ArrayList<>();
					filters.add(new GenotypeQualityFilter(30));
					filters.add(new ReadDepthFilter(10));
					VCFVariant v = new VCFVariant(ln, filters, false);

					Feature f = new Feature();

					f.setChromosome(Chromosome.parseChr(v.getChr()));
					f.setStart(v.getPos());
					f.setStop(v.getPos());

					if (region.overlaps(f)) {
						double maf = v.getMAF();
						if (maf > 0 && maf < 1) {
							int size = (int) Math.ceil((maf * 2) * dotsize);
							System.out.println(maf + "\t" + size);
							int relativeStart = f.getStart() - region.getStart();
							if (relativeStart < 0) {
								relativeStart = 0;
							}

							double percStart = (double) relativeStart / regionSize;
							int pixelStart = (int) Math.ceil(percStart * nrPixels);

							if (size < minDotSize) {
								size = minDotSize;
							}
							int halfdsize = size / 2;

							if (ctr.get(f) == null || ctr.get(f) == 1) {
								g2d.setColor(new Color(255, 0, 0, 128));
							} else {
								g2d.setColor(grey);
							}

							g2d.fillOval(margin + pixelStart - halfdsize, starty - halfdsize, size, size);
						}
					}
				}
				ln = tf.readLine();
			}
			tf.close();
		}

		// plotVariantsUniqueIneachDataset coordinates
		starty += betweenmargin + 10;

		g2d.setColor(Color.BLACK);

		g2d.drawLine(margin, starty, margin + nrPixels, starty);

		int unit = (int) Math.ceil(determineUnit(regionSize));

		for (int i = region.getStart(); i < region.getStop(); i++) {
			if (i % unit == 0) {
				int relativeStart = i - region.getStart();
				double percStart = (double) relativeStart / regionSize;
				int pixelStart = (int) Math.ceil(percStart * nrPixels);
				g2d.drawLine(margin + pixelStart, starty - 5, margin + pixelStart, starty + 5);
				g2d.drawString("" + i, margin + pixelStart, starty + 10);
			}
		}

		starty += betweenmargin + 10;

		int startx = margin;
		int between = 5;
		int c = 0;
		int nextX = startx;
		for (double maf = 0; maf < 1; maf += 0.1) {
			int size = (int) Math.ceil((maf) * dotsize);

			g2d.setColor(new Color(255, 0, 0, 128));
			int halfdsize = size / 2;
			g2d.fillOval(margin + nextX - halfdsize, starty - halfdsize, size, size);

			nextX += (between + size);

			c++;
		}
		close();
	}


	public void plotImputationRSquared(
			String[] datasetFiles,
			String[] datasetNames,
			String[] inputBeforeImputation,
			String sequencedRegionFile,
			Feature region,
			String gtf) throws IOException {

		TextFile tf1 = new TextFile(sequencedRegionFile, TextFile.R);

		HashSet<Feature> sequencedRegions = new HashSet<Feature>();

		String[] elems = tf1.readLineElems(TextFile.tab);
		while (elems != null) {
			Feature f = new Feature();
			if (elems.length > 2) {
				f.setChromosome(Chromosome.parseChr(elems[0]));
				f.setStart(Integer.parseInt(elems[1]));
				f.setStop(Integer.parseInt(elems[2]));
				sequencedRegions.add(f);
			}
			elems = tf1.readLineElems(TextFile.tab);
		}
		tf1.close();

		int regionSize = region.getStop() - region.getStart();
		int nrPixels = figureWidth - (2 * margin);

		System.out.println(nrPixels);

		drawOverlappingGenes(gtf, region, regionSize, nrPixels);

		g2d.setColor(new Color(208, 83, 77));

		// plot sequenced regions
		Font defaultfont = g2d.getFont();
		g2d.setFont(new Font("default", Font.BOLD, 16));

		g2d.drawString("Targeted regions in sequencing", margin, margin - 20);

		g2d.setFont(defaultfont);

		for (Feature f : sequencedRegions) {
			if (f.overlaps(region)) {
				int start = f.getStart();
				int featurewidth = f.getStop() - start;
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

				int pixelStart = (int) Math.ceil(percStart * nrPixels);
				int pixelStop = (int) Math.ceil(percStop * nrPixels);

				System.out.println(f.toString());
				System.out.println(pixelStart + "\t" + pixelStop);
				int y1 = margin;
				int y2 = margin + 10;
				g2d.fillRect(margin + pixelStart, y1, pixelStop - pixelStart, 10);
			}
		}

		int betweenmargin = 100;
		int starty = margin + betweenmargin + 10;
		int plotSize = 500 / datasetFiles.length;
		int minDotSize = 2;


		starty += betweenmargin + plotSize;


		Color grey = new Color(0, 0, 0, 28);
		g2d.setColor(grey);

		int plotStart = starty;
		int q = 0;

		Color[] colors = new Color[datasetFiles.length];
		for (int i = 0; i < colors.length; i++) {
			switch (i) {
				case 0:
					colors[i] = new Color(70, 67, 58);
					break;
				case 1:
					colors[i] = new Color(174, 164, 140);
					break;
				case 2:
					colors[i] = new Color(208, 83, 77);
					break;
				case 3:
					colors[i] = new Color(98, 182, 177);
					break;
				case 4:
					colors[i] = new Color(116, 156, 80);
					break;
				case 5:
					colors[i] = new Color(156, 67, 109);
					break;
				case 6:
					colors[i] = new Color(81, 65, 156);
					break;

			}
		}


		double unit = 0.1;
		Font originalfont = g2d.getFont();
		g2d.setFont(new Font("default", Font.BOLD, 16));
		FontMetrics metrics = g2d.getFontMetrics(g2d.getFont());

//		int starty =

		for (int dataset = 0; dataset < datasetFiles.length; dataset++) {
			int plotStarty = starty + (dataset * plotSize) + (dataset * betweenmargin);
			if (datasetFiles[dataset] != null) {

				g2d.setColor(new Color(98, 182, 177));

				// plot input before imputation (if any)
				if (inputBeforeImputation != null && inputBeforeImputation[dataset] != null) {
					TextFile tf = new TextFile(inputBeforeImputation[dataset], TextFile.R);
					elems = tf.readLineElems(Strings.whitespace);
					while (elems != null) {
						if (elems.length > 1 && !elems[0].startsWith("#")) {
							Chromosome chr = Chromosome.parseChr(elems[0]);
							Integer pos = Integer.parseInt(elems[1]);
//					System.out.println(chr.toString() + "\t" + pos);
							Feature f = new Feature();
							f.setChromosome(chr);
							f.setStart(pos);
							f.setStop(pos);
							if (region.overlaps(f)) {
								int relativeStart = pos - region.getStart();
								double percStart = (double) relativeStart / regionSize;
								int pixelStart = margin + (int) Math.ceil(percStart * nrPixels);
								int pixelY = 20;
								int dotsize = 4;
								g2d.fillOval(pixelStart - (dotsize / 2), plotStarty - plotSize - pixelY - (dotsize / 2), dotsize, dotsize);
							}
						}
						elems = tf.readLineElems(Strings.whitespace);
					}
					tf.close();
				}

				g2d.setColor(colors[dataset]);
				TextFile tf = new TextFile(datasetFiles[dataset], TextFile.R);

				elems = tf.readLineElems(Strings.whitespace);
				while (elems != null) {
					if (elems.length > 1 && !elems[0].startsWith("#")) {
						Chromosome chr = Chromosome.parseChr(elems[0]);
						Integer pos = Integer.parseInt(elems[1]);
//					System.out.println(chr.toString() + "\t" + pos);
						Feature f = new Feature();
						f.setChromosome(chr);
						f.setStart(pos);
						f.setStop(pos);
						if (region.overlaps(f)) {
							String name = elems[2];
							String[] infoElems = elems[7].split(";");
							for (int e = 0; e < infoElems.length; e++) {
								if (infoElems[e].startsWith("AR2")) {
									String[] infoElemsElems = infoElems[e].split("=");
									Double ar2 = Double.parseDouble(infoElemsElems[1]);
									// plot
									int relativeStart = pos - region.getStart();
									double percStart = (double) relativeStart / regionSize;
									int pixelStart = margin + (int) Math.ceil(percStart * nrPixels);

									int pixelY = (int) Math.ceil(ar2 * plotSize);

									int dotsize = 2 + (int) Math.ceil(ar2 * 10);
									g2d.fillOval(pixelStart - (dotsize / 2), plotStarty - pixelY - (dotsize / 2), dotsize, dotsize);
								}
							}

						}
					}
					elems = tf.readLineElems(Strings.whitespace);
				}
				int adv = metrics.stringWidth(datasetNames[dataset]);
				int hgt = metrics.getHeight();
				g2d.drawString(datasetNames[dataset], margin, plotStarty - plotSize - 40);
				tf.close();
			}
		}


		g2d.setFont(originalfont);

		// draw red line near 5E-8)


		double yperc = 0.8;
		int pixelY = (int) Math.ceil(yperc * plotSize);
		g2d.setColor(new Color(208, 83, 77));
		g2d.drawLine(margin, starty - pixelY, margin + nrPixels, starty - pixelY);

		// plotVariantsUniqueIneachDataset coordinates
//		starty += 5;

		g2d.setColor(new Color(70, 67, 58));

		// y-axis

		for (int dataset = 0; dataset < datasetFiles.length; dataset++) {
			// determine unit

			int plotStarty = starty + (dataset * plotSize) + (dataset * betweenmargin);
			g2d.drawLine(margin - 10, plotStarty, margin - 10, plotStarty - plotSize);
			double steps = 0.1;

			String pattern = "###,###,###.##";
			DecimalFormat decimalFormat = new DecimalFormat(pattern);

			for (double i = 0; i < 1 + steps; i += steps) {
				if (i <= 1) {
					int plusY = (int) Math.ceil(((double) i / 1) * plotSize);
					g2d.drawLine(margin - 13, plotStarty - plusY, margin - 7, plotStarty - plusY);
					String formattedStr = decimalFormat.format(i);
					int adv = metrics.stringWidth(formattedStr);
					int hgt = metrics.getHeight();
					Dimension size = new Dimension(adv + 10, hgt + 10);
					g2d.drawString(formattedStr, margin - (int) Math.ceil(size.getWidth()) - 10, plotStarty - plusY + 5);
				}
			}

			// x-axis

			g2d.drawLine(margin - 5, plotStarty + 5, margin + nrPixels + 5, plotStarty + 5);

			int xunit = (int) Math.ceil(determineUnit(regionSize));
			while (regionSize / xunit < 10 && xunit > 1) {
				xunit /= 2;
			}
			while (regionSize / xunit > 10) {
				xunit *= 2;
			}


			for (int i = region.getStart(); i < region.getStop(); i++) {
				if (i % xunit == 0) {
					int relativeStart = i - region.getStart();
					double percStart = (double) relativeStart / regionSize;
					int pixelStart = (int) Math.ceil(percStart * nrPixels);
					g2d.drawLine(margin + pixelStart, plotStarty, margin + pixelStart, plotStarty + 10);
					String formattedString = decimalFormat.format(i);
					int adv = metrics.stringWidth(formattedString);
					int hgt = metrics.getHeight();

					g2d.drawString(formattedString, margin + pixelStart - (int) Math.ceil((double) adv / 2), plotStarty + 5 + 20);
				}
			}
		}


		close();

	}

	private HashSet<Feature> loadFeaturesFromVCF(String vcfFile) throws IOException {
		VCFGenotypeData genotypeData = new VCFGenotypeData(vcfFile);
		int ctr = 0;
		System.out.println("Reading: " + vcfFile);
		HashSet<Feature> output = new HashSet<Feature>();
		while (genotypeData.hasNext()) {
			VCFVariant v = genotypeData.next();
			Feature f = new Feature();
			f.setChromosome(Chromosome.parseChr(v.getChr()));
			f.setStart(v.getPos());
			f.setStop(v.getPos());


			ctr++;
			if (ctr % 100 == 0) {
				System.out.println(ctr + " variants processed");
			}

			output.add(f);
		}
		genotypeData.close();
		return output;
	}

	private Pair<ArrayList<Feature>, ArrayList<Double>> readImputedVCF(String vcfFile, double mafthreshold) throws IOException {
		VCFGenotypeData genotypeDataAfterImputation = new VCFGenotypeData(vcfFile);

		ArrayList<Feature> variants = new ArrayList<Feature>();
		ArrayList<Double> rsquared = new ArrayList<Double>();

		int ctr = 0;
		while (genotypeDataAfterImputation.hasNext()) {
			VCFVariant v = genotypeDataAfterImputation.next();

			Feature f = new Feature();
			f.setChromosome(Chromosome.parseChr(v.getChr()));
			f.setStart(v.getPos());
			f.setStop(v.getPos());
			ctr++;
			if (ctr % 100 == 0) {
				System.out.println(ctr + " variants processed");
			}


			double maf = v.getMAF();
			if (maf > mafthreshold && maf < 1) {

				Double ar2 = v.getImputationQualityScore();
				if (ar2 != null) {

					variants.add(f);
					rsquared.add(ar2);
				}
			}

		}
		return new Pair<ArrayList<Feature>, ArrayList<Double>>(variants, rsquared);
	}

	private static GTFAnnotation annot = null;

	private void drawOverlappingGenes(String gtf, Feature region, int regionSize, int nrPixels) throws IOException {
		if (annot == null) {
			annot = new GTFAnnotation(gtf);
		}
		TreeSet<Gene> genes = annot.getGeneTree();
		Gene geneStart = new Gene("", region.getChromosome(), Strand.POS, region.getStart(), region.getStart());
		Gene geneStop = new Gene("", region.getChromosome(), Strand.POS, region.getStop(), region.getStop());
		SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);

		int ylevel = 0;

		HashMap<Gene, Integer> ylevels = new HashMap<Gene, Integer>();

		HashSet<String> geneNames = new HashSet<String>();

		for (Gene g : overlappingGenes) {
			int ctr = 2;
			if (geneNames.contains(g.getName())) {
				String name = g.getName() + "_" + ctr;
				while (geneNames.contains(name)) {
					name = g.getName() + "_" + ctr;
					ctr++;
				}
				g.setName(name);
			}
			geneNames.add(g.getName());
		}


		// get gene plot positions


		HashSet<Gene> genesOverlappingAtYLevel = new HashSet<Gene>();
		genesOverlappingAtYLevel.addAll(overlappingGenes);
		int geneXmargin = 200;
		while (!genesOverlappingAtYLevel.isEmpty()) {

//			System.out.println("start with: " + genesOverlappingAtYLevel.size());
			ArrayList<Gene> allgenesAtThisYLevel = new ArrayList<Gene>();
			for (Gene g : genesOverlappingAtYLevel) {
				allgenesAtThisYLevel.add(g);
			}
			Collections.sort(allgenesAtThisYLevel, new FeatureComparator(false));
//			System.out.println(allgenesAtThisYLevel.size() + " after sorting..");


			HashSet<Gene> visitedGenes = new HashSet<Gene>();
			for (int i = 0; i < allgenesAtThisYLevel.size(); i++) {
				Gene gene1 = allgenesAtThisYLevel.get(i);
				Pair<Integer, Integer> pixelPair = getGenePlotPositions(gene1, region, regionSize, nrPixels);
				int pixelStart = pixelPair.getLeft();
				int pixelStop = pixelPair.getRight();
				Feature f1 = new Feature();
				f1.setChromosome(gene1.getChromosome());

				f1.setStart(pixelStart);


				String name = gene1.getName();
				FontMetrics metrics = g2d.getFontMetrics();
				int strWidth = metrics.stringWidth(name);
				f1.setStop(pixelStop + strWidth + 100);
//				System.out.println(ylevel + " - " + gene1.getName() + " iterating " + allgenesAtThisYLevel.size());


				if (!visitedGenes.contains(gene1)) {
					ylevels.put(gene1, ylevel);
					boolean overlap = true;

					for (int j = i + 1; (j < allgenesAtThisYLevel.size() && overlap); j++) {
						Gene gene2 = allgenesAtThisYLevel.get(j);
						Feature f2 = new Feature();
						f2.setChromosome(gene2.getChromosome());
						Pair<Integer, Integer> pixelPair2 = getGenePlotPositions(gene2, region, regionSize, nrPixels);
						int pixelStart2 = pixelPair2.getLeft();
						int pixelStop2 = pixelPair2.getRight();
						f2.setStart(pixelStart2);

						String name2 = gene2.getName();
						strWidth = metrics.stringWidth(name2);
						f2.setStop(pixelStop2 + strWidth + 100);

						if (f1.overlaps(f2)) {
//							System.out.println("overlap with " + gene2.getName());
							visitedGenes.add(gene2);
						} else {
							overlap = false;
						}
						j++;

					}
				}

			}
			genesOverlappingAtYLevel = visitedGenes;


//			System.out.println(genesOverlappingAtYLevel.size() + " - " + ylevel + "-" + ylevels.size());
			ylevel++;
		}
		// System.exit(-1);
//		System.out.println(ylevel);


		g2d.setColor(new Color(70, 67, 58));
		int marginbetweengenes = 5;
		int geneheight = 10;

		Font defaultfont = g2d.getFont();
		g2d.setFont(new Font("default", Font.BOLD, 16));
		int y0 = margin - 100 - 20;
		int x0 = margin;
		g2d.drawString("Genes", x0, y0);
		g2d.setFont(defaultfont);

		for (Gene g : overlappingGenes) {
			ylevel = ylevels.get(g);

			int y1 = margin - 100 + ((geneheight + marginbetweengenes) * ylevel);
			// System.out.println(g.toString() + "\t" + ylevel + "\t" + y1);

			ArrayList<Transcript> transcripts = g.getTranscripts();
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
					if (relativeStop > region.getStop()) {
						relativeStop = regionSize;
					}

					double percStart = (double) relativeStart / regionSize;
					double percStop = (double) relativeStop / regionSize;

					int pixelStart = (int) Math.ceil(percStart * nrPixels);
					int pixelStop = (int) Math.ceil(percStop * nrPixels);

//					System.out.println(f.toString());
//					System.out.println(pixelStart + "\t" + pixelStop);

					int exonwidth = pixelStop - pixelStart;
					if (pixelStop - pixelStart <= 0) {
						exonwidth = 1;
					}

					g2d.fillRect(margin + pixelStart, y1, exonwidth, geneheight);
				}
			}

			Pair<Integer, Integer> pixelPair = getGenePlotPositions(g, region, regionSize, nrPixels);
			int pixelStart = pixelPair.getLeft();
			int pixelStop = pixelPair.getRight();

//			System.out.println(g.toString());
//			System.out.println(pixelStart + "\t" + pixelStop);


			g2d.drawLine(margin + pixelStart, y1 + 5, margin + pixelStop, y1 + 5);
			g2d.drawString(g.getGeneId(), margin + pixelStop + 5, y1 + 10);

		}

	}

	private Pair<Integer, Integer> getGenePlotPositions(Gene g, Feature region, int regionSize, int nrPixels) {
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

		int pixelStart = (int) Math.ceil(percStart * nrPixels);
		int pixelStop = (int) Math.ceil(percStop * nrPixels);


		return new Pair<Integer, Integer>(pixelStart, pixelStop);
	}


//	public void plotAssocPvalue(String gtf,
//								String[] datasetFiles,
//								String[] datasetNames,
//								String sequencedRegionFile,
//								Feature region,
//								double mafthreshold,
//								boolean logtransform) throws IOException {
//
//		TextFile tf1 = new TextFile(sequencedRegionFile, TextFile.R);
//
//		HashSet<Feature> sequencedRegions = new HashSet<Feature>();
//
//		String[] elems = tf1.readLineElems(TextFile.tab);
//		while (elems != null) {
//			Feature f = new Feature();
//			if (elems.length > 2) {
//				f.setChromosome(Chromosome.parseChr(elems[0]));
//				f.setStart(Integer.parseInt(elems[1]));
//				f.setStop(Integer.parseInt(elems[2]));
//				sequencedRegions.add(f);
//			}
//			elems = tf1.readLineElems(TextFile.tab);
//		}
//		tf1.close();
//
//		int regionSize = region.getStop() - region.getStart();
//		int nrPixelsX = figureWidth - (2 * margin);
//
//		System.out.println(nrPixelsX);
//
//		drawOverlappingGenes(gtf, region, regionSize, nrPixelsX);
//
//		// draw sequenced regions
//		g2d.setColor(new Color(208, 83, 77));
//
//
//		// plot sequenced regions
//		Font defaultfont = g2d.getFont();
//		g2d.setFont(new Font("default", Font.BOLD, 16));
//
//		g2d.drawString("Targeted regions in sequencing", margin, margin - 20);
//
//		g2d.setFont(defaultfont);
//		for (Feature f : sequencedRegions) {
//			if (f.overlaps(region)) {
//				int start = f.getStart();
//				int featurewidth = f.getStop() - start;
//				int relativeStart = start - region.getStart();
//				if (relativeStart < 0) {
//					featurewidth -= Math.abs(relativeStart);
//					relativeStart = 0;
//				}
//
//				int relativeStop = relativeStart + featurewidth;
//				if (relativeStop > region.getStop()) {
//					relativeStop = regionSize;
//				}
//
//				double percStart = (double) relativeStart / regionSize;
//				double percStop = (double) relativeStop / regionSize;
//
//				int pixelStart = (int) Math.ceil(percStart * nrPixelsX);
//				int pixelStop = (int) Math.ceil(percStop * nrPixelsX);
//
////				System.out.println(f.toString());
////				System.out.println(pixelStart + "\t" + pixelStop);
//				int y1 = margin;
//				int y2 = margin + 10;
//				g2d.fillRect(margin + pixelStart, y1, pixelStop - pixelStart, 10);
//			}
//		}
//
//		int betweenmargin = 100;
//		int starty = margin + betweenmargin + 10;
//
//		int nrPixelsY = figureHeight - (2 * margin);
//
//		int plotHeight = (int) Math.floor((double) (nrPixelsY - 300) / datasetFiles.length);
//		int plotWidth = (int) Math.floor((double) nrPixelsX / datasetFiles.length);
//
//		int minDotSize = 2;
//
//		starty += betweenmargin + 100;
//
//
//		Color grey = new Color(70, 67, 58);
//
//
//		int plotStart = starty;
//		int q = 0;
//
//		// plot the actual p-values..
//
//		// generate colors
//		double maxPval = -Double.MAX_VALUE;
//		ArrayList<ArrayList<Pair<Integer, Double>>> allPValues = new ArrayList<ArrayList<Pair<Integer, Double>>>();
//		for (int dataset = 0; dataset < datasetFiles.length; dataset++) {
//			String file = datasetFiles[dataset];
//			if (Gpio.exists(file)) {
//				Pair<HashSet<String>, ArrayList<Pair<Integer, Double>>> pair = readTabFile(file, region);
//				ArrayList<Pair<Integer, Double>> pvals = pair.getRight();
//				for (Pair<Integer, Double> p : pvals) {
//					if (p.getRight() > maxPval) {
//						maxPval = p.getRight();
//					}
//				}
//				allPValues.add(pvals);
//			} else {
//				allPValues.add(null);
//			}
//		}
//
//		Color[] colors = new Color[datasetFiles.length];
//		for (int i = 0; i < colors.length; i++) {
//			switch (i) {
//				case 0:
//					colors[i] = new Color(70, 67, 58);
//					break;
//				case 1:
//					colors[i] = new Color(174, 164, 140);
//					break;
//				case 2:
//					colors[i] = new Color(208, 83, 77);
//					break;
//				case 3:
//					colors[i] = new Color(98, 182, 177);
//					break;
//				case 4:
//					colors[i] = new Color(116, 156, 80);
//					break;
//
//			}
//		}
//
//
//		double unit = determineUnit(maxPval);
//		double remainder = maxPval % unit;
//		System.out.println(maxPval + " maxp");
//		System.out.println(unit + " unit");
//		System.out.println(remainder + " remainder");
//
//		maxPval += (unit - remainder); // round off using unit
//		System.out.println(maxPval + " maxp");
//
//		Font originalfont = g2d.getFont();
//		g2d.setFont(new Font("default", Font.BOLD, 16));
//		FontMetrics metrics = g2d.getFontMetrics(g2d.getFont());
//
//		for (int z = allPValues.size() - 1; z > -1; z--) {
//			int plotStarty = starty + (z * plotHeight) + (z * betweenmargin);
//			ArrayList<Pair<Integer, Double>> toPlot = allPValues.get(z);
//			if (toPlot != null) {
//				g2d.setColor(colors[z]);
//
//				for (Pair<Integer, Double> p : toPlot) {
//					Integer pos = p.getLeft();
//					Double pval = p.getRight();
//
//					// x-coord
//					int relativeStart = pos - region.getStart();
//					double percStart = (double) relativeStart / regionSize;
//					int pixelStart = margin + (int) Math.ceil(percStart * nrPixelsX);
//
//					// y-coord
//
//					double yperc = pval / maxPval;
//					int pixelY = (int) Math.ceil(yperc * plotHeight);
//
//					int dotsize = 2 + (int) Math.ceil(yperc * 10);
//					if (z == 0) {
////						dotsize = 8 + (int) Math.ceil(yperc * 10);
//					}
//
//					g2d.fillOval(pixelStart - (dotsize / 2), plotStarty - pixelY - (dotsize / 2), dotsize, dotsize);
//
//
//				}
//
//
//				int adv = metrics.stringWidth(datasetNames[z]);
//				int hgt = metrics.getHeight();
//				Dimension size = new Dimension(adv + 10, hgt + 10);
//				g2d.drawString(datasetNames[z], margin, plotStarty - plotHeight - 20);
//			}
//
//
//		}
//
//
//		// plotVariantsUniqueIneachDataset coordinates
//		starty += 5;
//
//
//		// determine unit
//
//
//		double steps = maxPval / 10;
//
//		String pattern = "###,###,###.##";
//		DecimalFormat decimalFormat = new DecimalFormat(pattern);
//
//
//		starty += 5;
//		for (int z = allPValues.size() - 1; z > -1; z--) {
//			int plotStarty = starty + (z * plotHeight) + (z * betweenmargin);
//
//
//			// draw red line near 5E-8)
//			double gwas = -Math.log10(5E-8);
//			if (maxPval >= gwas) {
//				double yperc = gwas / maxPval;
//				int pixelY = (int) Math.ceil(yperc * plotHeight);
//				g2d.setColor(new Color(208, 83, 77));
//				g2d.drawLine(margin, plotStarty - pixelY, margin + nrPixelsX, plotStarty - pixelY);
//			}
//
//			g2d.setFont(originalfont);
//			g2d.setColor(new Color(70, 67, 58));
//			// y-axis
//			g2d.drawLine(margin - 10, plotStarty, margin - 10, plotStarty - plotHeight);
//
//			for (double i = 0; i < maxPval + steps; i += steps) {
//				if (i <= maxPval) {
//					int plusY = (int) Math.ceil(((double) i / maxPval) * plotHeight);
//					g2d.drawLine(margin - 13, plotStarty - plusY, margin - 7, plotStarty - plusY);
//					String formattedStr = decimalFormat.format(i);
//					int adv = metrics.stringWidth(formattedStr);
//					int hgt = metrics.getHeight();
//					Dimension size = new Dimension(adv + 10, hgt + 10);
//					g2d.drawString(formattedStr, margin - (int) size.getWidth() - 10, plotStarty - plusY + 5);
//				}
//
//			}
//// x-axis
//			g2d.drawLine(margin - 5, plotStarty, margin + nrPixelsX + 5, plotStarty);
//
//			int xunit = (int) Math.ceil(determineUnit(regionSize));
//			while (regionSize / xunit < 10 && xunit > 1) {
//				xunit /= 2;
//			}
//			while (regionSize / xunit > 10) {
//				xunit *= 2;
//			}
//
//
//			for (int i = region.getStart(); i < region.getStop(); i++) {
//				if (i % xunit == 0) {
//					int relativeStart = i - region.getStart();
//					double percStart = (double) relativeStart / regionSize;
//					int pixelStart = (int) Math.ceil(percStart * nrPixelsX);
//					g2d.drawLine(margin + pixelStart, plotStarty - 5, margin + pixelStart, plotStarty + 5);
//					String formattedString = decimalFormat.format(i);
//					int adv = metrics.stringWidth(formattedString);
//					int hgt = metrics.getHeight();
//
//					g2d.drawString(formattedString, margin + pixelStart - (adv / 2), plotStarty + 20);
//				}
//			}
//
//
//		}
//
//
//		starty += betweenmargin + 10;
//
//
//		close();
//
//	}

//	public Pair<HashSet<String>, ArrayList<Pair<Integer, Double>>> readTabFile(String pvaluefile, Feature region) throws IOException {
//		HashSet<String> variantHash = new HashSet<String>();
//		TextFile textfile = new TextFile(pvaluefile, TextFile.R);
//
//		ArrayList<Pair<Integer, Double>> pvaluesPerPosition = new ArrayList<Pair<Integer, Double>>();
//
//		// skip header
//		textfile.readLine();
//		String[] elems = textfile.readLineElems(TextFile.tab);
//		int pvalctr = 0;
//		while (elems != null) {
//			// Marker	Chr	Position	PValue	Odds Ratio
//			if (pvaluefile.endsWith(".tab")) {
//				Chromosome chr = Chromosome.parseChr(elems[1]);
//				if (region.getChromosome().equals(chr)) {
//					String variant = elems[0] + "_" + elems[1] + "_" + elems[2];
//					Feature f2 = new Feature();
//					f2.setChromosome(chr);
//					f2.setStart(Integer.parseInt(elems[2]));
//					f2.setStop(Integer.parseInt(elems[2]));
//					if (f2.overlaps(region)) {
//						Integer position = Integer.parseInt(elems[2]);
//						try {
//							Double pval = Double.parseDouble(elems[elems.length - 2]); // need to check position...
//							double log10 = -Math.log10(pval);
//							variantHash.add(variant);
//							pvaluesPerPosition.add(new Pair<Integer, Double>(position, log10));
//						} catch (NumberFormatException e) {
//
//						}
//						pvalctr++;
//					}
//
//				}
//			} else {
//				Chromosome chr = Chromosome.parseChr(elems[1]);
//				if (region.getChromosome().equals(chr)) {
//					String variant = elems[0] + "_" + elems[1] + "_" + elems[2];
//					Feature f2 = new Feature();
//					f2.setChromosome(chr);
//					f2.setStart(Integer.parseInt(elems[2]));
//					f2.setStop(Integer.parseInt(elems[2]));
//					if (f2.overlaps(region)) {
//						Integer position = Integer.parseInt(elems[2]);
//						Double pval = Double.parseDouble(elems[elems.length - 1]); // need to check position...
//						variantHash.add(variant);
//						pvaluesPerPosition.add(new Pair<Integer, Double>(position, pval));
//						pvalctr++;
//					}
//
//				}
//			}
//			elems = textfile.readLineElems(TextFile.tab);
//		}
//		textfile.close();
//
//		System.out.println(pvalctr + " pvalues for " + pvaluesPerPosition.size() + " positions from file: " + pvaluefile);
//
//
//		return new Pair<HashSet<String>, ArrayList<Pair<Integer, Double>>>(variantHash, pvaluesPerPosition);
//	}


	public void plot2(String[] variantFiles, String[] variantFileNames, String sequencedRegionFile, String regionFile) throws IOException {
		TextFile tf1 = new TextFile(sequencedRegionFile, TextFile.R);

		HashSet<Feature> sequencedRegions = new HashSet<Feature>();

		String[] elems = tf1.readLineElems(TextFile.tab);
		while (elems != null) {
			Feature f = new Feature();
			if (elems.length > 2) {
				f.setChromosome(Chromosome.parseChr(elems[0]));
				f.setStart(Integer.parseInt(elems[1]));
				f.setStop(Integer.parseInt(elems[2]));
				sequencedRegions.add(f);
			}
			elems = tf1.readLineElems(TextFile.tab);
		}
		tf1.close();

		TextFile tf2 = new TextFile(regionFile, TextFile.R);
		elems = tf2.readLineElems(TextFile.tab);
		Feature impregion = new Feature();
		impregion.setChromosome(Chromosome.parseChr(elems[0]));
		impregion.setStart(Integer.parseInt(elems[1]));
		impregion.setStop(Integer.parseInt(elems[2]));
		tf2.close();

		ArrayList<Feature> overlappingRegions = new ArrayList<Feature>();

		for (Feature f : sequencedRegions) {
			if (f.overlaps(impregion)) {
				overlappingRegions.add(f);
			}
		}

		int nrPixels = figureWidth - (2 * margin);

		System.out.println(nrPixels);
		System.out.println(overlappingRegions.size() + " overlapping regions");
		int pixelsBetween = 5;
		int pixelsPerRegion = (nrPixels - (overlappingRegions.size() - 1) * pixelsBetween) / overlappingRegions.size();
		System.out.println(pixelsPerRegion + " pixels per region");

		int betweenmargin = 25;
		int origStarty = margin + betweenmargin + 10;
		int dotsize = 24;
		int minDotSize = 2;
		int origstartx = margin;

		HashMap<Feature, Integer> ctr = new HashMap<Feature, Integer>();
		for (int i = 0; i < variantFiles.length; i++) {
			String variantFile = variantFiles[i];
			TextFile tf = new TextFile(variantFile, TextFile.R);
			String ln = tf.readLine();
			while (ln != null) {

				if (!ln.startsWith("#")) {
					ArrayList<VCFGenotypeFilter> filters = new ArrayList<>();
					filters.add(new GenotypeQualityFilter(30));
					filters.add(new ReadDepthFilter(10));
					VCFVariant v = new VCFVariant(ln, filters, false);


					Feature f = new Feature();

					f.setChromosome(Chromosome.parseChr(v.getChr()));
					f.setStart(v.getPos());
					f.setStop(v.getPos());
					Integer ct = ctr.get(f);
					if (ct == null) {
						ct = 0;
					}
					ct++;
					ctr.put(f, ct);
				}
				ln = tf.readLine();
			}
			tf.close();

		}
		System.out.println(ctr.size() + " variants loaded");

		int fctr = 0;
		for (Feature region : overlappingRegions) {
			int starty = origStarty;

			int startx = origstartx + (fctr * pixelsBetween) + (fctr * pixelsPerRegion);

			int regionSize = region.getStop() - region.getStart();

			for (int i = 0; i < variantFiles.length; i++) {
				starty += betweenmargin + dotsize;
				String variantFile = variantFiles[i];
				String variantFileName = variantFileNames[i];

				TextFile tf = new TextFile(variantFile, TextFile.R);

				String ln = tf.readLine();

				Color grey = new Color(0, 0, 0, 64);
				g2d.setColor(grey);

				while (ln != null) {

					if (!ln.startsWith("#")) {
						ArrayList<VCFGenotypeFilter> filters = new ArrayList<>();
						filters.add(new GenotypeQualityFilter(30));
						filters.add(new ReadDepthFilter(10));
						VCFVariant v = new VCFVariant(ln, filters, false);


						Feature f = new Feature();

						f.setChromosome(Chromosome.parseChr(v.getChr()));
						f.setStart(v.getPos());
						f.setStop(v.getPos());

						if (region.overlaps(f)) {
							double maf = v.getMAF();
							if (maf > 0 && maf < 1) {
								int size = (int) Math.ceil((maf * 2) * 10);
								int relativeStart = f.getStart() - region.getStart();
								if (relativeStart < 0) {
									relativeStart = 0;
								}

								double percStart = (double) relativeStart / regionSize;
								int pixelStart = (int) Math.ceil(percStart * nrPixels);

								int halfdotsize = dotsize / 2;
								if (size < minDotSize) {
									size = minDotSize;
								}
								int halfdsize = size / 2;

								if (ctr.get(f) == null || ctr.get(f) == 1) {
									g2d.setColor(new Color(255, 0, 0, 128));
								} else {
									g2d.setColor(grey);
								}
								g2d.fillOval(startx + pixelStart - halfdsize, starty - halfdsize, size, size);
							}
						}
					}
					ln = tf.readLine();
				}
				tf.close();
			}
		}

		close();
	}
}
