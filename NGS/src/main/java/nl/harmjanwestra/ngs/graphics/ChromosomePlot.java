package nl.harmjanwestra.ngs.graphics;

import com.lowagie.text.DocumentException;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.Strand;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import umcg.genetica.io.text.TextFile;

import java.awt.*;
import java.awt.geom.Arc2D;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Created by hwestra on 6/6/15.
 */
public class ChromosomePlot extends DefaultGraphics {

	private int margin;

	public ChromosomePlot(String name, int pageWidth, int pageHeight) throws FileNotFoundException, DocumentException {
		super(name, pageWidth, pageHeight);
		figureWidth = pageWidth;
	}

	public void setMargin(int margin) {
		this.margin = margin;
	}

	private ArrayList<Feature> getCytobands(String cytobandfile, Chromosome chr) throws IOException {

		ArrayList<Feature> output = new ArrayList<Feature>();
		TextFile tf = new TextFile(cytobandfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

			Chromosome c = Chromosome.parseChr(elems[0]);
			if (c.equals(chr)) {
				Feature f = new Feature();
				f.setStart(Integer.parseInt(elems[1]));
				f.setStop(Integer.parseInt(elems[2]));
				f.setName(elems[4]);
				output.add(f);
			}

			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return output;
	}

	private ArrayList<Feature> readRegionFile(String f) throws IOException {
		TextFile tf = new TextFile(f, TextFile.R);

		System.out.println("reading regions: " + f);
		ArrayList<Feature> output = new ArrayList<Feature>();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 2) {
				Chromosome chr = Chromosome.parseChr(elems[0]);
				Feature feat = new Feature();
				feat.setChromosome(chr);
				feat.setStart(Integer.parseInt(elems[1]));
				feat.setStop(Integer.parseInt(elems[2]));
				output.add(feat);
			}
			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();
		return output;
	}

	public void plot(String cytobandfile, String[] gffLocusFiles, boolean onlySuggestedLoci, String sequencedRegionFile) throws IOException {

		// plot the chromosomes
		int maxChrHeight = 300;
		int chrWidth = 30;
		int nrRows = 3;

		Chromosome[] allChr = Chromosome.values();

		int ctr = 0;
		int maxSize = 0;
		int nrChr = 0;
		for (Chromosome chr : allChr) {
			if (chr.getLength() > 1) {
				nrChr++;
			}
			if (chr.getLength() > maxSize) {
				maxSize = chr.getLength();

			}
		}

		ArrayList<Feature> sequencedRegions = null;
		if (sequencedRegionFile != null) {
			sequencedRegions = readRegionFile(sequencedRegionFile);
		}

		ArrayList<ArrayList<Feature>> associatedLoci = new ArrayList<ArrayList<Feature>>();
		for (String gff : gffLocusFiles) {
			associatedLoci.add(readGFF(gff, onlySuggestedLoci));
		}

		int nrPerRow = nrChr / nrRows;
		Color[] colors = new Color[]{new Color(70, 67, 58), new Color(174, 164, 140)};

		int betweenChrMarginX = 250;
		int betweenChrMarginY = 100;

		g2d.setFont(new Font("Helvetica", Font.PLAIN, 20));


		double pixelsPerBp = (double) maxChrHeight / maxSize;
		System.out.println(pixelsPerBp);
		for (Chromosome chr : allChr) {
			// get cytobands


			if (chr.getLength() > 1) {


				int row = (int) Math.floor(ctr / nrPerRow);
				int colNr = ctr % nrPerRow;

				// get the centromere
				ArrayList<Feature> allBands = getCytobands(cytobandfile, chr);
				// draw a circle
				int starty = margin + (maxChrHeight * row) + (betweenChrMarginY * row);
				int startx = margin + (betweenChrMarginX * colNr) + (chrWidth * colNr);
				int chrHeight = (int) Math.ceil(pixelsPerBp * chr.getLength());
				System.out.println(chr.getName() + "\t" + row + "\t" + colNr + "\t" + startx + "\t" + starty + "\t" + chrHeight);

				boolean centromereseen = false;

				int bottomcentromerestop = 0;
				for (Feature f : allBands) {
					int bandStart = starty + (int) Math.ceil(f.getStart() * pixelsPerBp);
					int bandStop = starty + (int) Math.ceil(f.getStop() * pixelsPerBp);
					String name = f.getName();
					boolean draw = false;

					if (name.equals("gneg")) {

					} else if (name.equals("gpos25")) {
						g2d.setColor(new Color(70, 67, 58, 100));
						draw = true;
					} else if (name.equals("gpos50")) {
						g2d.setColor(new Color(70, 67, 58, 75));
						draw = true;
					} else if (name.equals("gpos75")) {
						g2d.setColor(new Color(70, 67, 58, 50));
						draw = true;
					} else if (name.equals("gpos100")) {

						g2d.setColor(new Color(70, 67, 58, 25));
						draw = true;
					} else if (name.equals("acen")) {
						int[] xpoints = new int[3];
						int[] ypoints = new int[3];
						g2d.setColor(colors[0]);
						xpoints[0] = startx;
						xpoints[1] = startx + chrWidth / 2;
						xpoints[2] = startx + chrWidth;

						if (!centromereseen) {
							g2d.drawLine(startx, starty, startx, bandStart);
							g2d.drawLine(startx + chrWidth, starty, startx + chrWidth, bandStart);
							ypoints[0] = bandStart;
							ypoints[1] = bandStop;
							ypoints[2] = bandStart;
							centromereseen = true;
						} else {
							bottomcentromerestop = bandStop;
							ypoints[0] = bandStop;
							ypoints[1] = bandStart;
							ypoints[2] = bandStop;
						}
						Polygon p = new Polygon(xpoints, ypoints, 3);

						g2d.fillPolygon(p);

					} else {
						System.out.println("unknown dinges: " + name);
					}

					if (draw) {

//					System.out.println(f.getStart() + "\t" + f.getStop() + "\t" + bandStart + "\t" + bandStop);
						g2d.fillRect(startx, bandStart, chrWidth, (bandStop - bandStart));

					}


				}
				ctr++;

				g2d.setColor(colors[1]);
				g2d.setFont(new Font("Helvetica", Font.BOLD, 20));
				FontMetrics metrics = g2d.getFontMetrics();
				g2d.drawString(chr.getName(), startx + (chrWidth / 2) - (metrics.stringWidth(chr.getName()) / 2), starty - 40);

				g2d.setColor(colors[0]);
				g2d.setFont(new Font("Helvetica", Font.PLAIN, 20));

				g2d.drawLine(startx, bottomcentromerestop, startx, bottomcentromerestop + (chrHeight - (bottomcentromerestop - starty)));
				g2d.drawLine(startx + chrWidth, bottomcentromerestop, startx + chrWidth, bottomcentromerestop + (chrHeight - (bottomcentromerestop - starty)));

				g2d.draw(new Arc2D.Double(startx, starty - (chrWidth / 2), chrWidth, chrWidth, 0, 180, Arc2D.OPEN));

				g2d.draw(new Arc2D.Double(startx, starty + chrHeight - (chrWidth / 2), chrWidth, chrWidth, 180, 180, Arc2D.OPEN));


				// print the loci...
				g2d.setColor(colors[0]);
				for (int dataset = 0; dataset < associatedLoci.size(); dataset++) {

					ArrayList<Feature> lociOnChr = new ArrayList<Feature>();
					ArrayList<Feature> lociForDs = associatedLoci.get(dataset);
					for (Feature feature : lociForDs) {

						if (feature.getChromosome() == null) {
							System.err.println(feature.getName() + " null chr");
						}
						if (feature.getChromosome().equals(chr)) {
							if (sequencedRegions != null) {
								boolean includeregion = false;
								for (Feature seqregion : sequencedRegions) {
									if (feature.overlaps(seqregion)) {
										includeregion = true;
									}
								}
								if (includeregion) {
									lociOnChr.add(feature);
								}
							} else {
								lociOnChr.add(feature);
							}
						}

					}

					// assign pixel locations
					for (int i = 0; i < lociOnChr.size(); i++) {
						Feature feature = lociOnChr.get(i);
						// starty + (int) Math.ceil(f.getStart() * pixelsPerBp)

						int geneY1 = starty + (int) Math.ceil(feature.getStart() * pixelsPerBp);
						System.out.println(feature.getName() + "\t" + feature.getStart() + "\t" + geneY1);
						int geneY2 = starty + (int) Math.ceil(feature.getStop() * pixelsPerBp);
						feature.setStart(geneY1);
						feature.setStop(geneY2);

					}

					Collections.sort(lociOnChr, new FeatureComparator(false));
					// plot the gene names
					int lastGeneY = starty;
					int geneMarginY = 10;

					metrics = g2d.getFontMetrics();
					int fontheight = metrics.getHeight();
					for (int i = 0; i < lociOnChr.size(); i++) {
						Feature f2 = lociOnChr.get(i);
						int start = f2.getStart();
						System.out.println(f2.getName() + "\t" + f2.getStart() + "\t" + start);
						if (start <= lastGeneY + geneMarginY + fontheight) {
							start = lastGeneY + geneMarginY + fontheight;
						}

						lastGeneY = start;
						if (dataset == 0) {
							int geneStartX = startx + chrWidth + chrWidth + 10;
							System.out.println(f2.getStart() + "\t" + start);
							g2d.drawLine(startx + chrWidth + 5, f2.getStart() - 5, startx + chrWidth + 30, start - 5);

							g2d.drawString(f2.getName(), geneStartX, start);
						} else if (dataset == 1) {
							int geneStartX = startx - chrWidth - 10;
							System.out.println(f2.getStart() + "\t" + start);
							int geneStrwidth = metrics.stringWidth(f2.getName());
							g2d.drawLine(startx - 5, f2.getStart() - 5, startx - 30, start - 5);

							g2d.drawString(f2.getName(), geneStartX - geneStrwidth, start);
						}
					}


				}
			}
		}

		close();
	}

	private ArrayList<Feature> readGFF(String gff, boolean onlySuggestedLoci) throws IOException {
		TextFile tf = new TextFile(gff, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		ArrayList<Feature> output = new ArrayList<Feature>();
		while (elems != null) {
			if (elems[0].startsWith("#")) {

			} else {

				Chromosome chr = Chromosome.parseChr(elems[0]);
				Integer start = Integer.parseInt(elems[3]);
				Integer stop = Integer.parseInt(elems[4]);
				String info = elems[8];
				String[] infoElems = info.split(";");
				boolean isCandidate = false;
				Feature f = new Feature();
				f.setChromosome(chr);
				f.setStrand(Strand.NEG);
				f.setStart(start);
				f.setStart(stop);
				for (int i = 0; i < infoElems.length; i++) {
					if (infoElems[i].startsWith("is_candidate")) {
						String[] data = infoElems[i].split("=");
						if (data[1].equals("1")) {
							isCandidate = true;
						}
					} else if (infoElems[i].startsWith("Name")) {
						String[] data = infoElems[i].split("=");
						f.setName(data[1]);
					}
				}

				if (onlySuggestedLoci) {
					if (isCandidate) {
						output.add(f);
					}
				} else {
					output.add(f);
				}


			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return output;
	}
}
