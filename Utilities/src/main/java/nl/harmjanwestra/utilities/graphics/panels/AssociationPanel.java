package nl.harmjanwestra.utilities.graphics.panels;

import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.themes.DefaultTheme;
import umcg.genetica.containers.Pair;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 10/27/15.
 */
public class AssociationPanel extends Panel {

	Feature region;
	HashSet<Feature> sequencedRegions;
	ArrayList<ArrayList<Pair<Integer, Double>>> allPValues;
	String[] datasetNames;
	private ArrayList<Pair<Integer, Double>> ld;
	private boolean[][] markDifferentColor;

	public AssociationPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	public void setDataSingleDs(Feature region,
								HashSet<Feature> sequencedRegions,
								ArrayList<Pair<Integer, Double>> allPValues,
								String datasetName) {
		this.region = region;
		this.sequencedRegions = sequencedRegions;
		ArrayList<ArrayList<Pair<Integer, Double>>> tmp = new ArrayList<ArrayList<Pair<Integer, Double>>>();
		tmp.add(allPValues);
		this.allPValues = tmp;
		this.datasetNames = new String[]{datasetName};
	}

	public void setData(Feature region, HashSet<Feature> sequencedRegions, ArrayList<ArrayList<Pair<Integer, Double>>> allPValues, String[] datasetNames) {
		this.region = region;
		this.sequencedRegions = sequencedRegions;
		this.allPValues = allPValues;
		this.datasetNames = datasetNames;
	}

	public void setMaxPVal(double d) {
		maxPval = d;
	}

	double maxPval = Double.NaN;

	boolean plotGWASSignificance = true;

	public void setPlotGWASSignificance(boolean b) {
		plotGWASSignificance = b;
	}

	public double getMaxP() {
		return maxPval;
	}

	public String getTitle() {
		return title;
	}

	@Override
	public void draw(DefaultGraphics g) {

		Graphics2D g2d = g.getG2d();

		int figureWidth = width;

		int figureHeight = height;
		int nrPixelsY = figureHeight - (2 * marginY);
		int nrDatasets = allPValues.size();

		int minDotSize = 2;
		int regionSize = region.getStop() - region.getStart();
		int nrPixelsX = figureWidth - (2 * marginX);

		int plotStarty = y0 + marginY + nrPixelsY;

		System.out.println(nrPixelsX);

		// drawOverlappingGenes(gtf, region, regionSize, nrPixelsX);

		// draw sequenced regions
		g2d.setColor(new Color(208, 83, 77));
		Color highlight = new Color(208, 83, 77);


		// plot sequenced regions
		Font defaultfont = g2d.getFont();
		g2d.setFont(new Font("default", Font.BOLD, 16));

		g2d.drawString("Targeted regions in sequencing", marginX, marginY - 20);

		g2d.setFont(defaultfont);
		if (sequencedRegions != null) {
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

					int pixelStart = x0 + marginX + (int) Math.ceil(percStart * nrPixelsX);
					int pixelStop = x0 + marginX + (int) Math.ceil(percStop * nrPixelsX);

					int y1 = marginY + y0 - 20;

					int boxwidth = pixelStop - pixelStart;
					if (boxwidth <= 0) {
						boxwidth = 1;
					}

					g2d.fillRect(pixelStart, y1, boxwidth, 10);
				}
			}
		}

		// determine max Y

		if (Double.isNaN(maxPval)) {

//			System.out.println("Determining max P");
//			System.exit(-1);
			maxPval = -Double.MAX_VALUE;
			for (int d = 0; d < allPValues.size(); d++) {
				ArrayList<Pair<Integer, Double>> pvals = allPValues.get(d);
				for (Pair<Integer, Double> p : pvals) {
					if (p.getRight() > maxPval) {
						maxPval = p.getRight();
					}
				}
			}
		}


		Color[] colors = new Color[nrDatasets];
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

			}
		}


		Range range = new Range(0, 0, 0, 0);
		double unit = range.determineUnit(maxPval);
		double remainder = maxPval % unit;
		System.out.println(maxPval + " maxp");
		System.out.println(unit + " unit");
		System.out.println(remainder + " remainder");

		maxPval += (unit - remainder); // round off using unit
		System.out.println(maxPval + " maxp");

		Font originalfont = g2d.getFont();
		g2d.setFont(new Font("default", Font.BOLD, 12));
		FontMetrics metrics = g2d.getFontMetrics(g2d.getFont());

		if (ld != null) {
			DefaultTheme theme = new DefaultTheme();

			Color currentColor = g2d.getColor();
			g2d.setColor(theme.getLightGrey());

			for (int q = 0; q < ld.size(); q++) {
				Pair<Integer, Double> d = ld.get(q);
				Integer pos = d.getLeft();
				Double val = d.getRight();
				// x-coord
				int relativeStart = pos - region.getStart();
				double percStart = (double) relativeStart / regionSize;
				int pixelStartX = x0 + marginX + (int) Math.ceil(percStart * nrPixelsX);

				// y-coord


				int pixelY = (int) Math.ceil(val * nrPixelsY);

				int dotsize = 2;
				g2d.fillOval(pixelStartX - (dotsize / 2), plotStarty - pixelY - (dotsize / 2), dotsize, dotsize);// x-coord
			}
			g2d.setColor(currentColor);
		}


		// plot the values
		for (int z = allPValues.size() - 1; z > -1; z--) {
			ArrayList<Pair<Integer, Double>> toPlot = allPValues.get(z);
			boolean[] mark = null;
			if (markDifferentColor != null) {
				mark = markDifferentColor[z];
			}
			if (toPlot != null) {
				g2d.setColor(colors[z]);

				for (int v = 0; v < toPlot.size(); v++) {
					Pair<Integer, Double> p = toPlot.get(v);
					Integer pos = p.getLeft();
					Double pval = p.getRight();

					// x-coord
					int relativeStart = pos - region.getStart();
					double percStart = (double) relativeStart / regionSize;
					int pixelStartX = x0 + marginX + (int) Math.ceil(percStart * nrPixelsX);

					// y-coord

					double yperc = pval / maxPval;
					int pixelY = (int) Math.ceil(yperc * nrPixelsY);

					int dotsize = 2 + (int) Math.ceil(yperc * 10);
					if (z == 0) {
//						dotsize = 8 + (int) Math.ceil(yperc * 10);
					}

					if (mark != null && mark[v]) {
						Color col = g2d.getColor();
						g2d.setColor(highlight);
						g2d.fillOval(pixelStartX - (dotsize / 2), plotStarty - pixelY - (dotsize / 2), dotsize, dotsize);
						g2d.setColor(col);
					} else {
						g2d.fillOval(pixelStartX - (dotsize / 2), plotStarty - pixelY - (dotsize / 2), dotsize, dotsize);
					}


				}

				int adv = metrics.stringWidth(datasetNames[z]);
				int hgt = metrics.getHeight();
				Dimension size = new Dimension(adv + 10, hgt + 10);
				g2d.drawString(datasetNames[z], x0 + marginX, y0 + marginY - 30);
			}
		}

		// determine unit
		double steps = maxPval / 10;

		// draw red line near 5E-8)
		if (plotGWASSignificance) {
			double gwas = -Math.log10(5E-8);
			if (maxPval >= gwas) {
				double yperc = gwas / maxPval;
				int pixelY = (int) Math.ceil(yperc * nrPixelsY);
				g2d.setColor(new Color(208, 83, 77));
				g2d.drawLine(x0 + marginX, plotStarty - pixelY, x0 + marginX + nrPixelsX, plotStarty - pixelY);
			}
		}

		String pattern = "###,###,###.##";
		DecimalFormat decimalFormat = new DecimalFormat(pattern);


		g2d.setFont(originalfont);
		g2d.setColor(new Color(70, 67, 58));
		// y-axis
		g2d.drawLine(x0 + marginX - 10,
				plotStarty,
				x0 + marginX - 10,
				y0 + marginY);

		// tick lines
		for (double i = 0; i < maxPval + steps; i += steps) {
			if (i <= maxPval) {
				int plusY = (int) Math.ceil(((double) i / maxPval) * nrPixelsY);
				g2d.drawLine(x0 + marginX - 13, plotStarty - plusY, x0 + marginX - 7, plotStarty - plusY);
				String formattedStr = decimalFormat.format(i);
				int adv = metrics.stringWidth(formattedStr);
				int hgt = metrics.getHeight();
				Dimension size = new Dimension(adv + 10, hgt + 10);
				g2d.drawString(formattedStr, x0 + marginX - (int) size.getWidth() - 10, plotStarty - plusY + 5);
			}
		}


// x-axis
		g2d.drawLine(x0 + marginX - 5, plotStarty, x0 + marginX + nrPixelsX + 5, plotStarty);

		int xunit = (int) Math.ceil(range.determineUnit(regionSize));
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
				int pixelStart = (int) Math.ceil(percStart * nrPixelsX);
				g2d.drawLine(x0 + marginX + pixelStart, plotStarty - 5, x0 + marginX + pixelStart, plotStarty + 5);
				String formattedString = decimalFormat.format(i);
				int adv = metrics.stringWidth(formattedString);
				int hgt = metrics.getHeight();

				g2d.drawString(formattedString, x0 + marginX + pixelStart - (adv / 2), plotStarty + 20);
			}
		}


	}

	public void setLdInfo(ArrayList<Pair<Integer, Double>> ldsqr) {
		this.ld = ldsqr;
	}

	public void setMarkDifferentColor(boolean[] markDifferentColor) {
		this.markDifferentColor = new boolean[1][0];
		this.markDifferentColor[0] = markDifferentColor;
	}
}
