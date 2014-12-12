package nl.harmjanwestra.ngs.graphics;

import com.lowagie.text.DocumentException;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import umcg.genetica.util.Primitives;

import java.awt.*;
import java.io.FileNotFoundException;

/**
 * Created by hwestra on 12/11/14.
 */
public class LocusCoveragePlot extends DefaultGraphics {


	private int plotIndividualHeight;
	private int plotIndividualWidth;
	private int betweenPlotMargin;
	private int pageMargin;
	private String[] rowLabels;
	private int[][] coverageMapPosStr;
	private int[][] coverageMapNegStr;
	private Feature feature;
	private int featureMargin;
	private int maxCoverage = -1;

	public LocusCoveragePlot(String name, int pageWidth, int pageHeight) throws FileNotFoundException, DocumentException {
		super(name, pageWidth, pageHeight);
	}


	public void setMargin(int pageMargin) {
		this.pageMargin = pageMargin;
	}

	public void setBetweenPlotMargin(int betweenPlotMargin) {
		this.betweenPlotMargin = betweenPlotMargin;
	}

	public void setPlotIndividualWidth(int plotIndividualWidth) {
		this.plotIndividualWidth = plotIndividualWidth;
	}

	public void setPlotIndividualHeight(int plotIndividualHeight) {
		this.plotIndividualHeight = plotIndividualHeight;
	}

	public void setRowLabels(String[] rowLabels) {
		this.rowLabels = rowLabels;
	}

	public void setCoverageData(int[] coverageMapPosStrand, int[] coverageMapNegStrand) {
		this.coverageMapPosStr = new int[1][];
		this.coverageMapNegStr = new int[1][];
		this.coverageMapPosStr[0] = coverageMapPosStrand;
		this.coverageMapNegStr[0] = coverageMapNegStrand;
	}

	public void setCoverageData(int[][] coverageMapPosStrand, int[][] coverageMapNegStrand) {
		this.coverageMapPosStr = coverageMapPosStrand;
		this.coverageMapNegStr = coverageMapNegStrand;
	}

	public void setFeature(Feature f) {
		this.feature = f;
	}

	public void setFeatureMargin(int featureMargin) {
		this.featureMargin = featureMargin;
	}

	public void setMaxCoverage(int maxCoverage) {
		this.maxCoverage = maxCoverage;
	}

	public void draw() {

		int nrSamples = coverageMapPosStr.length;
		int x0 = pageMargin;

		// determine scaling
		int featureSize = feature.getStop() - feature.getStart();
		int featureSizePlusMargin = featureSize + (2 * featureMargin);
		double marginPercentOfTotalWindow = featureMargin / featureSizePlusMargin;
		int nrMarginPixels = (int) Math.floor(marginPercentOfTotalWindow * plotIndividualWidth);

		for (int row = 0; row < nrSamples; row++) {

			int y0 = pageMargin + (row * betweenPlotMargin) + (row * plotIndividualHeight);

			// plot a grey line halfway (if strandedness)
			if (coverageMapNegStr != null) {
				int halfway = y0 + (plotIndividualHeight / 2);
				g2d.setColor(Color.gray);
				g2d.setStroke(line);
				g2d.drawLine(x0, halfway, x0 + plotIndividualWidth, halfway);
			}

			if (rowLabels != null) {
				g2d.setColor(new Color(0, 0, 0));
				g2d.setFont(SMALL_FONT);
				g2d.drawString(rowLabels[row], x0, y0 - 10);
			}

			// indicate the buffer zone
			g2d.setStroke(dashed);
			g2d.setColor(Color.lightGray);
			g2d.drawLine(x0 + nrMarginPixels, y0, x0 + nrMarginPixels, y0 + plotIndividualHeight);
			g2d.drawLine(x0 + plotIndividualWidth + nrMarginPixels, y0, x0 + plotIndividualWidth + nrMarginPixels, y0 + plotIndividualHeight);

			// plot the bins for this sample
			int[] negStr = coverageMapNegStr[row];
			int[] posStr = coverageMapPosStr[row];

			// determine max if max has not been set (although may now vary between samples)
			if (maxCoverage == -1) {
				maxCoverage = Primitives.max(posStr);
				if (coverageMapNegStr != null) {
					int maxNegStr = Primitives.max(negStr);
					if (maxNegStr > maxCoverage) {
						maxCoverage = maxNegStr;
					}
				}
			}

			for (int bin = 0; bin < posStr.length; bin++) {
				double binPerc = (double) bin / posStr.length;
				int binX0 = (int) Math.floor(binPerc * plotIndividualWidth);

				if (coverageMapNegStr != null) {
					// plot both str
					int maxPlotPixels = plotIndividualHeight / 2;

					double percentOfMaxPos = posStr[bin] / maxCoverage;
					double percentOfMaxNeg = negStr[bin] / maxCoverage;

					int pixelsPos = (int) Math.floor(percentOfMaxPos * maxPlotPixels);
					int pixelsNeg = (int) Math.floor(percentOfMaxNeg * maxPlotPixels);
					if (pixelsPos > maxPlotPixels) {
						pixelsPos = maxPlotPixels;
					}
					if (pixelsNeg > maxPlotPixels) {
						pixelsNeg = maxPlotPixels;
					}
					int halfway = y0 + (plotIndividualHeight / 2);

					g2d.setColor(new Color(255, 128, 0));
					g2d.fillOval(binX0 - 1, halfway - pixelsPos - 1, 2, 2);

					g2d.setColor(new Color(0, 128, 255));
					g2d.fillOval(binX0 - 1, halfway + pixelsNeg - 1, 2, 2);

				} else {
					// only plot positive str
					int maxPlotPixels = plotIndividualHeight;
					double percentOfMaxPos = posStr[bin] / maxCoverage;
					int pixelsPos = (int) Math.floor(percentOfMaxPos * maxPlotPixels);
					if (pixelsPos > maxPlotPixels) {
						pixelsPos = maxPlotPixels;
					}
					int halfway = y0 + (plotIndividualHeight);
					g2d.setColor(new Color(255, 128, 0));
					g2d.fillOval(binX0 - 1, halfway - pixelsPos - 1, 2, 2);
				}
			}
		}
	}
}
