package nl.harmjanwestra.ngs.graphics;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.features.Exon;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.features.Transcript;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.graphics.themes.DefaultTheme;
import nl.harmjanwestra.utilities.legacy.genetica.util.Primitives;

import java.awt.*;
import java.io.FileNotFoundException;
import java.util.ArrayList;

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
	private Gene[] genes;

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

		// draw genes
		g2d.setColor(new Color(0, 0, 0));
		int genesheight = 0;
		if (genes != null) {
			int nrGenes = genes.length;

			int y0 = pageMargin;

			for (int i = 0; i < nrGenes; i++) {

				int y1 = y0 + (i * 15);
				Gene g = genes[i];
				ArrayList<Transcript> transcripts = g.getTranscripts();
				g2d.drawString(g.getName(), x0, y1);
				double pixelPerBP = (double) plotIndividualWidth / (feature.getStop() - feature.getStart());
				for (Transcript t : transcripts) {
					ArrayList<Exon> exons = t.getExons();
					int tsRelativeStart = t.getStart() - feature.getStart();
					if (tsRelativeStart < 0) {
						tsRelativeStart = 0;
					}
					int tsRelativeStop = t.getStop() - feature.getStop();
					if (tsRelativeStart > 0) {
						tsRelativeStop = 0;
					}

					g2d.drawLine(x0 + tsRelativeStart, y1 + 5, x0 + plotIndividualWidth + tsRelativeStop, y1 + 5);
					for (Exon e : exons) {
						int start = e.getStart();
						int stop = e.getStop();
						if (start < feature.getStart()) {
							start = feature.getStart();
						}
						if (stop > feature.getStop()) {
							stop = feature.getStop();
						}

						if (start <= feature.getStop()) {
							int nrBpForExon = stop - start;
							int nrPixelsWidth = (int) Math.ceil(nrBpForExon * pixelPerBP);
							int relativeBpStart = start - feature.getStart();
							int nrPixelsStart = (int) Math.ceil(relativeBpStart * pixelPerBP);
							g2d.fillRect(x0 + nrPixelsStart, y1, nrPixelsWidth, 10);
						}
					}
				}

			}
			genesheight = (genes.length * 25);

		}

		DefaultTheme theme = new DefaultTheme();

		// determine scaling
		for (int row = 0; row < nrSamples; row++) {


			int y0 = pageMargin + (row * betweenPlotMargin) + (row * plotIndividualHeight) + genesheight;


			// plotVariantsUniqueIneachDataset a grey line halfway (if strandedness)
			if (coverageMapNegStr != null) {
				int halfway = y0 + (plotIndividualHeight / 2);
				g2d.setColor(Color.gray);
				g2d.setStroke(theme.getStroke());
				g2d.drawLine(x0, halfway, x0 + plotIndividualWidth, halfway);
			}


			// indicate the buffer zone
			g2d.setStroke(theme.getStrokeDashed());
			g2d.setColor(Color.lightGray);
			g2d.drawLine(x0, y0, x0, y0 + plotIndividualHeight);
			g2d.drawLine(x0 + plotIndividualWidth, y0, x0 + plotIndividualWidth, y0 + plotIndividualHeight);

			// plotVariantsUniqueIneachDataset the bins for this sample
			int[] negStr = null;
			if (coverageMapNegStr != null) {
				negStr = coverageMapNegStr[row];
			}

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

			if (rowLabels != null) {
				g2d.setColor(new Color(0, 0, 0));
				g2d.setFont(theme.getSmallFont());
				g2d.drawString(rowLabels[row] + " max coverage: " + maxCoverage, x0, y0 - 10);
			}


			for (int bin = 0; bin < posStr.length; bin++) {
				double binPerc = (double) bin / posStr.length;
				int binX0 = x0 + (int) Math.floor(binPerc * plotIndividualWidth);

				if (coverageMapNegStr != null) {
					// plotVariantsUniqueIneachDataset both str
					int maxPlotPixels = plotIndividualHeight / 2;

					double percentOfMaxPos = (double) posStr[bin] / maxCoverage;
					double percentOfMaxNeg = (double) negStr[bin] / maxCoverage;

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
					// only plotVariantsUniqueIneachDataset positive str
					int maxPlotPixels = plotIndividualHeight;
					double percentOfMaxPos = (double) posStr[bin] / maxCoverage;
					int pixelsPos = (int) Math.floor(percentOfMaxPos * maxPlotPixels);
					if (pixelsPos > maxPlotPixels) {
						pixelsPos = maxPlotPixels;
					}
					int halfway = y0 + (plotIndividualHeight);
					g2d.setColor(new Color(255, 128, 0));
					g2d.fillRect(binX0, halfway - pixelsPos, 1, pixelsPos);
				}
			}
		}
	}

	public void setGenes(Gene[] genes) {
		this.genes = genes;
	}
}
