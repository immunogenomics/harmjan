package nl.harmjanwestra.utilities.graphics.panels;

import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.themes.DefaultTheme;
import nl.harmjanwestra.utilities.graphics.themes.Theme;

import java.awt.*;
import java.text.DecimalFormat;

/**
 * Created by hwestra on 7/16/15.
 */
public class HistogramPanel extends Panel {

	private double[][] histogram;
	private String[] datasetLabels;
	private String xAxisLabel;
	private String yAxisLabel;
	private Range dataRange;
	private Theme theme = new DefaultTheme();
	private LOG xAxisLog = LOG.NONE;
	private LOG yAxisLog = LOG.NONE;

	public HistogramPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	public void setData(int[] dataset) {
		double[][] tmpData = new double[1][dataset.length];

		for (int j = 0; j < dataset.length; j++) {
			tmpData[0][j] = dataset[j];
		}
		this.histogram = tmpData;
	}

	public void setAxisLog(LOG xAxis, LOG yAxis) {
		this.xAxisLog = xAxis;
		this.yAxisLog = yAxis;

	}

	public enum PLOTTYPE {
		BAR,
		POLY
	}

	public enum LOG {
		NONE,
		TWO,
		TEN
	}

	PLOTTYPE[] types = null;

	public void setData(int[][] histogram) {

		double[][] tmpData = new double[histogram.length][histogram[0].length];
		for (int i = 0; i < histogram.length; i++) {
			for (int j = 0; j < histogram[i].length; j++) {
				tmpData[i][j] = histogram[i][j];
			}
		}
		this.histogram = tmpData;
	}

	public void setPlotTypes(PLOTTYPE[] types) {
		this.types = types;
	}


	public void setData(double[] dataset) {
		double[][] tmpData = new double[1][dataset.length];

		for (int j = 0; j < dataset.length; j++) {
			tmpData[0][j] = dataset[j];
		}
		this.histogram = tmpData;
	}

	public void setData(double[][] histogram) {
		this.histogram = histogram;
	}

	public void setRange(Range range) {
		this.dataRange = range;
	}

	public void setLabels(String xAxis, String yAxis) {
		this.xAxisLabel = xAxis;
		this.yAxisLabel = yAxis;
	}

	public void setDatasetLabels(String[] datasetLabels) {
		this.datasetLabels = datasetLabels;
	}

	public void draw(DefaultGraphics g) {

		if (dataRange == null) {
			dataRange = new Range(histogram);
		}


		Range plotRange = new Range(dataRange.getMinX(), dataRange.getMinY(), dataRange.getMaxX(), dataRange.getMaxY());
		plotRange.round();

		Graphics2D g2d = g.getG2d();


		// plot the data...
		int nrPixelsMaxX = width - (2 * marginX);
		int nrPixelsMaxY = height - (2 * marginY);
		for (int i = 0; i < histogram.length; i++) {
			double[] toPlot = histogram[i];
			PLOTTYPE type = PLOTTYPE.BAR;

			g2d.setColor(theme.getColor(i));

			if (types != null && types[i] != null) {
				type = types[i];
			}


			double plotRangeX = plotRange.getRangeX();
			int nrPixelsPerBar = (int) Math.floor(nrPixelsMaxX / plotRangeX);
			if (nrPixelsPerBar < 1) {
				nrPixelsPerBar = 1;
			}

			double dataStepsPerBar = dataRange.getRangeX() / toPlot.length;

			for (int j = 0; j < toPlot.length; j++) {
				// get current X step
				// determine position of currentX relative to plotRange
				if (type.equals(PLOTTYPE.BAR)) {
					double percY1 = plotRange.getRelativePositionY(toPlot[j]);
					int yHeight1 = (int) Math.ceil(percY1 * nrPixelsMaxY);
					int yPos1 = y0 + marginY + nrPixelsMaxY - yHeight1;
					double currentX1 = dataRange.getMinX() + (j * dataStepsPerBar);
					double percX1 = plotRange.getRelativePositionX(currentX1);
					int pixelX1 = x0 + marginX + (int) Math.ceil(percX1 * nrPixelsMaxX);
					g2d.fillRect(pixelX1, yPos1, nrPixelsPerBar, yHeight1);
				} else {

					if (j > 0) {
						double currentX1 = dataRange.getMinX() + ((j - 1) * dataStepsPerBar);
						double currentX2 = dataRange.getMinX() + (j * dataStepsPerBar);

						double percX1 = plotRange.getRelativePositionX(currentX1);
						int pixelX1 = x0 + marginX + (int) Math.ceil(percX1 * nrPixelsMaxX) + (int) (Math.floor(nrPixelsPerBar / 2));
						double percX2 = plotRange.getRelativePositionX(currentX2);
						int pixelX2 = x0 + marginX + (int) Math.ceil(percX2 * nrPixelsMaxX) + (int) (Math.floor(nrPixelsPerBar / 2));

						double percY1 = 1 - plotRange.getRelativePositionY(toPlot[j - 1]);
						int yHeight1 = y0 + marginY + (int) Math.ceil(percY1 * nrPixelsMaxY);
						double percY2 = 1 - plotRange.getRelativePositionY(toPlot[j]);
						int yHeight2 = y0 + marginY + (int) Math.ceil(percY2 * nrPixelsMaxY);
						g2d.setStroke(theme.getThickStroke());
						g2d.drawLine(pixelX1, yHeight1, pixelX2, yHeight2);
						g2d.setStroke(theme.getStroke());
					}
				}
			}
		}

		// draw axes
		// find out where axes intersect
		// leave the x/y intersection for now
//			double yx0Perc = plotRange.getRelativePositionX(0d);
//			double xy0Perc = plotRange.getRelativePositionY(0d);

		g2d.setColor(theme.getLightGrey());
		// plot y-axis
		double tickUnitY = plotRange.getRangeY() / 10;
		String pattern = "###,###,###.##";
		DecimalFormat decimalFormat = new DecimalFormat(pattern);

		g2d.setFont(theme.getSmallFont());
		FontMetrics metrics = g2d.getFontMetrics(g2d.getFont());

		int xPosYAxis = x0 + marginX - 10;
		int yPosYAxis = y0 + marginY;
		g2d.drawLine(xPosYAxis, yPosYAxis, xPosYAxis, yPosYAxis + nrPixelsMaxY);

		int maxlen = 0;
		for (double y = plotRange.getMinY(); y < plotRange.getMaxY() + (tickUnitY / 2); y += tickUnitY) {
			double yPerc = plotRange.getRelativePositionY(y);

			int ypos = y0 + marginY + (int) Math.ceil((1 - yPerc) * nrPixelsMaxY);
			int startx = xPosYAxis - 5;
			int stopx = startx + 10;
			g2d.drawLine(startx, ypos, stopx, ypos);
			String formattedStr = decimalFormat.format(y);
			int adv = metrics.stringWidth(formattedStr);
			if (adv > maxlen) {
				maxlen = adv;
			}
			g2d.setFont(theme.getSmallFont());
			g2d.drawString(formattedStr, startx - adv - 5, ypos);


		}

		if (yAxisLabel != null) {
			g2d.setFont(theme.getLargeFont());
			// determine middle of axis
			int middle = yPosYAxis + (nrPixelsMaxY / 2);
			int lengthOfAxisStr = metrics.stringWidth(yAxisLabel);
			int halfLength = lengthOfAxisStr / 2;
			int drawy = middle + halfLength;
			int drawx = xPosYAxis - maxlen - 20;

			drawRotate(g2d, drawx, drawy, -90, yAxisLabel);
		}

		// plot x-axis
		int yPosXAxis = y0 + marginY + nrPixelsMaxY + 10;

		int xPosXAxis = x0 + marginX;
		g2d.drawLine(xPosXAxis, yPosXAxis, xPosXAxis + nrPixelsMaxX, yPosXAxis);
		double tickUnitX = plotRange.getRangeX() / 10;

		for (double x = plotRange.getMinX(); x < plotRange.getMaxX() + (tickUnitX / 2); x += tickUnitX) {
			double xPerc = plotRange.getRelativePositionX(x);
			int xpos = xPosXAxis + (int) Math.ceil(xPerc * nrPixelsMaxX);
			int starty = yPosXAxis - 5;
			int stopy = starty + 10;
			String formattedStr = decimalFormat.format(x);
			g2d.drawLine(xpos, starty, xpos, stopy);
			int adv = metrics.stringWidth(formattedStr);
			g2d.setFont(theme.getSmallFont());
			g2d.drawString(formattedStr, xpos - (adv / 2), stopy + 10);

		}

		if (xAxisLabel != null) {
			g2d.setFont(theme.getLargeFont());
			// determine middle of axis
			int middle = xPosXAxis + (nrPixelsMaxX / 2);
			int lengthOfAxisStr = metrics.stringWidth(xAxisLabel);
			int halfLength = lengthOfAxisStr / 2;
			int drawx = middle - halfLength;
			int drawy = y0 + marginY + nrPixelsMaxY + 10 + (metrics.getHeight() * 2) + 10;
			g2d.drawString(xAxisLabel, drawx, drawy);
		}

		// draw title
		if (title != null) {
			int titlePosX = x0 + marginX;
			int titlePosY = y0 + marginY - theme.getLargeFont().getSize() - 5;
			g2d.setColor(theme.getDarkGrey());
			g2d.drawString(title, titlePosX, titlePosY);
		}

		// draw legend
		if (datasetLabels != null) {

		}

	}



}
