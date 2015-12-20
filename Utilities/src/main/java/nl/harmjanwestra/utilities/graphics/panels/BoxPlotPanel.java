package nl.harmjanwestra.utilities.graphics.panels;

import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.themes.DefaultTheme;
import nl.harmjanwestra.utilities.graphics.themes.Theme;
import umcg.genetica.io.text.TextFile;

import java.awt.*;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Arrays;

/**
 * Created by hwestra on 9/12/15.
 */
public class BoxPlotPanel extends Panel {

	private double[][] data;
	private Range range;
	private boolean drawDataPoints;
	private boolean useMeanAndSd;
	private String outputIQRS;
	private String[] labels;

	@Override
	public void draw(DefaultGraphics g) {

		double[][] iqrs = new double[data.length][0];
		double max = -Double.MAX_VALUE;
		double min = Double.MAX_VALUE;
		for (int i = 0; i < data.length; i++) {
			iqrs[i] = collectStats(data[i]);
			double localmin = iqrs[i][0];
			double localmax = iqrs[i][iqrs[i].length - 1];
//			System.out.println(localmin + "\t" + localmax);
			if (localmin < min) {
				min = localmin;
			}
			if (localmax > max) {
				max = localmax;
			}
		}


		if (range == null) {
			range = new Range(0, min, 0, max);
		}


		int pixelsY = height - (2 * marginY);


		// plot the data
		int marginBetweenBoxes = 10;

		int nrMargins = data.length - 1;
		int marginwidthTotal = nrMargins * marginBetweenBoxes;
		int widthPerBox = (int) Math.ceil((double) (width - marginwidthTotal) / data.length);
		int starty = y0 + marginY;


		Theme theme = new DefaultTheme();

		Color c = theme.getDarkGrey();
		Graphics2D g2d = g.getG2d();
		int halfbox = widthPerBox / 2;

		try {
			TextFile iqrout = null;
			if (outputIQRS != null) {
				iqrout = new TextFile(outputIQRS, TextFile.W);
			}


			for (int i = 0; i < data.length; i++) {

				int startX = x0 + marginX + (i * marginBetweenBoxes) + (i * widthPerBox);


				// determine IQR
				double m1 = iqrs[i][0];
				double q1 = iqrs[i][1];
				double q2 = iqrs[i][2];
				double q3 = iqrs[i][3];
				double m2 = iqrs[i][4];

				if (drawDataPoints) {
					for (int j = 0; j < data[i].length; j++) {
						// determine y Position

						double d = data[i][j];
						if (d > range.getMaxY()) {
							d = range.getMaxY();
						}
						if (d < range.getMinY()) {
							d = range.getMinY();
						}

						// determine where this point falls in the range
						double yperc = range.getRelativePositionY(d);
						int relY = (int) Math.ceil(yperc * pixelsY);

						int plotY = starty + pixelsY - relY;

						// add some jitter to the x-direction
						// should be dependent on density

						int direction = 1;
						if (Math.random() > 0.5) {
							direction = -1;
						}
						int jittersize = (int) Math.ceil(Math.random() * Math.random() * halfbox);
						int jitter = jittersize * direction;
						int plotX = startX + halfbox + jitter;
						g2d.setColor(c);
						if (d > q1 && d < q3) {
							// falls within IQR

							// set opacity to 30%;

							Color c2 = new Color(70, 67, 58, 128);
							g2d.setColor(c2);

						} else {
							g2d.setColor(c);
						}

						g2d.fillOval(plotX - 2, plotY - 2, 4, 4);

					}
				} else {
					// just draw the box plot
					g2d.setColor(theme.getLightGrey());
					g2d.setStroke(theme.getStroke());
					// draw line from min to q1

					boolean clippingbottom = false;
					boolean clippingtop = false;

					if (m1 < range.getMinY()) {
						m1 = range.getMinY();
						clippingbottom = true;
					}

					if (m2 > range.getMaxY()) {
						m2 = range.getMaxY();
						clippingtop = true;
					}


					double m1y = range.getRelativePositionY(m1);
					double q1y = range.getRelativePositionY(q1);

					double m2y = range.getRelativePositionY(m2);
					double q3y = range.getRelativePositionY(q3);


					int m1yPx = (int) Math.ceil(m1y * pixelsY);
					int q1yPx = (int) Math.ceil(q1y * pixelsY);


					int m2yPx = (int) Math.ceil(m2y * pixelsY);
					int q3yPx = (int) Math.ceil(q3y * pixelsY);

					int plotYm1y = starty + pixelsY - m1yPx;
					int plotYq1y = starty + pixelsY - q1yPx;


					if (clippingbottom) {
						g2d.setStroke(theme.getStrokeDashed());
						g2d.drawLine(startX + halfbox - halfbox, starty + pixelsY, startX + halfbox + halfbox, starty + pixelsY);
					}
					g2d.setStroke(theme.getStroke());

					g2d.drawLine(startX + halfbox, plotYm1y, startX + halfbox, plotYq1y);


					int plotYm2y = starty + pixelsY - m2yPx;
					int plotYq3y = starty + pixelsY - q3yPx;

					if (clippingtop) {
						g2d.setStroke(theme.getStrokeDashed());
						g2d.drawLine(startX + halfbox - halfbox, starty, startX + halfbox + halfbox, starty);
					}
					g2d.setStroke(theme.getStroke());
					g2d.drawLine(startX + halfbox, plotYm2y, startX + halfbox, plotYq3y);


					// draw the median
					double q2y = range.getRelativePositionY(q2);
					int q2yPx = (int) Math.ceil(q2y * pixelsY);
					g2d.setColor(theme.getColor(1));
					int plotYq2y = starty + pixelsY - q2yPx;
					g2d.fillOval(startX + halfbox - 3, plotYq2y - 3, 6, 6);
					g2d.setColor(theme.getColor(0));

					// draw a lightgrey box
					g2d.setColor(theme.getLightGrey());
					//g2d.drawRect(startX, plotYq3y, widthPerBox, plotYq1y - plotYq3y);
					if (outputIQRS != null) {
						iqrout.writeln(m1 + "\t" + q1 + "\t" + q2 + "\t" + q3 + "\t" + m2
								+ "\t" + plotYm1y + "\t" + plotYq1y + "\t" + plotYq2y + "\t" + plotYq3y + "\t" + plotYm2y);
					}
				}


				// plot an axis


			}
			if (outputIQRS != null) {
				iqrout.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		g2d.setColor(theme.getDarkGrey());

		// plot y-axis
		double tickUnitY = range.getRangeY() / 10;
		String pattern = "###,###,###.##";
		DecimalFormat decimalFormat = new DecimalFormat(pattern);

		g2d.setFont(theme.getSmallFont());
		FontMetrics metrics = g2d.getFontMetrics(g2d.getFont());

		int xPosYAxis = x0 + marginX - 10;
		int yPosYAxis = y0 + marginY;
		g2d.drawLine(xPosYAxis, yPosYAxis, xPosYAxis, yPosYAxis + pixelsY);

		int maxlen = 0;
		for (double y = range.getMinY(); y < range.getMaxY() + (tickUnitY / 2); y += tickUnitY) {
			double yPerc = range.getRelativePositionY(y);

			int ypos = y0 + marginY + (int) Math.ceil((1 - yPerc) * pixelsY);
			int startx = xPosYAxis - 5;
			int stopx = xPosYAxis;
			g2d.drawLine(startx, ypos, stopx, ypos);
			String formattedStr = decimalFormat.format(y);
			int adv = metrics.stringWidth(formattedStr);
			if (adv > maxlen) {
				maxlen = adv;
			}
			g2d.setFont(theme.getSmallFont());
			g2d.drawString(formattedStr, startx - adv - 5, ypos);
		}


		// draw an x-axis

		// plot x-axis
		int yPosXAxis = y0 + marginY + pixelsY + 10;

		int xPosXAxis = x0 + marginX;
		int nrPixelsX = width - (2 * marginX);
		g2d.drawLine(xPosXAxis, yPosXAxis, xPosXAxis + nrPixelsX, yPosXAxis);

		int halfmargin = marginBetweenBoxes / 2;
		for (int q = 1; q < data.length; q++) {
			int startX = x0 + marginX + (q * marginBetweenBoxes) + (q * widthPerBox) - halfmargin;
			g2d.drawLine(startX, yPosXAxis, startX, yPosXAxis + 5);
		}

		if (labels != null) {
			g2d.setColor(theme.getDarkGrey());
			g2d.setFont(theme.getSmallFont());
			metrics = g2d.getFontMetrics(g2d.getFont());
			int fontheight = metrics.getHeight();

			int y = yPosXAxis + 5;
			for (int q = 0; q < labels.length; q++) {
				int startX = x0 + marginX + (q * marginBetweenBoxes) + (q * widthPerBox) + halfbox;

				String str = labels[q];
				int widthOfStr = metrics.stringWidth(str);
				drawRotate(g2d, startX + (fontheight / 2), y + widthOfStr, -90, str);
			}
		}


	}

	public void drawRotate(Graphics2D g2d, double x, double y, int angle, String text) {
		g2d.translate((float) x, (float) y);
		g2d.rotate(Math.toRadians(angle));
		g2d.drawString(text, 0, 0);
		g2d.rotate(-Math.toRadians(angle));
		g2d.translate(-(float) x, -(float) y);
	}

	public void setRange(Range range) {
		this.range = range;
	}

	public BoxPlotPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	public void setData(double[][] data) {
		this.data = data;
	}

	private double[] collectStats(double[] dataset) {

		double min = Double.MAX_VALUE;
		double max = -Double.MAX_VALUE;
		double firstquartile = 0;
		double thirdquartile = 0;
		double median = 0;

		double[] datacopy = new double[dataset.length];
		System.arraycopy(dataset, 0, datacopy, 0, dataset.length);

		boolean even = false;
		if (dataset.length % 2 == 0) {
			even = true;
		}

		Arrays.sort(datacopy);
		if (even) {
			int middle = (int) Math.ceil((double) dataset.length / 2);
			double mid1 = datacopy[middle - 1];
			double mid2 = datacopy[middle];
			median = (mid1 + mid2) / 2;
			int firstqpos = (int) Math.floor((double) dataset.length / 4);
			int thirdqpos = (int) Math.floor((double) dataset.length * .75);
			firstquartile = datacopy[firstqpos];
			thirdquartile = datacopy[thirdqpos];


		} else {
			int middle = (int) Math.floor((double) dataset.length / 2);
			int firstqpos = (int) Math.floor((double) dataset.length / 4);
			int thirdqpos = middle + firstqpos;
			median = datacopy[middle];
			firstquartile = datacopy[firstqpos];
			thirdquartile = datacopy[thirdqpos];
		}


		min = datacopy[0];
		max = datacopy[datacopy.length - 1];

		return new double[]{
				min, firstquartile, median, thirdquartile, max
		};
	}

	public void setDrawDataPoints(boolean drawDataPoints) {
		this.drawDataPoints = drawDataPoints;
	}


	public void setUseMeanAndSd(boolean useMeanAndSd) {
		this.useMeanAndSd = useMeanAndSd;
	}

	public void setOutputIQRS(String outputIQRS) {
		this.outputIQRS = outputIQRS;
	}

	public void setLabels(String[] labels) {
		this.labels = labels;
	}
}
