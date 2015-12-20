package nl.harmjanwestra.utilities.graphics.panels;

import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.graphics.Range;

import java.awt.*;
import java.text.DecimalFormat;

/**
 * Created by hwestra on 10/15/15.
 */
public class ScatterplotPanel extends Panel {

	double[] x;
	double[] y;

	private Range dataRange;
	private String xAxisLabel;
	private String yAxisLabel;

	public ScatterplotPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	public void setData(double[] x, double[] y) {
		this.x = x;
		this.y = y;
	}

	public void setLabels(String xAxis, String yAxis) {
		this.xAxisLabel = xAxis;
		this.yAxisLabel = yAxis;
	}

	@Override
	public void draw(DefaultGraphics g) {


		Graphics2D g2d = g.getG2d();

		// determine range
		if (dataRange == null) {
			dataRange = new Range(x, y);
			dataRange.round();
		}

		Range plotRange = dataRange;
		// plot the points

		int nrPixelsMaxX = width - (2 * marginX);
		int nrPixelsMaxY = height - (2 * marginY);

		g2d.setColor(theme.getDarkGrey());

		for (int i = 0; i < x.length; i++) {
			double xval = x[i];
			double yval = y[i];

			double xperc = dataRange.getRelativePositionX(xval);
			double yperc = dataRange.getRelativePositionY(yval);


			int pixelX = x0 + marginX + (int) Math.ceil(nrPixelsMaxX * xperc);
			int pixelY = y0 + marginY + (int) Math.ceil(nrPixelsMaxY - (nrPixelsMaxY * yperc));


			g2d.fillOval(pixelX - 1, pixelY - 1, 2, 2);
		}

		// axis labels
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

		// X axis
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

	}
}
