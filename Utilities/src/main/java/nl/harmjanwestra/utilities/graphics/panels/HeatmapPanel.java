package nl.harmjanwestra.utilities.graphics.panels;

import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.graphics.Range;

import java.awt.*;
import java.text.DecimalFormat;

/**
 * Created by hwestra on 9/12/15.
 */
public class HeatmapPanel extends Panel {
	private double[][] data;
	private Range range;
	private String[][] labels;
	private String[] rowLabels;
	private String[] colLabels;
	private MODE plotMode;

	public void setRange(Range range) {
		this.range = range;
	}

	public void setPlotMode(MODE plotMode) {
		this.plotMode = plotMode;
	}

	public enum MODE {
		UPPER,
		LOWER,
		FULL
	}

	public HeatmapPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	public void setData(double[][] data) {
		this.data = data;
	}

	public void setData(double[][] data, String[] rowLabels, String[] colLabels) {
		this.data = data;
		this.rowLabels = rowLabels;
		this.colLabels = colLabels;
	}

	@Override
	public void draw(DefaultGraphics g) {

		if (range == null) {
			// determine min and max
			range = determineRange(data);
		}

		Graphics2D g2d = g.getG2d();

		// plot boxes
		int plotWidthX = width - (2 * marginX);
		int plotWidthY = height - (2 * marginY);
		int boxHeight = plotWidthY / data.length;
		int boxWidth = plotWidthX / data[0].length;

		int startY = marginY + y0;
		int startX = marginX + x0;

		if (plotMode.equals(MODE.FULL)) {
			for (int i = 0; i < data.length; i++) {

				int dy = (i * boxHeight) + startY;

				for (int j = 0; j < data[i].length; j++) {
					int dx = startX + (j * boxWidth);
					double v = data[i][j];
					double yPerc = range.getRelativePositionY(v);

					// generate a color
					int op = (int) Math.ceil(253 * yPerc);

					if (op > 255) {
						System.out.println("ERRORRR: " + yPerc + "\t" + op);
						System.exit(-1);
					}
					Color c = new Color(0, 128, 255, op);

					g2d.setColor(c);

					g2d.fillRect(dx, dy, boxWidth, boxHeight);

				}
			}
		} else if (plotMode.equals(MODE.UPPER)) {
			for (int i = 0; i < data.length; i++) {

				int dy = (i * boxHeight) + startY;

				for (int j = i + 1; j < data[i].length; j++) {
					int dx = startX + (j * boxWidth);
					double v = data[i][j];
					double yPerc = range.getRelativePositionY(v);

					// generate a color
					int op = (int) Math.ceil(253 * yPerc);

					if (op > 255) {
						System.out.println("ERRORRR: " + yPerc + "\t" + op);
						System.exit(-1);
					}
					Color c = new Color(0, 128, 255, op);

					g2d.setColor(c);

					g2d.fillRect(dx, dy, boxWidth, boxHeight);

				}
			}
		} else if (plotMode.equals(MODE.LOWER)) {
			for (int i = 0; i < data.length; i++) {

				int dy = (i * boxHeight) + startY;

				for (int j = 0; j < i + 1; j++) {
					int dx = startX + (j * boxWidth);
					double v = data[i][j];
					double yPerc = range.getRelativePositionY(v);

					// generate a color
					int op = (int) Math.ceil(253 * yPerc);

					if (op > 255) {
						System.out.println("ERRORRR: " + yPerc + "\t" + op);
						System.exit(-1);
					}
					Color c = new Color(0, 128, 255, op);

					g2d.setColor(c);

					g2d.fillRect(dx, dy, boxWidth, boxHeight);

				}
			}
		}


		// draw a small legend
		double rangeTicks = range.getRangeY() / 5;
		int tickNr = 0;
		String pattern = "###,###,###.##";
		DecimalFormat decimalFormat = new DecimalFormat(pattern);
		for (double i = range.getMinY(); i < range.getMaxY() + rangeTicks; i += rangeTicks) {

			int dy = y0 - (boxHeight * tickNr);

			double yPerc = range.getRelativePositionY(i);
			int op = (int) Math.ceil(253 * yPerc);


			Color c = new Color(0, 128, 255, op);

			g2d.setColor(c);

			int dx = x0 + marginX - 50;
			g2d.fillRect(dx, dy, boxWidth, boxHeight);

			g2d.setColor(theme.getDarkGrey());
			g2d.drawString(decimalFormat.format(i), dx + boxWidth + 5, dy);
			tickNr++;


		}

		// draw the labels, if any..
		if (colLabels != null || rowLabels != null) {
			g2d.setColor(theme.getDarkGrey());
			g2d.setFont(theme.getSmallFont());
			FontMetrics metrics = g2d.getFontMetrics(g2d.getFont());
			int fontheight = metrics.getHeight();
			if (colLabels != null) {
				int y = startY - 5;
				for (int col = 0; col < colLabels.length; col++) {
					int x = startX + (col * boxWidth);
					String str = colLabels[col];
					int widthOfStr = metrics.stringWidth(str);
					drawRotate(g2d, x + fontheight, y, -90, str);
				}
			}

			if (rowLabels != null) {
				int x = startX - 5;
				for (int row = 0; row < colLabels.length; row++) {
					int y = startY + (row * boxHeight) + fontheight;
					String str = rowLabels[row];
					int widthOfStr = metrics.stringWidth(str);
					g2d.drawString(str, x - widthOfStr, y);
				}
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

	private Range determineRange(double[][] data) {
		double max = -Double.MAX_VALUE;
		double min = Double.MAX_VALUE;
		for (int i = 0; i < data.length; i++) {
			for (int j = 0; j < data[i].length; j++) {
				if (data[i][j] > max) {
					max = data[i][j];
				}
				if (data[i][j] < min) {
					min = data[i][j];
				}
			}
		}

		return new Range(0, min, 0, max);
	}

}
