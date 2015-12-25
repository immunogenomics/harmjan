package nl.harmjanwestra.harmonics.graphics;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.harmonics.smoothing.AverageSmoothingFunction;
import nl.harmjanwestra.harmonics.smoothing.Smoother;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import umcg.genetica.util.Primitives;

import java.awt.*;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DecimalFormat;

/**
 * Created by hwestra on 6/29/15.
 */
public class HistogramPlot extends DefaultGraphics {

	public HistogramPlot(String outputFileName, int width, int height) throws FileNotFoundException, DocumentException {
		super(outputFileName, width, height);
	}

	public enum PLOTTYPE {
		BAR,
		POLY
	}

	PLOTTYPE[] types = null;
	double[][] data = null;

	public void setData(double[][] data, PLOTTYPE[] types) {
		this.types = types;
		this.data = data;
	}

	int margin = 0;

	public void setMargin(int pixels) {
		this.margin = pixels;
	}

	public void draw() throws IOException {

		double maxY = 0;
		for (int i = 0; i < data.length; i++) {
			double maxd = Primitives.max(data[i]);
			if (maxd > maxY) {
				maxY = maxd;
			}
		}

		Color grey = new Color(70, 67, 58);
		g2d.setColor(grey);

		Color[] colors = new Color[data.length];
		for (int i = 0; i < colors.length; i++) {
			switch (i) {
				case 0:
					colors[i] = new Color(70, 67, 58);
					break;
				case 1:
					colors[i] = new Color(208, 83, 77);
					break;
				case 2:
					colors[i] = new Color(174, 164, 140);
					break;
				case 3:
					colors[i] = new Color(98, 182, 177);
					break;
				case 4:
					colors[i] = new Color(116, 156, 80);
					break;

			}
		}

		int plotheight = figureHeight - (2 * margin);
		int plotwidth = figureWidth - (2 * margin);
//		System.out.println(figureWidth);

		// scale y
		double unitY = determineUnit(maxY);
		double remainder = maxY % unitY;
//		System.out.println(maxY + " maxp");
//		System.out.println(unitY + " unit");
//		System.out.println(remainder + " remainder");

		maxY += (unitY - remainder); // round off using unit
//		System.out.println(maxY + " maxp");

		int maxX = 0;
		for (int i = 0; i < data.length; i++) {
			double[] toPlot = data[i];
			PLOTTYPE type = PLOTTYPE.BAR;

			g2d.setColor(colors[i]);

			if (types != null && types[i] != null) {
				type = types[i];
			}

			if (toPlot.length > maxX) {
				maxX = toPlot.length;
			}

			int widthPerBar = 1;
//			System.out.println(toPlot.length + "\t" + plotwidth);
			if (toPlot.length > plotwidth) {
				// reduce resolution
				Smoother smoother = new Smoother(new AverageSmoothingFunction());
				toPlot = smoother.reduce(toPlot, plotwidth);
			} else {

				widthPerBar = (int) Math.floor(plotwidth / toPlot.length);

			}
			Stroke defaultStroke = g2d.getStroke();


			for (int j = 0; j < toPlot.length; j++) {


				if (type.equals(PLOTTYPE.POLY)) {
					// TODO: cubicsplines
					// solve with linsegments for now

					if (j > 0) {
						double yPerc1 = 1 - (toPlot[j - 1] / maxY);
						double yPerc2 = 1 - (toPlot[j] / maxY);
						int ypos1 = margin + (int) Math.ceil(yPerc1 * plotheight);
						int ypos2 = margin + (int) Math.ceil(yPerc2 * plotheight);

						int xpos1 = margin + ((j - 1) * widthPerBar) + (int) (Math.floor(widthPerBar / 2));
						int xpos2 = margin + (j * widthPerBar) + (int) (Math.floor(widthPerBar / 2));
						g2d.setStroke(new BasicStroke(6f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
						g2d.drawLine(xpos1, ypos1, xpos2, ypos2);
					}
				} else {
					g2d.setStroke(defaultStroke);
					double yPerc = toPlot[j] / maxY;
					int barheight = (int) Math.ceil(yPerc * plotheight);
					int xpos = margin + (j * widthPerBar);
					g2d.fillRect(xpos, margin + plotheight - barheight, widthPerBar, barheight);
				}
			}
			g2d.setStroke(defaultStroke);
		}

		g2d.setColor(grey);

		// draw y-axis
		String pattern = "###,###,###.##";
		DecimalFormat decimalFormat = new DecimalFormat(pattern);
		double steps = maxY / 10;

		int yAxisXPos = margin - 10;
		int yAxisStart = margin;
		int yAxisStop = margin + plotheight;
		g2d.drawLine(yAxisXPos, yAxisStart, yAxisXPos, yAxisStop);
		FontMetrics metrics = g2d.getFontMetrics(g2d.getFont());
		for (double i = 0; i < maxY + (steps / 2); i += steps) {
			double yPerc = i / maxY;

			int ypos = margin + (int) Math.ceil((1 - yPerc) * plotheight);
			int startx = yAxisXPos - 5;
			int stopx = yAxisXPos + 5;
			g2d.drawLine(startx, ypos, stopx, ypos);
			String formattedStr = decimalFormat.format(i);
			int adv = metrics.stringWidth(formattedStr);
			int hgt = metrics.getHeight();
			Dimension size = new Dimension(adv + 10, hgt + 10);
			g2d.drawString(decimalFormat.format(i), startx - adv - 5, ypos);

		}


		// draw x-axis
		int xAxisYPos = margin + plotheight + 10;
		int xAxisXStart = margin;
		int xAxisXStop = margin + plotwidth;
		g2d.drawLine(xAxisXStart, xAxisYPos, xAxisXStop, xAxisYPos);

		int numberXTicks = 10;
		int stepsX = maxX / numberXTicks;
		int pixelsPerStep = plotwidth / numberXTicks;
		for (int i = 0; i < numberXTicks + 1; i++) {
			int xStep = i * stepsX;
			int xpos = margin + (pixelsPerStep * i);
			g2d.drawLine(xpos, xAxisYPos - 5, xpos, xAxisYPos + 5);
			int adv = metrics.stringWidth("" + xStep);
			g2d.drawString("" + xStep, xpos - (adv / 2), xAxisYPos + 20);

		}
		super.close();

	}


}
