package nl.harmjanwestra.utilities.graphics.panels;

import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.themes.DefaultTheme;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.math.Goniometry;

import java.awt.*;
import java.awt.geom.Arc2D;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by hwestra on 9/12/16.
 */
public class CircularHeatmapPanel extends Panel {


	private double[][][] data; // format [dataset][bin][subbin]
	private String[] rownames;
	private String[] colnames;
	private ArrayList<Triple<Integer, Integer, String>> groups;

	public CircularHeatmapPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}


	public void setData(String[] rowNames, String[] columnNames, double[][][] data) {
		this.rownames = rowNames;
		this.colnames = columnNames;
		this.data = data;
	}

	public void setData(String[] rowNames, String[] columnNames, double[][] data) {
		this.rownames = rowNames;
		this.colnames = columnNames;

		this.data = new double[data.length][data[0].length][1];
		for (int r = 0; r < data.length; r++) {
			for (int c = 0; c < data[r].length; c++) {
				this.data[r][c][0] = data[r][c];
			}
		}
	}

	@Override
	public void draw(DefaultGraphics g) {

		Range rangeValues = new Range(data);

		int nrRows = data.length + 1;
//		int nrCols = data[0].length + 1;

		int plotWidth = width - (2 * marginX);
		int plotHeight = height - (2 * marginY);

		int maxWidth = plotWidth;
		if (plotHeight < plotWidth) {
			maxWidth = plotHeight;
		}


		Graphics2D g2d = g.getG2d();


		double widthPerDataset = maxWidth / nrRows;

		double startX = marginX + x0;
		double startY = marginY + y0;

		DefaultTheme theme = new DefaultTheme();

		System.out.println(marginX);
		System.out.println(marginY);
		System.out.println(x0);
		System.out.println(y0);
		System.out.println(maxWidth);

		g2d.setStroke(theme.getThickStroke());

		double degreesForRowNames = 45d;

		double angleOffSet = 90d - (degreesForRowNames / 2);
		double degreesPerSegment = (360d - degreesForRowNames) / data[0].length;

		double degreesPerGroup = 4;
		HashMap<Integer, Integer> colToGroup = null;
		if (groups != null) {
			// this is really really dumb, but what the hack.
			colToGroup = new HashMap<Integer, Integer>();
			System.out.println("found " + groups.size() + " groups");

			double groupdeg = (groups.size() - 1) * degreesPerGroup;
			degreesPerSegment = (360d - degreesForRowNames - groupdeg) / data[0].length;

			for (int group = 0; group < groups.size(); group++) {
				int s = groups.get(group).getLeft();
				int e = groups.get(group).getMiddle();
				for (int q = s; q < e; q++) {
					colToGroup.put(q, group);
				}
			}
		}

		// draw values per dataset.
		for (int dataset = 0; dataset < data.length; dataset++) {
			// plot a white circle first
			g2d.setColor(Color.white);
			int dsWidth = (int) Math.floor(maxWidth - (widthPerDataset * dataset));
			System.out.println(dataset + "\t" + dsWidth);
			double remainder = (maxWidth - dsWidth) / 2;
			int x0 = (int) Math.floor(startX + remainder);
			int y0 = (int) Math.floor(startY + remainder);
			g2d.fillOval(x0, y0, dsWidth, dsWidth);

			Color color = theme.getColor(dataset);

			for (int column = 0; column < data[dataset].length; column++) {

				int group = 0;
				if (colToGroup != null) {
					group = colToGroup.get(column);
				}

				double angle0 = angleOffSet - (degreesPerSegment * column) - (group * degreesPerGroup);
				double degreesPerSubCol = degreesPerSegment / data[dataset][column].length;
				System.out.println("colvsangle 0:\t" + column + "\t" + group + "\t" + angle0 + "\t" + degreesPerSegment + "\t" + degreesPerGroup);
				for (int subcol = 0; subcol < data[dataset][column].length; subcol++) {


					double angle1 = (angle0 - (degreesPerSubCol * subcol));
//					System.out.println("colvsangle 1:\t" + column + "\t" + angle1 + "\t" + degreesPerSubCol);


					double value = data[dataset][column][subcol];


					double perc = rangeValues.getRelativePositionY(value);

					double alpha = 255 * perc;
					if (alpha < 0) {
						alpha = 0;
					} else if (alpha > 255) {
						alpha = 255;
					}

					Color color2 = new Color(color.getRed(), color.getGreen(), color.getBlue(), (int) Math.floor(alpha));

					if (value >= 0) {
						g2d.setColor(color2);

						Arc2D arc = new Arc2D.Double(x0, y0, dsWidth, dsWidth, angle1, -degreesPerSubCol, Arc2D.PIE);
						g2d.fill(arc);

						if (data[dataset].length < 200) {

							g2d.setColor(new Color(220, 220, 220));
							g2d.draw(arc);
						}

					}
				}

			}
		}

//		// draw sector lines
//		g2d.setColor(theme.getDarkGrey());
//		for (int column = 0; column < data[0].length; column++) {
//			double angle0 = degreesPerSegment * column;
//
//			int originX = (int) Math.floor(startX + (maxWidth / 2));
//			int originY = (int) Math.floor(startY + (maxWidth / 2));
//			Pair<Integer, Integer> xy0 = Goniometry.calcPosOnCircle((maxWidth / 2), originX, originY, angle0);
//
//			g2d.drawLine(originX, originY, xy0.getLeft(), xy0.getRight());
//
//		}

		// fill inner circle with white

		int dsWidth = (int) Math.floor(maxWidth - (widthPerDataset * data.length));
		double remainder = (maxWidth - dsWidth) / 2;
		int x0 = (int) Math.floor(startX + remainder);
		int y0 = (int) Math.floor(startY + remainder);

		// first draw the outlines for the groups, if any
		if (groups != null) {
			g2d.setColor(theme.getDarkerColor(theme.getDarkGrey(), 0.5));
			g2d.setStroke(theme.getStroke());
			for (int group = 0; group < groups.size(); group++) {
				int col0 = groups.get(group).getLeft() - 1;
				int col1 = groups.get(group).getMiddle() - 1;
				int diff = col1 - col0;

				double angle0 = angleOffSet - (degreesPerSegment * col0) - (group * degreesPerGroup);

//				double degreesPerSubCol = degreesPerSegment / diff;
				double deg = degreesPerSegment; // degreesPerSubCol * diff;
				double angle1 = angle0 - deg;
//
//				double angle0 = degreesPerSegment * col0;
//				angle0 -= (degreesPerGroup / 2);
//				double deg = degreesPerSegment * diff;
//				angle0 -= (degreesForRowNames) + (degreesForRowNames / 2) - (degreesPerSegment / 2);
//				angle0 += (group * degreesPerGroup);
//
//				angle0 += 0.55;

				int dsWidth2 = (int) Math.floor(maxWidth - (widthPerDataset * 0));
				double remainder2 = (maxWidth - dsWidth2) / 2;
				int x02 = (int) Math.floor(startX + remainder2);
				int y02 = (int) Math.floor(startY + remainder2);
				Arc2D arc2 = new Arc2D.Double(x02, y02, maxWidth, maxWidth, angle1, -deg, Arc2D.PIE);
				g2d.draw(arc2);

				Arc2D arc = new Arc2D.Double(x0, y0, dsWidth, dsWidth, angle1, -deg, Arc2D.PIE);
//				g2d.setColor(new Color(220, 220, 220));
				g2d.draw(arc);

			}

		}

		g2d.setColor(Color.white);
		g2d.fillOval(x0 + 1, y0 + 1, dsWidth - 2, dsWidth - 2);


		// draw dataset names
		g2d.setFont(theme.getLargeFont());
		FontMetrics metrics = g2d.getFontMetrics(theme.getLargeFont());

		g2d.setColor(theme.getDarkGrey());

		for (int d = 0; d < rownames.length; d++) {

			int strlen = metrics.stringWidth(rownames[d]) / 2;

			int midpointX = (int) Math.floor(startX + (maxWidth / 2)) - strlen;

			int dsWidth1 = (int) Math.floor(maxWidth - (widthPerDataset * d));
			int dsWidth2 = (int) Math.floor(maxWidth - (widthPerDataset * (d + 1)));
			int remainder1 = (maxWidth - dsWidth1) / 2;
			int remainder2 = (maxWidth - dsWidth2) / 2;


			int midpointY1 = (int) Math.floor(startY + (remainder1));
			int midpointY2 = (int) Math.floor(startY + (remainder2));
			int actualMidPointY1 = (midpointY1 + midpointY2) / 2;
			g2d.drawString(rownames[d], midpointX, actualMidPointY1);

		}

		// draw column names
		for (int column = 0; column < data[0].length; column++) {

			double angle0 = degreesPerSegment * column;
			angle0 -= (degreesForRowNames) + (degreesForRowNames / 2) - (degreesPerSegment / 2);

			int group = 0;
			if (colToGroup != null) {
				group = colToGroup.get(column);
			}

			angle0 += (group * degreesPerGroup);

			// use the outer edge
			double originX = startX + (maxWidth / 2);
			double originY = startX + (maxWidth / 2);


//			g2d.drawString(colnames[column], xy0.getLeft(), xy0.getRight());


			// g2d.drawLine(originXi, originYi, xy0.getLeft(), xy0.getRight());


//			if (angle0 > 90) {
//
//
////				angle0 = degreesPerSegment * column;
////				angle0 -= (degreesForRowNames); // + (degreesForRowNames / 2) ;
////				// angle0 += (degreesPerSegment);
////				angle0 += (group * degreesPerGroup) - (degreesPerGroup / 2);
//				double radius = (double) (maxWidth + 30) / 2;
//				int strlen = metrics.stringWidth(colnames[column]);
//				radius += strlen;
//				Pair<Integer, Integer> xy0 = Goniometry.calcPosOnCircle(radius, originX, originY, angle0);
//				angle0 += 180;
//				drawRotate(g2d, xy0.getLeft(), xy0.getRight(), angle0, colnames[column]);
//			} else {
				double radius = (double) (maxWidth + 30) / 2;
				Pair<Integer, Integer> xy0 = Goniometry.calcPosOnCircle(radius, originX, originY, angle0);
				drawRotate(g2d, xy0.getLeft(), xy0.getRight(), angle0, colnames[column]);
//			}
			System.out.println(column + "\t" + colnames[column] + "\t" + angle0);


		}


	}

	public void setGroups(ArrayList<Triple<Integer, Integer, String>> groups) {
		this.groups = groups;
	}
}
