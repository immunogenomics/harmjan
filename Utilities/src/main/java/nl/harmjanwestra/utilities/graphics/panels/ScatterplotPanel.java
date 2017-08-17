package nl.harmjanwestra.utilities.graphics.panels;

import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.Correlation;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.Regression;
import nl.harmjanwestra.utilities.legacy.genetica.util.Primitives;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.List;

/**
 * Created by hwestra on 10/15/15.
 */
public class ScatterplotPanel extends Panel {
	
	double[][] x; // format [dataset1][point1] [dataset1][point2] [dataset1][pointetc]
	double[][] y;
	
	private Range dataRange;
	private String xAxisLabel;
	private String yAxisLabel;
	private String[] datasetLabels;
	private boolean plotLinearRegression = false;
	
	public ScatterplotPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}
	
	public void setData(double[][] x, double[][] y) {
		this.x = x;
		this.y = y;
	}
	
	public void setData(double[] x, double[] y, int[] rowGroups) {
		// determine number of unique row groups
		HashSet<Integer> uniquegroups = new HashSet<Integer>();
		
		HashMap<Integer, Integer> groupValueCounter = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> groupIndex = new HashMap<Integer, Integer>();
		
		int ctr = 0;
		for (int i = 0; i < rowGroups.length; i++) {
			int numValuesInGroup = 0;
			if (groupValueCounter.containsKey(rowGroups[i])) {
				numValuesInGroup = groupValueCounter.get(rowGroups[i]);
			} else {
				groupIndex.put(rowGroups[i], ctr);
				ctr++;
			}
			numValuesInGroup++;
			groupValueCounter.put(rowGroups[i], numValuesInGroup);
		}
		
		double[][] x1 = new double[groupValueCounter.size()][];
		double[][] y1 = new double[groupValueCounter.size()][];
		
		for (int i = 0; i < x.length; i++) {
		
		}
	}
	
	public void setData(double[] x, double[] y) {
		this.x = new double[1][];
		this.x[0] = x;
		this.y = new double[1][];
		this.y[0] = y;
	}
	
	public void setLabels(String xAxis, String yAxis) {
		this.xAxisLabel = xAxis;
		this.yAxisLabel = yAxis;
	}
	
	public void setDataRange(Range dataRange) {
		this.dataRange = dataRange;
	}
	
	@Override
	public void draw(DefaultGraphics g) {
		
		
		Graphics2D g2d = g.getG2d();
		
		// determine range
		if (dataRange == null) {
			// determine max and min using all data
			
			dataRange = new Range(Double.MAX_VALUE,
					Double.MAX_VALUE,
					-Double.MAX_VALUE,
					-Double.MAX_VALUE);
			
			for (int i = 0; i < x.length; i++) {
				Range r = new Range(x[i], y[i]);
				r.round();
				if (r.getMaxX() > dataRange.getMaxX()) {
					dataRange = new Range(dataRange.getMinX(),
							dataRange.getMinY(),
							r.getMaxX(),
							dataRange.getMaxY());
				}
				if (r.getMaxY() > dataRange.getMaxY()) {
					dataRange = new Range(dataRange.getMinX(),
							dataRange.getMinY(),
							dataRange.getMaxX(),
							r.getMaxY());
				}
				if (r.getMinX() < dataRange.getMinX()) {
					dataRange = new Range(r.getMinX(),
							dataRange.getMinY(),
							dataRange.getMaxX(),
							dataRange.getMaxY());
				}
				if (r.getMinY() < dataRange.getMinY()) {
					dataRange = new Range(dataRange.getMinX(),
							r.getMinY(),
							dataRange.getMaxX(),
							dataRange.getMaxY());
				}
			}
			
			dataRange.round();
			
		}
		
		Range plotRange = dataRange;
		// plot the points
		
		int nrPixelsMaxX = width - (2 * marginX);
		int nrPixelsMaxY = height - (2 * marginY);
		
		// axis labels
		// plot y-axis
		g2d.setColor(theme.getDarkGrey());
		double tickUnitY = plotRange.getRangeY() / 10;
		String pattern = "###,###,###.##";
		DecimalFormat decimalFormat = new DecimalFormat(pattern);
		
		g2d.setFont(theme.getMediumFont());
		FontMetrics metrics = g2d.getFontMetrics(g2d.getFont());
		
		int xPosYAxis = x0 + marginX - 10;
		int yPosYAxis = y0 + marginY;
		g2d.drawLine(xPosYAxis, yPosYAxis, xPosYAxis, yPosYAxis + nrPixelsMaxY);
		
		int maxlen = 0;
		for (double y = plotRange.getMinY(); y < plotRange.getMaxY() + (tickUnitY / 2); y += tickUnitY) {
			double yPerc = plotRange.getRelativePositionY(y);
			
			int ypos = y0 + marginY + (int) Math.ceil((1 - yPerc) * nrPixelsMaxY);
			int startx = xPosYAxis - 5;
			int stopx = startx + 5;
			g2d.drawLine(startx, ypos, stopx, ypos);
			String formattedStr = decimalFormat.format(y);
			int adv = metrics.stringWidth(formattedStr);
			if (adv > maxlen) {
				maxlen = adv;
			}
			g2d.setFont(theme.getMediumFont());
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
			int starty = yPosXAxis;
			int stopy = starty + 5;
			String formattedStr = decimalFormat.format(x);
			g2d.drawLine(xpos, starty, xpos, stopy);
			int adv = metrics.stringWidth(formattedStr);
			g2d.setFont(theme.getMediumFont());
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
		
		// plot the points
		for (int i = 0; i < x.length; i++) {
			g2d.setColor(theme.getColor(i));
			for (int j = 0; j < x[i].length; j++) {
				double xval = x[i][j];
				double yval = y[i][j];
				
				if (!Double.isNaN(xval) && !Double.isNaN(yval)) {
					double xperc = dataRange.getRelativePositionX(xval);
					double yperc = dataRange.getRelativePositionY(yval);
					int pixelX = x0 + marginX + (int) Math.ceil(nrPixelsMaxX * xperc);
					int pixelY = y0 + marginY + (int) Math.ceil(nrPixelsMaxY - (nrPixelsMaxY * yperc));
					
					g2d.fillOval(pixelX - 2, pixelY - 2, 4, 4);
					
				} else {
					// plot some markers somewhere?
				}
			}
			
			
		}
		
		if (plotLinearRegression) {
			for (int i = 0; i < x.length; i++) {
				g2d.setColor(theme.getColor(i));
				double[] xval = x[i];
				double[] yval = y[i];
				
				Pair<double[], double[]> filtered = removeNulls(xval, yval);
				xval = filtered.getLeft();
				yval = filtered.getRight();
				
				// y = a + bx
				double[] params = Regression.getLinearRegressionCoefficients(xval, yval);
				double beta = params[0];
				double alpha = params[1];
				
				double minX = dataRange.getMinX();
				double maxX = dataRange.getMaxX();
				
				double y1 = alpha + (minX * beta);
				double y2 = alpha + (maxX * beta);
				
				
				double x1perc = dataRange.getRelativePositionX(minX);
				double y1perc = dataRange.getRelativePositionY(y1);
				double x2perc = dataRange.getRelativePositionX(maxX);
				double y2perc = dataRange.getRelativePositionY(y2);
				
				
				int pixelX1 = x0 + marginX + (int) Math.ceil(nrPixelsMaxX * x1perc);
				int pixelY1 = y0 + marginY + (int) Math.ceil(nrPixelsMaxY - (nrPixelsMaxY * y1perc));
				int pixelX2 = x0 + marginX + (int) Math.ceil(nrPixelsMaxX * x2perc);
				int pixelY2 = y0 + marginY + (int) Math.ceil(nrPixelsMaxY - (nrPixelsMaxY * y2perc));
				
				g2d.drawLine(pixelX1, pixelY1, pixelX2, pixelY2);
				
				// for shits and giggles, calculate the R-squared using the correlation coefficient
				double corr = Correlation.correlate(xval, yval);
				double rsq = corr * corr;
				
				DecimalFormat format = new DecimalFormat("#.##");
				
				g2d.drawString("" + format.format(rsq), pixelX2, pixelY2);
				
			}
		}
		
		if (datasetLabels != null) {
			// plot a legend
			
			int maxStrWidth = 0;
			for (int i = 0; i < datasetLabels.length; i++) {
				String label = datasetLabels[i];
				int strwdth = getStringWidth(g2d, label);
				if (strwdth > maxStrWidth) {
					maxStrWidth = strwdth;
				}
			}
			
			for (int i = 0; i < datasetLabels.length; i++) {
				g2d.setColor(theme.getColor(i));
				int pixelX = x0 + marginX + nrPixelsMaxX - maxStrWidth - 15;
				int pixelY = y0 + marginY;
				
				g2d.fillRect(pixelX, pixelY + (i * 15), 10, 10);
				g2d.drawString(datasetLabels[i], pixelX + 15, pixelY + (i * 15) + 10);
				
			}
			
		}
		
		
	}
	
	private Pair<double[], double[]> removeNulls(double[] xval, double[] yval) {
		boolean[] isnull = new boolean[xval.length];
		int ctrnotnull = 0;
		for (int i = 0; i < xval.length; i++) {
			if (Double.isNaN(xval[i]) || Double.isNaN(yval[i])) {
				isnull[i] = true;
			} else {
				ctrnotnull++;
			}
		}
		
		double[] outx = new double[ctrnotnull];
		double[] outy = new double[ctrnotnull];
		ctrnotnull = 0;
		for (int i = 0; i < xval.length; i++) {
			if (Double.isNaN(xval[i]) || Double.isNaN(yval[i])) {
			
			} else {
				outx[ctrnotnull] = xval[i];
				outy[ctrnotnull] = yval[i];
				ctrnotnull++;
			}
		}
		return new Pair<double[], double[]>(outx, outy);
	}
	
	
	public void setDatasetLabels(String[] datasetLabels) {
		this.datasetLabels = datasetLabels;
	}
}
