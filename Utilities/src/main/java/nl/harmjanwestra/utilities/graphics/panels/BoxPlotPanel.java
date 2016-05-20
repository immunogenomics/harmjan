package nl.harmjanwestra.utilities.graphics.panels;

import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.themes.DefaultTheme;
import nl.harmjanwestra.utilities.graphics.themes.Theme;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.util.Primitives;

import java.awt.*;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by hwestra on 9/12/15.
 */
public class BoxPlotPanel extends Panel {

	private double[][][] data;
	private Range range;
	private boolean drawDataPoints;
	private boolean useMeanAndSd;
	private String outputIQRS;
	private String[] labels;
	private boolean useTukeysDefault;

	public BoxPlotPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	@Override
	public void draw(DefaultGraphics g) {

		int nrDatasets = data.length;
		int nrBins = data[0].length;

		double[][][] iqrs = new double[data.length][data[0].length][0];
		double max = -Double.MAX_VALUE;
		double min = Double.MAX_VALUE;
		for (int ds = 0; ds < nrDatasets; ds++) {
			for (int bin = 0; bin < nrBins; bin++) {
				iqrs[ds][bin] = collectStats(data[ds][bin]);
				double localmin = iqrs[ds][bin][0];
				double localmax = iqrs[ds][bin][4];
//			System.out.println(localmin + "\t" + localmax);
				if (localmin < min) {
					min = localmin;
				}
				if (localmax > max) {
					max = localmax;
				}
			}
		}


		if (range == null) {
			range = new Range(0, min, 0, max);
		}


		int pixelsY = height - (2 * marginY);
		int pixelsX = width - (2 - marginX);


		// plot the data
		int marginBetweenBins = 20;
		int marginBetweenDatasets = 2;

		int totalNrOfBins = nrDatasets * nrBins;
		int spaceRequiredByBetweenDatasetMargins = (totalNrOfBins - 1) * marginBetweenDatasets;
		int spaceRequiredByBetweenBinMargins = (nrBins - 1) * marginBetweenBins;
		int totalSpaceLeftAfterRemovalOfBins = pixelsX - spaceRequiredByBetweenBinMargins - spaceRequiredByBetweenDatasetMargins;
		int widthPerDatasetBin = totalSpaceLeftAfterRemovalOfBins / totalNrOfBins;

		int starty = y0 + marginY;

		Theme theme = new DefaultTheme();

		Color c = theme.getDarkGrey();
		Graphics2D g2d = g.getG2d();
		int halfWidthOfDataset = widthPerDatasetBin / 2;

		try {
			TextFile iqrout = null;
			if (outputIQRS != null) {
				iqrout = new TextFile(outputIQRS, TextFile.W);
			}


			for (int bin = 0; bin < nrBins; bin++) {
				for (int ds = 0; ds < nrDatasets; ds++) {
					int startX = x0 + marginX + (bin * marginBetweenBins) + (bin * nrDatasets * widthPerDatasetBin) + ((nrDatasets - 1) * marginBetweenDatasets);
					startX += (ds * widthPerDatasetBin) + (ds * marginBetweenDatasets);

					// retrieve IQRS
					double datamin = iqrs[ds][bin][0];
					double dataq1 = iqrs[ds][bin][1];
					double datamedian = iqrs[ds][bin][2];
					double dataq3 = iqrs[ds][bin][3];
					double datamax = iqrs[ds][bin][4];
					double mean = iqrs[ds][bin][5];
					double sd = iqrs[ds][bin][6];

					double[] datasetDataForBin = data[ds][bin];

					if (datasetDataForBin.length > 2) {
						if (drawDataPoints) {

							Pair<double[][], double[]> densityData = determineDensity(datasetDataForBin, range);

							double[][] values = densityData.getLeft();
							double[] frequencies = densityData.getRight();

							for (int j = 0; j < frequencies.length; j++) {
								double[] dataInFreqBin = values[j];
								double maxPercJitter = frequencies[j];
								for (int v = 0; v < dataInFreqBin.length; v++) {
									double d = dataInFreqBin[v];

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

									int jittersize = (int) Math.ceil(Math.random() * Math.random() * (maxPercJitter * halfWidthOfDataset));
									int jitter = jittersize * direction;
									int plotX = startX + halfWidthOfDataset + jitter;
									g2d.setColor(c);
									if (d > dataq1 && d < dataq3) {
										// falls within IQR
										// set opacity to 30%;
										Color c2 = new Color(70, 67, 58, 128);
										g2d.setColor(c2);
									} else {
										g2d.setColor(c);
									}
									g2d.fillOval(plotX - 2, plotY - 2, 4, 4);
								}
							}
						} else {
							// just draw the box plot
							g2d.setColor(theme.getLightGrey());
							g2d.setStroke(theme.getStroke());
							// draw line from min to q1


							if (useMeanAndSd) {

							} else if (useTukeysDefault) {

								double iqr = (dataq3 - dataq1);
								double iqrmax = dataq3 + (1.5 * iqr);
								double iqrmin = dataq1 - (1.5 * iqr);

								// plot the outliers
								g2d.setColor(theme.getLightGrey());
								for (int i = 0; i < datasetDataForBin.length; i++) {
									double d = datasetDataForBin[i];
									if (d > iqrmax || d < iqrmin) {
										double pos = range.getRelativePositionY(d);
										int outlierY = (int) Math.ceil(pos * pixelsY);
										int outlierYPx = starty + pixelsY - outlierY;
										g2d.fillOval(startX + halfWidthOfDataset - 2, outlierYPx - 2, 4, 4);
									}
								}

								g2d.setColor(theme.getColor(ds));

								boolean clippingbottom = false;
								boolean clippingtop = false;

								if (iqrmin < range.getMinY()) {
									iqrmin = range.getMinY();
									clippingbottom = true;
								}

								if (iqrmax > range.getMaxY()) {
									iqrmax = range.getMaxY();
									clippingtop = true;
								}


								double iqrminy = range.getRelativePositionY(iqrmin);
								double q1y = range.getRelativePositionY(dataq1);
								double m2y = range.getRelativePositionY(datamedian);
								double q3y = range.getRelativePositionY(dataq3);
								double iqrmaxy = range.getRelativePositionY(iqrmax);

								int iqrMinYPx = (int) Math.ceil(iqrminy * pixelsY);
								int q1yPx = (int) Math.ceil(q1y * pixelsY);
								int m2yPx = (int) Math.ceil(m2y * pixelsY);
								int q3yPx = (int) Math.ceil(q3y * pixelsY);
								int iqrMaxYPx = (int) Math.ceil(iqrmaxy * pixelsY);

								int plotYiqrMinY = starty + pixelsY - iqrMinYPx;
								int plotYq1yPx = starty + pixelsY - q1yPx;
								int plotYm2yPx = starty + pixelsY - m2yPx;
								int plotYq3yPx = starty + pixelsY - q3yPx;
								int plotYiqrMaxY = starty + pixelsY - iqrMaxYPx;

//								if (clippingbottom) {
//									g2d.setStroke(theme.getStrokeDashed());
//									g2d.drawLine(startX + halfWidthOfDataset - halfWidthOfDataset, starty + pixelsY, startX + halfWidthOfDataset + halfWidthOfDataset, starty + pixelsY);
//								}
//								g2d.setStroke(theme.getStroke());
//
//								// draw horizontal leg of the whisker
//								g2d.drawLine(startX + halfWidthOfDataset, plotYiqrMinY, startX + halfWidthOfDataset, plotYiqrMinY);
//
//								if (clippingtop) {
//									g2d.setStroke(theme.getStrokeDashed());
//									g2d.drawLine(startX + halfWidthOfDataset - halfWidthOfDataset, starty, startX + halfWidthOfDataset + halfWidthOfDataset, starty);
//								}

								// draw vertical parts of leg
								g2d.setStroke(theme.getStroke());
								g2d.drawLine(startX + halfWidthOfDataset, plotYiqrMinY, startX + halfWidthOfDataset, plotYq1yPx);
								g2d.drawLine(startX + halfWidthOfDataset, plotYq3yPx, startX + halfWidthOfDataset, plotYiqrMaxY);

								// draw the median
								g2d.drawLine(startX, plotYm2yPx, startX + widthPerDatasetBin, plotYm2yPx);
								g2d.fillOval(startX + halfWidthOfDataset - 3, plotYm2yPx - 3, 6, 6);


								// draw a lightgrey box

								g2d.drawRect(startX, plotYq3yPx, widthPerDatasetBin, plotYq1yPx - plotYq3yPx);

								// System.out.println(startX + "\t" + plotYq1y + "\t" + widthPerDataset + "\t" + (plotYq1y - plotYq3y));

							} else {

								boolean clippingbottom = false;
								boolean clippingtop = false;

								if (datamin < range.getMinY()) {
									datamin = range.getMinY();
									clippingbottom = true;
								}

								if (datamax > range.getMaxY()) {
									datamax = range.getMaxY();
									clippingtop = true;
								}

								double m1y = range.getRelativePositionY(datamin);
								double q1y = range.getRelativePositionY(dataq1);
								double m2y = range.getRelativePositionY(datamax);
								double q3y = range.getRelativePositionY(dataq3);

								int m1yPx = (int) Math.ceil(m1y * pixelsY);
								int q1yPx = (int) Math.ceil(q1y * pixelsY);


								int m2yPx = (int) Math.ceil(m2y * pixelsY);
								int q3yPx = (int) Math.ceil(q3y * pixelsY);

								int plotYm1y = starty + pixelsY - m1yPx;
								int plotYq1y = starty + pixelsY - q1yPx;


								if (clippingbottom) {
									g2d.setStroke(theme.getStrokeDashed());
									g2d.drawLine(startX + halfWidthOfDataset - halfWidthOfDataset, starty + pixelsY, startX + halfWidthOfDataset + halfWidthOfDataset, starty + pixelsY);
								}
								g2d.setStroke(theme.getStroke());

								g2d.drawLine(startX + halfWidthOfDataset, plotYm1y, startX + halfWidthOfDataset, plotYq1y);


								int plotYm2y = starty + pixelsY - m2yPx;
								int plotYq3y = starty + pixelsY - q3yPx;

								if (clippingtop) {
									g2d.setStroke(theme.getStrokeDashed());
									g2d.drawLine(startX + halfWidthOfDataset - halfWidthOfDataset, starty, startX + halfWidthOfDataset + halfWidthOfDataset, starty);
								}
								g2d.setStroke(theme.getStroke());
								g2d.drawLine(startX + halfWidthOfDataset, plotYm2y, startX + halfWidthOfDataset, plotYq3y);


								// draw the median
								double q2y = range.getRelativePositionY(datamedian);
								int q2yPx = (int) Math.ceil(q2y * pixelsY);
								g2d.setColor(theme.getColor(1));
								int plotYq2y = starty + pixelsY - q2yPx;
								g2d.fillOval(startX + halfWidthOfDataset - 3, plotYq2y - 3, 6, 6);
								g2d.setColor(theme.getColor(0));

								// draw a lightgrey box
								g2d.setColor(theme.getLightGrey());

								g2d.drawRect(startX, plotYq1y, widthPerDatasetBin, plotYq3y - plotYq1y);
							}


//						//g2d.drawRect(startX, plotYq3y, widthPerBox, plotYq1y - plotYq3y);
//						if (outputIQRS != null) {
//							iqrout.writeln(datamin + "\t" + dataq1 + "\t" + datamedian + "\t" + dataq3 + "\t" + datamax
//									+ "\t" + plotYm1y + "\t" + plotYq1y + "\t" + plotYq2y + "\t" + plotYq3y + "\t" + plotYm2y);
//						}
						}
					}


					// plot an axis
				}
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
		for (int bin = 0; bin < nrBins; bin++) {
			int startX = x0 + marginX + (bin * marginBetweenBins) + (bin * nrDatasets * widthPerDatasetBin) + ((nrDatasets - 1) * marginBetweenDatasets);

			int totalBinWidth = (nrDatasets * widthPerDatasetBin) + ((nrDatasets - 1) * marginBetweenDatasets);

			startX += totalBinWidth / 2;
			g2d.drawLine(startX, yPosXAxis, startX, yPosXAxis + 5);
		}

		if (labels != null) {
			g2d.setColor(theme.getDarkGrey());
			g2d.setFont(theme.getSmallFont());
			metrics = g2d.getFontMetrics(g2d.getFont());
			int fontheight = metrics.getHeight();

			int y = yPosXAxis + 5;
			for (int bin = 0; bin < labels.length; bin++) {
				int startX = x0 + marginX + (bin * marginBetweenBins) + (bin * nrDatasets * widthPerDatasetBin) + ((nrDatasets - 1) * marginBetweenDatasets);
				int totalBinWidth = (nrDatasets * widthPerDatasetBin) + ((nrDatasets - 1) * marginBetweenDatasets);
				startX += totalBinWidth / 2;

				String str = labels[bin];
				int widthOfStr = metrics.stringWidth(str);
				drawRotate(g2d, startX + (fontheight / 2), y + widthOfStr + 10, -90, str);
			}
		}


	}

	private Pair<double[][], double[]> determineDensity(double[] doubles, Range yrange) {
		ArrayList<ArrayList<Double>> binobj = new ArrayList<ArrayList<Double>>();
		int nrBins = 100;
		for (int i = 0; i < nrBins; i++) {
			binobj.add(new ArrayList<Double>());
		}
		for (int i = 0; i < doubles.length; i++) {
			double d = doubles[i];
			double perc = yrange.getRelativePositionY(d);
			double bin = perc * nrBins;
			int ibin = (int) Math.ceil(bin);
			if (ibin < 0) {
				ibin = 0;
			} else if (ibin >= nrBins) {
				ibin = nrBins - 1;
			}
			binobj.get(ibin).add(d);
		}

		// make 100 bins
		double[] freqs = new double[nrBins];
		double[][] binnedData = new double[nrBins][];

		// determine number of values in each bin
		for (int i = 0; i < nrBins; i++) {
			binnedData[i] = Primitives.toPrimitiveArr(binobj.get(i).toArray(new Double[0]));
			double percOfTotalData = (double) binobj.get(i).size() / doubles.length;
			freqs[i] = percOfTotalData;
		}

		return new Pair<double[][], double[]>(binnedData, freqs);

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

	public void setData(double[][] data) {
		this.data = new double[1][data.length][];
		for (int i = 0; i < data.length; i++) {
			this.data[0][i] = data[i];
		}
	}

	private double[] collectStats(double[] dataset) {

		double min = Double.MAX_VALUE;
		double max = -Double.MAX_VALUE;
		double firstquartile = 0;
		double thirdquartile = 0;
		double median = 0;

		double[] datacopy = new double[dataset.length];
		System.arraycopy(dataset, 0, datacopy, 0, dataset.length);

		double mean = Descriptives.mean(datacopy);
		double sd = Math.sqrt(Descriptives.variance(datacopy));

		boolean even = false;
		if (dataset.length % 2 == 0) {
			even = true;
		}

		if (datacopy.length < 2) {
			if (datacopy.length == 1) {
				return new double[]{
						datacopy[0], datacopy[0], datacopy[0], datacopy[0], datacopy[0], datacopy[0], 0
				};
			} else {
				return new double[]{
						0, 0, 0, 0, 0, 0, 0
				};
			}

		} else {
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
					min, firstquartile, median, thirdquartile, max, mean, sd
			};
		}


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

	public void setData(double[][][] bins) {
		data = bins;
	}

	public void useTukeysDefault(boolean b) {
		useTukeysDefault = b;
	}
}
