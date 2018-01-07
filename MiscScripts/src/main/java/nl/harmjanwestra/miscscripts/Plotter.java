package nl.harmjanwestra.miscscripts;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.Panel;
import nl.harmjanwestra.utilities.graphics.themes.DefaultTheme;
import nl.harmjanwestra.utilities.graphics.themes.Theme;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.Descriptives;

import java.awt.*;
import java.awt.font.FontRenderContext;
import java.awt.font.GlyphVector;
import java.awt.font.LineMetrics;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.ArrayList;

public class Plotter {
	
	public static void main(String[] args) {
		Plotter p = new Plotter();
		
		try {
			p.run("D:\\Sync\\OneDrive\\Postdoc\\2017-01-SNPSeq\\v3\\data\\");
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String dir) throws IOException, DocumentException {
		File folder = new File(dir);
		File[] lsof = folder.listFiles(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				if (name.endsWith(".txt")) {
					return true;
				} else {
					return false;
				}
			}
		});
		System.out.println(lsof.length);
		Grid grid = new Grid(150, 100, lsof.length, 4, 100, 50);
		for (int f = 0; f < lsof.length; f++) {
			
			TextFile tf = new TextFile(lsof[f], TextFile.R);
			
			ArrayList<double[]> data = new ArrayList<double[]>();
			ArrayList<String> dsnames = new ArrayList<>();
			String[] header = tf.readLineElems(TextFile.tab);
			double max = Double.NaN;
			if (header[0].length() > 0) {
				max = Double.parseDouble(header[0]);
			}
			
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				double[] lndt = new double[elems.length - 1];
				for (int q = 1; q < elems.length; q++) {
					lndt[q - 1] = Double.parseDouble(elems[q]);
				}
				data.add(lndt);
				dsnames.add(elems[0]);
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			
			
			DotPlotPanel p = new DotPlotPanel(1, 1);
			p.setTitle(lsof[f].getName());
			p.setData(data);
			p.setMax(max);
			p.setDsnames(dsnames);
			grid.addPanel(p);
			
			
		}
		grid.draw(dir + "plots.png", DefaultGraphics.Output.PNG);
		
	}
	
	
}

class DotPlotPanel extends Panel {
	
	private ArrayList<double[]> data;
	private double max;
	private ArrayList<String> dsnames;
	
	public DotPlotPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}
	
	@Override
	public void draw(DefaultGraphics g) {
		
		
		// draw some x axis
		Graphics2D g2d = g.getG2d();
		Theme theme = new DefaultTheme();
		setTheme(theme);
		
		g2d.setColor(theme.getDarkGrey());
		g2d.setStroke(theme.getStroke());
		
		
		int innermargin = 10;
		
		g2d.drawLine(x0, y0 + innermargin, x0, y0 + height - innermargin); // y axis
		
		g2d.drawLine(x0 - 3, y0 + innermargin, x0, y0 + innermargin); // y axis
		g2d.drawLine(x0 - 3, y0 + height - innermargin, x0, y0 + height - innermargin); // y axis
		
		
		g2d.drawLine(x0 + innermargin, y0 + height, x0 + width - innermargin, y0 + height); // x- axis
//		g2d.drawLine(x0 + innermargin, y0 + height + 3, x0 + innermargin, y0 + height); // x- axis
//		g2d.drawLine(x0 + width - innermargin, y0 + height + 3, x0 + width - innermargin, y0 + height); // x- axis
		
		
		// figure out number of datasets to plot
		int spacingbetweendatasets = 10;
		int widthperdataset = 10;
		
		
		// figure out a scale
		if (Double.isNaN(max)) {
			max = -1;
			for (double[] ds : data) {
				// get stdev and mean
				
				for (double d : ds) {
					if (d > max) {
						max = d;
					}
				}
			}
		}
		
		
		int innerheight = height - (2 * innermargin);
		int innerwidth = width - (2 * innermargin);
		int widthperds = innerwidth / data.size();
		
		// round up
		double log = Math.ceil(-Math.log10(max)); // get the nearest full log10
		double dint = max * Math.pow(10, log); // convert to integer space
		double roundup = Math.ceil(dint); // ceil
		max = roundup * Math.pow(10, -log); // convert back down
		
		
		// draw file name
		g2d.drawString(title, x0, y0);
		
		// draw min and max
		DecimalFormat format = new DecimalFormat("#.####");
		String maxStr = format.format(max);
		int strwidthmax = super.getStringWidth(g2d, maxStr);
		FontMetrics fm = g2d.getFontMetrics(g2d.getFont());
		int halffontheight = (fm.getHeight()) / 4;
		g2d.drawString(maxStr, x0 - 3 - 3 - strwidthmax, y0 + innermargin + halffontheight);
		strwidthmax = super.getStringWidth(g2d, "0");
		g2d.drawString("0", x0 - 3 - 3 - strwidthmax, y0 + height - innermargin + halffontheight);
		
		
		int dsctr = 0;
		
		for (double[] ds : data) {
			// get stdev and mean
			double mean = Descriptives.mean(ds);
			double sd = Math.sqrt(Descriptives.variance(ds));
			
			int halfwidth = widthperds / 2;
			int xoffset = 0;
			
			int xpos = x0 + innermargin + (widthperds * dsctr) + halfwidth;
			g2d.setColor(new Color(0, 0, 0, 64));
			
			
			
			for (double d : ds) {
				// conver to pixel space
				double perc = d / max;
				int nrPixels = (int) Math.floor(perc * innerheight);
				int ypos = y0 + innermargin + innerheight - nrPixels;
				g2d.fillOval(xpos - 4 - xoffset, ypos - 4, 8, 8);
			}
			
			
			// draw mean and sd
			g2d.setColor(new Color(255, 0, 0));
			double perc = mean / max;
			int nrPixels = (int) Math.floor(perc * innerheight);
			double percsd = sd / max;
			int nrPixelssdev = (int) Math.floor(percsd * innerheight);
			
			
			int ypos = y0 + innermargin + innerheight - nrPixels;
			g2d.setStroke(new BasicStroke(2.5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
			
			g2d.drawLine(xpos - 5 + xoffset, ypos, xpos + 5 + xoffset, ypos); // mean
			g2d.drawLine(xpos + xoffset, ypos - nrPixelssdev, xpos + xoffset, ypos + nrPixelssdev); // sd
			
			
			g2d.setColor(new Color(0, 0, 0, 255));
			g2d.setStroke(theme.getStroke());
			g2d.drawLine(xpos, y0 + height, xpos, y0 + height + 3); // line at the bottom of the x-axis
			
			if (dsnames != null) {
				
				int strwidth = super.getStringWidth(g2d, dsnames.get(dsctr));
				int halfwidthstr = strwidth / 2;
				g2d.drawString(dsnames.get(dsctr), xpos - halfwidthstr, y0 + height + 17);
				
				
			}
			
			// draw mean and s.d.
			dsctr++;
		}
		
		
	}
	
	public void setData(ArrayList<double[]> data) {
		this.data = data;
	}
	
	public void setMax(double max) {
		this.max = max;
	}
	
	public void setDsnames(ArrayList<String> dsnames) {
		this.dsnames = dsnames;
	}
}