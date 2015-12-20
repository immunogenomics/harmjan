package nl.harmjanwestra.utilities.graphics.panels;

import nl.harmjanwestra.utilities.features.BedGraphFeature;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;

import java.awt.*;
import java.io.File;
import java.util.ArrayList;

/**
 * Created by hwestra on 12/3/15.
 */
public class GraphAnnotationPanel extends Panel {
	public GraphAnnotationPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	Feature region;
	ArrayList<ArrayList<BedGraphFeature>> data;
	String[] filenames;

	double maxHeight = 100;
	int trackheight = 50;

	public void setData(Feature region, ArrayList<ArrayList<BedGraphFeature>> data, String[] bfiles) {
		this.region = region;
		this.data = data;
		this.filenames = bfiles;

	}


	@Override
	public void draw(DefaultGraphics g) {

		Graphics2D g2d = g.getG2d();

		int figureWidth = width;
		int regionSize = region.getStop() - region.getStart();
		int nrPixelsX = figureWidth - (2 * marginX);


		int marginBetween = 5;


		Color defaultLightGrey = new Color(175, 175, 175);
		Color defaultColor = new Color(90, 90, 90);
		Color highlight = new Color(208, 83, 77);
		g2d.setColor(defaultColor);


		for (int i = 0; i < data.size(); i++) {


			int trackYpos = marginY + y0 + (i * trackheight) + (i * marginBetween);

			int startX = x0 + marginX;
			g2d.setColor(defaultLightGrey);
			g2d.drawLine(startX, trackYpos + trackheight, startX + width, trackYpos + trackheight);
			g2d.setColor(defaultColor);


			// get subset within region
			ArrayList<BedGraphFeature> t = data.get(i);

			for (int q = 0; q < t.size(); q++) {
				BedGraphFeature f = t.get(q);
				int start = f.getStart();
				int stop = f.getStop();
				if (start > region.getStart() && stop < region.getStop()) {
					int featurewidth = f.getStop() - start;
					int relativeStart = start - region.getStart();
					if (relativeStart < 0) {
						featurewidth -= Math.abs(relativeStart);
						relativeStart = 0;
					}

					int relativeStop = relativeStart + featurewidth;
					if (relativeStop > regionSize) {
						relativeStop = regionSize;
					}

					double percStart = (double) relativeStart / regionSize;
					double percStop = (double) relativeStop / regionSize;

					int pixelStart = x0 + marginX + (int) Math.ceil(percStart * nrPixelsX);
					int pixelStop = x0 + marginX + (int) Math.ceil(percStop * nrPixelsX);

					int y1 = trackYpos;

					int boxwidth = pixelStop - pixelStart;
					if (boxwidth <= 0) {
						boxwidth = 1;
					}

					double v = f.getValue();
					int boxHeight = 0;
					if (v > 5) {
						System.out.println("gotone");
					}
					if (v > maxHeight) {
						v = maxHeight;
					}


					boxHeight = (int) ((Math.ceil((v / maxHeight) * trackheight)));

					g2d.fillRect(pixelStart, y1 + trackheight - boxHeight, boxwidth, boxHeight);
				}

			}

			g2d.setColor(defaultColor);

			// plot the name of the annotation
			String name = filenames[i];
			File file = new File(name);
			String filename = file.getName();
			g2d.setFont(new Font("default", Font.BOLD, 8));

			int pixelStart = x0 + marginX + 5 + nrPixelsX;
			g2d.drawString(filename, pixelStart, trackYpos + trackheight);


		}

	}


}
