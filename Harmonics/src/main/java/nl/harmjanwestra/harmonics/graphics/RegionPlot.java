package nl.harmjanwestra.harmonics.graphics;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;

import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Created by hwestra on 6/25/15.
 */
public class RegionPlot extends DefaultGraphics {


	int margin;
	private int[][] data;
	private Feature region;

	private int betweenPlotMargin = 50;
	private int geneHeight;

	public RegionPlot(String filename, int width, int height) throws FileNotFoundException, DocumentException {
		super(filename, width, height);

	}

	public void setMargin(int margin) {
		this.margin = margin;
	}

	public void setGeneHeight(int geneheight) {
		this.geneHeight = geneheight;
	}

	public void setData(int[][] data, Feature region, String genefile) {
		this.data = data;
		this.region = region;
	}

	public void draw() throws IOException {

		int nrPlots = data.length;

		int heightminmargin = figureHeight - betweenPlotMargin - geneHeight - (margin * 2);
		int heightPerPlot = (heightminmargin / data.length) - (data.length * betweenPlotMargin);
		int widthPerPlot = figureWidth - (margin * 2);




		close();
	}


}
