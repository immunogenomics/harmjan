package nl.harmjanwestra.miscscripts;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.CircularHeatmapPanel;

import java.io.IOException;

/**
 * Created by hwestra on 9/12/16.
 */
public class Circularheatmapplottest {

	public static void main(String[] args) {
		Grid grid = new Grid(1000, 1000, 1, 1, 100, 100);
		CircularHeatmapPanel panel = new CircularHeatmapPanel(1, 1);
		try {


			String[] columnnames = new String[16];
			double[][][] data = new double[3][16][];
			for (int d = 0; d < data.length; d++) {

				for (int i = 0; i < data[d].length; i++) {
					data[d][i] = new double[1];
					for (int q = 0; q < data[d][i].length; q++) {
						data[d][i][q] = 1;
					}
					columnnames[i] = "column-" + i;
				}
			}


			data[0][0] = new double[3];
			data[0][0][0] = 1;
			data[0][0][1] = 0.75;
			data[0][0][2] = 0.5;
			data[0][1][0] = 0.75;
			data[0][2][0] = 0.5;


			String[] rownames = new String[]{"disease1", "disease2", "disease3"};

			panel.setData(rownames, columnnames, data);

			grid.addPanel(panel);
			grid.draw("/Data/tmp/circplot.pdf");
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
	}

}
