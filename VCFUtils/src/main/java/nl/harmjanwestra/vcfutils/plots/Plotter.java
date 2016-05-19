package nl.harmjanwestra.vcfutils.plots;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.BoxPlotPanel;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import umcg.genetica.io.text.TextFile;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 5/18/16.
 */
public class Plotter {


	public void plotCorr() {
		String[] files = new String[]{};
		
		String out = "";
		int width = 1000;
		int height = 1000;

		try {
			// plot 1: x-axis nr of variants, y-axis correlation,
			ArrayList<ArrayList<Double>> vals = new ArrayList<ArrayList<Double>>();
			int maxSize = 0;
			for (int i = 0; i < files.length; i++) {
				String file = files[i];
				ArrayList<Double> corvals = new ArrayList<>();
				TextFile tf = new TextFile(file, TextFile.R);
				String[] elems = tf.readLineElems(TextFile.tab);
				while (elems != null) {
					double val = Double.parseDouble(elems[1]);
					corvals.add(val);
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();
				if (corvals.size() > maxSize) {
					maxSize = corvals.size();
				}
				vals.add(corvals);
			}

			double[][] x = new double[files.length][maxSize];
			double[][] y = new double[files.length][maxSize];
			for (int ds = 0; ds < files.length; ds++) {
				for (int i = 0; i < maxSize; i++) {
					x[ds][i] = i;
				}
				ArrayList<Double> corvals = vals.get(ds);
				for (int i = 0; i < corvals.size(); i++) {
					y[ds][i] = corvals.get(i);
				}
				for (int i = corvals.size(); i < maxSize; i++) {
					y[ds][i] = Double.NaN;
				}
			}

			vals = null;
			Grid grid = new Grid(1000, 1000, 1, 1, 0, 0);
			ScatterplotPanel panel = new ScatterplotPanel(1, 1);
			panel.setData(x, y);
			grid.addPanel(panel);
			grid.draw(out);

			// plot 2: x-axis maf, y-axis correlation (boxplot)
			ArrayList<ArrayList<ArrayList<Double>>> bins = new ArrayList<ArrayList<ArrayList<Double>>>();
			double[] lowerthreshold = new double[]{
					0,
					0.005,
					0.01,
					0.02,
					0.03,
					0.04,
					0.05,
					0.1,
					0.2,
					0.3,
					0.4
			};
			double[] upperthreshold = new double[]{
					0.005,
					0.01,
					0.02,
					0.03,
					0.04,
					0.05,
					0.1,
					0.2,
					0.3,
					0.4,
					0.5
			};

			for (int i = 0; i < files.length; i++) {
				String file = files[i];

				TextFile tf = new TextFile(file, TextFile.R);
				String[] elems = tf.readLineElems(TextFile.tab);
				ArrayList<ArrayList<Double>> filebins = new ArrayList<>();
				for (int bin = 0; bin < bins.size(); bin++) {
					filebins.add(new ArrayList<Double>());
				}
				while (elems != null) {
					double val = Double.parseDouble(elems[1]);
					double maf = Double.parseDouble(elems[2]);

					for (int bin = 0; bin < upperthreshold.length; bin++) {
						if (maf >= lowerthreshold[bin] && maf <= upperthreshold[bin]) {
							filebins.get(bin).add(val);
						}
					}
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();
				bins.add(filebins);
			}

			BoxPlotPanel panel1 = new BoxPlotPanel(1, 1);
			grid = new Grid(1000, 1000, 1, 1, 0, 0);
			grid.addPanel(panel1);
			grid.draw(out);
			bins = null;

			// plot 3: maf vs maf
			grid = new Grid(1000, 1000, 1, files.length, 0, 0);
			for (int i = 0; i < files.length; i++) {
				String file = files[i];

				ArrayList<Double> fy = new ArrayList<>();
				ArrayList<Double> fx = new ArrayList<>();
				int ctr = 0;
				TextFile tf = new TextFile(file, TextFile.R);
				String[] elems = tf.readLineElems(TextFile.tab);
				while (elems != null) {
					double maf1 = Double.parseDouble(elems[1]);
					double maf2 = Double.parseDouble(elems[1]);
					fx.add(maf1);
					fy.add(maf2);
					ctr++;
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();
				ScatterplotPanel panel2 = new ScatterplotPanel(1, 1);

				grid.addPanel(panel2);
			}
			grid.draw(out);

			// plot 4: imputation qual vs beta

			grid = new Grid(1000, 1000, 1, files.length, 0, 0);
			for (int i = 0; i < files.length; i++) {
				String file = files[i];
				int ctr = 0;
				TextFile tf = new TextFile(file, TextFile.R);
				String[] elems = tf.readLineElems(TextFile.tab);
				ArrayList<Double> fy = new ArrayList<>();
				ArrayList<Double> fx = new ArrayList<>();
				while (elems != null) {
					double beta = Double.parseDouble(elems[1]);
					double impqual = Double.parseDouble(elems[1]);
					fx.add(beta);
					fy.add(impqual);
					ctr++;
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();

				ScatterplotPanel panel3 = new ScatterplotPanel(1, 1);
				grid.addPanel(panel3);
			}
			grid.draw(out);

			// plot 5: imputation qual vs correlation
			grid = new Grid(1000, 1000, 1, files.length, 0, 0);
			for (int i = 0; i < files.length; i++) {
				String file = files[i];
				int ctr = 0;
				TextFile tf = new TextFile(file, TextFile.R);
				String[] elems = tf.readLineElems(TextFile.tab);
				ArrayList<Double> fy = new ArrayList<>();
				ArrayList<Double> fx = new ArrayList<>();
				while (elems != null) {
					double impqual = Double.parseDouble(elems[1]);
					double correlation = Double.parseDouble(elems[1]);
					fx.add(impqual);
					fy.add(correlation);
					ctr++;
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();

				ScatterplotPanel panel3 = new ScatterplotPanel(1, 1);
				grid.addPanel(panel3);
			}
			grid.draw(out);


		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
