package nl.harmjanwestra.vcfutils;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.Panel;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import nl.harmjanwestra.utilities.graphics.panels.SpacerPanel;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

/**
 * Created by Harm-Jan on 12/24/15.
 */
public class CorrelationResultPlotter {

	public void run(String summaryFile, String output) throws IOException, DocumentException {

		CorrelationResultCombiner c = new CorrelationResultCombiner();
		ArrayList<CorrelationResult> data = c.getData(summaryFile);
		Grid grid = new Grid(400, 400, 2, 4, 100, 100);

		// make panel with maf1 vs maf2
		grid.addPanel(makePanel(data, CorrelationResult.TYPE.maf1, CorrelationResult.TYPE.maf2, "Maf1", "Maf2"), 0, 0);
		grid.addPanel(makePanel(data, CorrelationResult.TYPE.maf1, CorrelationResult.TYPE.rsqb1, "Maf1", "Beagle R-square"), 0, 1);
		grid.addPanel(makePanel(data, CorrelationResult.TYPE.maf1, CorrelationResult.TYPE.rsqp, "Maf1", "Pearson R-square"), 0, 2);
		grid.addPanel(new SpacerPanel(1, 1), 0, 3);
		grid.addPanel(makePanel(data, CorrelationResult.TYPE.maf2, CorrelationResult.TYPE.rsqb1, "Maf2", "Beagle R-square"), 1, 1);
		grid.addPanel(makePanel(data, CorrelationResult.TYPE.maf2, CorrelationResult.TYPE.rsqp, "Maf2", "Pearson R-square"), 1, 2);
		grid.addPanel(makePanel(data, CorrelationResult.TYPE.rsqb1, CorrelationResult.TYPE.rsqp, "Beagle R-square", "Pearson R-square"), 1, 3);
		grid.addPanel(new SpacerPanel(1, 1), 1, 0);


		grid.draw(output);

	}

	public void run(String listoffiles, String nameslist, String outfile) throws IOException, DocumentException {

		CorrelationResultCombiner c = new CorrelationResultCombiner();

		String[] files = listoffiles.split(",");
		String[] names = nameslist.split(",");

		ArrayList<ArrayList<CorrelationResult>> alldata = new ArrayList<>();
		HashSet<String> allVariants = new HashSet<String>();
		for (int i = 0; i < files.length; i++) {
			System.out.println("Reading: " + files[i]);
			HashSet<String> variantsInDs = new HashSet<String>();
			ArrayList<CorrelationResult> data = c.getData(files[i]);
			alldata.add(data);
			for (CorrelationResult r : data) {

				if (variantsInDs.contains(r.variant)) {
					System.out.println("already found variant: " + r.variant);
				}
				variantsInDs.add(r.variant);

				allVariants.add(r.variant);
			}
		}

		System.out.println(allVariants.size() + " variants in total.");

		Grid grid = new Grid(400, 300, 1, 1, 100, 100);

		ScatterplotPanel panel = new ScatterplotPanel(1, 1);

		double[][] x = new double[files.length][allVariants.size()];
		double[][] y = new double[files.length][allVariants.size()];

		Range dataRange = new Range(0, 0, 1, 1);
		int ctr = 0;
		String[] refNames = new String[names.length];
		for (int i = 0; i < files.length; i++) {

			ArrayList<CorrelationResult> data = alldata.get(i);
			ArrayList<Double> d = new ArrayList<>();

			for (CorrelationResult r : data) {
				Double p = r.rsqPearson;
				if (p != null && !Double.isInfinite(p) && !Double.isNaN(p)) {
					d.add(p);
				}
			}

			Collections.sort(d, Collections.reverseOrder());
			int nrAboveThreshold = 0;
			for (int v = 0; v < d.size(); v++) {
				if (v > 0) {
					x[ctr][v] = (double) v / allVariants.size();
					y[ctr][v] = d.get(v);

				} else {
					x[ctr][v] = 0;
					y[ctr][v] = d.get(v);
				}
			}
			for (int v = d.size(); v < y[0].length; v++) {
				y[ctr][v] = Double.NaN;
				x[ctr][v] = Double.NaN;
			}
			refNames[ctr] = names[i] + " (" + d.size() + ")";
			ctr++;
		}


		panel.setData(x, y);
		panel.setDataRange(dataRange);
		panel.setTitle("Correlation with genotyped variants (" + allVariants.size() + " total variants)");
		panel.setLabels("Variants", "Pearson R-squared");

		panel.setDatasetLabels(refNames);

		grid.addPanel(panel);


		grid.draw(outfile);
		System.out.println("Done. Plot is here: " + outfile);


	}


	private Panel makePanel(ArrayList<CorrelationResult> data, CorrelationResult.TYPE t1, CorrelationResult.TYPE t2, String xLabel, String yLabel) {
		ScatterplotPanel panel = new ScatterplotPanel(1, 1);
		Pair<double[], double[]> set = collect(data, t1, t2);
		panel.setData(set.getLeft(), set.getRight());
		panel.setLabels(xLabel, yLabel);
		return panel;
	}

	private Pair<double[], double[]> collect(ArrayList<CorrelationResult> data, CorrelationResult.TYPE t1, CorrelationResult.TYPE t2) {
		double[] x = new double[data.size()];
		double[] y = new double[data.size()];
		for (int i = 0; i < data.size(); i++) {
			CorrelationResult result = data.get(i);
			x[i] = result.get(t1);
			x[i] = result.get(t2);
		}
		return new Pair<double[], double[]>(x, y);
	}

}
