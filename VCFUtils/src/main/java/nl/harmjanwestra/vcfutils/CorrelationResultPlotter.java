package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.Panel;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import umcg.genetica.containers.Pair;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by Harm-Jan on 12/24/15.
 */
public class CorrelationResultPlotter {

	public void run(String summaryFile, String output) throws IOException {

		CorrelationResultCombiner c = new CorrelationResultCombiner();
		ArrayList<CorrelationResult> data = c.getData(summaryFile);
		Grid grid = new Grid(400, 400, 1, 6 , 100, 100);

		// make panel with maf1 vs maf2

		grid.addPanel(makePanel(data, CorrelationResult.TYPE.maf1, CorrelationResult.TYPE.maf2, "Maf1", "Maf2"));
		grid.addPanel(makePanel(data, CorrelationResult.TYPE.maf1, CorrelationResult.TYPE.rsqb, "Maf1", "Beagle R-square"));
		grid.addPanel(makePanel(data, CorrelationResult.TYPE.maf1, CorrelationResult.TYPE.rsqp, "Maf1", "Pearson R-square"));

		grid.addPanel(makePanel(data, CorrelationResult.TYPE.maf2, CorrelationResult.TYPE.rsqb, "Maf2", "Beagle R-square"));
		grid.addPanel(makePanel(data, CorrelationResult.TYPE.maf2, CorrelationResult.TYPE.rsqp, "Maf2", "Pearson R-square"));
		grid.addPanel(makePanel(data, CorrelationResult.TYPE.rsqb, CorrelationResult.TYPE.rsqp, "Beagle R-square", "Pearson R-square"));

		grid.draw(output);



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
