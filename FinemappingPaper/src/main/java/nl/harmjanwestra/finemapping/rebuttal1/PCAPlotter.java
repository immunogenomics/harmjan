package nl.harmjanwestra.finemapping.rebuttal1;

import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.math.matrix2.DoubleMatrixDataset;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;

public class PCAPlotter {
	
	
	public static void main(String[] args) {
	
	}
	
	public void draw(String f, String pcamatrix) throws IOException {
		// read sample annotation
		
		TextFile tf = new TextFile(f, TextFile.R);
		String[] elems = tf.readLineElems(Strings.tab);
		while (elems != null) {
			String sample = elems[0];
			String dataset = elems[1];
			elems = tf.readLineElems(Strings.tab);
		}
		tf.close();
		
		
		DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(pcamatrix);
		
		Grid grid = new Grid(100, 100, 10, 10, 25, 25);
		
		// put each row in a group
		
		for (int i = 0; i < ds.columns(); i++) {
			double[] x = ds.getMatrix().viewColumn(i).toArray();
			for (int j = i + 1; j < ds.columns(); j++) {
				double[] y = ds.getMatrix().viewColumn(i).toArray();
				
				// matrix-ify?
				
				// panel!
				ScatterplotPanel panel = new ScatterplotPanel(1, 1);
				panel.setData(x, y);
				panel.setLabels("PCA " + i, "PCA " + j);
				
				
			}
		}
		
	}
}
