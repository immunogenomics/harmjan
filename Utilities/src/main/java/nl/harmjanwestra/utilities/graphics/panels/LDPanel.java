package nl.harmjanwestra.utilities.graphics.panels;

import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 * Created by hwestra on 5/5/16.
 */
public class LDPanel extends Panel {
	private DoubleMatrixDataset<String, String> ds;

	public LDPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}


	public void setData(DoubleMatrixDataset<String, String> ds) {
		this.ds = ds;
	}

	@Override
	public void draw(DefaultGraphics g) {

		int cols = ds.columns();
		int rows = ds.rows();

		int defaultBlockSize = 10;
		double pixelsPerBlock = (double) width / defaultBlockSize;
		if (pixelsPerBlock < 1) {
			pixelsPerBlock = 1;
		}


	}
}
