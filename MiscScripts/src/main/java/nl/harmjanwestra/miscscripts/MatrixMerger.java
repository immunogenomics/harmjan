package nl.harmjanwestra.miscscripts;


import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.math.matrix.DoubleMatrixDataset;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Created by hwestra on 3/15/16.
 */
public class MatrixMerger {

	public static void main(String[] args) {
		MatrixMerger m = new MatrixMerger();
		try {
			m.rewriteMDSMatrix();
			m.merge();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void rewriteMDSMatrix() throws IOException {
		String f = "/Data/tmp/2016-03-15/157_pca.eigenvec";
		String f2 = "/Data/tmp/2016-03-15/157_pca.eigenvec.rewriteMDSMatrix.txt";
		String f2trans = "/Data/tmp/2016-03-15/157_pca.eigenvec.rewriteMDSMatrix_t.txt";
		TextFile tf = new TextFile(f, TextFile.R);
		TextFile tfout = new TextFile(f2, TextFile.W);

		String[] elems = tf.readLineElems(Strings.whitespace);

		String header = "-";
		for (int i = 2; i < elems.length; i++) {
			header += "\tMDS" + (i - 2);
		}
		tfout.writeln(header);

		while (elems != null) {
			tfout.writeln(Strings.concat(elems, Strings.tab, 1, elems.length));
			elems = tf.readLineElems(Strings.whitespace);
		}

		tfout.close();
		tf.close();
		DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<>(f2);
		ds.transposeDataset();
		ds.save(f2trans);
	}


	public void merge() throws IOException {

		String f1 = "/Data/tmp/2016-03-15/eqtl_157samples_pca2_pca_covariates.filtered.txt";
		String f2 = "/Data/tmp/2016-03-15/157_pca.eigenvec.rewriteMDSMatrix_t.txt";
		String out = "/Data/tmp/2016-03-15/157_covariatesMerged.txt";


		// pcs on rows, samples on columns
		DoubleMatrixDataset<String, String> ds1 = new DoubleMatrixDataset<>(f1);
		DoubleMatrixDataset<String, String> ds2 = new DoubleMatrixDataset<>(f2);

		HashMap<String, Integer> sharedSamples = new HashMap<String, Integer>();
		ArrayList<String> sampleArr = new ArrayList<String>();
		List<String> samples1 = ds1.colObjects;
		int ctr = 0;
		for (String sample : samples1) {
			if (ds2.hashCols.containsKey(sample)) {
				sharedSamples.put(sample, ctr);
				sampleArr.add(sample);
				ctr++;
			}
		}
		double[][] matout = new double[ds1.nrRows + ds2.nrRows][sharedSamples.size()];

		double[][] data1 = ds1.rawData;
		for (int col = 0; col < ds1.nrCols; col++) {
			String sample = ds1.colObjects.get(col);
			Integer id = sharedSamples.get(sample);
			if (id != null) {
				for (int row = 0; row < ds1.nrRows; row++) {
					matout[row][id] = data1[row][col];
				}
			}
		}

		double[][] data2 = ds2.rawData;
		for (int col = 0; col < ds2.nrCols; col++) {
			String sample = ds2.colObjects.get(col);
			Integer id = sharedSamples.get(sample);
			if (id != null) {
				for (int row = 0; row < ds2.nrRows; row++) {
					matout[row + ds1.nrRows][id] = data2[row][col];
				}
			}
		}

		ArrayList<String> newRows = new ArrayList<String>();
		newRows.addAll(ds1.rowObjects);
		newRows.addAll(ds2.rowObjects);

		DoubleMatrixDataset<String, String> ds3 = new DoubleMatrixDataset<>();
		ds3.setRawData(matout);
		ds3.rowObjects = newRows;
		ds3.colObjects = sampleArr;
		ds3.recalculateHashMaps();
		ds3.save(out);

	}


}
