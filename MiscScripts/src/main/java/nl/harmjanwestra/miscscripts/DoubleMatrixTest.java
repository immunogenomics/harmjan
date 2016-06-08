package nl.harmjanwestra.miscscripts;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;

/**
 * Created by hwestra on 6/8/16.
 */
public class DoubleMatrixTest {

	public static void main(String[] args) {
		double[][] matrixA = new double[3][3];
		double[][] matrixB = new double[3][3];
		int ctr = 0;
		for (int i = 0; i < matrixA.length; i++) {
			for (int j = 0; j < matrixA.length; j++) {
				matrixA[i][j] = ctr;
				matrixB[i][j] = ctr;
				ctr++;
			}
		}


		DoubleMatrix2D mA = new DenseDoubleMatrix2D(matrixA);
		DoubleMatrix2D mB = new DenseDoubleMatrix2D(matrixB);



		DoubleMatrix2D mC = DoubleFactory2D.dense.appendColumns(mA,mB);
		double[][] output = mC.toArray();
		System.out.println(matrixA[0][1]);
	}

}
