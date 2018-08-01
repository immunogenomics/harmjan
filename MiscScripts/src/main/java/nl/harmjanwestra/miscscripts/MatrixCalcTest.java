package nl.harmjanwestra.miscscripts;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;

/**
 * Created by hwestra on 7/8/16.
 */
public class MatrixCalcTest {

	public static void main(String[] args) {
		DoubleMatrix2D mat = new DenseDoubleMatrix2D(4, 3);
		mat.setQuick(0, 0, 0);
		mat.setQuick(0, 1, 0);
		mat.setQuick(0, 2, 0);
		mat.setQuick(1, 0, 0);
		mat.setQuick(1, 1, 0);
		mat.setQuick(1, 2, 1);
		mat.setQuick(2, 0, 1);
		mat.setQuick(2, 1, 0);
		mat.setQuick(2, 2, 1);
		mat.setQuick(3, 0, 1);
		mat.setQuick(3, 1, 1);
		mat.setQuick(3, 2, 0);

		DoubleMatrix2D mat2 = new DenseDoubleMatrix2D(4, 2);
		mat2.setQuick(0, 0, 0);
		mat2.setQuick(0, 1, 0);
		mat2.setQuick(1, 0, 0);
		mat2.setQuick(1, 1, 0);

		cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra d= new DenseDoubleAlgebra();
		DoubleMatrix2D mat3 = d.mult(mat,mat2);
		for (int i = 0; i < mat3.rows(); i++) {
			String ln = "";
			for (int j = 0; j < mat3.columns(); j++) {
				ln += mat3.get(i, j);
			}
			System.out.println(ln);
		}
	}

}
