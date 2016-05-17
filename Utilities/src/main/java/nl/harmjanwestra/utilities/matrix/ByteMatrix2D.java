package nl.harmjanwestra.utilities.matrix;

import cern.colt.matrix.AbstractMatrix2D;

/**
 * Created by hwestra on 5/16/16.
 */
public class ByteMatrix2D extends AbstractMatrix2D {

	private byte[] matrix;

	public ByteMatrix2D(byte[][] matrix) {
		rows = matrix.length;
		columns = matrix[0].length;
		this.matrix = new byte[rows * columns];
		for (int i = 0; i < rows; i++) {
			int start = i * columns;
			System.arraycopy(matrix[i], 0, this.matrix, start, columns);
		}
	}

	public ByteMatrix2D(int rows, int columns) {
		this.rows = rows;
		this.columns = columns;
		this.matrix = new byte[rows * columns];
	}

	public byte getQuick(int i, int j) {
		int index = i * columns + j;
		return matrix[index];
	}

	public void setQuick(int i, int j, byte b) {
		int index = i * columns + j;
		matrix[index] = b;
	}

	public byte[][] toArray() {
		byte[][] output = new byte[rows][columns];
		for (int i = 0; i < rows; i++) {
			int start = i * columns;
			System.arraycopy(this.matrix, start, output[i], 0, columns);
		}
		return output;
	}


}
