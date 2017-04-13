package nl.harmjanwestra.ngs;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.math.matrix.DoubleMatrixDataset;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

/**
 * Created by hwestra on 2/17/16.
 */
public class MatrixFilter {

	public static void main(String[] args) {
		try {

			MatrixFilter f = new MatrixFilter();
			f.filterRows("/Data/tmp/nixup/mixups-matrix.txt.gz", "/Data/tmp/nixup/mixups-matrix-filtered.txt", "/Data/tmp/nixup/selectsamples.txt");
			f.transpose("/Data/tmp/nixup/mixups-matrix-filtered.txt", "/Data/tmp/nixup/mixups-matrix-filtered-transposed.txt");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void transpose(String s, String s1) throws IOException {
		DoubleMatrixDataset<String, String> d = new DoubleMatrixDataset<String, String>(s);
		d.transposeDataset();
		d.save(s1);

	}

	private HashSet<String> loadInclude(String filteron) throws IOException {
		TextFile tf = new TextFile(filteron, TextFile.R);
		HashSet<String> include = new HashSet<String>();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			include.add(elems[0]);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return include;
	}

	public void filterRows(String file, String out, String filteron) throws IOException {
		HashSet<String> include = loadInclude(filteron);

		TextFile matin = new TextFile(file, TextFile.R);
		TextFile matout = new TextFile(out, TextFile.W);
		matout.writeln(matin.readLine());
		String[] elems = matin.readLineElems(TextFile.tab);
		while (elems != null) {
			if (include.contains(elems[0])) {
				matout.writeln(Strings.concat(elems, Strings.tab));
			}
			elems = matin.readLineElems(TextFile.tab);
		}
		matout.close();
		matin.close();

	}

	public void filterCols(String file, String out, String filteron) throws IOException {
		HashSet<String> include = loadInclude(filteron);

		TextFile matin = new TextFile(file, TextFile.R);
		TextFile matout = new TextFile(out, TextFile.W);
		String[] elems = matin.readLineElems(TextFile.tab);
		boolean[] includecol = new boolean[elems.length];
		includecol[0] = true;

		for (int i = 1; i < elems.length; i++) {
			if (include.contains(elems[i])) {
				includecol[i] = true;
			}
		}

		matout.writeln(Strings.concat(elems, includecol, Strings.tab));

		elems = matin.readLineElems(TextFile.tab);
		while (elems != null) {

			matout.writeln(Strings.concat(elems, includecol, Strings.tab));

			elems = matin.readLineElems(TextFile.tab);
		}
		matout.close();
		matin.close();

	}
}
