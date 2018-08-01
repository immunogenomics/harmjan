package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.math.LogisticRegressionOptimized;
import nl.harmjanwestra.utilities.math.LogisticRegressionResult;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 7/14/17.
 */
public class LRTestTest {

	public static void main(String[] args) {


		try {
			String xfile = "/Data/tmp/2017-07-14/x.txt";
			String yfile = "/Data/tmp/2017-07-14/y.txt";

			ArrayList<double[]> datay = new ArrayList<>();
			ArrayList<double[]> datax = new ArrayList<>();

			TextFile tf = new TextFile(yfile, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length >= 2) {
					double[] tmp = new double[elems.length - 1];
					for (int i = 1; i < elems.length; i++) {
						tmp[i - 1] = Double.parseDouble(elems[i]);
					}
					datay.add(tmp);
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();


			tf = new TextFile(xfile, TextFile.R);
			elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length >= 2) {
					double[] tmp = new double[elems.length - 1];
					for (int i = 1; i < elems.length; i++) {
						tmp[i - 1] = Double.parseDouble(elems[i]);
					}
					datax.add(tmp);
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

			// normal analysis

			double[] y1 = new double[datay.size()];
			int size = (datay.get(0).length + 1);
			System.out.println(size);

			double[][] y2 = new double[datay.size()][];
			double[][] x1 = new double[datax.size()][2];
			for (int i = 0; i < datay.size(); i++) {
				y1[i] = datay.get(i)[0];

				double[] tmp = datay.get(i);
				y2[i] = tmp;
//				for (int j = 0; j < tmp.length; j++) {
//					y2[i][j] = tmp[j];
//				}
//				y2[i][y2[i].length - 1] = 0;

				x1[i][0] = 1;
				x1[i][1] = datax.get(i)[0];
			}

			System.out.println("Binomial");
			LogisticRegressionOptimized r = new LogisticRegressionOptimized();
			r.debug = true;

			LogisticRegressionResult result = r.binomial(y1, x1);
			printResult(result);
			System.exit(-1);
			System.out.println();
			System.out.println("Multinomial");

			LogisticRegressionResult result2 = r.multinomial(y2, x1);
			printResult(result2);
//

		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public static void printResult(LogisticRegressionResult result2) {

		double[][] beta = result2.getBeta();
		double[][] se = result2.getStderrs();

		for (int i = 0; i < beta.length; i++) {
			for (int j = 0; j < beta[i].length; j++) {
				System.out.println("disease: " + i + "\tcov: " + j + "\t" + beta[i][j] + "\t" + se[i][j]);
			}
		}
	}
}
