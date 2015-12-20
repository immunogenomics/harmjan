package hms.hwestra.harmonics;

import cern.jet.random.NegativeBinomial;
import cern.jet.random.engine.DRand;
import hms.hwestra.harmonics.graphics.HistogramPlot;
import umontreal.iro.lecuyer.probdist.NegativeBinomialDist;

/**
 * Created by hwestra on 6/29/15.
 */
public class Tests {

	public static void main(String[] args) {


		try {

			int nrsamples = 1000000;
			int mean = 10;
			double sigma = 3;

			NegativeBinomial d = new NegativeBinomial(2, 0.3, new DRand());

			int[] obs = new int[nrsamples];
			int sum = 0;
			for (int i = 0; i < nrsamples; i++) {
				obs[i] = d.nextInt();

			}

			int maxCt = JSci.maths.ArrayMath.max(obs);


			String outfilenae = "/Data/tmp/HMTest/plot.pdf";
			int width = 1000;
			int height = 1000;

			int margin = 100;

			double[][] histograms = new double[2][maxCt + 1];

			HistogramPlot.PLOTTYPE[] types = new HistogramPlot.PLOTTYPE[2];
			types[0] = HistogramPlot.PLOTTYPE.BAR;
			types[1] = HistogramPlot.PLOTTYPE.POLY;

			NegativeBinomialDist umontreal = NegativeBinomialDist.getInstanceFromMLE(obs, nrsamples);

			hms.hwestra.harmonics.math.NegativeBinomial hjnb = new hms.hwestra.harmonics.math.NegativeBinomial(obs);
//			double[] params = hjnb.estimateParameters();
//			System.out.println("HJ1: " + params[0]);
//			System.out.println("HJ2: " + params[1]);

			System.out.println("N: " + umontreal.getN() + "\tP: " + umontreal.getP());

			for (int i = 0; i < nrsamples; i++) {
				histograms[0][obs[i]]++;
			}

			for (int i = 0; i < histograms[1].length; i++) {
				sum += histograms[0][i];
			}

			for (int i = 0; i < histograms[1].length; i++) {
				histograms[0][i] /= sum;
				histograms[1][i] = umontreal.prob(i);
			}


			HistogramPlot plot = new HistogramPlot(outfilenae, width, height);
			plot.setMargin(margin);
			plot.setData(histograms, types);

			plot.draw();


		} catch (Exception e) {

			e.printStackTrace();
		}
	}
}
