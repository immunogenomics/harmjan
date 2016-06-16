package nl.harmjanwestra.utilities.vcf;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by hwestra on 6/6/16.
 * Adapted from: SNPSummaryComputation.cpp
 * <p>
 * //          Copyright Gavin Band 2008 - 2012.
 * // Distributed under the Boost Software License, Version 1.0.
 * //    (See accompanying file LICENSE_1_0.txt or copy at
 * //          http://www.boost.org/LICENSE_1_0.txt)
 */
public class VCFImputationQualScoreImpute {


	public double getInfo() {
		return info;
	}

	public double getImpinfo() {
		return impinfo;
	}

	public static void main(String[] args) {
		String f = "/Data/tmp/2016-06-06/test-sorted.vcf";
		VCFImputationQualScoreImpute c = new VCFImputationQualScoreImpute();
		try {
			int ctr = 0;
			TextFile tf = new TextFile(f, TextFile.R);
			String ln = tf.readLine();
			while (ln != null) {
				if (!ln.startsWith("#")) {
					VCFVariant var = new VCFVariant(ln, VCFVariant.PARSE.ALL);

					c.computeAutosomal(var);

				}

				if (ctr > 10) {
					ln = null;
				} else {
					ln = tf.readLine();
				}
				ctr++;
			}


			tf.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	double info = 0;
	double impinfo = 0;

	// this only works for biallelic variants
	public void computeAutosomal(VCFVariant variant) {
		DoubleMatrix2D probs = variant.getGenotypeProbabilies(); //  samples | ALLELES

//		double theta_mle = 0; // sum gprob(AB)+ 2 * gprob(BB) / 2 * sum(gprob)
		double theta_mle = (probs.viewColumn(1).zSum() + 2.0 * probs.viewColumn(2).zSum()) / (2.0 * probs.zSum());
		double[] imputeFallBackDist = new double[3];
		double[] fallBackDist = new double[3];
		fallBackDist[0] = (1 - theta_mle) * (1 - theta_mle);
		fallBackDist[1] = 2 * theta_mle * (1 - theta_mle);
		fallBackDist[2] = theta_mle * theta_mle;


		info = 1 - compExpVar(probs, fallBackDist) / (2 * theta_mle * (1 - theta_mle));
		impinfo = 1 - compExpVar(probs, imputeFallBackDist) / (2 * theta_mle * (1 - theta_mle));

		if (Double.isNaN(info)) {
			info = 0;
		}

		if (Double.isNaN(impinfo)) {
			impinfo = 0;
		}

//		VCFImputationQualScoreBeagle vbs = new VCFImputationQualScoreBeagle(variant, true);
//		System.out.println(variant.toString() + "\t" + info + "\t" + impinfo + "\t" + vbs.allelicR2() + "\t" + vbs.doseR2());


	}

	final double[] levels = new double[]{0, 1, 2};
	final double[] levelsSquared = new double[]{0, 1, 4};

	private double compExpVar(DoubleMatrix2D probs, double[] fallBackDist) {

		double result = 0;
		for (int i = 0; i < probs.rows(); i++) { // probabilities.rows() == nrSamples
			double c = probs.viewRow(i).zSum();
			double[] indprobs = new double[3];
			indprobs[0] = probs.getQuick(i, 0) + (1 - c) * fallBackDist[0];
			indprobs[1] = probs.getQuick(i, 1) + (1 - c) * fallBackDist[1];
			indprobs[2] = probs.getQuick(i, 2) + (1 - c) * fallBackDist[2];
			double v = compVar(indprobs);
			result += v;
		}

		double output = result / probs.rows();
		if (Double.isNaN(output)) {
			return 1;
		} else {
			return output;
		}
	}


	private double compVar(double[] probs) {


		double mean = 0;
		for (int i = 0; i < levels.length; i++) {
			mean += probs[i] * levels[i];
		}


		double variance = 0;
		for (int i = 0; i < levels.length; i++) {
			variance += probs[i] * levelsSquared[i];
		}

		variance -= (mean * mean);


		// std::cerr << "compute_variance: levels = " << levels.transpose() << ", levels_squared = " << levels_squared.transpose() << ", probs = " << probs.transpose() << ", variance = " << variance << ".\n" ;
		return variance;

	}

}
