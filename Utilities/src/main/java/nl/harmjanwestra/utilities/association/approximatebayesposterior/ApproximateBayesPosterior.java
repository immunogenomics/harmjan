package nl.harmjanwestra.utilities.association.approximatebayesposterior;

import JSci.maths.statistics.OutOfRangeException;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.AssociationResultPValuePair;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;
import java.util.Collections;

/**
 * Created by hwestra on 10/1/15.
 */
public class ApproximateBayesPosterior {

	NormalDistribution standardNormal = new NormalDistribution();

	/**
	 * Determine the expected variance for the Bayesian prior distribution
	 *
	 * @param thetaprior expected effect size (thetaprior = log(thetaexpected))
	 * @param e          probability for the expected effect size (e = (1-eta))
	 * @return the expected variance
	 */
	public double calculateW(double thetaprior, double e) {
		double w = thetaprior / standardNormal.inverseCumulativeProbability(1 - (e / 2));
		return (w * w);
	}


	public double calculateBFDP(double thetahat, double se, double p0, double alpha, double w) {
		double thetasd = Math.sqrt(se);

		NormalDistribution priorH0 = new NormalDistribution(0, w);
		NormalDistribution priorH1 = new NormalDistribution(0, se + w);

		double abfNormal = priorH0.probability(thetahat) / priorH1.probability(thetahat);

		double zsquared = thetahat / thetasd;
		double r = w / (se + w);
		double abf = (1 / (1 - r)) * Math.exp(-(zsquared / 2) * r);

		double bfdp = (abf * p0) / ((abf * p0) + 1);

		// determine confidence interval for theta
		LogNormalDistribution logNormalDistribution = new LogNormalDistribution((r * thetahat), (r * se));
		double posteriorTheta = logNormalDistribution.inverseCumulativeProbability(thetahat);
		double rthetahat = (r * thetahat);

		//double interval =
		double upperbound = Math.exp(rthetahat + standardNormal.inverseCumulativeProbability(alpha / 2) * Math.sqrt(r * se));
		double lowerbound = Math.exp(rthetahat - standardNormal.inverseCumulativeProbability(alpha / 2) * Math.sqrt(r * se));


		return abf;
	}

	public boolean isNoteWorthy(double cAlpha, double cBeta, double pr) {
		double ratio = cBeta / cAlpha;
		return (pr < ratio / (1 + ratio));

	}


	public ArrayList<AssociationResult> createCredibleSet(ArrayList<AssociationResult> associationResults, double bayesthreshold) {

		ArrayList<AssociationResultPValuePair> pairs = new ArrayList<AssociationResultPValuePair>();
		for (AssociationResult r : associationResults) {
			AssociationResultPValuePair p = new AssociationResultPValuePair(r, r.getPosterior(), false);
			if (!Double.isInfinite(p.getP()) && !Double.isNaN(p.getP())) {
				pairs.add(p);
			}
		}

		ArrayList<AssociationResult> credibleSet = new ArrayList<AssociationResult>();
		if (!pairs.isEmpty()) {
			Collections.sort(pairs);

			double sum = 0;
			int ctr = 0;
			while (sum < bayesthreshold && ctr < pairs.size()) {
				AssociationResultPValuePair p = pairs.get(ctr);
				double abf = p.getP();
				sum += abf;
				credibleSet.add(p.getAssociationResult());
				ctr++;
			}

			if (sum < bayesthreshold && ctr == pairs.size()) {
				System.out.println("Error when determining credible set: sigmaPosterior == " + sum + " after adding " + ctr + " out of " + pairs.size());
			}
		}

		return credibleSet;

	}

	public void calculatePosterior(ArrayList<AssociationResult> assocResults) {
		double prior = 1.5;
		double sum = 0;
		for (int i = 0; i < assocResults.size(); i++) {
			AssociationResult r = assocResults.get(i);
//			Double beta = Math.log(Math.abs(gwas.getMiddle()));
			double totalABF = 0;

			for (int b = 0; b < r.getBeta().length; b++) {
				double beta = Math.abs(r.getBeta()[b]);
				double se = r.getSe()[b];

				if (!Double.isNaN(beta) && !Double.isNaN(se)) {


					double abf = abf(beta, se, prior);
					totalABF += abf;


					if (!Double.isNaN(abf) && !Double.isInfinite(abf)) {
						sum += abf;
					}

				}
			}
			r.setBf(totalABF);
		}

		for (int i = 0; i < assocResults.size(); i++) {
			AssociationResult r = assocResults.get(i);
			double bf = r.getBf();
			if (!Double.isNaN(bf) && !Double.isInfinite(bf)) {
				double posterior = bf / sum;
				if (Double.isNaN(posterior) || Double.isInfinite(posterior)) {
					r.setPosterior(0);
				} else {
					r.setPosterior(posterior);
				}
			}
		}

		// abf
	}

	public static void main(String[] args) {
		ApproximateBayesPosterior p = new ApproximateBayesPosterior();
		System.out.println(p.abf(1.316, 0.10, 2));
	}


	private double abf(double beta, double se, double prior) {
		double variance = se * se;
		double theta = beta; // Math.log(beta);
		try {
			double nullVariance = (Math.log(prior) / 1.96) * (Math.log(prior) / 1.96); // 0.04

			JSci.maths.statistics.NormalDistribution n1 = new JSci.maths.statistics.NormalDistribution(0, variance);
			JSci.maths.statistics.NormalDistribution n2 = new JSci.maths.statistics.NormalDistribution(0, variance + nullVariance);

			double p0 = n1.probability(theta);
			double p1 = n2.probability(theta);
			double abf = p1 / p0;
			return abf;


		} catch (OutOfRangeException e) {
			System.out.println("Distribution out of range: b: " + beta + "\tv: " + variance);
		}
		/*
		double z = beta / se;
						double r2 = nullVariance / (variance + nullVariance);
						double otherbf = 1 / Math.sqrt(1 - r2) * Math.exp(-((z * z) / 2) * r2);
		 */

		return 0;
	}

	private double abf2(double beta, double se, double prior) {
		double variance = se * se;
		double theta = beta; // Math.log(beta);
		double nullVariance = (Math.log(prior) / 1.96) * (Math.log(prior) / 1.96); // 0.04
		double r = variance / (variance + nullVariance);
		double z = theta / se;
		double abf = (1 / (Math.sqrt(1 - r))) * Math.exp(-(z * z) / 2 * r);
		return abf;

	}

}
