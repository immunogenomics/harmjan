package nl.harmjanwestra.utilities.association.approximatebayesposterior;

import JSci.maths.statistics.OutOfRangeException;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.AssociationResultPValuePair;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import umcg.genetica.containers.Triple;

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
			AssociationResultPValuePair p = new AssociationResultPValuePair(r, r.getAbf(), false);
			if (!Double.isInfinite(p.getP()) && !Double.isNaN(p.getP())) {
				pairs.add(p);
			}
		}

		ArrayList<AssociationResult> credibleSet = new ArrayList<AssociationResult>();
		if (!pairs.isEmpty()) {
			Collections.sort(pairs);

			double sum = 0;
			int ctr = 0;
			while (sum < bayesthreshold) {
				AssociationResultPValuePair p = pairs.get(ctr);
				double abf = p.getP();
				sum += abf;
				credibleSet.add(p.getAssociationResult());
				ctr++;
			}
		}

		return credibleSet;

	}

	public void calculateABF(ArrayList<AssociationResult> assocResults) {
		double nullVariance = (Math.log(1.5) / 1.96) * (Math.log(1.5) / 1.96); // 0.4
		ArrayList<Triple<Integer, Double, Double>> output = new ArrayList<Triple<Integer, Double, Double>>();
		double[] abfs = new double[assocResults.size()];
		double sum = 0;
		for (int i = 0; i < assocResults.size(); i++) {
			AssociationResult r = assocResults.get(i);
//			Double beta = Math.log(Math.abs(assoc.getMiddle()));
			double beta = Math.abs(r.getBeta());
			double se = r.getSe();


			if (!Double.isNaN(beta) && !Double.isNaN(se)) {
				double variance = se * se;
				try {
					JSci.maths.statistics.NormalDistribution n1 = new JSci.maths.statistics.NormalDistribution(0, variance);
					JSci.maths.statistics.NormalDistribution n2 = new JSci.maths.statistics.NormalDistribution(0, variance + nullVariance);
					double p0 = n1.probability(beta);
					double p1 = n2.probability(beta);

					double abf = p1 / p0;
					r.setBf(abf);
					if (!Double.isNaN(abf) && !Double.isInfinite(abf)) {
						abfs[i] = abf;
						sum += abf;
					}
				} catch (OutOfRangeException e) {
					System.out.println("Distribution out of range: b: " + beta + "\tv: " + variance);
				}

			}

		}


		for (int i = 0; i < assocResults.size(); i++) {
			if (abfs[i] != 0d) {
				double relativeABF = abfs[i] / sum;

				AssociationResult r = assocResults.get(i);
				r.setAbf(relativeABF);
			}
		}
	}

}
