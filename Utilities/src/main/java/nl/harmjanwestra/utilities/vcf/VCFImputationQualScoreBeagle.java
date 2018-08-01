/*
 * Copyright (C) 2014 Brian L. Browning
 *
 * This file is part of Beagle
 *
 * Beagle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Beagle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package nl.harmjanwestra.utilities.vcf;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.text.DecimalFormat;
import java.util.Arrays;

/**
 * Adapted by Harm-Jan on 06/02/16.
 * <p>
 * Modified from GprobsStatistics.java in the Beagle imputation software package.
 * https://faculty.washington.edu/browning/beagle/beagle.html
 * Copyright (C) 2014 Brian L. Browning
 * <p>
 * This path is part of Beagle
 */
public class VCFImputationQualScoreBeagle {

	private static final double MIN_R2_DEN = 1e-8;

	private final VCFVariant marker;
	private final int nSamples;
	private final double[] alleleFreq;

	private double sumCall = 0;
	private double sumSquareCall = 0;
	private double sumExpected = 0;
	private double sumExpectedSquare = 0;
	private double sumSquareExpected = 0;
	private double sumCallExpected = 0;

	/**
	 * Constructs a new {@code GprobsStatistics} instance from the
	 * specified scaled genotype probabilities.
	 *
	 * @throws IndexOutOfBoundsException if
	 *                                   {@code marker < 0 || marker >= gv.nMarkers()}
	 * @throws NullPointerException      if {@code gv == null}
	 */
	public VCFImputationQualScoreBeagle(VCFVariant variant, boolean useAlleleProbs) {
		int nAlleles = variant.getAlleles().length;
		this.marker = variant;
		this.nSamples = variant.getNrSamples();
		this.alleleFreq = new double[nAlleles];
		double[] alProbs = new double[nAlleles];
		double[] gtProbs = new double[3];


		for (int j = 0; j < this.nSamples; ++j) {
			setProbs(j, gtProbs, alProbs, useAlleleProbs);
			for (int a = 0; a < nAlleles; ++a) {
				alleleFreq[a] += alProbs[a];
			}
			int call = maxIndex(gtProbs);
			double exp = (gtProbs[1] + 2 * gtProbs[2]);
			double expSquare = (gtProbs[1] + 4 * gtProbs[2]);
			sumCall += call;
			sumSquareCall += call * call;
			sumExpected += exp;
			sumExpectedSquare += expSquare;
			sumSquareExpected += (exp * exp);
			sumCallExpected += (call * exp);

		}


		double sum = sum(alleleFreq);
		if (variant.getId().equals("rs113022009")) {
			System.out.println("got it!");
		}
		divideBy(alleleFreq, sum);
	}

	private static String format(DecimalFormat df, double d) {
		if (Double.isNaN(d)) {
			return "NaN";
		} else if (d == Double.POSITIVE_INFINITY) {
			return "Infinity";
		} else if (d == Double.NEGATIVE_INFINITY) {
			return "-Infinity";
		} else {
			return df.format(d);
		}
	}

	private void setProbs(int sample,
						  double[] gtProbs, double[] alProbs, boolean useAlleleProbs) {
		Arrays.fill(gtProbs, 0.0f);
		Arrays.fill(alProbs, 0.0f);
		if (useAlleleProbs) {
			int gt = 0;
			DoubleMatrix2D probs = marker.getGenotypeProbabilies();
			for (int a2 = 0; a2 < alProbs.length; ++a2) {
				for (int a1 = 0; a1 <= a2; ++a1) {
					double gprob = probs.get(sample, gt++); // gv.value(marker, sample, gt++);
					alProbs[a1] += gprob;
					alProbs[a2] += gprob;
					if (a2 == 0) {
						gtProbs[0] += gprob;
					} else if (a1 == 0) {
						gtProbs[1] += gprob;
					} else {
						gtProbs[2] += gprob;
					}
				}
			}
			double sum = sum(gtProbs);
			divideBy(gtProbs, sum);
			divideBy(alProbs, 2 * sum);
		} else {
			// use genotype dosages in stead
//
			int gt = 0;

			/*
			double[] alProbs = new double[nAlleles];
		double[] gtProbs = new double[3];
			 */

			DoubleMatrix2D probs = marker.getGenotypeAllelesAsMatrix2D();
			byte a1 = (byte) probs.getQuick(sample, 0);
			byte a2 = (byte) probs.getQuick(sample, 1);

			alProbs[a1] += 0.5;
			alProbs[a2] += 0.5;
			gtProbs[a1 + a2] = 1;
//
//			Arrays.fill(gtProbs, 0.0f);
//			Arrays.fill(alProbs, 0.0f);
//			int gt = 0;
//			for (int a2 = 0; a2 < alProbs.length; ++a2) {
//				for (int a1 = 0; a1 <= a2; ++a1) {
//					float gprob = gv.value(marker, sample, gt++);
//					alProbs[a1] += gprob;
//					alProbs[a2] += gprob;
//					if (a2 == 0) {
//						gtProbs[0] += gprob;
//					} else if (a1 == 0) {
//						gtProbs[1] += gprob;
//					} else {
//						gtProbs[2] += gprob;
//					}
//				}
//			}
//			double sum = sum(gtProbs);
//			divideBy(gtProbs, sum);
//			divideBy(alProbs, 2 * sum);
		}

	}

	private int maxIndex(double[] fa) {
		int maxIndex = 0;
		for (int j = 1; j < fa.length; ++j) {
			if (fa[j] > fa[maxIndex]) {
				maxIndex = j;
			}
		}
		return maxIndex;
	}

	private double sum(double[] fa) {
		double sum = 0.0f;
		for (double f : fa) {
			sum += f;
		}
		return sum;
	}

	private void divideBy(double[] fa, double divisor) {
		for (int j = 0; j < fa.length; ++j) {
			fa[j] /= divisor;
		}
	}

	/**
	 * Returns the marker.
	 *
	 * @return the marker.
	 */
	public VCFVariant marker() {
		return marker;
	}

	/**
	 * Returns an array of length {@code this.marker().nAlleles()} whose
	 * {@code j}-th element is the estimated sample frequency of allele
	 * {@code j}.
	 *
	 * @return an array of length {@code this.marker().nAlleles()} whose
	 * {@code j}-th element is the estimated sample frequency of allele
	 * {@code j}
	 */
	public double[] alleleFreq() {
		return alleleFreq.clone();
	}

	/**
	 * Returns the estimated squared correlation between the most probable
	 * ALT allele dose and the true ALT allele dose.
	 * Returns 0 if the marker is monomorphic or if most probable ALT
	 * allele dose is monomorphic.
	 *
	 * @return the estimated squared correlation between the most likely
	 * allele dose and the true allele dose
	 */
	public double allelicR2() {
		double f = 1.0f / nSamples;
		double cov = sumCallExpected - (sumCall * sumExpected * f);
		double varBest = sumSquareCall - (sumCall * sumCall * f);
		double varExp = sumExpectedSquare - (sumExpected * sumExpected * f);
		double den = varBest * varExp;
		return (den < MIN_R2_DEN) ? 0.0f : (cov * cov / den);
	}

	/**
	 * Returns the estimated squared correlation between the estimated
	 * ALT allele dose and the true ALT allele dose.  Returns 0 if the
	 * marker is monomorphic.
	 *
	 * @return the estimated squared correlation between the estimated
	 * ALT allele dose and the true ALT allele dose
	 */
	public double doseR2() {
		double f = 1.0f / (double) nSamples;
		double num = sumSquareExpected - (sumExpected * sumExpected * f);
		double den = sumExpectedSquare - (sumExpected * sumExpected * f);
		if (num < 0.0) {
			num = 0.0;
		}
		return (den < MIN_R2_DEN) ? 0.0f : (num / den);
	}

	/**
	 * Returns the estimated squared correlation between the estimated
	 * ALT allele dose and the true ALT allele dose where the variance of
	 * the true ALT allele dose is estimated from the estimated
	 * ALT allele frequency. Returns 0 if the marker is monomorphic.
	 *
	 * @return the estimated squared correlation between the estimated
	 * ALT allele dose and the true ALT allele dose
	 */
	public double hweDoseR2() {
		double f = 1.0f / nSamples;
		double altFreq = sumExpected / (2.0f * nSamples);
		double num = (sumSquareExpected - (sumExpected * sumExpected * f)) / nSamples;
		double den = 2.0f * altFreq * (1.0f - altFreq);
		if (num < 0.0) {
			num = 0.0;
		}
		return (den < MIN_R2_DEN) ? 0.0f : (num / den);
	}

	/**
	 * Returns a string representation of {@code this}.  The exact
	 * details of the representation are unspecified and subject to change.
	 *
	 * @return a string representation of {@code this}
	 */
	@Override
	public String toString() {
		DecimalFormat df = new DecimalFormat("0.####");
		StringBuilder sb = new StringBuilder(80);
		sb.append(marker);
		sb.append(Strings.tab);
		for (int j = 0; j < alleleFreq.length; ++j) {
			sb.append((j == 0) ? "AF=" : Strings.comma);
			sb.append(alleleFreq[j]);
		}
		sb.append(Strings.tab);
		sb.append("AR2=");
		sb.append(format(df, allelicR2()));
		sb.append(Strings.tab);
		sb.append("DR2=");
		sb.append(format(df, doseR2()));
		sb.append(Strings.tab);
		sb.append("HDR2=");
		sb.append(format(df, hweDoseR2()));
		return sb.toString();
	}

}
