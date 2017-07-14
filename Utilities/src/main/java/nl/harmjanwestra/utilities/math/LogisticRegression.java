package nl.harmjanwestra.utilities.math;

/**
 * Created by hwestra on 10/26/15.
 * Derived from Scott. A. Czepiel's implementation (http://czep.net/)
 */

import cern.jet.stat.Gamma;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.ChiSquare;

public class LogisticRegression {

	int max_iter = 50;
	double EPSILON = 1E-6;

	public LogisticRegressionResult univariate(double[][] y, double[][] x) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("Error in logistic regression: length of y does not match length of x");
		}
		int J = 2;
		int N = x.length;
		int K = x[0].length;
		double[] n = new double[y.length];
		for (int i = 0; i < y.length; i++) {
			n[i] = 1;
		}
		return mlelr(J, N, K, n, x, y);
	}

	public LogisticRegressionResult univariate(double[] y, double[][] x) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("Error in logistic regression: length of y does not match length of x");
		}
		double[][] yreg = new double[x.length][2];
		double[] n = new double[y.length];
		for (int i = 0; i < y.length; i++) {
			if (y[i] > 1 || y[i] < 0) {
				throw new IllegalArgumentException("Error in univariate logistic regression; unexpected categorical value: " + y[i] + ". Expected 0 or 1");
			}
			yreg[i][(int) y[i]] = 1;
			n[i] = 1;
		}

		int J = 2;
		int N = x.length;
		int K = x[0].length;

		return mlelr(J, N, K, n, x, yreg);

	}

	public LogisticRegressionResult mlelr(int J,
										  int N,
										  int K,
										  double[] n,
										  double[][] X,
										  double[][] Y) {

		int i, j, k;

		int iter = 0;
		boolean converged = false;

		int kJMinusOne = K * (J - 1);
		// start values for beta can be determined by linear regression of log(pij/piJ) on design matrix x
//		System.out.println("Nr betas: " + kJMinusOne);
		double[] beta = new double[kJMinusOne];

		double[] beta0 = new double[kJMinusOne];
		//double[] beta_inf = new double[kJMinusOne];
		double[][] xtwx = new double[kJMinusOne][kJMinusOne];


		double[] loglike = new double[1];
		double[] deviance = new double[1];
		double loglike0 = 0;

		while (iter < max_iter && !converged) {

			for (k = 0; k < kJMinusOne; k++) {
				beta0[k] = beta[k];
			}

			// ! nr function needs error handling
			nr(X, Y, n, J, N, K, beta0, beta, xtwx, loglike, deviance);

			if (loglike[0] < loglike0 && iter > 0) {
				// backtracking code
				// run subiterations to halve the distance to prior iteration
				// until difference in log_like increases (which is when the model has converged)
				// remember: ml is about maximizing loglike
			}

//			// test for infinity of beta
//			for (k = 0; k < kJMinusOne; k++) {
//				if (beta_inf[k] != 0) {
//					beta[k] = beta_inf[k];
//				} else {
//					//  Math.sqrt(xtwx[k][k]) contains the variance of beta[k]
//					double betak = beta[k];
//					double absbeta = Math.abs(betak);
//					if (absbeta > (5d / xrange[k]) && Math.sqrt(xtwx[k][k]) >= (3 * absbeta)) {
//						beta_inf[k] = betak;
//					}
//				}
//			}

			// test for convergence
			converged = true;
			for (k = 0; k < kJMinusOne; k++) {
				if (Math.abs(beta[k] - beta0[k]) > EPSILON * Math.abs(beta0[k])) {
					converged = false;
					break;
				}
			}

			if (iter == 0) {
				loglike0 = loglike[0];
			}
			iter++;
		}


		double[] sigprms = new double[kJMinusOne];
		double[] stderrs = new double[kJMinusOne];
		double[] wald = new double[kJMinusOne];
		AssociationResult result = new AssociationResult();
		if (converged) {
			// chi2 tests for significance
//			double chi1 = 2 * (loglike[0] - loglike0);
//			int df1 = (K * (J - 1)) - J - 1;
//			double chiTest = 1 - ChiSquare.getP(df1, chi1);
//
//			double chi2 = deviance[0];
//			int df2 = (N * (J - 1)) - (K * (J - 1));
//			double chiTest2 = 1 - ChiSquare.getP(df2, chi2);


			for (i = 0; i < kJMinusOne; i++) {
				if (xtwx[i][i] > 0) {
					stderrs[i] = Math.sqrt(xtwx[i][i]);
					wald[i] = Math.pow(beta[i] / stderrs[i], 2);
					sigprms[i] = 1 - ChiSquare.getP(1, wald[i]);
				} else {
					sigprms[i] = -1;
				}
			}
			// VariantID	N	MAF	DevianceNull	DfNull	DevianceGeno	DfAlt	Beta(Genotype)	SE(Genotype)	OR	OR-Hi	OR-Lo	Pval	#NAME?
			return null;
			//return new LogisticRegressionResult(beta, stderrs, sigprms, deviance[0], loglike, loglike0);
		} else {
			return null;
		}
	}


	private int nr(
			double[][] X,
			double[][] Y,
			double[] n,
			int J,
			int N,
			int K,
			double[] beta0,
			double[] beta1,
			double[][] xtwx,
			double[] loglike,
			double[] deviance

	) {


		int kJMinusOne = (K * (J - 1));
		double[][] pi = new double[N][J];
		double[] g = new double[kJMinusOne];
		double[][] H = new double[kJMinusOne][kJMinusOne];
		double[] numer = new double[J];

		double denom, q1, w1, w2, sum1;
		double devtmp;

		int i, j, k, jj, kk, kprime, jprime;

		loglike[0] = 0;
		deviance[0] = 0;

		int jMinusOne = J - 1;
		for (i = 0; i < N; i++) {
			denom = 1d;
			jj = 0;
			for (j = 0; j < jMinusOne; j++) {
				sum1 = 0;
				for (k = 0; k < K; k++) {
					sum1 += X[i][k] * beta0[jj++];
				}
				numer[j] = Math.exp(sum1);
				denom += numer[j];

			}

			/* calculate predicted probabilities */
			for (j = 0; j < J - 1; j++) {
				pi[i][j] = numer[j] / denom;
			}

        /* omitted category */
			pi[i][j] = 1.0d / denom;

        /* increment log likelihood */
			loglike[0] += Gamma.logGamma(n[i] + 1);
			for (j = 0; j < J; j++) {
				double yij = Y[i][j];
				loglike[0] = loglike[0] - Gamma.logGamma(yij + 1) + (yij * Math.log(pi[i][j]));
			}

			/* increment deviance */
			for (j = 0; j < J; j++) {
				double yij = Y[i][j];
				if (yij > 0) {
					devtmp = 2 * yij * Math.log(yij / (n[i] * pi[i][j]));
				} else {
					devtmp = 0;
				}
				deviance[0] += devtmp;
			}

			double ni = n[i];
			/* increment first and second derivatives */
			for (j = 0, jj = 0; j < J - 1; j++) {

				double yij = Y[i][j];
				double pij = pi[i][j];
			/* terms for first derivative, see Eq. 32 */
				q1 = yij - (ni * pij);

            /* terms for second derivative, see Eq. 37 */
				w1 = ni * pij * (1 - pij);

				for (k = 0; k < K; k++) {

					double xik = X[i][k];
				/* first derivative term in Eq. 23 */
					g[jj] += q1 * xik;

                /* increment the current pop's contribution to the 2nd derivative */

                /* jprime = j (see Eq. 37) */

					kk = jj - 1;
					for (kprime = k; kprime < K; kprime++) {
						kk += 1;
						H[jj][kk] += w1 * xik * X[i][kprime];
						H[kk][jj] = H[jj][kk];
					}

                /* jprime != j (see Eq. 37) */

					for (jprime = j + 1; jprime < J - 1; jprime++) {
						w2 = -ni * pij * pi[i][jprime];
						for (kprime = 0; kprime < K; kprime++) {
							kk += 1;
							H[jj][kk] += w2 * xik * X[i][kprime];
							H[kk][jj] = H[jj][kk];
						}
					}
					jj++;
				}
			}

		}

		/* compute xtwx * beta0 + x(y-mu) (see Eq. 40) */
		for (i = 0; i < kJMinusOne; i++) {
			sum1 = 0;
			for (j = 0; j < kJMinusOne; j++) {
				sum1 += H[i][j] * beta0[j];
			}
			g[i] += sum1;
		}

		/* invert xtwx */
		if (cholesky(H, kJMinusOne)) return -1;
		if (backsub(H, kJMinusOne)) return -2;
		if (trimult(H, xtwx, kJMinusOne)) return -3;

		/* solve for new betas */
		for (i = 0; i < kJMinusOne; i++) {
			sum1 = 0;
			for (j = 0; j < kJMinusOne; j++) {
				sum1 += xtwx[i][j] * g[j];
			}
			beta1[i] = sum1;
		}

		return 0;
	}

	private boolean trimult(double[][] in, double[][] out, int order) {
		int i, j, k, m;
		double sum;

		for (i = 0; i < order; i++) {
			for (j = 0; j < order; j++) {
				sum = 0;
				if (i > j)
					m = i;
				else
					m = j;
				for (k = m; k < order; k++)
					sum += in[i][k] * in[j][k];
				out[i][j] = sum;
			}
		}
		return false;
	}

	private boolean backsub(double[][] x, int order) {
		int i, j, k;
		double sum;

		if (x[0][0] == 0) return true; // throw new ArithmeticException("Problem with backsubstitution: x[0][0] == 0d");

		x[0][0] = 1 / x[0][0];
		for (i = 1; i < order; i++) {
			if (x[i][i] == 0)
				return true; // throw new ArithmeticException("Problem with backsubstitution: x[i][i] == 0d");
			x[i][i] = 1 / x[i][i];
			for (j = 0; j < i; j++) {
				sum = 0;
				for (k = j; k < i; k++)
					sum += x[j][k] * x[k][i];
				x[j][i] = -sum * x[i][i];
			}
		}
		return false;
	}

	private boolean cholesky(double[][] x, int order) {
		int i, j, k;
		double sum;

		for (i = 0; i < order; i++) {
			sum = 0;
			for (j = 0; j < i; j++) {
				sum += x[j][i] * x[j][i];
			}
			if (sum >= x[i][i]) {
				return true; // throw new ArithmeticException("Problem with Cholesky operation: sum>x[i][i]");
			}
			x[i][i] = Math.sqrt(x[i][i] - sum);
			for (j = i + 1; j < order; j++) {
				sum = 0;
				for (k = 0; k < i; k++) {
					sum += x[k][i] * x[k][j];
				}
				x[i][j] = (x[i][j] - sum) / x[i][i];
			}
		}

		return false;
	}


}
