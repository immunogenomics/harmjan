package nl.harmjanwestra.utilities.math;

/**
 * Created by hwestra on 10/26/15.
 * Derived from Scott. A. Czepiel's implementation (http://czep.net/)
 */

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.jet.stat.Gamma;
import umcg.genetica.math.stats.ChiSquare;

public class LogisticRegressionOptimized {

	int max_iter = 50;
	double EPSILON = 1E-3;

	public LogisticRegressionOptimized(){

	}

	public LogisticRegressionOptimized(int maxiter) {
		max_iter = maxiter;
	}

	public LogisticRegressionOptimized(int maxiter, double epsilon) {
		max_iter = maxiter;
		EPSILON = epsilon;
	}

	public LogisticRegressionResult univariate(DoubleMatrix2D y, DoubleMatrix2D x) {
		if (x.rows() != y.rows()) {
			throw new IllegalArgumentException("Error in logistic regression: length of y does not match length of x");
		}
		int J = 2;
		int N = x.rows();
		int K = x.columns();
		double[] n = new double[y.rows()];
		for (int i = 0; i < y.rows(); i++) {
			n[i] = 1;
		}
		return mlelr(J, N, K, n, x, y);
	}

	public LogisticRegressionResult univariate(double[][] y, double[][] x) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("Error in logistic regression: length of y does not match length of x");
		}
		DoubleMatrix2D xtmp = new DenseDoubleMatrix2D(x);
		DoubleMatrix2D ytmp = new DenseDoubleMatrix2D(y);

		return univariate(ytmp, xtmp);
	}

	public LogisticRegressionResult univariate(double[] y, DoubleMatrix2D x) {
		if (x.rows() != y.length) {
			throw new IllegalArgumentException("Error in logistic regression: length of y does not match length of x");
		}
		DoubleMatrix2D yreg = new DenseDoubleMatrix2D(x.rows(), 2);
		for (int i = 0; i < y.length; i++) {
			if (y[i] > 1 || y[i] < 0) {
				throw new IllegalArgumentException("Error in univariate logistic regression; unexpected categorical value: " + y[i] + ". Expected 0 or 1");
			}
			yreg.setQuick(i, (int) y[i], 1);
		}

		return univariate(yreg, x);
	}

	public LogisticRegressionResult univariate(double[] y, double[][] x) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("Error in logistic regression: length of y does not match length of x");
		}
		DoubleMatrix2D xtmp = new DenseDoubleMatrix2D(x);
		DoubleMatrix2D yreg = new DenseDoubleMatrix2D(x.length, 2);

		for (int i = 0; i < y.length; i++) {
			if (y[i] > 1 || y[i] < 0) {
				throw new IllegalArgumentException("Error in univariate logistic regression; unexpected categorical value: " + y[i] + ". Expected 0 or 1");
			}
			yreg.setQuick(i, (int) y[i], 1);

		}
		return univariate(yreg, xtmp);
	}

	public LogisticRegressionResult mlelr(int J,
										  int N,
										  int K,
										  double[] n,
										  DoubleMatrix2D X,
										  DoubleMatrix2D Y) {

		int i, j, k;

		int iter = 0;
		boolean converged = false;

		int kJMinusOne = K * (J - 1);
		// start values for beta can be determined by linear regression of log(pij/piJ) on design matrix x
//		System.out.println("Nr betas: " + kJMinusOne);
		double[] beta = new double[kJMinusOne];

		double[] beta0 = new double[kJMinusOne];
//		double[] diff = new double[kJMinusOne];
//		boolean[] diffb = new boolean[kJMinusOne];
		//double[] beta_inf = new double[kJMinusOne];
		DoubleMatrix2D xtwx = new DenseDoubleMatrix2D(kJMinusOne, kJMinusOne);


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
				double beta0k = beta0[k];
//				diff[k] = Math.abs(beta[k] - beta0k);
//				diffb[k] = true;
				if (Math.abs(beta[k] - beta0k) > EPSILON * Math.abs(beta0k)) {
					converged = false;
//					diffb[k] = false;
					break;
				}
			}

			if (iter == 0) {
				loglike0 = loglike[0];
			}
			iter++;
		}

//		for (int q = 0; q < diffb.length; q++) {
//
//			double[] bliep = new double[X.rows()];
//			for (int row = 0; row < bliep.length; row++) {
//				bliep[row] = X.getQuick(row, q);
//			}
//
//			System.out.println(Descriptives.variance(bliep) + "\t" + diffb[q] + "\t" + diff[q]);
//
//
//		}

		double[] sigprms = new double[kJMinusOne];
		double[] stderrs = new double[kJMinusOne];
		double[] wald = new double[kJMinusOne];
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
				double xtwxii = xtwx.get(i, i);
				if (xtwxii > 0) {
					stderrs[i] = Math.sqrt(xtwxii);
					wald[i] = Math.pow(beta[i] / stderrs[i], 2);
					sigprms[i] = 1 - ChiSquare.getP(1, wald[i]);
				} else {
					sigprms[i] = -1;
				}
			}
			// VariantID	N	MAF	DevianceNull	DfNull	DevianceGeno	DfAlt	Beta(Genotype)	SE(Genotype)	OR	OR-Hi	OR-Lo	Pval	#NAME?
			return new LogisticRegressionResult(beta, stderrs, sigprms, deviance[0], loglike, loglike0);
		} else {
			return null;
		}
	}

	private double[] g;
	private double[] numer;
	private DoubleMatrix2D pi;
	private DoubleMatrix2D H;

	private int nr(
			DoubleMatrix2D X,
			DoubleMatrix2D Y,
			double[] n,
			int J,
			int N,
			int K,
			double[] beta0,
			double[] beta1,
			DoubleMatrix2D xtwx,
			double[] loglike,
			double[] deviance
	) {
		int kJMinusOne = (K * (J - 1));

		// if this is the first iteration, initialize (and save some GC time in the following iterations)
		if (pi != null) {
			// check the size
			if (pi.rows() != N || pi.columns() != J) {
				pi = new DenseDoubleMatrix2D(N, J);
			}
			if (H.rows() != kJMinusOne || H.columns() != kJMinusOne) {
				H = new DenseDoubleMatrix2D(kJMinusOne, kJMinusOne);
			}
			if (g.length != kJMinusOne) {
				g = new double[kJMinusOne];
			}
			if (numer.length != J) {
				numer = new double[J];
			}
		}

		if (pi == null) {
			pi = new DenseDoubleMatrix2D(N, J);
			H = new DenseDoubleMatrix2D(kJMinusOne, kJMinusOne);
			g = new double[kJMinusOne];
			numer = new double[J];
		} else {
			// set everything to null
			for (int i = 0; i < kJMinusOne; i++) {
				for (int j = 0; j < kJMinusOne; j++) {
					pi.setQuick(i, j, 0);
					H.setQuick(i, j, 0);
				}
				g[i] = 0;

			}
			for (int i = 0; i < numer.length; i++) {
				numer[i] = 0d;
			}
		}


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
					sum1 += X.get(i, k) * beta0[jj++];
				}
				numer[j] = Math.exp(sum1);
				denom += numer[j];

			}

			/* calculate predicted probabilities */
			for (j = 0; j < J - 1; j++) {
				pi.setQuick(i, j, numer[j] / denom);
			}

        /* omitted category */
			pi.setQuick(i, j, 1.0d / denom);

        /* increment log likelihood */
			loglike[0] += Gamma.logGamma(n[i] + 1);
			for (j = 0; j < J; j++) {
				double pij = pi.getQuick(i, j);
				double yij = Y.getQuick(i, j);
				loglike[0] = loglike[0] - Gamma.logGamma(yij + 1) + (yij * Math.log(pij));
//			}
//
//			/* increment deviance */
//			for (j = 0; j < J; j++) {
//				double yij = Y.get(i, j);
				if (yij > 0) {
					devtmp = 2 * yij * Math.log(yij / (n[i] * pij));
				} else {
					devtmp = 0;
				}
				deviance[0] += devtmp;
			}

			double ni = n[i];
			/* increment first and second derivatives */
			for (j = 0, jj = 0; j < J - 1; j++) {

				double yij = Y.getQuick(i, j);
				double pij = pi.getQuick(i, j);
			/* terms for first derivative, see Eq. 32 */
				q1 = yij - (ni * pij);

            /* terms for second derivative, see Eq. 37 */
				w1 = ni * pij * (1 - pij);

				for (k = 0; k < K; k++) {

					double xik = X.getQuick(i, k);
				/* first derivative term in Eq. 23 */
					g[jj] += q1 * xik;

                /* increment the current pop's contribution to the 2nd derivative */

                /* jprime = j (see Eq. 37) */
					kk = jj - 1;
					for (kprime = k; kprime < K; kprime++) {
						kk += 1;
						double hjjkk = H.getQuick(jj, kk);
						hjjkk += w1 * xik * X.getQuick(i, kprime);
						H.setQuick(jj, kk, hjjkk);
						H.setQuick(kk, jj, hjjkk);
					}

                /* jprime != j (see Eq. 37) */

					for (jprime = j + 1; jprime < J - 1; jprime++) {
						w2 = -ni * pij * pi.getQuick(i, jprime);
						for (kprime = 0; kprime < K; kprime++) {
							kk += 1;
							double hjjkk = H.getQuick(jj, kk);
							hjjkk += w2 * xik * X.getQuick(i, kprime);
							H.setQuick(jj, kk, hjjkk);
							H.setQuick(kk, jj, hjjkk);
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
				sum1 += H.getQuick(i, j) * beta0[j];
			}
			g[i] += sum1;
		}

		/* invert xtwx */
		if (cholesky(H)) return -1;
		if (backsub(H)) return -2;
		if (trimult(H, xtwx)) return -3;

		/* solve for new betas */
		for (i = 0; i < kJMinusOne; i++) {
			sum1 = 0;
			for (j = 0; j < kJMinusOne; j++) {
				sum1 += xtwx.getQuick(i, j) * g[j];
			}
			beta1[i] = sum1;
		}

		return 0;
	}

	private boolean trimult(DoubleMatrix2D in, DoubleMatrix2D out) {
		int i, j, k, m;
		double sum;
		int order = in.rows();

		for (i = 0; i < order; i++) {
			for (j = 0; j < order; j++) {
				sum = 0;
				if (i > j) {
					m = i;
				} else {
					m = j;
				}
				for (k = m; k < order; k++) {
					sum += in.getQuick(i, k) * in.getQuick(j, k);
				}
				out.setQuick(i, j, sum);
			}
		}
		return false;
	}

	private boolean backsub(DoubleMatrix2D x) {
		int i, j, k;
		double sum;
		int order = x.rows();

		if (x.get(0, 0) == 0) {
			return true; // throw new ArithmeticException("Problem with backsubstitution: x[0][0] == 0d");
		}

		x.setQuick(0, 0, (1d / x.getQuick(0, 0)));
		for (i = 1; i < order; i++) {
			if (x.getQuick(i, i) == 0) {
				return true; // throw new ArithmeticException("Problem with backsubstitution: x[i][i] == 0d");
			}
			x.setQuick(i, i, (1d / x.getQuick(i, i)));
			for (j = 0; j < i; j++) {
				sum = 0;
				for (k = j; k < i; k++) {
					sum += x.getQuick(j, k) * x.getQuick(k, i);
				}
				x.setQuick(j, i, (-sum * x.getQuick(i, i)));
			}
		}
		return false;
	}

	private boolean cholesky(DoubleMatrix2D x) {
		int i, j, k;
		double sum;
		int order = x.rows();

		for (i = 0; i < order; i++) {
			sum = 0;
			double xii = x.getQuick(i, i);
			for (j = 0; j < i; j++) {
				double xji = x.getQuick(j, i);
				sum += xji * xji;
			}
			if (sum >= xii) {
				return true; // throw new ArithmeticException("Problem with Cholesky operation: sum>x[i][i]");
			}
			x.setQuick(i, i, Math.sqrt(xii - sum));
			for (j = i + 1; j < order; j++) {
				sum = 0;
				for (k = 0; k < i; k++) {
					sum += x.getQuick(k, i) * x.getQuick(k, j);
				}
				x.setQuick(i, j, (x.getQuick(i, j) - sum) / x.getQuick(i, i));
			}
		}
		return false;
	}


}
