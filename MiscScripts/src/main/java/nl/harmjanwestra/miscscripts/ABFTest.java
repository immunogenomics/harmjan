package nl.harmjanwestra.miscscripts;

import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * Created by hwestra on 11/1/15.
 */
public class ABFTest {
	public static void main(String[] args) {
		double s1 = 0.1 * 0.1; // s.e.
		double s2 = (0.1 * 0.1) + (0.35 * 0.35); // s.e. + prior
		double i = (Math.log(2) / 1.96) * (Math.log(2) / 1.96);
		NormalDistribution d1 = new NormalDistribution(0, s1);
		JSci.maths.statistics.NormalDistribution dist = new JSci.maths.statistics.NormalDistribution(0, s1);
		JSci.maths.statistics.NormalDistribution dist2 = new JSci.maths.statistics.NormalDistribution(0, s2);
		NormalDistribution d2 = new NormalDistribution(0, s2);
		double theta = Math.log(1.316);
//		theta = 1.316;
		double p0 = d1.probability(theta);
		double p1 = d2.probability(theta);

		double p2 = dist.probability(theta);
		double p3 = dist2.probability(theta);

		System.out.println(p0);
		System.out.println(p1);
		System.out.println((p0 / p1));
		System.out.println();
		System.out.println(p2);
		System.out.println(p3);
		System.out.println((p2 / p3));
		System.out.println();
		System.out.println(i);
		System.out.println((0.35 * 0.35));
		System.out.println(theta);
		i = (Math.log(2) / 1.5) * (Math.log(2) / 1.5);
		System.out.println(Math.log(1.5));

	}
}
