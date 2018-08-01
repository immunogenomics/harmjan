package nl.harmjanwestra.vcfutils;

/**
 * Created by Harm-Jan on 12/24/15.
 */
public class CorrelationResult {


	double maf1;
	double maf2;
	double rsqPearson;
	double rsqBeagle1;
	double rsqBeagle2;
	int nrSamples;
	String variant;

	public double get(TYPE t) {
		switch (t) {
			case maf1:
				return maf1;
			case maf2:
				return maf2;
			case rsqb1:
				return rsqBeagle1;
			case rsqb2:
				return rsqBeagle2;
			case rsqp:
				return rsqPearson;
			case nrsa:
				return nrSamples;
		}
		return 0;
	}

	enum TYPE {
		maf1,
		maf2,
		rsqp,
		rsqb1,
		rsqb2,
		nrsa
	}
}
