package nl.harmjanwestra.finemappingtools.gwas;

import cern.colt.matrix.tbit.BitVector;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by hwestra on 11/2/16.
 */
public class HaplotypeCombiner {
	private final ArrayList<BitVector> inputHaplotypes;
	ArrayList<List<BitVector>> results = new ArrayList<List<BitVector>>();
	private List<BitVector> output = new ArrayList<BitVector>();

	public HaplotypeCombiner(final ArrayList<BitVector> haps) {
		inputHaplotypes = haps;
	}

	public void combine() {
		combine(0);
	}

	private void combine(int start) {
		for (int i = start; i < inputHaplotypes.size(); ++i) {
			output.add(inputHaplotypes.get(i));
			List<BitVector> tmpOutput = new ArrayList<BitVector>();
			tmpOutput.addAll(output);
			results.add(tmpOutput);
			if (i < inputHaplotypes.size()) {
				combine(i + 1);
			}
			output = output.subList(0, output.size() - 1);
		}
	}

	public ArrayList<List<BitVector>> getResults() {
		return results;
	}
}