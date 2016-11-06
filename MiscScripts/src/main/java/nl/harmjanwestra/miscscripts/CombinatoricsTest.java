package nl.harmjanwestra.miscscripts;

import org.apache.commons.math3.util.CombinatoricsUtils;
import umcg.genetica.text.Strings;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Created by Harm-Jan on 06/23/16.
 */
public class CombinatoricsTest {


	public static void main(String[] args) {

		Iterator<int[]> combos = CombinatoricsUtils.combinationsIterator(2, 2);

		while (combos.hasNext()) {
			int[] combo = combos.next();
			System.out.println(Strings.concat(combo, Strings.tab));
		}

		int numBits = 5;

		int val = 0;
		int[] values = new int[]{0, 1};
		values[0] = 0;
		values[1] = 1;

		for (int i = 1; i < numBits; i++) {
			int[] moreValues = new int[values.length * 2];
			int start = (int) Math.pow(2, i);
			for (int j = 0; j < values.length; j++) {
				moreValues[j * 2] = values[j] << 1;
				moreValues[j * 2 + 1] = values[j] << 1 | 1;
			}
			values = moreValues;
		}

		//print the values
		for (int value : values) {

//			BitSet set = Bits.convert(value);


			System.out.println(Integer.toBinaryString(value));
		}


//		CombinatoricsTest t = new CombinatoricsTest();
//		t.test();
	}


	public void test() {
		String[] letters = new String[]{"B", "C", "D", "E"};

		Combinations c = new Combinations(letters);
		c.combine();
		ArrayList<List<String>> results = c.getResults();
		for (List<String> s : results) {
			System.out.println(Strings.concat(s, Strings.tab));
		}
	}

	public class Combinations {
		private final String[] inputstring;
		ArrayList<List<String>> results = new ArrayList<List<String>>();
		private List<String> output = new ArrayList<String>();

		public Combinations(final String[] str) {
			inputstring = str;
			System.out.println("The input string  is  : " + Strings.concat(inputstring, Strings.tab));
		}

		public void combine() {
			combine(0);
		}

		private void combine(int start) {
			for (int i = start; i < inputstring.length; ++i) {
				output.add(inputstring[i]);
				List<String> tmpOutput = new ArrayList<String>();
				tmpOutput.addAll(output);
				results.add(tmpOutput);
				if (i < inputstring.length) {
					combine(i + 1);
				}
				output = output.subList(0, output.size() - 1);
			}
		}

		public ArrayList<List<String>> getResults() {
			return results;
		}
	}
}

