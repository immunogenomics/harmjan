package nl.harmjanwestra.miscscripts;

import umcg.genetica.text.Strings;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by Harm-Jan on 06/23/16.
 */
public class CombinatoricsTest {


	public static void main(String[] args) {

		CombinatoricsTest t = new CombinatoricsTest();
		t.test();
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
