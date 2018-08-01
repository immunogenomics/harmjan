package nl.harmjanwestra.utilities.sets;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 11/30/15.
 */
public class StringSets {

	public static HashMap<String, Integer> index(String[] s) {
		int c = 0;
		HashMap<String, Integer> index = new HashMap<String, Integer>();
		for (String s1 : s) {
			index.put(s1, c);
			c++;
		}
		return index;
	}

	public static HashSet<String> intersect(ArrayList<String> s1, ArrayList<String> s2) {
		HashSet<String> intersect = new HashSet<String>();
		HashSet<String> set1 = new HashSet<String>();
		set1.addAll(s1);
		for (String s : s2) {
			if (set1.contains(s)) {
				intersect.add(s);
			}
		}
		return intersect;
	}


}
