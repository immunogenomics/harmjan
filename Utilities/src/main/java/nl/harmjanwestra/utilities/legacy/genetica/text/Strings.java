/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.utilities.legacy.genetica.text;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author harmjan
 */
public class Strings {

	public static final Pattern tab = Pattern.compile("\t");
	public static final Pattern whitespace = Pattern.compile("\\s");
	public static final Pattern comma = Pattern.compile(",");

	public static final Pattern semicolon = Pattern.compile(";");
	public static final Pattern equalssign = Pattern.compile("=");
	public static final Pattern colon = Pattern.compile(":");
	public static final Pattern pipe = Pattern.compile("\\|");
	public static final Pattern forwardslash = Pattern.compile("/");
	public static final Pattern backwardslash = Pattern.compile("\\\\");
	public static final Pattern dot = Pattern.compile("\\.");
	public static final Pattern space = Pattern.compile(" ");
	public static final Pattern dash = Pattern.compile("-");

	public static String concat(String[] s, Pattern t) {

		return concat(s, t, null, null);
	}

	// this variant of Pattern.split uses primitive arrays instead of an array list
	// should be slower (since iterating the same string twice
	// but may prevent some unwanted out of memory messages when parsing lots of long strings
	// also, use the string cache when possible.
	public static String[] split(CharSequence input, int limit, Pattern p) {
		int index = 0;
		boolean matchLimited = limit > 0;

		Matcher m = p.matcher(input);

		// count number of matches
		int nrMatch = 0;
		while (m.find()) {
			if (!matchLimited || nrMatch < limit - 1) {
				if (index == 0 && index == m.start() && m.start() == m.end()) {
					// no empty leading substring included for zero-width match
					// at the beginning of the input char sequence.
					continue;
				}
//				String match = input.subSequence(index, m.start()).toString();
//				matchList.add(match);
				nrMatch++;
				index = m.end();
			} else if (nrMatch == limit - 1) { // last one
//				String match = input.subSequence(index,
//						input.length()).toString();
//				matchList.add(match);
				nrMatch++;
				index = m.end();
			}
		}

		String[] output = new String[nrMatch + 1];
		nrMatch = 0;
		index = 0;
		m.reset();

		// Add segments before each match found
		while (m.find()) {
			if (!matchLimited || nrMatch < limit - 1) {
				if (index == 0 && index == m.start() && m.start() == m.end()) {
					// no empty leading substring included for zero-width match
					// at the beginning of the input char sequence.
					continue;
				}
				String match = new String(input.subSequence(index, m.start()).toString()).intern();
				output[nrMatch] = match;
				nrMatch++;
				index = m.end();
			} else if (nrMatch == limit - 1) { // last one
				String match = new String(input.subSequence(index, input.length()).toString()).intern();
				output[nrMatch] = match;
				nrMatch++;
				index = m.end();
			}
		}

		// If no match was found, return this
		if (index == 0)
			return new String[]{input.toString()};

		// Add remaining segment
		if (!matchLimited || nrMatch < limit)
			output[nrMatch] = new String(input.subSequence(index, input.length()).toString()).intern();

		// Construct result
		int resultSize = output.length;
		if (limit == 0)
			while (resultSize > 0 && output[resultSize - 1].equals(""))
				resultSize--;
		String[] result = new String[resultSize];
		for (int i = 0; i < output.length && i < resultSize; i++) {
			result[i] = output[i];
		}
		return result;
	}


	public static String concat(String[] s, Pattern t, String replaceNull) {
		return concat(s, t, null, replaceNull);
	}

	public static String concat(String[] s, Pattern t, boolean[] includeElem, String replaceNull) {

		if (s == null) {
			return null;
		}

		int approximateFinalStrLen = 0;
		for (int i = 0; i < s.length; i++) {
			if (includeElem != null && includeElem[i]) {
				approximateFinalStrLen += s[i].length();
			}
		}
		approximateFinalStrLen += s.length;

		StringBuilder output = new StringBuilder(approximateFinalStrLen);
		for (int i = 0; i < s.length; i++) {
			if (includeElem == null || includeElem[i]) {
				if (s[i] == null) {
					s[i] = replaceNull;
				}
				if (i == 0) {
					output.append(s[i]);
				} else {
					output.append(t.toString()).append(s[i]);
				}
			}
		}
		return output.toString();
	}

	public static String concat(Object[] s, Pattern t) {

		StringBuilder output = new StringBuilder();
		for (int i = 0; i < s.length; i++) {
			if (i == 0) {
				output.append(s[i].toString());
			} else {
				output.append(t.toString()).append(s[i].toString());
			}
		}
		return output.toString();
	}


	public static String concat(double[] s, Pattern t) {
		return concat(s, t, null);
	}

	public static String concat(double[] s, Pattern t, DecimalFormat f) {

		StringBuilder output = new StringBuilder();
		for (int i = 0; i < s.length; i++) {
			if (i == 0) {
				if (f == null) {
					output.append(s[i]);
				} else {
					output.append(f.format(s[i]));
				}
			} else {
				if (f == null) {
					output.append(t.toString()).append(s[i]);
				} else {
					output.append(t.toString()).append(f.format(s[i]));
				}
			}
		}
		return output.toString();
	}

	public static String concat(double[] s, DecimalFormat f, Pattern t) {

		StringBuilder output = new StringBuilder();
		for (int i = 0; i < s.length; i++) {
			if (i == 0) {
				output.append(f.format(s[i]));
			} else {
				output.append(t.toString()).append(f.format(s[i]));
			}
		}
		return output.toString();
	}

	public static String concat(float[] s, DecimalFormat f, Pattern t) {

		StringBuilder output = new StringBuilder();
		for (int i = 0; i < s.length; i++) {
			if (i == 0) {
				output.append(f.format(s[i]));
			} else {
				output.append(t.toString()).append(f.format(s[i]));
			}
		}
		return output.toString();
	}

	public static String concat(int[] s, Pattern t) {

		StringBuilder output = new StringBuilder();
		for (int i = 0; i < s.length; i++) {
			if (i == 0) {
				output.append(s[i]);
			} else {
				output.append(t.toString()).append(s[i]);
			}
		}
		return output.toString();
	}

	public static String concat(List<String> s, Pattern t) {
		String[] data = s.toArray(new String[0]);
		return concat(data, t);
	}

	public static String concat(String[] s, Pattern t, int start, int end) {
		String[] data = new String[end - start];
		for (int i = start; i < end; i++) {
			data[i - start] = s[i];
		}
		return concat(data, t);
	}

	public static String[] split(String in) {
		List<String> list = new ArrayList<String>();
		StringBuilder sb = new StringBuilder();
		int len = in.length();
		int i = 0;
		char c;
		while (i < len) {
			c = in.charAt(i);
			if (c == '\t' || i == len) {
				list.add(sb.toString());
				sb.delete(0, len - 1);
			} else {
				sb.append(c);
			}
			i++;
		}


		return list.toArray(new String[0]);
	}

	public static String reverse(String str) {
		return new StringBuffer(str).reverse().toString();
	}

	public static String concat(String[] elems, boolean[] includeelem, Pattern tab) {

		return concat(elems, tab, includeelem, null);
	}

	public static String concat(double[] d, Pattern t, int start, int end) {
		String[] data = new String[end - start];
		for (int i = start; i < end; i++) {
			data[i - start] = "" + d[i];
		}
		return concat(data, t);
	}

	public static String concat(int[] d, Pattern t, int start, int end) {
		String[] data = new String[end - start];
		for (int i = start; i < end; i++) {
			data[i - start] = "" + d[i];
		}
		return concat(data, t);
	}

	public static String repeat(String toRepeat, int times) {

		if (toRepeat == null) {
			toRepeat = "";
		}

		if (times < 0) {
			times = 0;
		}

		final int length = toRepeat.length();
		final int total = length * times;
		final char[] src = toRepeat.toCharArray();
		char[] dst = new char[total];

		for (int i = 0; i < total; i += length) {
			System.arraycopy(src, 0, dst, i, length);
		}

		return String.copyValueOf(dst);

	}
}
