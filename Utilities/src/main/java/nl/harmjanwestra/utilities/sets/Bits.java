package nl.harmjanwestra.utilities.sets;

import java.util.BitSet;

/**
 * Created by hwestra on 11/3/16.
 */
public class Bits {

	public static BitSet convert(long value) {
		BitSet bits = new BitSet();
		int index = 0;
		while (value != 0L) {
			if (value % 2L != 0) {
				bits.set(index);
			}
			++index;
			value = value >>> 1;
		}
		return bits;
	}

	public static long convert(BitSet bits) {
		long value = 0L;
		for (int i = 0; i < bits.length(); ++i) {
			value += bits.get(i) ? (1L << i) : 0L;
		}
		return value;
	}

	public static int[] valuesForBitCombinations(int numBits) {
//		int val = 0;
		int[] values = new int[]{0, 1};
		values[0] = 0;
		values[1] = 1;

		for (int i = 1; i < numBits; i++) {
			int[] moreValues = new int[values.length * 2];
//			int start = (int) Math.pow(2, i);
			for (int j = 0; j < values.length; j++) {
				moreValues[j * 2] = values[j] << 1;
				moreValues[j * 2 + 1] = values[j] << 1 | 1;
			}
			values = moreValues;
		}
		return values;
	}

	public static boolean[] convertAsBoolean(long value, int numBits) {
		BitSet bits = convert(value);
		boolean[] setB = new boolean[numBits];
		for (int q = 0; q < numBits; q++) {
			setB[(numBits - 1) - q] = bits.get(q);
		}
		return setB;
	}
}