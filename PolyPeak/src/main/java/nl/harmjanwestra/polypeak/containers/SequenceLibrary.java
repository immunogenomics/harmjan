package nl.harmjanwestra.polypeak.containers;


import java.util.HashMap;

/**
 * Created by hwestra on 1/7/15.
 */
public class SequenceLibrary {
	private static final HashMap<String, Sequence> sequenceHash = new HashMap<String, Sequence>();

	public static boolean checkSequence(String sequenceName, int sequenceLength) {
		Sequence s = new Sequence(sequenceName, sequenceLength);
		if (sequenceHash.containsKey(s.getName())) {
			// this only checks the name
			Sequence other = sequenceHash.get(s.getName());
			return s.compareLength(other);

		} else {
			// add to library
			sequenceHash.put(s.getName(), s);
			System.out.println("Adding sequence: " + sequenceName + " as " + s.getName() + " with length " + sequenceLength);
			return true;
		}
	}

	public static Sequence getSequence(Sequence s) {
		return sequenceHash.get(s.getName());
	}

	public static void setSequence(Sequence sequence) {
		sequenceHash.put(sequence.getName(), sequence);
	}

	public static HashMap<String, Sequence> getSequences() {
		return sequenceHash;
	}
}
