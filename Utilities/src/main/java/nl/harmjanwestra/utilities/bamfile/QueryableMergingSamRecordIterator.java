package nl.harmjanwestra.utilities.bamfile;

import htsjdk.samtools.SAMRecord;
import nl.harmjanwestra.utilities.features.Feature;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * Created by hwestra on 7/14/15.
 */
public class QueryableMergingSamRecordIterator implements Iterator<SAMRecord> {

	ArrayList<BamFileReader> readers;
	Feature f;

	SAMRecord next = null;

	Iterator<SAMRecord> currentIterator = null;
	int currentFileIndex = 0;

	public QueryableMergingSamRecordIterator(ArrayList<BamFileReader> readers, Feature f) {
		this.f = f;
		this.readers = readers;
		next();
	}

	@Override
	public boolean hasNext() {
		return next != null;
	}

	@Override
	public SAMRecord next() {

		SAMRecord current = next;
		// TODO: get chromosome name from somewhere.
		if (currentIterator == null) {

			String chrName = readers.get(currentFileIndex).matchChromosomeName(f);
			if (chrName != null) {
				currentIterator = readers.get(currentFileIndex).query(chrName, f.getStart(), f.getStop(), false);
				System.out.println("path index: " + currentFileIndex);
			}
		}

		if (currentIterator.hasNext()) {
			next = currentIterator.next();
		} else if (currentFileIndex + 1 < readers.size()) {
			currentFileIndex++;
			String chrName = readers.get(currentFileIndex).matchChromosomeName(f);
			if (chrName != null) {
				currentIterator = readers.get(currentFileIndex).query(chrName, f.getStart(), f.getStop(), false);
			}
			System.out.println("path index: " + currentFileIndex);
		} else {
			next = null;
		}


		return current;
	}

	@Override
	public void remove() {

	}


}
