package nl.harmjanwestra.utilities.bamfile;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.AggregateFilter;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.Strand;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

/**
 * Created by hwestra on 6/23/15.
 */
public class Coverage {

	private AggregateFilter filter = null;
	private int baseqQualThreshold = 20;
	private HashMap<String, Integer> readGroupMap = null;
	private int[][][] coverageStranded;

	public void setReadGroupMap(HashMap<String, Integer> readGroupMap) {
		this.readGroupMap = readGroupMap;
	}

	public void setBaseqQualThreshold(int b) {
		this.baseqQualThreshold = b;
	}

	public void setFilter(AggregateFilter f) {
		this.filter = f;
	}

	public void calculate(Iterator<SAMRecord> it, Feature chromosomeWindow) {
		calculate(it, chromosomeWindow, 0, 0);
	}

	public void calculate(Iterator<SAMRecord> it, Feature chromosomeWindow, int posShift, int negShift) {
		int chromosomeWindowStart = chromosomeWindow.getStart();
		int chromosomeWindowEnd = chromosomeWindow.getStop();
		int chromosomeWindowSize = chromosomeWindowEnd - chromosomeWindowStart;

		coverageStranded = null;
		if (readGroupMap == null) {
			coverageStranded = new int[2][1][chromosomeWindowSize]; // [strand][ind][bp]
		} else {
			coverageStranded = new int[2][readGroupMap.size()][chromosomeWindowSize];
		}


		while (it.hasNext()) {
			SAMRecord record = it.next();

			Feature f = new Feature();
			f.setChromosome(Chromosome.parseChr(record.getReferenceName()));
			f.setStart(record.getAlignmentStart());
			f.setStop(record.getAlignmentEnd());
			if (filter == null || !filter.filterOut(record)) {

				if (chromosomeWindow.overlaps(f)) {

					String refname = record.getReferenceName();
//					int len = record.getReadLength();
//					int sta = record.getAlignmentStart();
//					int sto = record.getAlignmentEnd();
//					int insert = record.getInferredInsertSize();


					byte[] baseQual = record.getBaseQualities();
					byte[] bases = record.getReadBases();

					Cigar cigar = record.getCigar();
					List<CigarElement> cigarElements = cigar.getCigarElements();

					if (bases == null || bases.length == 0 || baseQual == null || baseQual.length == 0) {
						System.err.println("No bases for read: " + record.getReadUnmappedFlag() + "\t" + record.toString() + "\t" + bases.length + "\t" + baseQual.length);
					} else {


						int shift = posShift;
						Strand strand = Strand.POS;
						if (record.getReadNegativeStrandFlag()) {
							strand = Strand.NEG;
							shift = negShift;
						}

						int mapPos = record.getAlignmentStart();
						int windowRelativePosition = mapPos - chromosomeWindowStart + shift;
						int readPosition = 0;

						Integer sampleId = 0;
						int[] tmpCoverage = null;
						if (readGroupMap != null) {
							SAMReadGroupRecord rg = record.getReadGroup();
							sampleId = readGroupMap.get(rg);
						}

						tmpCoverage = coverageStranded[strand.getNumber()][sampleId];


						for (CigarElement e : cigarElements) {
							int cigarElementLength = e.getLength();

							switch (e.getOperator()) {
								case H:
									break;
								case P:
									break;
								case S: // soft clip
									readPosition += cigarElementLength;
									break;
								case N: // ref skip
									windowRelativePosition += cigarElementLength;
									break;
								case D: // deletion
									windowRelativePosition += cigarElementLength;
									break;
								case I: // insertion
									windowRelativePosition += cigarElementLength;
									break;
								case M:
								case EQ:
								case X:
									int endPosition = readPosition + cigarElementLength;
									for (int pos = readPosition; pos < endPosition; pos++) {
										boolean properbase = false;
										byte base = bases[pos];

										if (windowRelativePosition >= 0 && windowRelativePosition < chromosomeWindowSize) { // the read could overlap the leftmost edge of this window
											if (base == 78 || base == 110) {
												// N

											} else if (readPosition < baseQual.length && baseQual[readPosition] > baseqQualThreshold) { //    -- for each base pos: check whether basequal > 50
												//    -- determine number of A/T/C/G/N bases
												if (base == 65 || base == 97) {

													properbase = true;
												} else if (base == 67 || base == 99) {

													properbase = true;
												} else if (base == 71 || base == 103) {

													properbase = true;
												} else if (base == 84 || base == 116) { // extend to capture U?

													properbase = true;
												}
											} else if (baseQual[readPosition] > 0) {
												System.err.println("Unknown base found! " + base);
											}
										}

										if (properbase) {
											tmpCoverage[windowRelativePosition]++;

										}
										windowRelativePosition++;
									} // if pos < readposition
									readPosition += cigarElementLength;
									break;
								default:
									System.err.println("Unknown CIGAR operator found: " + e.getOperator().toString());
									System.err.println("In read: " + record.toString());
									break;
							} // switch operator
						} // for each cigar element
						coverageStranded[strand.getNumber()][sampleId] = tmpCoverage;
					}

				}

			}
		}
	}


	public int[][][] getCoverageStranded() {
		return coverageStranded;
	}


	public int[][] getCoverageUnStranded() {
		// rewrite to unstranded
		int[][] unstranded = new int[coverageStranded[0].length][coverageStranded[0][0].length]; //[ind][bp]


		/// [strand][ind][bp]
		for (int strand = 0; strand < 2; strand++) {
			for (int i = 0; i < unstranded.length; i++) {
				for (int j = 0; j < unstranded[0].length; j++) {
					unstranded[i][j] += coverageStranded[strand][i][j];
				}
			}
		}

		return unstranded;
	}


	public int[] getCoverageCombineReadroups() {
		int[] rgMerged = new int[coverageStranded[0][0].length]; // [bp]
		/// [strand][ind][bp]
		for (int strand = 0; strand < 2; strand++) {
			for (int i = 0; i < coverageStranded[strand].length; i++) {
				for (int j = 0; j < coverageStranded[strand][i].length; j++) {
					rgMerged[i] += coverageStranded[strand][i][j];
				}
			}
		}
		return rgMerged;
	}


}
