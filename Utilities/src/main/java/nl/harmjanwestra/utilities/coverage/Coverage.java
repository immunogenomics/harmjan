package nl.harmjanwestra.utilities.coverage;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.AggregateFilter;
import nl.harmjanwestra.utilities.bamfile.MinimalRead;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.enums.Strand;
import umcg.genetica.containers.Pair;

import java.util.*;

/**
 * Created by hwestra on 6/23/15.
 */
public class Coverage {

	private AggregateFilter filter = null;
	private int baseqQualThreshold = 10;
	private HashMap<String, Integer> readGroupMap = null;
	private int[][][] coverageStranded;
	private ArrayList<HashMap<String, MinimalRead>> readPositionsPerSample = null;
	private boolean countReads;
	private boolean countDuplicates;
	private long nrReads1 = 0;
	private long nrReads2 = 0;
	private long nrDups = 0;


	private boolean determineInsertSize;
	private int[] insertSize;
	private int[] fragmentsizeAfterFilter;
	private long nrReadsPassingFilter;
	private long nrDupsPassingFilter;
	private int maxInsertSize;
	private boolean countMultiMapping;
	private boolean countSupplementaryReads;
	private boolean countFragments;


	private long nrFragments;
	private long nrFragmentsDup;
	private long nrFragmentsPassingFilter;


	public AggregateFilter getFilter() {
		return filter;
	}

	public int getBaseqQualThreshold() {
		return baseqQualThreshold;
	}

	public HashMap<String, Integer> getReadGroupMap() {
		return readGroupMap;
	}

	public void setCoverageStranded(int[][][] coverageStranded) {
		this.coverageStranded = coverageStranded;
	}

	public ArrayList<HashMap<String, MinimalRead>> getReadPositionsPerSample() {
		return readPositionsPerSample;
	}

	public void setReadPositionsPerSample(ArrayList<HashMap<String, MinimalRead>> readPositionsPerSample) {
		this.readPositionsPerSample = readPositionsPerSample;
	}

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
		readPositionsPerSample = new ArrayList<HashMap<String, MinimalRead>>();

		coverageStranded = null;
		if (readGroupMap == null) {
			coverageStranded = new int[2][1][chromosomeWindowSize]; // [strand][ind][bp]
			readPositionsPerSample.add(new HashMap<String, MinimalRead>(1000));
		} else {
			coverageStranded = new int[2][readGroupMap.size()][chromosomeWindowSize]; // [strand][ind][bp]
			for (int i = 0; i < readGroupMap.size(); i++) {
				readPositionsPerSample.add(new HashMap<String, MinimalRead>(1000));
			}
		}

		nrFragments = 0;
		nrFragmentsDup = 0;
		nrFragmentsPassingFilter = 0;

		nrReads1 = 0;
		nrReads2 = 0;
		nrDups = 0;
		nrReadsPassingFilter = 0;
		nrDupsPassingFilter = 0;



		if (determineInsertSize) {
			insertSize = new int[maxInsertSize];
		}


		while (it.hasNext()) {
			SAMRecord record = it.next();

			Feature f = new Feature();
			f.setChromosome(Chromosome.parseChr(record.getReferenceName()));
			f.setStart(record.getAlignmentStart());
			f.setStop(record.getAlignmentEnd());

			if (chromosomeWindow.overlaps(f)) {
				if (!record.isSecondaryOrSupplementary()) {
					if (record.getFirstOfPairFlag()) {
						nrReads1++;
					} else {
						nrReads2++;
					}
					if (countDuplicates && record.getDuplicateReadFlag()) {
						nrDups++;
					}
					if (countFragments) {
						if (record.getFirstOfPairFlag()) {
							nrFragments++;
						}
					}
				}

				if (record.getMappingQuality() > 0 && !record.isSecondaryOrSupplementary()) {



					if (filter == null || !filter.filterOut(record)) {

						if (countReads) {
							nrReadsPassingFilter++;
						}

						if (record.getFirstOfPairFlag()) {
							nrFragmentsPassingFilter++;
						}


						if (countDuplicates && record.getDuplicateReadFlag()) {
							nrDupsPassingFilter++;
						}

						if (determineInsertSize) {
							if (insertSize == null) {
								insertSize = new int[maxInsertSize];
							}
							if (record.getFirstOfPairFlag()) {
								int insertSize = record.getInferredInsertSize();
								if (insertSize > maxInsertSize - 1) {
									insertSize = maxInsertSize - 1;
								}
								if (insertSize > 0) {
									this.insertSize[insertSize]++;
								}
							}
						}

						byte[] baseQual = record.getBaseQualities();
						byte[] bases = record.getReadBases();

						Cigar cigar = record.getCigar();
						List<CigarElement> cigarElements = cigar.getCigarElements();

						boolean includeReadInPermutations = false;

						if (bases == null || bases.length == 0 || baseQual == null || baseQual.length == 0) {
							// System.err.println("No bases for readAsTrack: " + record.getReadUnmappedFlag() + "\t" + record.toString() + "\t" + bases.length + "\t" + baseQual.length);
						} else {

							Strand strand = Strand.POS;
							int mapPos = record.getAlignmentStart();
							int windowRelativePosition = mapPos - chromosomeWindowStart;
							if (record.getReadNegativeStrandFlag()) {
								strand = Strand.NEG;
								windowRelativePosition += negShift;
							} else {
								windowRelativePosition += posShift;
							}

							int readPosition = 0;

							Integer sampleId = 0;
							int[] tmpCoverage = null;
							String sampleName = null;
							if (readGroupMap != null) {
								SAMReadGroupRecord rg = record.getReadGroup();
								if (rg == null) {
									// no readgroups set...
									String filename = record.getHeader().getAttribute("filename");
									rg = new SAMReadGroupRecord(filename);
									sampleName = filename;
								}
								sampleId = readGroupMap.get(rg.getId());
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

											if (windowRelativePosition >= 0 && windowRelativePosition < chromosomeWindowSize) { // the readAsTrack could overlap the leftmost edge of this window
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
												} else if (readPosition < baseQual.length && baseQual[readPosition] > 0) {
													// System.err.println("Low Qual base: " + base + ". Qual: " + baseQual[readPosition]);
												} else {
													System.err.println("Unknown base found: " + base);
												}
											}

											if (properbase) {
												tmpCoverage[windowRelativePosition]++;
												includeReadInPermutations = true;
											}
											windowRelativePosition++;
										} // if pos < readposition
										readPosition += cigarElementLength;
										break;
									default:
										System.err.println("Unknown CIGAR operator found: " + e.getOperator().toString());
										System.err.println("In readAsTrack: " + record.toString());
										break;
								} // switch operator
							} // for each cigar element
							coverageStranded[strand.getNumber()][sampleId] = tmpCoverage;
							if (includeReadInPermutations) {
								HashMap<String, MinimalRead> features = readPositionsPerSample.get(sampleId);
								String recordName = new String(record.getReadName()).intern();
								MinimalRead read = features.get(recordName);

								if (read == null) {
									read = new MinimalRead(recordName,
											chromosomeWindow.getChromosome(),
											record.getAlignmentStart(),
											record.getAlignmentEnd(),
											record.getCigar(),
											strand);
									features.put(recordName, read);
								} else {
									MinimalRead mate = new MinimalRead(recordName,
											chromosomeWindow.getChromosome(),
											record.getAlignmentStart(),
											record.getAlignmentEnd(),
											record.getCigar(),
											strand);
									if (read.getMate() != null) {
										System.err.println("WARNING: readAsTrack: " + read.getName() + " already has a mate!: " + read.getMate().getChr().toString() + "\t" + read.getMate().getStart() + "\t" + read.getMate().getEnd());
									}
									read.setMate(mate);
								}
							}
						}
					} else {
						// not passing filter...
					}
				}

			} // does not overlap window
		} // main while loop.
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

	public int[][] permuteReadPositionsPerSample(Feature window, int posShift, int negShift) {
		int windowStart = window.getStart();
		int windowEnd = window.getStop();
		int windowSize = windowEnd - windowStart;

		int[][] coverageOverall = new int[readPositionsPerSample.size()][0];

		for (int sample = 0; sample < readPositionsPerSample.size(); sample++) {
			HashMap<String, MinimalRead> f = readPositionsPerSample.get(sample);
			Set<Map.Entry<String, MinimalRead>> reads = f.entrySet();

			int[] coverage = new int[windowSize];

			for (Map.Entry<String, MinimalRead> readEntry : reads) {
				// randomize midpoint
				double direction = Math.random();
				if (direction >= 0.5) {
					direction = 1;
				} else {
					direction = -1;
				}

				double distance = direction * (Math.random() * windowSize);

				// make sure to also shift the mate in the same direction;
				MinimalRead read = readEntry.getValue();
				MinimalRead mate = read.getMate();

				int alignmentStart = read.getStart();


				Strand strand = read.getStrand();
				if (strand.equals(Strand.NEG)) {
					alignmentStart += distance - negShift;
				} else {
					alignmentStart += distance + posShift;
				}


				Cigar readCigar = read.getCigar();


				updateCoverageUsingCigar(readCigar,
						windowSize,
						alignmentStart,
						coverage);

				if (mate != null) {

					int alignmentStartMate = mate.getStart();

					Strand strandMate = mate.getStrand();
					if (strandMate.equals(Strand.NEG)) {
						alignmentStartMate += distance - negShift;
					} else {
						alignmentStartMate += distance + posShift;
					}

					Cigar mateCigar = mate.getCigar();
					updateCoverageUsingCigar(mateCigar,
							windowStart,
							alignmentStartMate,
							coverage);
				}
			}
			coverageOverall[sample] = coverage;
		}
		return coverageOverall;
	}

	private void updateCoverageUsingCigar(Cigar c, int chromosomeWindowStart, int start, int[] coverage) {
		int readPosition = 0;
		List<CigarElement> cigarElements = c.getCigarElements();

		int mapPos = start;

		int windowRelativePosition = mapPos - chromosomeWindowStart;
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
						coverage[windowRelativePosition]++;
						windowRelativePosition++;
					} // if pos < readposition
					readPosition += cigarElementLength;
					break;
				default:
					System.err.println("Unknown CIGAR operator found: " + e.getOperator().toString());
					break;
			} // switch operator
		} // for each cigar element
	}

	public Pair<Integer, Integer> getReadStartAndEnd(MinimalRead read, Feature window) {
		int readPosition = 0;
		Cigar c = read.getCigar();
		List<CigarElement> cigarElements = c.getCigarElements();
		int windowStart = window.getStart();
		int mapPos = read.getStart();

		int windowRelativePosition = mapPos - windowStart;
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

						windowRelativePosition++;
					} // if pos < readposition
					readPosition += cigarElementLength;
					break;
				default:
					System.err.println("Unknown CIGAR operator found: " + e.getOperator().toString());
					break;
			} // switch operator
		} // for each cigar element
		return new Pair<Integer, Integer>(mapPos, windowStart + windowRelativePosition);
	}

	public void setCountReads(boolean countReads) {
		this.countReads = countReads;
	}

	public void setCountDuplicates(boolean countDuplicates) {
		this.countDuplicates = countDuplicates;
	}

	public long getNrReads() {
		return nrReads1 + nrReads2;
	}

	public long getNrDups() {
		return nrDups;
	}

	public void setDetermineInsertSize(boolean determineInsertSize, int maxInsertSize) {
		this.determineInsertSize = determineInsertSize;
		this.maxInsertSize = maxInsertSize;
	}

	public int[] getInsertSize() {
		return insertSize;
	}

	public long getNrReadsPassingFilter() {
		return nrReadsPassingFilter;
	}

	public long getnrDupsPassingFilter() {
		return nrDupsPassingFilter;
	}

	public void setCountMultiMapping(boolean countMultiMapping) {
		this.countMultiMapping = countMultiMapping;
	}

	public void setCountSupplementaryReads(boolean countSupplementaryReads) {
		this.countSupplementaryReads = countSupplementaryReads;
	}

	public boolean isCountReads() {
		return countReads;
	}

	public boolean isCountDuplicates() {
		return countDuplicates;
	}

	public void setNrReads(long nrReads) {
		this.nrReads1 = nrReads;
	}

	public void setNrDups(long nrDups) {
		this.nrDups = nrDups;
	}

	public boolean isDetermineInsertSize() {
		return determineInsertSize;
	}

	public void setDetermineInsertSize(boolean determineInsertSize) {
		this.determineInsertSize = determineInsertSize;
	}

	public void setInsertSize(int[] insertSize) {
		this.insertSize = insertSize;
	}

	public int[] getFragmentsizeAfterFilter() {
		return fragmentsizeAfterFilter;
	}

	public void setFragmentsizeAfterFilter(int[] fragmentsizeAfterFilter) {
		this.fragmentsizeAfterFilter = fragmentsizeAfterFilter;
	}

	public void setNrReadsPassingFilter(long nrReadsPassingFilter) {
		this.nrReadsPassingFilter = nrReadsPassingFilter;
	}

	public long getNrDupsPassingFilter() {
		return nrDupsPassingFilter;
	}

	public void setNrDupsPassingFilter(long nrDupsPassingFilter) {
		this.nrDupsPassingFilter = nrDupsPassingFilter;
	}

	public int getMaxInsertSize() {
		return maxInsertSize;
	}

	public void setMaxInsertSize(int maxInsertSize) {
		this.maxInsertSize = maxInsertSize;
	}

	public boolean isCountMultiMapping() {
		return countMultiMapping;
	}

	public boolean isCountSupplementaryReads() {
		return countSupplementaryReads;
	}

	public boolean isCountFragments() {
		return countFragments;
	}

	public void setCountFragments(boolean countFragments) {
		this.countFragments = countFragments;
	}

	public long getNrFragments() {
		return nrFragments;
	}

	public void setNrFragments(long nrFragments) {
		this.nrFragments = nrFragments;
	}

	public long getNrFragmentsDup() {
		return nrFragmentsDup;
	}

	public void setNrFragmentsDup(long nrFragmentsDup) {
		this.nrFragmentsDup = nrFragmentsDup;
	}

	public long getNrFragmentsPassingFilter() {
		return nrFragmentsPassingFilter;
	}

	public void setNrFragmentsPassingFilter(long nrFragmentsPassingFilter) {
		this.nrFragmentsPassingFilter = nrFragmentsPassingFilter;
	}
}
