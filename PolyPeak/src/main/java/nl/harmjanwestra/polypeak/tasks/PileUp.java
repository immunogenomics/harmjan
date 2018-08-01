package nl.harmjanwestra.polypeak.tasks;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.filter.*;
import nl.harmjanwestra.polypeak.containers.Window;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.Descriptives;
import nl.harmjanwestra.utilities.legacy.genetica.util.Primitives;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by hwestra on 1/7/15.
 */
public class PileUp {

	private static AggregateFilter filter = null;

	public static void setFilter(AggregateFilter f) {
		filter = f;
	}


	public static void pileUp(SAMRecordIterator recordIterator, Window window) {
		pileUp(recordIterator, window, 0);
	}

	public static void pileUp(SAMRecordIterator recordIterator, Window window, int d) {
		int actualWindowStart = window.getStart();
		int nrRecordsPassingFilter = 0;
		int nrRecordsWindowTotal = 0;
		int maxInsertSize = 0;
		int nA = 0;
		int nT = 0;
		int nC = 0;
		int nG = 0;
		int nN = 0;
		int basequalthreshold = 0;

		int[][] tmpCoverage = null;
		if (window.isStranded()) {
			tmpCoverage = new int[2][window.getStop() - window.getStart()];
		} else {
			tmpCoverage = new int[1][window.getStop() - window.getStart()];
		}

		// TODO: filter for MAPQ, insert size, mate pair mapping, etc
		// TODO: setting that allows to only use the first readAsTrack of a pair (also useful if data is unpaired).

		ArrayList<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();
		filters.add(new DuplicateReadFilter());
		filters.add(new NotPrimaryAlignmentFilter());
		filters.add(new FailsVendorReadQualityFilter());
		filters.add(new AlignedFilter(true));
		filters.add(new SecondaryAlignmentFilter());
		filters.add(new WholeReadClippedFilter());


		AggregateFilter aggregateFilter = new AggregateFilter(filters);

		int windowSize = window.getStop() - window.getStart();

		ArrayList<Integer> insertSizes = new ArrayList<Integer>();
		ArrayList<Integer> readLengths = new ArrayList<Integer>();

		int nrRecordsPositive = 0;
		int nrRecordsNegative = 0;
		while (recordIterator.hasNext()) {
			SAMRecord record = recordIterator.next();

			int matchingBases = 0;
			// if (filter == null || !filter.filterOut(record)) {
			if (!aggregateFilter.filterOut(record)) {


				if (record.getFirstOfPairFlag()) {
					int insertSize = record.getInferredInsertSize();
					insertSizes.add(insertSize);
				}

				nrRecordsPassingFilter++;
				Cigar cigar = record.getCigar();
				List<CigarElement> cigarElements = cigar.getCigarElements();

				int mapPos = record.getAlignmentStart();
				int windowRelativePosition = mapPos - actualWindowStart;
				int strand = 0;

				if (window.isStranded()) {
					boolean negStrandFlag = record.getReadNegativeStrandFlag();
					if (negStrandFlag) {
						strand = 1;
						nrRecordsNegative++;
					} else {
						nrRecordsPositive++;
					}

				}
				int readPosition = 0;

				byte[] baseQual = record.getBaseQualities();
				byte[] bases = record.getReadBases();

				if (bases.length == 0) {
					System.err.println("No bases for readAsTrack: " + record.getReadUnmappedFlag() + "\t" + record.toString());
				} else {
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

									if (windowRelativePosition >= 0 && windowRelativePosition < windowSize) { // the readAsTrack could overlap the leftmost edge of this window
										if (base == 78 || base == 110) {
											nN++;
										} else if (baseQual[readPosition] > basequalthreshold) { //    -- for each base pos: check whether basequal > 50
											//    -- determine number of A/T/C/G/N bases
											if (base == 65 || base == 97) {
												nA++;
												properbase = true;
											} else if (base == 67 || base == 99) {
												nC++;
												properbase = true;
											} else if (base == 71 || base == 103) {
												nG++;
												properbase = true;
											} else if (base == 84 || base == 116) { // extend to capture U?
												nT++;
												properbase = true;
											}
										}
									}


									if (properbase) {
										matchingBases++;
										tmpCoverage[strand][windowRelativePosition]++;
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
				}


				readLengths.add(matchingBases);
			} // if filter != null || !filter.filterout
			nrRecordsWindowTotal++;
		}

		if (nrRecordsPassingFilter > 0) {

			window.setNrReadsPos(nrRecordsPositive);
			window.setNrReadsNeg(nrRecordsNegative);

			// System.out.println(window.toString() + "\tNr records passing filter: " + nrRecordsPassingFilter + " / " + nrRecordsWindowTotal + "\tMax: " + Primitives.max(tmpCoverage[0]) + " / " + Primitives.max(tmpCoverage[1]) + "  means: " + Descriptives.mean(tmpCoverage[0]) + " / " + Descriptives.mean(tmpCoverage[1]));

			window.setReadLengths(readLengths);
			window.setInsertSizes(insertSizes);

			window.setCoverage(tmpCoverage);

			window.setNrRecordsPassingFilter(nrRecordsPassingFilter);
			window.setNrRecordsTotal(nrRecordsWindowTotal);

			pileupWindowStatistics(window);
		}

	}

	private static void pileupWindowStatistics(Window pileupWindow) {

		int[][] coverage = pileupWindow.getCoverage();

		pileupWindow.setMaxCoverage(Math.max(Primitives.max(coverage[0]), Primitives.max(coverage[1])));

		// calculate mean + SD
		double meanPos = Descriptives.mean(coverage[0]);
		double meanNeg = Descriptives.mean(coverage[1]);
		double posVariance = Descriptives.variance(coverage[0], meanPos);
		double negVariance = Descriptives.variance(coverage[1], meanNeg);
		pileupWindow.setPosMean(meanPos);
		pileupWindow.setNegMean(meanNeg);
		pileupWindow.setPosVariance(posVariance);
		pileupWindow.setNegVariance(negVariance);
	}


}
