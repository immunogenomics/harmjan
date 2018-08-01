package nl.harmjanwestra.polypeak.tasks;

import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import nl.harmjanwestra.polypeak.comparators.WindowCoverageComparator;
import nl.harmjanwestra.polypeak.containers.Sample;
import nl.harmjanwestra.polypeak.containers.Sequence;
import nl.harmjanwestra.polypeak.containers.SequenceLibrary;
import nl.harmjanwestra.polypeak.containers.Window;
import nl.harmjanwestra.utilities.bamfile.BamFileReader;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.Correlation;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.Descriptives;
import nl.harmjanwestra.utilities.legacy.genetica.util.Primitives;
import nl.harmjanwestra.utilities.legacy.genetica.util.RankArray;
import nl.harmjanwestra.utilities.legacy.genetica.util.RunTimer;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.Callable;

/**
 * Created by hwestra on 1/8/15.
 * calculate statistics for each individual sample
 */
public class ModelBuildingTask implements Callable<Boolean> {


	private HashSet<Sequence> selectedChromosomes;

	public void ModelBuildingTask() {

	}


	public void setWindowSize(int windowSize) {
		this.windowSize = windowSize;
	}

	public void setWindowOverlap(int windowOverlap) {
		this.windowOverlap = windowOverlap;
	}

	public void setSample(Sample sample) {
		this.sample = sample;
	}

	public void setFinalWindowBufferSize(int finalWindowBufferSize) {
		this.finalWindowBufferSize = finalWindowBufferSize;
	}

	public void setTmpWindowBufferSize(int tmpWindowBufferSize) {
		this.tmpWindowBufferSize = tmpWindowBufferSize;
	}


	public void setOutputDirectory(String outputDirectory) {
		this.outputDirectory = outputDirectory;
	}

	int windowSize = 1000;
	int windowOverlap = 100;


	int[] coverageDistributionPosStr = new int[1000];
	int[] coverageDistributionNegStr = new int[1000];
	int[] insertSizeDistribution = new int[1000];
	int[] readLengthDistribution = new int[200];
	int nrReadsPosStrand = 0;
	int nrReadsNegStrand = 0;

	Sample sample = null;

	int finalWindowBufferSize = 1000;
	int tmpWindowBufferSize = 100;
	private boolean bufferHasOverflown = false;
	private int maxCoverageValueInBuffer = 0;

	int nrWindowsInTMPBuffer = 0;
	int nrWindowsInFinalBuffer = 0;

	String outputDirectory = "";

	DecimalFormat format = new DecimalFormat("#.##");

	@Override
	public Boolean call() {


		// TODO: check the regions list and blacklisted regions

		// get readers for sample
		try {
			BamFileReader[] readers = sample.getReaders();
			TextFile sampleOutput = new TextFile(outputDirectory + sample.getName() + "-coverage.txt", TextFile.W);
			String header = "File\tWindow\tWindowSize\tNrRecords\tNrRecordsPassingFilter\t" +
					"nrPosRecords\tmeanPosCoverage\tvarPosCoverage\tnrNegRecords\tmeanNegCoverage\tvarNegCoverage\t" +
					"maxCoverage\tmeanInsertSize\tvarInsertSize\tmeanReadLength\tvarReadLength";
			sampleOutput.writeln(header);
			for (BamFileReader reader : readers) {

				// build model for each reader...
				System.out.println("Processing sample: " + sample.getName() + "\t" + reader.getFile().getName());
				SAMSequenceDictionary sequenceDictionary = reader.getHeader().getSequenceDictionary();
				List<SAMSequenceRecord> sequenceRecords = sequenceDictionary.getSequences();

				// initialize two buffers
				Window[] tmpBuffer = new Window[tmpWindowBufferSize];
				Window[] finalBuffer = new Window[finalWindowBufferSize + tmpWindowBufferSize];
				nrWindowsInFinalBuffer = 0;
				nrWindowsInTMPBuffer = 0;
				bufferHasOverflown = false;
				// for each sequence


				for (SAMSequenceRecord sequenceRecord : sequenceRecords) {
					RunTimer timer = new RunTimer();
					timer.start();
					int sequenceLength = sequenceRecord.getSequenceLength();
					String sequenceName = sequenceRecord.getSequenceName();
					int position = 0;
					int windowCt = 0;
					long nrWindowsTotal = sequenceLength / windowSize;

					Sequence correspondingSequence = SequenceLibrary.getSequence(new Sequence(sequenceName, sequenceLength));

					if (selectedChromosomes == null || selectedChromosomes.contains(correspondingSequence)) {
						System.out.println("Sample: " + sample.getName() + "\tParsing sequence: " + sequenceName + "\tLength: " + sequenceLength + "\tWindows: " + nrWindowsTotal + ")");


						// divide up in regions
						while (position < sequenceLength) {

							int windowStart = position;
							int windowStop = position + windowSize;
							if (windowStart < 0) {
								windowStart = 0;
							}
							if (windowStop > sequenceLength) {
								windowStop = sequenceLength;
							}

							Window pileupWindow = new Window(windowStart, windowStop, correspondingSequence);
							pileupWindow.setStranded(true);
							if (windowCt % 10000 == 0) {
								double perc = ((double) windowCt / nrWindowsTotal) * 100;
								System.out.println("Sample: " + sample.getName() + "\t" + pileupWindow.toString() + "\t" + windowCt + "/" + nrWindowsTotal + "\t" + format.format(perc) + "%\tdeltaT: " + timer.getTimeDesc());
							}

							// get region from sample
							SAMRecordIterator iterator = reader.query(sequenceName, windowStart, windowStop, true);

							// determine pileup
							PileUp.pileUp(iterator, pileupWindow);
							iterator.close();

							// store top list of pile-ups
							if (pileupWindow.getCoverage() != null) {
								// determine coverage statistics for this window

								updateGlobalStats(pileupWindow);

								if (pileupWindow.getMaxCoverage() > maxCoverageValueInBuffer) {
									tmpBuffer[nrWindowsInTMPBuffer] = pileupWindow;
									nrWindowsInTMPBuffer++;
								}


								// create bedlike path:
								// - store pileup average + std. deviance
								// - store average insert size
								// - store average readAsTrack length

								if (nrWindowsInTMPBuffer == tmpWindowBufferSize) {
									updateFinalWindowBuffer(tmpBuffer, finalBuffer);
									nrWindowsInTMPBuffer = 0;
								}

							}
							if (pileupWindow.getNrRecordsPassingFilter() > 0) {

								sampleOutput.writeln(reader.getFile().getName()
												+ "\t" + pileupWindow.toString()
												+ "\t" + (pileupWindow.getStop() - pileupWindow.getStart())
												+ "\t" + pileupWindow.getNrRecordsTotal()
												+ "\t" + pileupWindow.getNrRecordsPassingFilter()
												+ "\t" + pileupWindow.getPosNrReads()
												+ "\t" + pileupWindow.getPosMean()
												+ "\t" + pileupWindow.getPosVariance()
												+ "\t" + pileupWindow.getNegNrReads()
												+ "\t" + pileupWindow.getNegMean()
												+ "\t" + pileupWindow.getNegVariance()
												+ "\t" + pileupWindow.getMaxCoverage()
												+ "\t" + pileupWindow.getInsertSizeMean()
												+ "\t" + pileupWindow.getInsertSizeVariance()
												+ "\t" + pileupWindow.getReadLengthMean()
												+ "\t" + pileupWindow.getReadLengthVariance()
								);
							}

							position += windowSize;
							windowCt++;
						}
						System.out.println("Sample: " + sample.getName() + "\tNr windows processed: " + windowCt + "\tdeltaT: " + timer.getTimeDesc());
					}
				}
				if (nrWindowsInTMPBuffer > 0) {
					updateFinalWindowBuffer(tmpBuffer, finalBuffer);
					nrWindowsInTMPBuffer = 0;
				}
				sampleOutput.close();


				int nrToOutput = finalWindowBufferSize;
				if (nrWindowsInFinalBuffer < finalWindowBufferSize) {
					nrToOutput = nrWindowsInFinalBuffer;
				}


				// store regions used for D-calculation
				TextFile tfout1 = new TextFile(outputDirectory + sample.getName() + "-selectedWindows.txt", TextFile.W);
				TextFile tfout2 = new TextFile(outputDirectory + sample.getName() + "-selectedWindowsDCorrelation.txt", TextFile.W);
				TextFile tfout3 = new TextFile(outputDirectory + sample.getName() + "-selectedWindowsDCorrelationMaxD.txt", TextFile.W);

				RankArray ra = new RankArray();
				for (int i = 0; i < nrToOutput; i++) {
					Window pileupWindow = finalBuffer[i];

					String outputln1 = reader.getFile().getName()
							+ "\t" + pileupWindow.toString()
							+ "\t" + (pileupWindow.getStop() - pileupWindow.getStart())
							+ "\t" + pileupWindow.getNrRecordsTotal()
							+ "\t" + pileupWindow.getNrRecordsPassingFilter()
							+ "\t" + pileupWindow.getPosNrReads()
							+ "\t" + pileupWindow.getPosMean()
							+ "\t" + pileupWindow.getPosVariance()
							+ "\t" + pileupWindow.getNegNrReads()
							+ "\t" + pileupWindow.getNegMean()
							+ "\t" + pileupWindow.getNegVariance()
							+ "\t" + pileupWindow.getMaxCoverage()
							+ "\t" + pileupWindow.getInsertSizeMean()
							+ "\t" + pileupWindow.getInsertSizeVariance()
							+ "\t" + pileupWindow.getReadLengthMean()
							+ "\t" + pileupWindow.getReadLengthVariance();

					String outputln2 = outputln1 + "\t-";
					outputln1 += "\t+";

					// add coverage info
					int[][] coverageInfo = pileupWindow.getCoverage();
					double[][] coverageDouble = new double[coverageInfo.length][coverageInfo[coverageInfo.length - 1].length];

					int windowSize = pileupWindow.getStop() - pileupWindow.getStart();


					for (int j = 0; j < coverageInfo.length; j++) {
						int[] strand = coverageInfo[j];


						String strandCoverage = "";
						for (int k = 0; k < windowSize; k++) {
							int coverage = strand[k];
							strandCoverage += "\t" + coverage;
							coverageDouble[j][k] = coverage;
						}
						if (j == 0) {
							outputln1 += strandCoverage;
						} else {
							outputln2 += strandCoverage;
						}

						// should we center + scale?
						// coverageDouble[j] = ra.rank(coverageDouble[j], true);
					}
					tfout1.writeln(outputln1);
					tfout1.writeln(outputln2);


					// determine mean rank of window
					double meanStr1 = Descriptives.mean(coverageDouble[0]);
					double meanStr2 = Descriptives.mean(coverageDouble[1]);
					int maxD = 0;
					double maxCorrelation = 0;

					String correlationOutput = reader.getFile().getName()
							+ "\t" + pileupWindow.getChr().getName()
							+ "\t" + pileupWindow.getStart()
							+ "\t" + pileupWindow.getStop();


					// perform cross correlation/covariance analysis
					// -- when overlapping windows: pick windows that have highest coverage.
					for (int j = 0; j < windowSize; j++) {
						if (j > 0) {
							shift(coverageDouble[1], 1);
						}
						double correlation = Correlation.correlate(meanStr1, meanStr2, coverageDouble[0], coverageDouble[1]);

						if (correlation > maxCorrelation) {
							maxD = j;
							maxCorrelation = correlation;
						}
						correlationOutput += "\t" + correlation;

					}
					tfout2.writeln(correlationOutput);
					tfout3.writeln(reader.getFile().getName()
							+ "\t" + pileupWindow.getChr().getName()
							+ "\t" + pileupWindow.getStart()
							+ "\t" + pileupWindow.getStop()
							+ "\t" + maxD + "\t" + maxCorrelation);

				}
				tfout1.close();
				tfout2.close();
				tfout3.close();


				// write readlength and insert size distributions
				//outputDirectory + sample.getName() + "-selectedWindows.txt", TextFile.W);

				TextFile readLengthOut = new TextFile(outputDirectory + sample.getName() + "-readLength.txt", TextFile.W);
				String rlHeader = "Length\tNr";
				readLengthOut.writeln(rlHeader);
				for (int i = 0; i < readLengthDistribution.length; i++) {
					readLengthOut.writeln(i + "\t" + readLengthDistribution[i]);
				}
				readLengthOut.close();

				TextFile insertSizeOut = new TextFile(outputDirectory + sample.getName() + "-insertSize.txt", TextFile.W);
				insertSizeOut.writeln(rlHeader);
				for (int i = 0; i < insertSizeDistribution.length; i++) {
					insertSizeOut.writeln(i + "\t" + insertSizeDistribution[i]);
				}
				insertSizeOut.close();

				TextFile coverageDistributionOut = new TextFile(outputDirectory + sample.getName() + "-coverageDistribution.txt", TextFile.W);
				String cvHeader = "Length\tNrPos\tnrNeg";
				coverageDistributionOut.writeln(cvHeader);

				for (int i = 0; i < coverageDistributionPosStr.length; i++) {
					coverageDistributionOut.writeln(i + "\t" + coverageDistributionPosStr[i] + "\t" + coverageDistributionNegStr[i]);
				}
				coverageDistributionOut.close();

			}

		} catch (IOException e) {
			e.printStackTrace();
		}
		return true;
	}

	private void updateGlobalStats(Window pileupWindow) {
		// get nr reads positive + negative strand
		nrReadsPosStrand += pileupWindow.getPosNrReads();
		nrReadsNegStrand += pileupWindow.getNegNrReads();
		// update coverage distribution


		// get insert size distribution
		updateInsertSizeDistribution(pileupWindow.getInsertSizes(), pileupWindow);
		updateReadLengthDistribution(pileupWindow.getReadLengths(), pileupWindow);
		updateCoverageDistribution(pileupWindow);

	}

	private void updateCoverageDistribution(Window pileupWindow) {
		int windowSize = pileupWindow.getStop() - pileupWindow.getStart();
		int[][] coverage = pileupWindow.getCoverage();
		for (int i = 0; i < windowSize; i++) {
			int coveragePos = coverage[0][i];
			int coverageNeg = coverage[1][i];

			if (coverageNeg >= coverageDistributionNegStr.length) {
				coverageNeg = coverageDistributionNegStr.length - 1;
			}
			if (coveragePos >= coverageDistributionPosStr.length) {
				coveragePos = coverageDistributionPosStr.length - 1;
			}
			coverageDistributionPosStr[coveragePos]++;
			coverageDistributionNegStr[coverageNeg]++;
		}
	}

	// TODO: filter negative insert sizes
	private void updateInsertSizeDistribution(ArrayList<Integer> insertSizes, Window pileupWindow) {
		int[] arr = Primitives.toPrimitiveArr(insertSizes.toArray(new Integer[0]));
		double mean = Descriptives.mean(arr);
		double variance = Descriptives.variance(arr, mean);

		pileupWindow.setInsertSizeMean(mean);
		pileupWindow.setInsertSizeVariance(variance);

		for (int i : arr) {
			if (i < 0) {
				i = 0;
				// System.err.println("WARNING: Negative insert size for window: " + pileupWindow.toString());
			}
			if (i >= insertSizeDistribution.length) {
				i = insertSizeDistribution.length - 1;
			}
			insertSizeDistribution[i]++;
		}

	}

	private void updateReadLengthDistribution(ArrayList<Integer> readLengths, Window pileupWindow) {
		int[] arr = Primitives.toPrimitiveArr(readLengths.toArray(new Integer[0]));
		double mean = Descriptives.mean(arr);
		double variance = Descriptives.variance(arr, mean);

		pileupWindow.setReadLengthMean(mean);
		pileupWindow.setReadLengthVariance(variance);

		for (int i : arr) {
			if (i < 0) {
				i = 0;
				System.err.println("WARNING: Negative readAsTrack length for window: " + pileupWindow.toString());
			}
			if (i >= readLengthDistribution.length) {
				i = readLengthDistribution.length - 1;
			}
			readLengthDistribution[i]++;
		}

	}

	// TODO: replace with arraycopy commands
	private void shift(double[] doubles, int shiftLen) {
		double[] firstValue = new double[shiftLen];
		for (int i = 0; i < shiftLen; i++) {
			firstValue[i] = doubles[i];
		}

		for (int i = shiftLen; i < doubles.length; i++) {
			doubles[i - shiftLen] = doubles[i];
		}

		for (int i = 0; i < shiftLen; i++) {
			doubles[doubles.length - shiftLen - i] = firstValue[i];
		}

	}


	private final WindowCoverageComparator windowCoverageComparator = new WindowCoverageComparator();

	private void updateFinalWindowBuffer(Window[] tmpBuffer, Window[] finalBuffer) {
		// keep a buffer with size of windowBuffer
		// - sort by:
		// -- height of pileup
		// -- number of covering reads


		if (bufferHasOverflown) {
			System.arraycopy(tmpBuffer, 0, finalBuffer, finalWindowBufferSize, nrWindowsInTMPBuffer);
			Arrays.sort(finalBuffer, windowCoverageComparator);
			maxCoverageValueInBuffer = finalBuffer[finalBuffer.length - 1].getMaxCoverage();
		} else {
			if (nrWindowsInFinalBuffer + nrWindowsInTMPBuffer > finalBuffer.length) {
				// do something with the surplus
				Window[] tmp = new Window[nrWindowsInFinalBuffer + nrWindowsInTMPBuffer];
				System.arraycopy(finalBuffer, 0, tmp, 0, nrWindowsInFinalBuffer);
				System.arraycopy(tmpBuffer, 0, tmp, nrWindowsInFinalBuffer, nrWindowsInTMPBuffer);
				Arrays.sort(finalBuffer, windowCoverageComparator);
				System.arraycopy(tmp, 0, tmp, 0, finalBuffer.length);
				nrWindowsInFinalBuffer = finalBuffer.length;
				bufferHasOverflown = true;

			} else {
				System.arraycopy(tmpBuffer, 0, finalBuffer, nrWindowsInFinalBuffer, nrWindowsInTMPBuffer);
				nrWindowsInFinalBuffer += nrWindowsInTMPBuffer;
			}

			if (nrWindowsInFinalBuffer == finalBuffer.length) {
				Arrays.sort(finalBuffer, windowCoverageComparator);
				maxCoverageValueInBuffer = finalBuffer[finalBuffer.length - 1].getMaxCoverage();
				bufferHasOverflown = true;
			}
		}
		// System.out.println("Updated buffer...." + nrWindowsInFinalBuffer + "\t" + nrWindowsInTMPBuffer + "\t" + maxCoverageValueInBuffer);
	}

	private void printWindows(Window[] finalBuffer) {
		for (int i = 0; i < finalBuffer.length; i++) {
			//System.out.println(i + "\t" + finalBuffer[i].getMaxCoverage());
		}

		// System.exit(0);
	}


	public void setSelectedChromosomes(HashSet<Sequence> selectedChromosomes) {
		this.selectedChromosomes = selectedChromosomes;
	}
}
