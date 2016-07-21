package nl.harmjanwestra.polypeak;

import htsjdk.samtools.filter.AggregateFilter;
import nl.harmjanwestra.polypeak.containers.Sample;
import nl.harmjanwestra.polypeak.containers.Sequence;
import nl.harmjanwestra.polypeak.containers.SequenceLibrary;
import nl.harmjanwestra.polypeak.tasks.ModelBuildingTask;
import nl.harmjanwestra.polypeak.tasks.PileUp;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Track;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.util.RunTimer;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Created by hwestra on 1/6/15.
 */
public class PolyPeak {

	private boolean checkSampleNamesInBAMFiles = true;

	private Sample[] samples;
	private HashMap<String, Sample> sampleNameToSample;
	private Track callRegions;
	private Track blackListedRegions;
	private AggregateFilter filter = null;

	private String outputDirectory = "";
	private int nrThreads;

	public static void main(String[] args) {
		PolyPeak p = new PolyPeak();

		if (args.length < 3) {
			System.out.println("usage: files.txt outdir nrthreads");
		} else {

			String sampleFileName = args[0]; //"/Data/ATAC-seq/GSE47753/dataWSampleNames.txt";
			String outputdirectory = args[1]; //"/Data/ATAC-seq/TestOutput/";
			int nrThreads = Integer.parseInt(args[2]);
			String regionListFile = null;
			String blackListFile = null;
			try {
				p.run(sampleFileName, regionListFile, blackListFile, outputdirectory, nrThreads);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		System.exit(0);
	}

	public void run(String sampleFileName, String regionListFile, String blackListFile, String outputDirectory, int nrThreads) throws IOException {

		RunTimer timer = new RunTimer();
		timer.start();

		this.nrThreads = nrThreads;
		this.outputDirectory = outputDirectory;

		// readAsTrack sample definition
		System.out.println("Initializing samples");
		initializeSamples(sampleFileName);
		System.out.println(samples.length + " samples found.");
		System.out.println("The data contains information on: " + SequenceLibrary.getSequences().size() + " sequences");

		// readAsTrack blacklisted regions
		if (blackListFile != null) {
			blackListedRegions = loadRegionsFromFile(blackListFile);
		} else {
			System.out.println("No blacklisted regions loaded.");
		}

		// initialize regions to call
		if (regionListFile != null) {
			callRegions = loadRegionsFromFile(regionListFile);
		} else {
			System.out.println("Running on all available regions.");
		}

		// TODO: initialize filter
		// TODO: initialize sequence library to add only sequences that we are interested in


		// TODO: initialize filters
		// determine whether PileUp should filter records:
		if (filter != null) {
			PileUp.setFilter(filter);
		}

		// determine D for each sample, for each reader
		System.out.println("Building model for samples individually");

		// get autosomes
		HashSet<Sequence> selectedChromosomes = new HashSet<Sequence>();
		for (int i = 1; i < 23; i++) {
			selectedChromosomes.add(new Sequence("" + i));
		}


		determineModelForEachSample(selectedChromosomes);

		// perform joint pileup
		performJointCalling();

		// close all sample BAM path readers
		for (Sample sample : samples) {
			sample.close();
		}

		System.out.println("Done processing. Processing time: " + timer.getTimeDesc());

	}


	private void determineModelForEachSample(HashSet<Sequence> selectedChromosomes) {
		// perform calculations for each sample on individual threads
		System.out.println("Opening threadpool for " + nrThreads + " threads.");
		ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
		CompletionService<Boolean> pool = new ExecutorCompletionService<Boolean>(threadPool);

		for (Sample s : samples) {
			ModelBuildingTask task = new ModelBuildingTask();
			task.setSample(s);
			task.setSelectedChromosomes(selectedChromosomes);
			task.setOutputDirectory(outputDirectory);
			pool.submit(task);
		}

		int returned = 0;
		while (returned < samples.length) {

			try {
				Boolean result = pool.take().get();
				if (result) {
					returned++;
				}
			} catch (Exception e) {
				e.printStackTrace();
			}

		}
		System.out.println("Closed down threads");

	}

	private void performJointCalling() {
		// initialize threadpool
		// perform calculations for each sample on individual threads

		// load windows over all samples in parallel from disk
		// pile-up the reads in another thread
		// we can only run one window at once
	}


	// TODO: this needs to be changed to something that can use flexible sequence names
	// this method currently uses the Chromosome object, which cannot account for weird contig names
	private Track loadRegionsFromFile(String regionListFile) throws IOException {
		BedFileReader reader = new BedFileReader();
		return null; //reader.readAsTrack(regionListFile, "");
	}

	private void initializeSamples(String sampleFileName) throws IOException {

		System.out.println("Parsing sample information: " + sampleFileName);
		TextFile tf = new TextFile(sampleFileName, TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, ArrayList<File>> sampleFiles = new HashMap<String, ArrayList<File>>();

		int ln = 1;
		HashMap<String, Boolean> controlStatus = new HashMap<String, Boolean>();
		HashSet<String> allSamples = new HashSet<String>();
		while (elems != null) {

			if (elems.length < 3) {
				System.err.println("Error: " + sampleFileName + " wrong format on line: " + ln);
				// throw new IllegalFormatException("Error: "+sampleFileName+ " does not have the correct format.");
			} else {
				if (elems.length >= 3) {
					String sampleName = elems[0].trim();
					String control = elems[1].trim();
					String fileName = elems[2].trim();

					if (sampleName.length() == 0) {
						System.err.println("Warning: sample name is empty on line: " + ln);
						// throw new IllegalFormatException("Error: "+sampleFileName+ " does not have the correct format.");
						break; // just break for now
					}

					if (fileName.length() == 0) {
						System.err.println("Warning: path name is empty on line: " + ln);
						// throw new IllegalFormatException("Error: "+sampleFileName+ " does not have the correct format.");
						break; // just break for now
					}

					boolean isControl = false;
					try {
						isControl = Boolean.parseBoolean(control);
					} catch (NumberFormatException e) {
						System.err.println("Expecting boolean at column 2 for line " + ln + " found " + control + " assuming control == false");
						// throw new IllegalFormatException("Error: "+sampleFileName+ " does not have the correct format.");
						break; // just break for now
					}

					File f = new File(fileName);
					if (!f.exists()) {
						System.err.println("Warning: " + fileName + " defined for sample " + sampleName + " does not exist.");
					} else {
						allSamples.add(sampleName);
						ArrayList<File> filesForSample = sampleFiles.get(sampleName);
						if (filesForSample == null) {
							filesForSample = new ArrayList<File>();
						}

						Boolean currentStatus = controlStatus.get(sampleName);
						if (currentStatus == null) {
							controlStatus.put(sampleName, isControl);
						} else {
							if (!controlStatus.equals(isControl)) {
								System.err.println("Sample " + sampleName + " control status: " + controlStatus + " but found " + isControl + " on line: " + ln);
							}
						}
						filesForSample.add(f);
						sampleFiles.put(sampleName, filesForSample);
					}
				}
			}
			elems = tf.readLineElems(TextFile.tab);
			ln++;
		}
		tf.close();

		for (String sample : allSamples) {
			if (!sampleFiles.containsKey(sample)) {
				System.err.println("Warning: " + sample + " is defined in sample definition, but suitable BAM files found.");
			}
		}

		// now initializing sample objects
		ArrayList<Sample> tmpSamples = new ArrayList<Sample>();
		Set<Map.Entry<String, ArrayList<File>>> sampleFileEntrySet = sampleFiles.entrySet();
		sampleNameToSample = new HashMap<String, Sample>();
		for (Map.Entry<String, ArrayList<File>> entry : sampleFileEntrySet) {
			String sampleName = entry.getKey();
			ArrayList<File> fileNames = entry.getValue();
			Sample tmpSample = new Sample(sampleName, fileNames.toArray(new File[0]), checkSampleNamesInBAMFiles, controlStatus.get(sampleName));
			if (tmpSample.hasReaders()) {
				if (sampleNameToSample.containsKey(sampleName)) {
					System.err.println("Error: something went wrong when loading samples");
				} else {
					sampleNameToSample.put(sampleName, tmpSample);
				}
				tmpSamples.add(tmpSample);
			}
		}
		samples = tmpSamples.toArray(new Sample[0]);
	}

}
