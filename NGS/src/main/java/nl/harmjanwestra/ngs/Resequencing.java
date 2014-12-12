/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.ngs;

import com.lowagie.text.DocumentException;
import htsjdk.samtools.*;
import htsjdk.samtools.filter.*;
import nl.harmjanwestra.ngs.containers.ReadGroup;
import nl.harmjanwestra.ngs.graphics.LocusCoveragePlot;
import nl.harmjanwestra.ngs.wrappers.MarkDups;
import nl.harmjanwestra.ngs.wrappers.SortBAM;
import nl.harmjanwestra.utilities.bamfile.BamFileReader;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import org.broadinstitute.gatk.engine.filters.MalformedReadFilter;
import org.broadinstitute.gatk.engine.filters.MappingQualityUnavailableFilter;
import org.broadinstitute.gatk.engine.filters.MateSameStrandFilter;
import org.broadinstitute.gatk.engine.filters.UnmappedReadFilter;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.HCMappingQualityFilter;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * @author hwestra
 */
public class Resequencing {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		// TODO code application logic here

		Resequencing s = new Resequencing(args);

	}

	private CLI cli;

	private Resequencing(String[] args) {

		String regions = args[0];
		String bamfile = args[1];
		String output = args[2];

		try {

//            determineCoverageForRegions(regions, bamfile, output);
		} catch (Exception e) {
			e.printStackTrace();
		}

//        if (args.length < 3) {
//            System.out.println("usage: input outdir stripoldqualityscores");
//        } else {
//            try {
//                splitPerChromosomeRemoveDuplicateReadsAndQualityScores(args[0], args[1], Boolean.parseBoolean(args[2]));
//
////        try {
////
////            args = new String[]{
////                "/Data/ATAC-seq/GSE47753/data2.txt",
////                "/Data/ATAC-seq/GSE47753/Stats"};
////
////            Resequencing.readsPerChromosomePerReadGroup(args);
////
////        } catch (Exception e) {
////            e.printStackTrace();
////        }
////        // dedup sort index
////        if (args.length < 2) {
////            System.out.println("Usage: in tmp");
////        } else {
////            String fileIn = args[0];
////
////            String tmpdir = args[1];
////
////            this.indexDedupSortIndex(fileIn, tmpdir);
////
////        }
//                //
//            /*
//                 // rewrite platform tags
//                 if (args.length < 3) {
//                 System.out.println("Usage: in out tmp");
//                 } else {
//                 String fileIn = args[0];
//                 String fileOut = args[1];
//                 String tmpdir = args[2];
//                 try {
//                 this.rewritePlatform(fileIn, fileOut, tmpdir);
//                 } catch (IOException ex) {
//                 Logger.getLogger(Resequencing.class.getName()).log(Level.SEVERE, null, ex);
//                 }
//                 }
//                 */
//            } catch (IOException ex) {
//                Logger.getLogger(Resequencing.class.getName()).log(Level.SEVERE, null, ex);
//            }
//        }
	}

	private void combine(String[] args) {
		cli = new CLI(args);

		// stuff for merger
		// get list of samples to use
		try {
			CoverageMatrix mat = new CoverageMatrix(cli.getCoverageMatrix());

			String indir = cli.getBamFileIndir();
			String outdir = cli.getBamFileOutdir();

			indir = Gpio.formatAsDirectory(indir);
			outdir = Gpio.formatAsDirectory(outdir);
			Gpio.createDir(outdir);
			String tmpdir = outdir + "tmp";
			tmpdir = Gpio.formatAsDirectory(tmpdir);
			Gpio.createDir(tmpdir);
			String suffix = cli.getBamSuffix();

			ArrayList<Pair<String, ArrayList<String>>> samples = mat.selectSamples(cli.getCoverageThresholdColumn(), cli.getCoverageThreshold(), outdir);
			System.out.println(samples.size() + " samples selected");
			TextFile outfile = new TextFile(outdir + "SamplesToUse.txt", TextFile.W);

			for (Pair<String, ArrayList<String>> p : samples) {

				String sampleName = p.getLeft();
				ArrayList<String> dups = p.getRight();

				if (dups.size() > 1) {
					// check whether files are actually there..
					int d = 0;
					ArrayList<File> finalFiles = new ArrayList<File>();
					System.out.println(sampleName + "\thas " + dups.size() + " duplicates");
					for (String f : dups) {
						if (Gpio.exists(indir + f + "." + suffix)) {
							d++;
							finalFiles.add(new File(indir + f + "." + suffix));
						} else {
							System.out.println(sampleName + "\tFile not found: " + indir + f + "." + suffix);
						}
					}
					System.out.println(sampleName + "\tFinal duplication size: " + finalFiles.size());
					if (finalFiles.size() > 1) {

						// merge bam files
						File outputFile = new File(outdir + sampleName + "_merged.bam");
						System.out.println(sampleName + "\tMerging BAM files");
						mergeBamFiles(sampleName, finalFiles, outputFile);
						// dedup

						File deduppedOut = new File(outdir + sampleName + "_merged.dedup.bam");
						File deduppedMetricsOut = new File(outdir + sampleName + "_merged.dedupmetrics.txt");
						System.out.println(sampleName + "\tDeduplicating alignments in BAM file" + outputFile.getAbsolutePath());
						System.out.println(sampleName + "\tWriting to: " + deduppedOut.getAbsolutePath());
						MarkDups md = new MarkDups(outputFile, deduppedOut, deduppedMetricsOut, new File(tmpdir));
// md.go();
						// sort file
						File deduppedSortOut = new File(outdir + sampleName + "_merged.dedup.sort.bam");
						SortBAM sb = new SortBAM(deduppedOut, deduppedSortOut, new File(tmpdir));

						// index
						indexBAM(deduppedSortOut);
						outfile.writeln(deduppedSortOut.getAbsolutePath());
					} else {
						// print filename to file
						if (finalFiles.size() == 1) {
							// check whether the file has been indexed..

							if (!Gpio.exists(finalFiles.get(0).getAbsolutePath() + ".bai")) {
								System.out.println(sampleName + "\tIndexing BAM: " + finalFiles.get(0).getAbsolutePath());
								indexBAM(finalFiles.get(0));
							}
							outfile.writeln(finalFiles.get(0).getAbsolutePath());
						}
					}

				} else {
					// check whether file is there,
					String f = dups.get(0);
					System.out.println(sampleName + "\thas no duplicates.");
					if (Gpio.exists(indir + f + "." + suffix)) {
						outfile.writeln(indir + f + "." + suffix);
						// check whether file has been indexed...
						if (!Gpio.exists(indir + f + "." + suffix + ".bai")) {
							System.out.println(sampleName + "\tIndexing BAM: " + indir + f + "." + suffix);
							indexBAM(new File(indir + f + "." + suffix));
						}
						System.out.println(sampleName + "\tFound file: " + indir + f + "." + suffix);
					} else {
						System.out.println(sampleName + "\tCould not find file: " + indir + f + "." + suffix);
					}
				}

			}
			outfile.close();

		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private void splitPerChromosomeRemoveDuplicateReadsAndQualityScores(String fileIn, String outDir, boolean stripOldQualityMetrics) throws IOException {
		outDir = Gpio.formatAsDirectory(outDir);
		Gpio.createDir(outDir);
		SAMFileReader reader = new SAMFileReader(new File(fileIn));
		SAMFileHeader header = reader.getFileHeader();

		SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();

		List<SAMSequenceRecord> sequences = sequenceDictionary.getSequences();
		List<SAMReadGroupRecord> var = header.getReadGroups();

		// create a new BAM file writer for each sequence
		HashMap<String, Integer> writerIndex = new HashMap<String, Integer>();
		ArrayList<SAMFileWriter> writers = new ArrayList<SAMFileWriter>();
		ArrayList<String> resultingBAMFiles = new ArrayList<String>();
		int samWriterIndex = 0;
		for (SAMSequenceRecord sequence : sequences) {
			ArrayList<SAMSequenceRecord> sequenceList = new ArrayList<SAMSequenceRecord>();
			sequenceList.add(sequence);
			SAMSequenceDictionary dict = new SAMSequenceDictionary(sequences);
			SAMFileHeader head = new SAMFileHeader();
			head.setSequenceDictionary(dict);
			head.setGroupOrder(header.getGroupOrder());
			head.setComments(header.getComments());
			head.setProgramRecords(header.getProgramRecords());
			head.setReadGroups(header.getReadGroups());
			head.setSortOrder(header.getSortOrder());
			head.setTextHeader(header.getTextHeader());
			head.setValidationErrors(header.getValidationErrors());

			String outfile = outDir + sequence.getSequenceName() + ".bam";
			resultingBAMFiles.add(outfile);
			System.out.println("Found sequence: " + sequence.getSequenceName() + ". Will output here: " + outfile);

			SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, new File(outfile));

			writers.add(outputSam);

			writerIndex.put(sequence.getSequenceName(), samWriterIndex);
			samWriterIndex++;
		}

		SAMRecordIterator iterator = reader.iterator();
		SAMFileWriter currentWriter = null;
		String currentReference = null;
		int readCtr = 0;
		int nrDups = 0;
		while (iterator.hasNext()) {
			SAMRecord record = iterator.next();
			String reference = record.getReferenceName();

			if (!record.getDuplicateReadFlag()) {
				if (currentReference == null || !reference.equals(currentReference)) {
					System.out.println("Switching reference: " + reference);
					Integer referenceIndex = writerIndex.get(reference);
					if (referenceIndex == null) {
						System.err.println("ERROR could not find writer for reference: " + reference);
					} else {
						currentWriter = writers.get(referenceIndex);
						currentReference = reference;
					}
				}

				// strip
				if (stripOldQualityMetrics) {
					record.setOriginalBaseQualities(null);
				}

				if (currentWriter != null) {
					currentWriter.addAlignment(record);
				} else {
					System.err.println("ERROR: there is no reference writer for record: " + record.getSAMString());
				}

			} else {
				nrDups++;
			}
			readCtr++;
			if (readCtr % 1000000 == 0) {
				System.out.println("Reads parsed: " + readCtr + "\t" + nrDups);
			}
		}

		System.out.println("Reads parsed total: " + readCtr);

		iterator.close();
		reader.close();

		for (SAMFileWriter w : writers) {
			w.close();

			// index resulting BAM files
		}

		for (String file : resultingBAMFiles) {
			indexBAM(new File(file));
		}

	}

	private void rewritePlatform(String fileIn, String fileOut, String tmp) throws IOException {

		SAMFileReader reader = new SAMFileReader(new File(fileIn));
		SAMFileHeader header = reader.getFileHeader();
		List<SAMReadGroupRecord> var = header.getReadGroups();
		System.out.println(var.size() + " readgroups found");
		for (SAMReadGroupRecord r : var) {
			System.out.println(String.format("Detected read group ID=%s PL=%s LB=%s SM=%s%n", r.getId(), r.getPlatform(), r.getLibrary(), r.getSample()));
		}

		List<SAMReadGroupRecord> varOut = new ArrayList<SAMReadGroupRecord>();

		for (SAMReadGroupRecord r : var) {
			final SAMReadGroupRecord rg = new SAMReadGroupRecord(r.getReadGroupId());
			rg.setLibrary(r.getLibrary());
			rg.setPlatform("ILLUMINA");
			rg.setSample(r.getSample());
			rg.setPlatformUnit(r.getPlatformUnit());
			System.out.println(String.format("\tCreated read group ID=%s PL=%s LB=%s SM=%s%n", rg.getId(), rg.getPlatform(), rg.getLibrary(), rg.getSample()));
			varOut.add(rg);
		}
		header.setReadGroups(varOut);

		String tmpdir = Gpio.formatAsDirectory(tmp);
		Gpio.createDir(tmp);
		String tmpfile = tmpdir + "output.bam";

		SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, new File(tmpfile));
		System.out.println("Writing to: " + tmpfile);
		SAMRecordIterator iterator = reader.iterator();
		while (iterator.hasNext()) {
			SAMRecord record = iterator.next();

			outputSam.addAlignment(record);

		}

		outputSam.close();
		reader.close();

		System.out.println("Sorgint BAM: " + tmpfile + " saving as: " + fileOut);
		SortBAM sb = new SortBAM(new File(tmpfile), new File(fileOut), new File(tmpdir));

		// index
		System.out.println("Indexing BAM: " + fileOut);
		indexBAM(new File(fileOut));

		new File(tmpfile).delete();
		new File(tmpdir).delete();
	}

	private void mergeBamFiles(String sample, ArrayList<File> finalFiles, File outfile) {
		ArrayList<SAMFileReader> readers = new ArrayList<SAMFileReader>();
		ArrayList<SAMFileHeader> headers = new ArrayList<SAMFileHeader>();
		for (File f : finalFiles) {
			System.out.println(sample + "\tIncluding file: " + f.getAbsolutePath());
			SAMFileReader tmpreader = new SAMFileReader(f);
			headers.add(tmpreader.getFileHeader());
			readers.add(tmpreader);
		}

		SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate, headers, true);
		MergingSamRecordIterator iterator = new MergingSamRecordIterator(headerMerger, readers, false);

		SAMFileHeader header = headerMerger.getMergedHeader();
		header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
		List<SAMReadGroupRecord> var = header.getReadGroups();
		System.out.println(var.size() + " readgroups found");
		for (SAMReadGroupRecord r : var) {
			System.out.println(String.format(sample + "\tDetected read group ID=%s PL=%s LB=%s SM=%s%n", r.getId(), r.getPlatform(), r.getLibrary(), r.getSample()));
		}

		String RGID = sample;
		String RGLB = var.get(0).getLibrary();
		String RGPL = var.get(0).getPlatform();
		String RGSM = sample;
		String RGPU = var.get(0).getPlatformUnit();

		final SAMReadGroupRecord rg = new SAMReadGroupRecord(RGID);
		rg.setLibrary(RGLB);
		rg.setPlatform(RGPL);
		rg.setSample(RGSM);
		rg.setPlatformUnit(RGPU);
		System.out.println(String.format(sample + "\tCreated read group ID=%s PL=%s LB=%s SM=%s%n", rg.getId(), rg.getPlatform(), rg.getLibrary(), rg.getSample()));
		header.setReadGroups(Arrays.asList(rg));

		SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(headerMerger.getMergedHeader(), false, outfile);
		System.out.println(sample + "\tWriting to: " + outfile);
		while (iterator.hasNext()) {
			SAMRecord record = iterator.next();
			record.setAttribute(SAMTag.RG.name(), RGID);
			outputSam.addAlignment(record);

		}

		outputSam.close();

		for (SAMFileReader r : readers) {
			r.close();
		}

		indexBAM(outfile);

	}

	private void indexBAM(File in) {
		SAMFileReader reader = new SAMFileReader(in);
		BAMIndexer.createIndex(reader, new File(in.getAbsolutePath() + ".bai"));
	}

	private void indexDedupSortIndex(String fileIn, String tmpdir) {

		String[] filesInDir = Gpio.getListOfFiles(fileIn);
		for (String f : filesInDir) {
			String path = fileIn + f;
			if (Gpio.isDir(path)) {
				String[] filesInDir2 = Gpio.getListOfFiles(path);
				for (String f2 : filesInDir2) {
					f2 = Gpio.formatAsDirectory(path) + f2;
					if (f2.endsWith(".bam")) {
						System.out.println("Running file:" + f2);
						System.out.println("-----------------------");
						System.out.println("");
						System.out.println("");
						File infile = new File(f2);
						File sortout = new File(f2.substring(0, f2.length() - 4) + "-sort.bam");

						if (!sortout.exists()) {
							SortBAM sb = new SortBAM(infile, sortout, new File(tmpdir));
							indexBAM(sortout);
						}

						// sort
						// dedup
						File dedupped = new File(f2.substring(0, f2.length() - 4) + "-sort-dedup.bam");
						File deduppedmetrics = new File(f2.substring(0, f2.length() - 4) + "-sort-dedupmetrics.txt");

						if (!dedupped.exists()) {
							MarkDups dupMark = new MarkDups(sortout, dedupped, deduppedmetrics, new File(tmpdir));
							indexBAM(dedupped);
							File sortout2 = new File(f2.substring(0, f2.length() - 4) + "-sort-dedup-sort.bam");
							SortBAM sb = new SortBAM(dedupped, sortout2, new File(tmpdir));
							indexBAM(sortout2);
							dedupped.delete();
						}
					}

				}
			}
		}
	}

	private static void mappingStatistics(String[] args) throws IOException {

		if (args.length < 1) {
			System.out.println("Usage: list.txt");
			System.exit(-1);
		}

		String list = args[1];
		TextFile tf = new TextFile(list, TextFile.R);
		String[] fileList = tf.readAsArray();
		tf.close();

		// open a wrapper around the HTSJDK
//        BamFileReader[] readerList = new BamFileReader[fileList.length];
		for (String fileName : fileList) {
			BamFileReader reader = new BamFileReader(new File(fileName), null, true);

			SAMRecordIterator iterator = reader.iterator();
			while (iterator.hasNext()) {

				SAMRecord record = iterator.next();

				Cigar c = record.getCigar();
//                int cigarElements = c.numCigarElements();
				List<CigarElement> cigarElements = c.getCigarElements();

				for (CigarElement e : cigarElements) {
					CigarOperator o = e.getOperator();
					o.consumesReadBases();
					o.consumesReferenceBases();

					int len = e.getLength();

				}

				c.getPaddedReferenceLength();
				c.getReadLength();
				c.getReferenceLength();

				List<AlignmentBlock> alignmentBlocks = record.getAlignmentBlocks();
				for (AlignmentBlock b : alignmentBlocks) {
					int l = b.getLength();
					int s1 = b.getReadStart();
					int s2 = b.getReferenceStart();
				}
				List<SAMRecord.SAMTagAndValue> attributes = record.getAttributes();
				for (SAMRecord.SAMTagAndValue a : attributes) {
					String tag = a.tag;
					Object val = a.value;
				}

				int start = record.getAlignmentStart();
				int stop = record.getAlignmentEnd();
				byte[] basequalities = record.getBaseQualities();

				boolean duplicateReadFlag = record.getDuplicateReadFlag();
				boolean firstOfPair = record.getFirstOfPairFlag();
				int alignmentFlag = record.getFlags();
				int inferredInsertSize = record.getInferredInsertSize();
				int mapQ = record.getMappingQuality();
				record.getMateAlignmentStart();
				record.getMateNegativeStrandFlag();
				record.getMateUnmappedFlag();
				byte[] originalbasequalities = record.getOriginalBaseQualities();

				record.getProperPairFlag();

				byte[] bases = record.getReadBases();
				record.getReadFailsVendorQualityCheckFlag();

				record.getReadLength();
				record.getReadNegativeStrandFlag();
				record.getReadPairedFlag();
//                readgroup = record.getReadGroup();

			}

			reader.close();
		}

		// coverage
	}

	private static void readsPerChromosomePerReadGroup(String[] args) throws IOException {
		if (args.length < 2) {
			System.out.println("Usage: list.txt output.txt");
			System.exit(-1);
		}

		String list = args[0];
		String output = args[1];
		TextFile tf = new TextFile(list, TextFile.R);
		String[] fileList = tf.readAsArray();
		tf.close();

		// open a wrapper around the HTSJDK
//        BamFileReader[] readerList = new BamFileReader[fileList.length];
		/*

         // sample
         // chr
         // 
        
         */
		HashMap<String, ReadGroup> strToSample = new HashMap<String, ReadGroup>();
		ArrayList<String> samples = new ArrayList<String>();
		for (String fileName : fileList) {
			int q = 0;
			System.out.println("Opening file: " + fileName);
			BamFileReader reader = new BamFileReader(new File(fileName), null, true);

			SAMRecordIterator iterator = reader.iterator();
			while (iterator.hasNext()) {

				SAMRecord record = iterator.next();
				String sample = fileName;
				String[] sampleNameElems = sample.split("/");
				sample = sampleNameElems[sampleNameElems.length - 1];
				Chromosome chr = Chromosome.parseChr(record.getReferenceName());
				SAMReadGroupRecord readgroup = record.getReadGroup();
				if (readgroup != null) {
					sample = readgroup.getSample();
				}

				ReadGroup s = strToSample.get(sample);
				if (s == null) {
					s = new ReadGroup();
					samples.add(sample);
					s.name = sample;
					strToSample.put(sample, s);
					System.out.println("found new Readgroup: " + sample);
				}

				if (record.getDuplicateReadFlag()) {
					if (record.getReadNegativeStrandFlag()) {
						s.dupsPerChr[0][chr.getNumber()]++;
					} else {
						s.dupsPerChr[1][chr.getNumber()]++;
					}
				}

				if (!record.getDuplicateReadFlag() && record.getProperPairFlag() && !record.getReadUnmappedFlag()) {
					// let's do a naïve approach first

					byte[] qual = record.getBaseQualities();
					for (int position = 0; position < qual.length; position++) {
						byte qualval = qual[position];
						if (qualval >= 60) {
							System.out.println(qualval);
							qualval = 59;
						}
						if (record.getReadNegativeStrandFlag()) {
							s.baseQualPerPosNegStrand[position][qualval]++;
						} else {
							s.baseQualPerPosPosStrand[position][qualval]++;
						}

					}

//                if (record.getProperPairFlag()) {
					if (record.getFirstOfPairFlag()) {
//                        if (!record.getMateUnmappedFlag()) {
						if (record.getReadNegativeStrandFlag()) {
							s.readsPairedPerChr[1][chr.getNumber()]++;
						} else {
							s.readsPairedPerChr[0][chr.getNumber()]++;
						}

					}
//                }

					if (record.getReadNegativeStrandFlag()) {
						s.readsPerChr[1][chr.getNumber()]++;
					} else {
						s.readsPerChr[0][chr.getNumber()]++;
					}

					int mapq = record.getMappingQuality();
					if (mapq > 99) {
						mapq = 0;
					}

					if (record.getReadNegativeStrandFlag()) {
						s.mapqPerChrNegStrand[chr.getNumber()][mapq]++;
					} else {
						s.mapqPerChrPosStrand[chr.getNumber()][mapq]++;
					}

					q++;
					if (q % 10000000 == 0) {
						System.out.println(q + " fragments processed ");
					}

				}

			}

			reader.close();

			System.out.println(strToSample.size() + " read groups found");

		}

		Chromosome[] allChr = Chromosome.values();

		TextFile outputf1 = new TextFile(output + "FragmentsPerChrPerSample.txt", TextFile.W);
		TextFile outputf2 = new TextFile(output + "DuplicatesPerChrPerSample.txt", TextFile.W);
		TextFile outputf3 = new TextFile(output + "ReadPairsPerChrPerSample.txt", TextFile.W);

		String header = "Chr";

		for (String sample : samples) {
			header += "\t" + sample;
		}

		outputf1.writeln(header);
		outputf2.writeln(header);
		outputf3.writeln(header);
		for (Chromosome chr : allChr) {
			String chrOutput = chr.getName();
			String chrOutput2 = chr.getName();
			String chrOutput3 = chr.getName();

			for (String sample : samples) {
				chrOutput += "\t" + strToSample.get(sample).readsPerChr[0][chr.getNumber()] + "; " + strToSample.get(sample).readsPerChr[1][chr.getNumber()];
				;
				chrOutput2 += "\t" + strToSample.get(sample).dupsPerChr[0][chr.getNumber()] + "; " + strToSample.get(sample).dupsPerChr[1][chr.getNumber()];
				chrOutput3 += "\t" + strToSample.get(sample).readsPairedPerChr[0][chr.getNumber()] + "; " + strToSample.get(sample).readsPairedPerChr[1][chr.getNumber()];
			}
			outputf1.writeln(chrOutput);
			outputf2.writeln(chrOutput2);
			outputf3.writeln(chrOutput3);
		}
		outputf1.close();
		outputf2.close();
		outputf3.close();

		// per sample stats
		for (String sample : samples) {
			String[] sampleNameElems = sample.split("/");
			ReadGroup s = strToSample.get(sample);
			TextFile outf = new TextFile(output + "-" + sampleNameElems[sampleNameElems.length - 1] + "-MapQPerChr.txt", TextFile.W);
			header = "Chr";
			for (int i = 0; i < 61; i++) {
				header += "\t" + i;
			}
			outf.writeln(header);
			for (Chromosome chr : allChr) {
				String chrOutput = chr.getName();
				for (int i = 0; i < 61; i++) {
					chrOutput += "\t" + s.mapqPerChrPosStrand[chr.getNumber()][i] + "; " + s.mapqPerChrNegStrand[chr.getNumber()][i];
				}
				outf.writeln(chrOutput);
			}
			outf.close();

			TextFile outf2 = new TextFile(output + "-" + sampleNameElems[sampleNameElems.length - 1] + "-BaseQualPerPos.txt", TextFile.W);
			header = "Pos";
			for (int i = 0; i < 51; i++) {
				header += "\t" + i; //+"\tAverage";
			}
			outf2.writeln(header);
			for (int pos = 0; pos < 50; pos++) {
				String outputStr = "" + pos;
				for (int qual = 0; qual < 60; qual++) {
					outputStr += "\t" + s.baseQualPerPosPosStrand[pos][qual] + "; " + s.baseQualPerPosNegStrand[pos][qual];
				}
				outf2.writeln(outputStr);
			}

			outf2.close();

			s.plotQualPerPos(output);
			s.plotMapQPerChr(output);

		}
	}

	private void determineCoverageForRegions(String regionsBedFile3, String bamfile, int readLength, String output) throws IOException {
		int doubleReadLength = readLength * 2;
		TextFile tf = new TextFile(regionsBedFile3, TextFile.R);
		ArrayList<Feature> features = new ArrayList<Feature>();
		String[] header = tf.readLineElems(TextFile.tab); // if any
		String[] ln = tf.readLineElems(TextFile.tab);
		int maxLen = 0;
		int minLen = Integer.MAX_VALUE;
		while (ln != null) {
			Feature f = new Feature();
			f.setChromosome(Chromosome.parseChr(ln[0]));
			f.setStart(Integer.parseInt(ln[1]));
			f.setStop(Integer.parseInt(ln[2]));

			int len = f.getStop() - f.getStart();
			if (len > maxLen) {
				maxLen = len;
			}
			if (len < minLen) {
				minLen = len;
			}
			features.add(f);
			ln = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		Collections.sort(features, new FeatureComparator(false));

		System.out.println("Number of features: " + features.size() + " min length: " + minLen + " max length: " + maxLen);

		// This is the aggregate filter that haplotypecaller uses.
		ArrayList<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();
		filters.add(new DuplicateReadFilter());
		filters.add(new NotPrimaryAlignmentFilter());
		filters.add(new FailsVendorReadQualityFilter());
		filters.add(new UnmappedReadFilter());
		filters.add(new MappingQualityUnavailableFilter());
		filters.add(new HCMappingQualityFilter());
		filters.add(new MalformedReadFilter());


		// add read mate should map to same strand
		filters.add(new MateSameStrandFilter());

		AggregateFilter aggregateFilter = new AggregateFilter(filters);
		BamFileReader reader = new BamFileReader(new File(bamfile), aggregateFilter, false);
		List<SAMReadGroupRecord> readGroups = reader.getReadGroups();

		// map readgroups to samplenames
		HashMap<SAMReadGroupRecord, Integer> readgroupMap = new HashMap<SAMReadGroupRecord, Integer>();
		String[] samples = new String[readGroups.size()];
		int rgctr = 0;
		for (SAMReadGroupRecord r : readGroups) {
			readgroupMap.put(r, rgctr);
			samples[rgctr] = r.getSample();
			rgctr++;
		}



		for (Feature f : features) {

			Chromosome featureChr = f.getChromosome();
			int featureStart = f.getStart();
			int featureStop = f.getStop();
			int featureSize = (featureStop - featureStart) + 2 * doubleReadLength;
			SAMRecordIterator iterator = reader.query("" + featureChr.getNumber(), featureStart - doubleReadLength, featureStop + doubleReadLength, false);

			int[][] coverageMapNegStrand = new int[readGroups.size()][featureSize];
			int[][] coverageMapPosStrand = new int[readGroups.size()][featureSize];

			int basequalthreshold = 50;

			while (iterator.hasNext()) {

				SAMRecord record = iterator.next();

				// record passes filter
				if (!aggregateFilter.filterOut(record)) {
					// do stuff to the record.

					// get the records readgroup
					SAMReadGroupRecord readGroup = record.getReadGroup();
					Integer rgId = readgroupMap.get(readGroup);
					if (rgId == null) {
						System.err.println("Readgroup not found: " + readGroup.getDescription());
					} else {
						// map it to the feature coverage map
						int[][] toMapTo = coverageMapPosStrand;
						if (record.getReadNegativeStrandFlag()) {
							toMapTo = coverageMapNegStrand;
						}

						//    -- if paired end, check whether other pair is present
						Cigar cigar = record.getCigar();
						List<CigarElement> cigarElements = cigar.getCigarElements();

						int mapPos = record.getAlignmentStart();

						int windowRelativePosition = mapPos - featureStart - doubleReadLength; // left most pos
						int readPosition = 0; // relative position in read/fragment

						byte[] baseQual = record.getBaseQualities();
						byte[] bases = record.getReadBases();

						// main loop
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
										byte base = bases[pos];
										if (windowRelativePosition >= 0) { // the read could overlap the leftmost edge of this window
											//    -- for each base pos: check basequal
											if (baseQual[readPosition] > basequalthreshold) {
												//    -- determine number of A/T/C/G/N bases
												if (base == 64 || base == 97) {
													toMapTo[rgId][windowRelativePosition]++;
												} else if (base == 67 || base == 99) {
													toMapTo[rgId][windowRelativePosition]++;
												} else if (base == 71 || base == 103) {
													toMapTo[rgId][windowRelativePosition]++;
												} else if (base == 84 || base == 116) {
													toMapTo[rgId][windowRelativePosition]++;
												}
											}
										}
									} // if pos < readposition
									readPosition += cigarElementLength;
									break;
								default:
									System.out.println("Unknown CIGAR operator found: " + e.getOperator().toString());
									System.out.println("In read: " + record.toString());
									break;
							} // switch operator
						} // for each cigarelement
					}
				}
			}

			// now plot the coverage per sample

			int pageMargin = 50;
			int betweenPlotMargin = 25;
			int plotIndividualWidth = maxLen;
			int plotIndividualHeight = 10;

			int pageHeight = (pageMargin * 2) + (betweenPlotMargin * (readGroups.size() - 1)) + (plotIndividualHeight * readGroups.size());
			int pageWidth = pageMargin * 2 + plotIndividualWidth;


			String name = output + "Locus-" + features.toString() + ".pdf";
			try {
				LocusCoveragePlot plot = new LocusCoveragePlot(name, pageWidth, pageHeight);

				plot.setMargin(pageMargin);
				plot.setBetweenPlotMargin(betweenPlotMargin);
				plot.setPlotIndividualWidth(plotIndividualWidth);
				plot.setPlotIndividualHeight(plotIndividualHeight);

				plot.setRowLabels(samples);
				plot.setCoverageData(coverageMapPosStrand, coverageMapNegStrand);
				plot.setFeature(f);
				plot.setFeatureMargin(doubleReadLength);
				// plot.setVariants();
				plot.setMaxCoverage(50);
				plot.draw();


				plot.close();
			} catch (DocumentException e) {
				e.printStackTrace();
			}


		}


		reader.close();

	}

}
