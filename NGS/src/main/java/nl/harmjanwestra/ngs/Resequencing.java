/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.ngs;

import com.lowagie.text.DocumentException;
import htsjdk.samtools.*;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.FailsVendorReadQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import nl.harmjanwestra.ngs.containers.ReadGroup;
import nl.harmjanwestra.ngs.graphics.LocusCoveragePlot;
import nl.harmjanwestra.ngs.wrappers.MarkDups;
import nl.harmjanwestra.ngs.wrappers.SortBAM;
import nl.harmjanwestra.utilities.bamfile.BamFileReader;
import nl.harmjanwestra.utilities.features.*;
import nl.harmjanwestra.utilities.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.vcf.VCFVariantType;
import org.broadinstitute.gatk.engine.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.gatk.engine.filters.MappingQualityUnavailableFilter;
import org.broadinstitute.gatk.engine.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.gatk.engine.filters.UnmappedReadFilter;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.HWE;
import umcg.genetica.text.Strings;

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

		// Resequencing s = new Resequencing(args);
		Resequencing s = new Resequencing(args);


	}

	private CLI cli;


	private Resequencing(String[] args) {


//		String indir = "/Data/Projects/2014-FR-Reseq/2014-12-18-HaplotypeCallerVariants-762samples/merged.vqsr.vcf";
//		String outdir = "/Data/Projects/2014-FR-Reseq/2014-12-18-HaplotypeCallerVariants-762samples/summary.txt";
//		try {
////
////			summarizeVCF(indir, outdir);
////
////			System.out.println("done");
////			System.exit(0);
//
//			String regions = args[0]; //"/Data/Projects/2014-Epipilot/atac/feat.txt";
//			String bam = args[1]; //"/Data/ATAC-seq/GSE47753/SRR891270/SRR891270-sort-dedup-sort.bam";
//			int readLength = 35;
//			String output = args[2]; // "/Data/Projects/2014-Epipilot/atac/";
//			String gtf = args[3];
//			determineCoverageForRegions(regions, bam, readLength, gtf, output);
//			if (args.length < 2) {
//				System.out.println("usage: indir outdir");
//			}
//
//        try {
//            vcfMerge(args[0], args[1]);
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//        indexDedupSortIndex(args[0], args[1]);

		if (args.length < 2) {
			// System.out.println("usage: input outdir stripoldqualityscores");
			System.out.println("usage: bamfile outfile");
		} else {
			try {
				String listFile = args[0];
				String outdir = args[1];
				String tmpdir = args[2];

				TextFile tf = new TextFile(listFile, TextFile.R);
				String[] filesIn = tf.readAsArray();
				tf.close();

				TextFile listOut = new TextFile(outdir + "SamplesRecoded.list", TextFile.W);

				for (String s : filesIn) {

					String[] bamfilelems = s.split("/");
					String bamfile = bamfilelems[bamfilelems.length - 1];
					String[] sampleelems = bamfile.split("\\.");
					String sample = sampleelems[0];
					String outfile = outdir + sample + ".bam";

					this.rewritePlatform(s, outfile, tmpdir);

					listOut.writeln(outfile);


				}


				listOut.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
//			try {
//
//				splitPerChromosomeRemoveDuplicateReadsAndQualityScores(args[0], args[1], Boolean.parseBoolean(args[2]));
////
//////        try {
//////
//////            args = new String[]{
//////                "/Data/ATAC-seq/GSE47753/data2.txt",
//////                "/Data/ATAC-seq/GSE47753/Stats"};
//////
//////            Resequencing.readsPerChromosomePerReadGroup(args);
//////
//////        } catch (Exception e) {
//////            e.printStackTrace();
//////        }
//////        // dedup sort index
//////        if (args.length < 2) {
//////            System.out.println("Usage: in tmp");
//////        } else {
//////            String fileIn = args[0];
//////
//////            String tmpdir = args[1];
//////
//////            this.indexDedupSortIndex(fileIn, tmpdir);
//////
//////        }
////                //
////            /*
////                 // rewrite platform tags
////                 if (args.length < 3) {
////                 System.out.println("Usage: in out tmp");
////                 } else {
////                 String fileIn = args[0];
////                 String fileOut = args[1];
////                 String tmpdir = args[2];
////                 try {
////                 this.rewritePlatform(fileIn, fileOut, tmpdir);
////                 } catch (IOException ex) {
////                 Logger.getLogger(Resequencing.class.getName()).log(Level.SEVERE, null, ex);
////                 }
////                 }
////                 */
////            } catch (IOException ex) {
////                Logger.getLogger(Resequencing.class.getName()).log(Level.SEVERE, null, ex);
////            }
////        }
//			} catch (IOException ex) {
//				Logger.getLogger(Resequencing.class.getName()).log(Level.SEVERE, null, ex);
//			}
		}
	}


	public void overlapFamWithSequencingIDs(String sampleList, String famFile) throws IOException {

	}

//	private void combine(String[] args) {
//		cli = new CLI(args);
//
//		// stuff for merger
//		// get list of samples to use
//		try {
//			CoverageMatrix mat = new CoverageMatrix(cli.getCoverageMatrix());
//
//			String indir = cli.getBamFileIndir();
//			String outdir = cli.getBamFileOutdir();
//
//			indir = Gpio.formatAsDirectory(indir);
//			outdir = Gpio.formatAsDirectory(outdir);
//			Gpio.createDir(outdir);
//			String tmpdir = outdir + "tmp";
//			tmpdir = Gpio.formatAsDirectory(tmpdir);
//			Gpio.createDir(tmpdir);
//			String suffix = cli.getBamSuffix();
//
//			ArrayList<Pair<String, ArrayList<String>>> samples = mat.selectSamples(cli.getCoverageThresholdColumn(), cli.getCoverageThreshold(), outdir);
//			System.out.println(samples.size() + " samples selected");
//			TextFile outfile = new TextFile(outdir + "SamplesToUse.txt", TextFile.W);
//
//			for (Pair<String, ArrayList<String>> p : samples) {
//
//				String sampleName = p.getLeft();
//				ArrayList<String> dups = p.getRight();
//
//				if (dups.size() > 1) {
//					// check whether files are actually there..
//					int d = 0;
//					ArrayList<File> finalFiles = new ArrayList<File>();
//					System.out.println(sampleName + "\thas " + dups.size() + " duplicates");
//					for (String f : dups) {
//						if (Gpio.exists(indir + f + "." + suffix)) {
//							d++;
//							finalFiles.add(new File(indir + f + "." + suffix));
//						} else {
//							System.out.println(sampleName + "\tFile not found: " + indir + f + "." + suffix);
//						}
//					}
//					System.out.println(sampleName + "\tFinal duplication size: " + finalFiles.size());
//					if (finalFiles.size() > 1) {
//
//						// merge bam files
//						File outputFile = new File(outdir + sampleName + "_merged.bam");
//						System.out.println(sampleName + "\tMerging BAM files");
//						mergeBamFiles(sampleName, finalFiles, outputFile);
//						// dedup
//
//						File deduppedOut = new File(outdir + sampleName + "_merged.dedup.bam");
//						File deduppedMetricsOut = new File(outdir + sampleName + "_merged.dedupmetrics.txt");
//						System.out.println(sampleName + "\tDeduplicating alignments in BAM file" + outputFile.getAbsolutePath());
//						System.out.println(sampleName + "\tWriting to: " + deduppedOut.getAbsolutePath());
//						MarkDups md = new MarkDups(outputFile, deduppedOut, deduppedMetricsOut, new File(tmpdir));
//// md.go();
//						// sort file
//						File deduppedSortOut = new File(outdir + sampleName + "_merged.dedup.sort.bam");
//						SortBAM sb = new SortBAM(deduppedOut, deduppedSortOut, new File(tmpdir));
//
//						// index
//						indexBAM(deduppedSortOut);
//						outfile.writeln(deduppedSortOut.getAbsolutePath());
//					} else {
//						// print filename to file
//						if (finalFiles.size() == 1) {
//							// check whether the file has been indexed..
//
//							if (!Gpio.exists(finalFiles.get(0).getAbsolutePath() + ".bai")) {
//								System.out.println(sampleName + "\tIndexing BAM: " + finalFiles.get(0).getAbsolutePath());
//								indexBAM(finalFiles.get(0));
//							}
//							outfile.writeln(finalFiles.get(0).getAbsolutePath());
//						}
//					}
//
//				} else {
//					// check whether file is there,
//					String f = dups.get(0);
//					System.out.println(sampleName + "\thas no duplicates.");
//					if (Gpio.exists(indir + f + "." + suffix)) {
//						outfile.writeln(indir + f + "." + suffix);
//						// check whether file has been indexed...
//						if (!Gpio.exists(indir + f + "." + suffix + ".bai")) {
//							System.out.println(sampleName + "\tIndexing BAM: " + indir + f + "." + suffix);
//							indexBAM(new File(indir + f + "." + suffix));
//						}
//						System.out.println(sampleName + "\tFound file: " + indir + f + "." + suffix);
//					} else {
//						System.out.println(sampleName + "\tCould not find file: " + indir + f + "." + suffix);
//					}
//				}
//
//			}
//			outfile.close();
//
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//
//	}

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
			SAMSequenceDictionary dict = new SAMSequenceDictionary(sequenceList);
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
		ArrayList<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();
		filters.add(new NotPrimaryAlignmentFilter());
		filters.add(new FailsVendorQualityCheckFilter());
		filters.add(new DuplicateReadFilter());
		filters.add(new UnmappedReadFilter());
		filters.add(new MappingQualityUnavailableFilter());


		AggregateFilter filter = new AggregateFilter(filters);
		while (iterator.hasNext()) {
			SAMRecord record = iterator.next();
			String reference = record.getReferenceName();

			if (!filter.filterOut(record)) {

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

	private void indexDedupSortIndex(String f2, String tmpdir) {

//        String[] filesInDir = Gpio.getListOfFiles(fileIn);
//        for (String f : filesInDir) {
//            String path = fileIn + f;
//            System.out.println("f: "+path);
//            if (Gpio.isDir(path)) {
//                String[] filesInDir2 = Gpio.getListOfFiles(path);
//                for (String f2 : filesInDir2) {
//                    f2 = Gpio.formatAsDirectory(path) + f2;
		System.out.println(f2);
		if (f2.endsWith(".bam") || f2.endsWith(".sam")) {
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
//
//                }
//            }
//        }
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
			BamFileReader reader = new BamFileReader(new File(fileName));

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
			BamFileReader reader = new BamFileReader(new File(fileName));

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

	private void determineCoverageForRegions(String regionsBedFile3, String bamfile, int readLength, String gtf, String output) throws IOException {
		int doubleReadLength = readLength * 2;
		TextFile tf = new TextFile(regionsBedFile3, TextFile.R);
		ArrayList<Feature> features = new ArrayList<Feature>();
		String[] header = tf.readLineElems(TextFile.tab); // if any
		String[] ln = tf.readLineElems(TextFile.tab);


		GTFAnnotation annot = new GTFAnnotation(gtf);
		TreeSet<Gene> annotationTree = annot.getGeneTree();


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
//        filters.add(new NotPrimaryAlignmentFilter());
		filters.add(new FailsVendorReadQualityFilter());
		filters.add(new UnmappedReadFilter());
		filters.add(new MappingQualityUnavailableFilter());
//        filters.add(new HCMappingQualityFilter());
		// filters.add(new MalformedReadFilter());

		// add read mate should map to same strand
//        filters.add(new MateSameStrandFilter());


		BamFileReader reader = new BamFileReader(new File(bamfile));
		List<SAMReadGroupRecord> readGroups = reader.getReadGroups();

		// map readgroups to samplenames
		HashMap<SAMReadGroupRecord, Integer> readgroupMap = new HashMap<SAMReadGroupRecord, Integer>();
		String[] samples = null;

		if (readGroups.size() > 0) {
			samples = new String[readGroups.size()];
			int rgctr = 0;
			for (SAMReadGroupRecord r : readGroups) {
				readgroupMap.put(r, rgctr);
				samples[rgctr] = r.getSample();
				rgctr++;
			}
		} else {
			System.err.println("WARNING: no readgroups found");
			samples = new String[1];
			samples[0] = bamfile;

		}

//        TextFile tf1 = new TextFile(output + "LocusCoveragePerBP.txt", TextFile.W);
//        TextFile tf2 = new TextFile(output + "LocusCoveragePerSample.txt", TextFile.W);
		for (Feature region : features) {

			Chromosome featureChr = region.getChromosome();
			int featureStart = region.getStart();
			int featureStop = region.getStop();
			int featureSize = featureStop - featureStart;

			SAMRecordIterator iterator = reader.query("chr" + featureChr.getNumber(), featureStart - doubleReadLength, featureStop + doubleReadLength, false);

			int[][] coverageMapNegStrand = new int[samples.length][featureSize];
//            int[][] coverageMapPosStrand = new int[samples.length][featureSize];

			int basequalthreshold = 0;
			int nrRecords = 0;
			int nrFilteredRecords = 0;
			while (iterator.hasNext()) {

				SAMRecord record = iterator.next();
				nrRecords++;
				// record passes filter
//				if (!f.filterOut(record)) {
				if (true) {
					nrFilteredRecords++;
					// do stuff to the record.

					// get the records readgroup
					Integer rgId = null;
					if (readGroups.size() > 1) {
						SAMReadGroupRecord readGroup = record.getReadGroup();
						rgId = readgroupMap.get(readGroup);
					} else {
						rgId = 0;
					}

					// map it to the feature coverage map
					int[][] toMapTo = coverageMapNegStrand;
//                    if (record.getReadNegativeStrandFlag()) {
//                        toMapTo = coverageMapNegStrand;
//                    }

					//    -- if paired end, check whether other pair is present
					Cigar cigar = record.getCigar();
					List<CigarElement> cigarElements = cigar.getCigarElements();

					int mapPos = record.getAlignmentStart();

					int windowRelativePosition = mapPos - featureStart; // left most pos
					int readPosition = 0; // relative position in read/fragment

					byte[] baseQual = record.getBaseQualities();
					byte[] bases = record.getReadBases();

					// main loop
//                    System.out.println(nrFilteredRecords + "\t" + record.getCigarString());
//                    System.out.println("mappos: " + mapPos);
//                    System.out.println("feats: " + featureStart);
//                    System.out.println("diff: " + (mapPos - featureStart));
//                    System.out.println("W: " + windowRelativePosition);
					for (CigarElement e : cigarElements) {
						int cigarElementLength = e.getLength();

						switch (e.getOperator()) {
							case H: // hard clip
								break;
							case P: // padding
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
//                                System.out.println("M: " + cigarElementLength);
//                                System.out.println("r: " + readPosition);
								int endPosition = readPosition + cigarElementLength;
//                                System.out.println("e: " + endPosition);
//                                System.out.println("w:" + windowRelativePosition);
								for (int pos = readPosition; pos < endPosition; pos++) {
									byte base = 0;
									if (bases.length > pos) {
										base = bases[pos];

										if (windowRelativePosition >= 0 && windowRelativePosition < featureSize) { // the read could overlap the leftmost edge of this window
											//    -- for each base pos: check basequal
//                                        if (windowRelativePosition > 685) {
////                                            System.out.println(baseQual[readPosition] + "\t" + base);
//
//                                        }
											if (baseQual[readPosition] > basequalthreshold) {

												//    -- determine number of A/T/C/G/N bases
												if (base == 65 || base == 97) {
													toMapTo[rgId][windowRelativePosition]++;
												} else if (base == 67 || base == 99) {
													toMapTo[rgId][windowRelativePosition]++;
												} else if (base == 71 || base == 103) {
													toMapTo[rgId][windowRelativePosition]++;
												} else if (base == 84 || base == 116) {
													toMapTo[rgId][windowRelativePosition]++;
												} else {
													System.err.println("unparsed base: " + base);
												}
											}
										}
										windowRelativePosition++;
									} else {
										System.err.println("Trying to access base: " + pos + " but length == " + bases.length);
										System.err.println(record.toString());
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
			iterator.close();

			System.out.println(region.getChromosome().getName() + "-" + region.getStart() + ":" + region.getStop());
			System.out.println("Fragments: " + nrRecords);
			System.out.println("Fragments filtered: " + nrFilteredRecords);

			// now plot the coverage per sample
			int pageMargin = 50;
			int betweenPlotMargin = 25;
			int plotIndividualWidth = featureSize;
			int plotIndividualHeight = featureSize / 4;


			NavigableSet<Gene> s = annot.getGeneTree().subSet(new Gene("", featureChr, Strand.POS, featureStart, featureStart), true,
					new Gene("", featureChr, Strand.NEG, featureStop, featureStop), true);
			Gene[] genes = s.toArray(new Gene[0]);


			int pageHeight = (pageMargin * 2) + (betweenPlotMargin * (samples.length - 1)) + (plotIndividualHeight * samples.length) + (genes.length * 25);
			int pageWidth = pageMargin * 2 + plotIndividualWidth;

			String name = output + "Locus-" + region.getChromosome().getName() + "-" + region.getStart() + ":" + region.getStop() + ".pdf";
			try {


				LocusCoveragePlot plot = new LocusCoveragePlot(name, pageWidth, pageHeight);

				plot.setMargin(pageMargin);
				plot.setBetweenPlotMargin(betweenPlotMargin);
				plot.setPlotIndividualWidth(plotIndividualWidth);
				plot.setPlotIndividualHeight(plotIndividualHeight);

				plot.setRowLabels(samples);
				plot.setCoverageData(coverageMapNegStrand, null);
				plot.setFeature(region);
				plot.setGenes(genes);


				// plot.setVariants();
				plot.setMaxCoverage(10);
				plot.draw();

				plot.close();
			} catch (DocumentException e) {
				e.printStackTrace();
			}

//            for (int i = 0; i < featureSize; i++) {
//                System.out.println(i + "\t" + coverageMapNegStrand[0][i]);
//            }
//            // write average coverage to text file
//            // locus coverage per bp
//            int[] coveragePerBp = new int[featureSize];
//            for (int pos = doubleReadLength; pos < doubleReadLength + featureSize; pos++) {
//                int sum1 = 0;
//                for (int s = 0; s < coverageMapPosStrand.length; s++) {
//                    sum1 += coverageMapNegStrand[s][pos] + coverageMapPosStrand[s][pos];
//                }
//
//                sum1 = (int) Math.ceil((double) sum1 / coverageMapNegStrand.length);
//                coveragePerBp[pos - doubleReadLength] = sum1;
//            }
//
//            int[] coveragePerLocus = new int[readGroups.size()];
		}

		reader.close();

	}

	private void vcfMerge(String dir, String outputdir) throws IOException {

		String[] files = new String[25];
		for (int i = 1; i < 23; i++) {
			files[i - 1] = dir + "/" + i + ".vcf";
		}
		files[22] = dir + "/X.vcf";
		files[23] = dir + "/Y.vcf";
		files[24] = dir + "/MT.vcf";

//        Arrays.sort(files, new AlphaNumericComparator());
		for (String f : files) {
			System.out.println(f);
		}

		System.out.println("Found " + files.length + " VCF files in your dir: " + dir);

		ArrayList<String> sampleNames = new ArrayList<String>();
		HashMap<String, Integer> sampleToId = new HashMap<String, Integer>();

		// index the samples
		HashSet<String> headerLines = new HashSet<String>();

		TextFile merged = new TextFile(outputdir + "merged.vcf", TextFile.W);

		int q = 0;
		for (String vcffile : files) {
			TextFile tf = new TextFile(vcffile, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);

			while (elems != null) {

				if (elems[0].startsWith("##")) {
					// superheader
					if (!elems[0].startsWith("##GATKCommandLine=")) { // don't need this.
						String ln = Strings.concat(elems, Strings.tab);
						if (!headerLines.contains(ln)) {
							merged.writeln(ln);
							headerLines.add(ln);
						}

					}
				} else if (elems[0].startsWith("#CHROM")) {
					// header
					// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1

					if (q != 0) {
						if (elems.length - 9 != sampleNames.size()) {
							System.err.println("ERROR: not same number of samples in vcf: " + vcffile);
							System.exit(-1);
						}
					}

					for (int i = 9; i < elems.length; i++) {
						String sample = elems[i];
						if (!sampleToId.containsKey(sample)) {
							if (q != 0) {
								System.out.println("Don't support non-shared samples at this time.");
								System.out.println("New sample detected: " + sample + " in file: " + vcffile);
								System.exit(-1);
							}
							sampleNames.add(sample);
							sampleToId.put(sample, sampleToId.size());
						}
					}
				} else {
					break;
				}
				elems = tf.readLineElems(TextFile.tab);
			}

			tf.close();

			q++;
			System.out.println("Total samples: " + sampleToId.size() + " after file: " + vcffile);
		}

		// write the rest of the header.
		String headerLn = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
		for (String sampleName : sampleNames) {
			headerLn += "\t" + sampleName;
		}
		merged.writeln(headerLn);

		for (String vcffile : files) {
			TextFile tf = new TextFile(vcffile, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);

			int[] sampleToNewSample = new int[sampleNames.size()];
			for (int i = 0; i < sampleToNewSample.length; i++) {
				sampleToNewSample[i] = -9;
			}

			HashSet<String> samples = new HashSet<String>();
			while (elems != null) {

				if (elems[0].startsWith("##")) {
					// superheader // skip it for now
				} else if (elems[0].startsWith("#CHROM")) {
					// header
					// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1
					for (int i = 9; i < elems.length; i++) {
						String sample = elems[i];

						// index
						Integer newIndex = sampleToId.get(sample);
						sampleToNewSample[newIndex] = i;
						samples.add(sample);
					}
				} else {
					if (samples.size() != sampleNames.size()) {
						System.err.println("ERROR: not same number of samples: ");
						System.exit(-1);
					}
					// must be something else.
					StringBuilder builder = new StringBuilder();
					for (int i = 0; i < 9; i++) {
						builder.append(elems[i]);
						if (i < 8) {
							builder.append("\t");
						}
					}
					for (int i : sampleToNewSample) {
						if (i < 0) {
							System.err.println("Error in file - missing sample:" + vcffile);
							System.exit(-1);
						}
						builder.append("\t");
						builder.append(elems[i]);
					}
					merged.writeln(builder.toString());
				}
				elems = tf.readLineElems(TextFile.tab);
			}

			tf.close();

		}
		merged.close();

	}

	public void summarizeVCF(String vcffile, String outputFileLoc) throws IOException {

		TextFile tf = new TextFile(vcffile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);

		String[] samples = null;

		HashMap<String, Integer> infoToColumn = new HashMap<String, Integer>();

		int novelSNPs = 0;
		int novelIndels = 0;
		int knownIndels = 0;
		int knownSNPs = 0;

		int nrVariants = 0;
		ArrayList<String> columnNames = new ArrayList<String>();

		while (elems != null) {

			if (elems[0].startsWith("##")) {
				// superheader

			} else if (elems[0].startsWith("#CHROM")) {
				// header
				// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1   SAMPLE2
			} else {
				// line
				String[] infoColumns = elems[7].split(";");
				for (String infoColumn : infoColumns) {
					// column name separated by  =
					String[] infocolumnelems = infoColumn.split("=");
					String infocolumnName = infocolumnelems[0];
					if (!infoToColumn.containsKey(infocolumnName)) {
						infoToColumn.put(infocolumnName, infoToColumn.size());
						columnNames.add(infocolumnName);
					}
				}

				nrVariants++;
			}
			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();


		System.out.println(nrVariants + " variants and  " + infoToColumn.size() + " pieces of info per variant");

		double[][] values = new double[nrVariants][infoToColumn.size()];
		boolean[] biAllelic = new boolean[nrVariants];
		boolean[] isKnown = new boolean[nrVariants];

		int variant = 0;
		tf.open();
		elems = tf.readLineElems(TextFile.tab);
		TextFile outputFileWriter = new TextFile(outputFileLoc, TextFile.W);
		String header = "Chrom\tPos\tId\tRef\tAlt\tQual\tFilter\tType\tisBiallelic\tCallRate\tHWEP\tMAF\tminorAllele\tnrAlleles\t" + Strings.concat(columnNames, Strings.tab);

		outputFileWriter.writeln(header);
		while (elems != null) {

			Chromosome chr = Chromosome.NA;
			if (elems[0].startsWith("##")) {
				// superheader

			} else if (elems[0].startsWith("#CHROM")) {
				// header
				// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1   SAMPLE2
				samples = new String[elems.length - 9];

				System.arraycopy(elems, 9, samples, 0, elems.length - 9);

			} else {
				// line
				String ref = elems[3];
				String alt = elems[4];
				chr = Chromosome.parseChr(elems[0]);
				String[] alternates = alt.split(",");
				String[] alleles = new String[1 + alternates.length];
				alleles[0] = ref.intern();
				for (int i = 0; i < alternates.length; i++) {
					alleles[i + 1] = alternates[i].intern();
				}


				VCFVariantType type = VCFVariantType.parseType(ref, alt);
				biAllelic[variant] = type.isBiallelic();
				String[] infoColumns = elems[7].split(";");
				for (String infoColumn : infoColumns) {
					// column name separated by  =
					String[] infocolumnelems = infoColumn.split("=");
					String infocolumnName = infocolumnelems[0];
					if (infocolumnelems.length == 1) {
						Integer id = infoToColumn.get(infocolumnName);
						values[variant][id] = 1d;

					} else {
						String value = infocolumnelems[1];

						Integer id = infoToColumn.get(infocolumnName);

						if (infocolumnName.equals("culprit")) {
							Integer id2 = infoToColumn.get(value);
							values[variant][id] = id2;
						} else {

							try {
								values[variant][id] = Double.parseDouble(value);
							} catch (NumberFormatException e) {
								String[] valueElems = value.split(",");
								if (valueElems.length > 1 && infocolumnName.equals("MLEAC") ||
										infocolumnName.equals("MLEAF") ||
										infocolumnName.equals("AC") ||
										infocolumnName.equals("AF")) {
									values[variant][id] = -1; // now we know that these are multi-allelic
								} else {
									System.out.println("Error: could not parse " + value + " for infoElem: " + infocolumnName + " and variant: " + elems[0] + "\t" + elems[1] + "\t" + elems[2]);
								}

							}

						}
					}
				}

				String[] format = elems[8].split(":");
				int gtCol = -1; // genotype
				int adCol = -1; // Allelic depths for the ref and alt alleles in the order listed
				int dpCol = -1; // Approximate read depth (reads with MQ=255 or with bad mates are filtered)
				int gqCol = -1; // Genotype Quality
				int plCol = -1; // Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification
				for (int c = 0; c < format.length; c++) {
					if (format[c].equals("GT")) {
						gtCol = c;
					} else if (format[c].equals("AD")) {
						adCol = c;
					} else if (format[c].equals("DP")) {
						dpCol = c;
					} else if (format[c].equals("GQ")) {
						gqCol = c;
					} else if (format[c].equals("PL")) {
						plCol = c;
					}
				}

				byte[] sampleGenotypes = new byte[elems.length - 9];

				byte[][] sampleAlleles = new byte[elems.length - 9][2];

				int nrHets = 0;
				int nrHomAA = 0;
				int nrHomBB = 0;

				int nrTotal = 0;
				int nrCalled = 0;
				for (int e = 9; e < elems.length; e++) {
					int samplePos = e - 9;
					nrTotal += 2;
					String sampleColumn = elems[e];
					if (sampleColumn.equals("./.")) {
						// not called
						sampleGenotypes[samplePos] = -1;
					} else {

						String[] sampleElems = sampleColumn.split(":");

						if (gtCol != -1) {
							String gt = sampleElems[gtCol];
							String[] gtElems = gt.split("/");
							String allele1 = gtElems[0];
							String allele2 = gtElems[1];
							byte allele1b = -1;
							byte allele2b = -1;
							if (allele1.equals(".")) {
								// missing
								sampleAlleles[samplePos][0] = -1;
								System.out.println("allele 1 mising");
							} else {
								nrCalled++;
								allele1b = Byte.parseByte(allele1);

							}
							if (allele2.equals(".")) {
								// missing
								sampleAlleles[samplePos][1] = -1;
								System.out.println("allele 2 mising");
							} else {
								nrCalled++;
								allele2b = Byte.parseByte(allele2);
							}


							if (type.isBiallelic()) {
								// this should change for chr X: there we should only count alleles for females
								if (allele1b != -1 && allele2b != -1) {
									if (allele1b == allele2b) {
										sampleGenotypes[samplePos] = allele1b; // homozygote AA
										if (allele1b == 1) {
											sampleGenotypes[samplePos]++; // homozygote BB
											nrHomBB++;
										} else {
											nrHomAA++;
										}
									} else {
										sampleGenotypes[samplePos] = 1; // heterozygote AB
										nrHets++;
									}

								}

							} else {
								// leave it for now..
								sampleGenotypes[samplePos] = -1;
							}

						}

						if (adCol != -1) {
							// Allelic depths for the ref and alt alleles in the order listed
							String ad = sampleElems[adCol];
							String[] adElems = ad.split(",");

						}

						if (dpCol != -1) {
							// Approximate read depth (reads with MQ=255 or with bad mates are filtered)
							String dp = sampleElems[dpCol];
							String[] dpElems = dp.split(",");
						}

						if (gqCol != -1) {
							// Genotype Quality
							String gq = sampleElems[gqCol];
							String[] gqElems = gq.split(",");
						}

						if (plCol != -1) {
							// Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification
							String pl = sampleElems[plCol];
							String[] plElems = pl.split(",");
						}


					}
				}

				// determine MAF, HWEP, CR
				double maf = 0;
				double hwep = 1;
				double cr = (double) nrCalled / nrTotal;
				String minorAllele = alleles[0];
				if (type.isBiallelic()) {
					int freqAlleleA = 2 * nrHomAA + nrHets; // genotypeFreq[0] + genotypeFreq[1];
					int freqAlleleB = 2 * nrHomBB + nrHets; // genotypeFreq[2] + genotypeFreq[1];

					maf = (double) (freqAlleleA) / (nrCalled);
					minorAllele = alleles[0];
					if (freqAlleleA > freqAlleleB) {
						minorAllele = alleles[1];
						maf = 1 - maf;
					}
					hwep = HWE.calculateExactHWEPValue(nrHets, nrHomAA, nrHomBB);

				}

				// write the info elems to disk
				outputFileWriter.append(elems[0] + "\t" +
						elems[1] + "\t" +
						elems[2] + "\t" +
						elems[3] + "\t" +
						elems[4] + "\t" +
						elems[5] + "\t" +
						elems[6]);
				outputFileWriter.append("\t");
				outputFileWriter.append(type.toString());
				outputFileWriter.append("\t");
				outputFileWriter.append("" + type.isBiallelic());
				outputFileWriter.append("\t");
				outputFileWriter.append("" + cr);
				outputFileWriter.append("\t");
				outputFileWriter.append("" + hwep);
				outputFileWriter.append("\t");
				outputFileWriter.append("" + maf);
				outputFileWriter.append("\t");
				outputFileWriter.append("" + minorAllele);
				outputFileWriter.append("\t");
				outputFileWriter.append("" + alleles.length);

				for (int i = 0; i < values[variant].length; i++) {
					outputFileWriter.append("\t");
					outputFileWriter.append("" + values[variant][i]);

				}
				outputFileWriter.append("\n");

				variant++;
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		outputFileWriter.close();

		//

	}

	public void determineReadsAndDuplicationsPerChr(String bamfile, String outputfile) throws IOException {

		int nrReadsUnmapped = 0;
		int nrWithMateUnmapped = 0;

		int total = 0;

		int[] nrMarkedDuplicatesPerChr = new int[27];
		int[] nrReadsMappedPerChr = new int[27];

		int[] nrFragmentsMappedPerChr = new int[27];
		int[] nrFragmentsMappedWithMatesOnOppositeStrandsPerChr = new int[27];
		int[] nrUniqueReadsWMapQ30PerChr = new int[27];

		BamFileReader reader = new BamFileReader(bamfile);
		SAMRecordIterator iterator = reader.iterator();

		int nrReadsProcessed = 0;
		SAMRecord record = iterator.next();
		while (iterator.hasNext()) {

			total++;

			if (record.getReadUnmappedFlag()) {
				nrReadsUnmapped++;
			} else {
				Chromosome c = Chromosome.parseChr(record.getReferenceName());
				boolean dup = record.getDuplicateReadFlag();
				int chr = c.getNumber();
				if (chr > 25) {
					chr = 26;
				}
				nrReadsMappedPerChr[chr]++;
				if (dup) {
					nrMarkedDuplicatesPerChr[chr]++;
				} else {

					int mapq = record.getMappingQuality();
					if (mapq > 30) {
						nrUniqueReadsWMapQ30PerChr[chr]++;
					}
				}
			}


			if (record.getFirstOfPairFlag()) {
				if (record.getMateUnmappedFlag() || record.getReadUnmappedFlag()) {
					// unpaired mate
					nrWithMateUnmapped++;
				} else {
					Chromosome c = Chromosome.parseChr(record.getReferenceName());
					Chromosome c2 = Chromosome.parseChr(record.getMateReferenceName());


					boolean negStr = record.getReadNegativeStrandFlag();
					boolean negStrMate = record.getMateNegativeStrandFlag();
					if (c.equals(c2)) {
						// both mate pairs on same chr

						int chr = c.getNumber();
						if (chr > 25) {
							chr = 26;
						}

						if ((negStr && !negStrMate) || (!negStr && negStrMate)) {
							nrFragmentsMappedWithMatesOnOppositeStrandsPerChr[chr]++;
						}
						nrFragmentsMappedPerChr[chr]++;
					}
				}
			}
			record = iterator.next();
			nrReadsProcessed++;
			if (nrReadsProcessed % 1000000 == 0) {
				System.out.println(nrReadsProcessed + " reads processed");
			}
		}

		TextFile outfile = new TextFile(outputfile, TextFile.W);

		outfile.writeln("Chr" +
				"\tChrSize" +
				"\tnrReadsMappedPerChr" +
				"\tnrReadsMarkedDuplicatesPerChr" +
				"\tnrUniqueReadsWMapQ30PerChr" +
				"\tnrFragmentsMappedPerChr" +
				"\tnrFragmentsMappedWithMatesOnOppositeStrandsPerChr");

		for (int i = 1; i < nrFragmentsMappedPerChr.length; i++) {
			Chromosome c = Chromosome.parseChr("" + i);
			if (i == 23) {
				c = Chromosome.parseChr("ChrX");
			} else if (i == 24) {
				c = Chromosome.parseChr("ChrY");
			} else if (i == 25) {
				c = Chromosome.parseChr("MT");
			}
			outfile.writeln(c.getName() +
							"\t" + c.getLength() +
							"\t" + nrReadsMappedPerChr[i] +
							"\t" + nrMarkedDuplicatesPerChr[i] +
							"\t" + nrUniqueReadsWMapQ30PerChr[i] +
							"\t" + nrFragmentsMappedPerChr[i] +
							"\t" + nrFragmentsMappedWithMatesOnOppositeStrandsPerChr[i]
			);
		}

		outfile.writeln("total nr reads:\t" + total);
		outfile.writeln("total nr reads unmapped:\t" + nrReadsUnmapped);
		outfile.writeln("total nr reads mate unmapped:\t" + nrWithMateUnmapped);
		outfile.close();
	}

}
