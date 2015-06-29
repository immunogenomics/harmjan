/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.ngs;

import htsjdk.samtools.*;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import nl.harmjanwestra.ngs.containers.ReadGroup;
import nl.harmjanwestra.ngs.wrappers.MarkDups;
import nl.harmjanwestra.ngs.wrappers.SortBAM;
import nl.harmjanwestra.utilities.bamfile.BamFileReader;
import nl.harmjanwestra.utilities.features.Chromosome;
import org.broadinstitute.gatk.engine.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.gatk.engine.filters.MappingQualityUnavailableFilter;
import org.broadinstitute.gatk.engine.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.gatk.engine.filters.UnmappedReadFilter;
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
		// Resequencing s = new Resequencing(args);
		Resequencing s = new Resequencing();

//		s.indexBAM(in);

//		s.mergeBamFiles(sample,files,outfile);
//		s.determineReadsAndDuplicationsPerChr(bamfile,outputfile);
//		s.overlapFamWithSequencingIDs(samplelist,famfile);

//		s.splitPerChromosomeRemoveDuplicateReadsAndQualityScores(filein,outdir,stripoldqualmetrics);
//			s.replaceReadGroupDedupAndSort(bamin, bamout, tmp, rg);
//		try {
//			s.readsPerChromosomePerReadGroup(args);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

		try {
//			String filein = args[0];
//			String tmp = args[1];
//			s.indexDedupSortIndex(filein,tmp);
			if(args.length <4){
				System.out.println("Usage: bamin bamout tmpdir readgroup");
			} else {
				String bamin = args[0];
				String bamout = args[1];
				String tmpdir = args[2];
				String rg = args[3];
				s.replaceReadGroupDedupAndSort(bamin, bamout, tmpdir, rg);
			}
//			s.rewritePlatform(filein, tmp);
		} catch (IOException e) {
			e.printStackTrace();
		}

//		s.indexDedupSortIndex(args[0],args[1]);
	}


	public void merge(String dirIn1, String dirIn2, String dirOut) throws IOException {
		String[] files1 = Gpio.getListOfFiles(dirIn1);
		String[] files2 = Gpio.getListOfFiles(dirIn2);

		HashSet<String> fileToInt = new HashSet<String>();
		for (int i = 0; i < files1.length; i++) {
			fileToInt.add(files1[i]);
		}

		for (String f : files2) {

			if (fileToInt.contains(f) && f.endsWith("sorted.bam")) {

				String fileout = dirOut + f;
				String[] nameElems = f.split("\\.");
				String name = nameElems[0];
				String fileIn1 = dirIn1 + f;
				String fileIn2 = dirIn2 + f;

				System.out.println("Merging sample: " + name);
				System.out.println(fileIn1);
				System.out.println(fileIn2);
				System.out.println(fileout);
				this.mergeBamFiles(name, new String[]{fileIn1, fileIn2}, null, false, fileout);
				System.out.println();


			}
		}

	}

//	private Resequencing(String[] args) {


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

//		if (args.length < 2) {
//			// System.out.println("usage: input outdir stripoldqualityscores");
//			System.out.println("usage: bamfile outfile");
//		} else {
//			try {
//				String listFile = args[0];
//				String outdir = args[1];
//				String tmpdir = args[2];
//
//				TextFile tf = new TextFile(listFile, TextFile.R);
//				String[] filesIn = tf.readAsArray();
//				tf.close();
//
//				TextFile listOut = new TextFile(outdir + "SamplesRecoded.list", TextFile.W);
//
//				for (String s : filesIn) {
//
//					String[] bamfilelems = s.split("/");
//					String bamfile = bamfilelems[bamfilelems.length - 1];
//					String[] sampleelems = bamfile.split("\\.");
//					String sample = sampleelems[0];
//					String outfile = outdir + sample + ".bam";
//
//					this.rewritePlatform(s, outfile, tmpdir);
//
//					listOut.writeln(outfile);
//
//
//				}
//
//
//				listOut.close();
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
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
//		}
//	}


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

	public String MD5(String md5) {
		try {
			java.security.MessageDigest md = java.security.MessageDigest.getInstance("MD5");
			byte[] array = md.digest(md5.getBytes());
			StringBuffer sb = new StringBuffer();
			for (int i = 0; i < array.length; ++i) {
				sb.append(Integer.toHexString((array[i] & 0xFF) | 0x100).substring(1, 3));
			}
			return sb.toString();
		} catch (java.security.NoSuchAlgorithmException e) {
		}
		return null;
	}

	private void rewritePlatform(String fileIn, String tmp) throws IOException {

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


		String randomStr = MD5(System.nanoTime() + "-bam-" + Math.random());
		String tmpfile = tmpdir + randomStr + ".bam";

		SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, new File(tmpfile));
		System.out.println("Writing to: " + tmpfile);
		SAMRecordIterator iterator = reader.iterator();
		while (iterator.hasNext()) {
			SAMRecord record = iterator.next();
			outputSam.addAlignment(record);
		}
		outputSam.close();
		reader.close();

		Gpio.moveFile(tmpfile, fileIn);

		// remove old index (if any)
		File f = new File(fileIn + ".bai");
		if (f.exists()) {
			f.delete();
		}
		indexBAM(new File(fileIn));

	}

	private void mergeBamFiles(String sample, String[] listOfFiles, String fileListStr, boolean filter, String outfileStr) throws IOException {
		ArrayList<SAMFileReader> readers = new ArrayList<SAMFileReader>();
		ArrayList<SAMFileHeader> headers = new ArrayList<SAMFileHeader>();
		File outfile = new File(outfileStr);
		String[] fileListArr = null;
		if (fileListStr != null && listOfFiles == null) {
			TextFile tf = new TextFile(fileListStr, TextFile.R);
			fileListArr = tf.readAsArray();
			tf.close();
		} else {
			fileListArr = listOfFiles;
		}

		File[] files = new File[fileListArr.length];
		for (int f = 0; f < files.length; f++) {
			files[f] = new File(fileListArr[f]);
		}

		for (File f : files) {
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
		String RGID = sample;
		String RGLB = null;
		String RGPL = null;
		String RGSM = sample;
		String RGPU = null;
		if (var.isEmpty()) {
			RGID = sample;
			RGPL = "Illumina";
			RGSM = sample;
		} else {
			for (SAMReadGroupRecord r : var) {
				System.out.println(String.format(sample + "\tDetected read group ID=%s PL=%s LB=%s SM=%s%n", r.getId(), r.getPlatform(), r.getLibrary(), r.getSample()));
			}
			RGID = sample;
			RGLB = var.get(0).getLibrary();
			if (RGLB != null && RGLB.length() == 0) {
				RGLB = null;
			}
			RGPL = var.get(0).getPlatform();
			if (RGPL != null && RGPL.length() == 0) {
				RGPL = null;
			}
			RGSM = sample;
			RGPU = var.get(0).getPlatformUnit();
			if (RGPU != null && RGPU.length() == 0) {
				RGPU = null;
			}
		}

		final SAMReadGroupRecord rg = new SAMReadGroupRecord(RGID);
		rg.setLibrary(RGLB);
		rg.setPlatform(RGPL);
		rg.setSample(RGSM);
		rg.setPlatformUnit(RGPU);
		System.out.println(String.format(sample + "\tCreated read group ID=%s PL=%s LB=%s SM=%s%n", rg.getId(), rg.getPlatform(), rg.getLibrary(), rg.getSample()));
		header.setReadGroups(Arrays.asList(rg));

		SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(headerMerger.getMergedHeader(), false, outfile);
		System.out.println(sample + "\tWriting to: " + outfile);
		int ctr = 0;
		int written = 0;
		int dups = 0;
		int mapq = 0;
		int second = 0;
		int supplementary = 0;
		int unpaired = 0;
		int unmapped = 0;
		int mateunmapped = 0;
		int improperpair = 0;

		while (iterator.hasNext()) {
			SAMRecord record = iterator.next();
			record.setAttribute(SAMTag.RG.name(), RGID);
			record.setAttribute(SAMTag.SM.name(), RGID);

			if (record.getDuplicateReadFlag()) {
				dups++;
			}
			if (record.getMappingQuality() == 0) {
				mapq++;
			}
			if (record.getReadUnmappedFlag()) {
				unmapped++;
			}
			if (record.getSupplementaryAlignmentFlag()) {
				second++;
			}
			if (record.getNotPrimaryAlignmentFlag()) {
				second++;
			}
			if (!record.getReadPairedFlag()) {
				unpaired++;
			}
			if (record.getMateUnmappedFlag()) {
				mateunmapped++;
			}
			if (!record.getProperPairFlag()) {
				improperpair++;
			}

			if (filter) {
				if (!record.getDuplicateReadFlag() &&
						record.getMappingQuality() > 0 &&
						!record.getSupplementaryAlignmentFlag() &&
						!record.getNotPrimaryAlignmentFlag() &&
						!record.getReadUnmappedFlag()) {
					if (record.getReadPairedFlag()) {
						if (!record.getMateUnmappedFlag() &&
								record.getProperPairFlag()) {
							outputSam.addAlignment(record);
							written++;
						}
					} else {
						outputSam.addAlignment(record);
						written++;
					}
				}
			} else {
				outputSam.addAlignment(record);
			}


			ctr++;
			if (ctr % 1000000 == 0) {
				System.out.println(ctr + "\trecords\t" + written + "\twritten\t" + ((double) written / ctr));
				System.out.println("dups\t" + dups);
				System.out.println("mapq\t" + mapq);
				System.out.println("secd\t" + second);
				System.out.println("supl\t" + supplementary);
				System.out.println("unma\t" + unmapped);
				System.out.println("unpa\t" + unpaired);
				System.out.println("maun\t" + mateunmapped);
				System.out.println("impr\t" + improperpair);


				System.out.println();
			}
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
			File sortout = new File(f2.substring(0, f2.length() - 4) + "-sorted.bam");

			if (!sortout.exists()) {
				SortBAM sb = new SortBAM(infile, sortout, new File(tmpdir));
				indexBAM(sortout);
			}

			// sort
			// dedup
			File dedupped = new File(f2.substring(0, f2.length() - 4) + "-sorted-dedup.bam");
			File deduppedmetrics = new File(f2.substring(0, f2.length() - 4) + "-sorted-dedupmetrics.txt");

			if (!dedupped.exists()) {
				MarkDups dupMark = new MarkDups(sortout, dedupped, deduppedmetrics, new File(tmpdir));
				indexBAM(dedupped);

			}

			sortout.delete();
			File sortoutindex = new File(f2.substring(0, f2.length() - 4) + "-sorted.bam.bai");
			sortoutindex.delete();
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

	public void readsPerChromosomePerReadGroup(String[] args) throws IOException {
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

				if (!record.getDuplicateReadFlag() && !record.getReadUnmappedFlag()) {
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
				if (chr.getNumber() > 0 && chr.getNumber() < strToSample.get(sample).readsPerChr[0].length) {
					chrOutput += "\t" + strToSample.get(sample).readsPerChr[0][chr.getNumber()] + "; " + strToSample.get(sample).readsPerChr[1][chr.getNumber()];

					chrOutput2 += "\t" + strToSample.get(sample).dupsPerChr[0][chr.getNumber()] + "; " + strToSample.get(sample).dupsPerChr[1][chr.getNumber()];
					chrOutput3 += "\t" + strToSample.get(sample).readsPairedPerChr[0][chr.getNumber()] + "; " + strToSample.get(sample).readsPairedPerChr[1][chr.getNumber()];
				}
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
				if (chr.getNumber() > 0 && chr.getNumber() < s.mapqPerChrPosStrand.length) {
					for (int i = 0; i < 61; i++) {
						chrOutput += "\t" + s.mapqPerChrPosStrand[chr.getNumber()][i] + "; " + s.mapqPerChrNegStrand[chr.getNumber()][i];
					}
					outf.writeln(chrOutput);
				}

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


			if (record.getFirstOfPairFlag() || !record.getReadPairedFlag()) {
				if (record.getReadPairedFlag() && record.getMateUnmappedFlag() || record.getReadUnmappedFlag()) {
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

	public void replaceReadGroupDedupAndSort(String bamIn, String bamOut, String tmpdir, String rg) throws IOException {
		SamReader reader = SamReaderFactory.makeDefault().open(new File(bamIn));
		SAMFileHeader header = reader.getFileHeader();
		SAMFileHeader newHeader = new SAMFileHeader();

		ArrayList<SAMReadGroupRecord> newReadGroups = new ArrayList<SAMReadGroupRecord>();
		SAMReadGroupRecord rgrecord = new SAMReadGroupRecord(rg);
		rgrecord.setPlatform("Illumina");
		rgrecord.setSample(rg);
		rgrecord.setAttribute("SM", rg);

		newReadGroups.add(rgrecord);
		newHeader.setReadGroups(newReadGroups);

		newHeader.setComments(header.getComments());
		newHeader.setProgramRecords(header.getProgramRecords());
		newHeader.setGroupOrder(header.getGroupOrder());
		newHeader.setSequenceDictionary(header.getSequenceDictionary());
		newHeader.setSortOrder(header.getSortOrder());
		newHeader.setTextHeader(header.getTextHeader());
		newHeader.setValidationErrors(header.getValidationErrors());

		SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(newHeader, false, new File(bamOut));
		SAMRecordIterator iterator = reader.iterator();
		while (iterator.hasNext()) {
			SAMRecord record = iterator.next();
			if (!record.isSecondaryOrSupplementary()) {
				record.setAttribute(SAMTag.RG.name(), rg);
				record.setAttribute(SAMTag.SM.name(), rg);
//				record.
				outputSam.addAlignment(record);
			}
		}
		iterator.close();
		outputSam.close();
		reader.close();


		indexDedupSortIndex(bamOut, tmpdir);

	}


}
