package nl.harmjanwestra.ngs;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import nl.harmjanwestra.utilities.bamfile.BamFileReader;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by hwestra on 6/3/15.
 */
public class QuickRPKM {

	public static void main(String[] args) {
		QuickRPKM q = new QuickRPKM();

		if (args.length < 3) {
			System.out.println("Usage: region bam out [mapq]");
			System.exit(0);
		}


		String targetregions = args[0];
		String bamfile = args[1];
		String outdir = args[2];
		Integer mapq = 10;
		if (args.length > 3) {
			mapq = Integer.parseInt(args[3]);
		}

		try {
			q.run(targetregions, bamfile, outdir, mapq);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void run(String targetregions, String bamfile, String outdir, int mapqthreshold) throws Exception {

		if (!Gpio.exists(outdir)) {
			Gpio.createDir(outdir);
		}

		Coverage cov = new Coverage();

		ArrayList<Feature> regions = cov.loadRegions(targetregions);


		BamFileReader reader = new BamFileReader(new File(bamfile));

		CoverageTask ct = new CoverageTask(bamfile, outdir, regions, true);
		ct.call();
		HashMap<String, String> chromosomeToSequence = ct.matchChromosomeNames(reader);

		// count total number of mapping reads..
		int totalReadsMapped = 0;
		SAMRecordIterator it = reader.iterator();
		int readcounter = 0;
		while (it.hasNext()) {

			SAMRecord record = it.next();

			readcounter++;
			if (readcounter % 100000 == 0) {
				System.out.println(readcounter + " reads processed.");
			}
			// quality checks..
			if (!record.getReadUnmappedFlag()) {
//				if (!record.getDuplicateReadFlag()) {
					if (record.getProperPairFlag()) {
						if (!record.isSecondaryOrSupplementary()) {
							if (!record.getNotPrimaryAlignmentFlag()) {
								if (record.getMappingQuality() >= mapqthreshold) { // according to Maria, 10 is OK.
									totalReadsMapped++;
								}
							}
						}
					}
//				}
			}
		}

		it.close();

		TextFile tf = new TextFile(outdir + "rpkms.txt", TextFile.W);

		tf.writeln("region\tnrReadsMappedInRegion\ttotalNrOfMappedReadsInBam\tRPKM");
		for (Feature f : regions) {
			Chromosome c = f.getChromosome();
			int start = f.getStart();
			int stop = f.getStop();
			int featuresize = Math.abs(stop - start);
			String chrName = chromosomeToSequence.get(c.getName());
			int nrReadsMappedInRegion = 0;
			if (chrName != null) {
				it = reader.query(chrName, start, stop, true);
				while (it.hasNext()) {

					SAMRecord record = it.next();

					// quality checks..
					if (!record.getReadUnmappedFlag()) {
//						if (!record.getDuplicateReadFlag()) {
							if (record.getProperPairFlag()) {
								if (!record.isSecondaryOrSupplementary()) {
									if (!record.getNotPrimaryAlignmentFlag()) {
										if (record.getMappingQuality() >= mapqthreshold) { // according to Maria, 10 is OK.
											nrReadsMappedInRegion++;
										}
									}
								}
							}
//						}
					}

				}
				it.close();
				// RPKM = (10^9 * C)/(N * L)
				double rpkm = (Math.pow(10, 9) * nrReadsMappedInRegion) / ((long) featuresize * (long) totalReadsMapped);

				tf.writeln(f.toString() + "\t" + nrReadsMappedInRegion + "\t" + totalReadsMapped + "\t" + featuresize + "\t" + rpkm + "\t" + (Math.pow(10, 9) * nrReadsMappedInRegion));


			}
		}
		tf.close();
		reader.close();

	}

}
