package nl.harmjanwestra.ngs;

import com.xeiam.xchart.BitmapEncoder;
import com.xeiam.xchart.Chart;
import nl.harmjanwestra.ngs.GenotypeFormats.PedAndMapFunctions;
import nl.harmjanwestra.ngs.GenotypeFormats.VCFFunctions;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.HWE;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

/**
 * Created by hwestra on 2/4/15.
 */
public class GenotypeTools {


	public void getAssocForImputedVariants(String assocout, String vcf, String filter, String imputedOnly) throws IOException {


		HashSet<Feature> vcfVariantHash = new HashSet<Feature>();
		TextFile filterf = new TextFile(filter, TextFile.R);

		String ln = filterf.readLine();
		while (ln != null) {
			String[] felems = ln.split("_");
			Feature f = new Feature();
			Chromosome chr = Chromosome.parseChr(felems[0]);
			int start = Integer.parseInt(felems[1]);
			f.setChromosome(chr);
			f.setStart(start);
			f.setStop(start);

			vcfVariantHash.add(f);
			ln = filterf.readLine();
		}


		TextFile tfin = new TextFile(assocout, TextFile.R);
		TextFile tfout = new TextFile(imputedOnly, TextFile.W);


		tfout.writeln(tfin.readLine());
		String[] elems = tfin.readLineElems(TextFile.tab);
		while (elems != null) {

			Chromosome chr = Chromosome.parseChr(elems[1]);
			int start = Integer.parseInt(elems[2]);
			Feature f = new Feature();
			f.setChromosome(chr);
			f.setStart(start);
			f.setStop(start);
			if (vcfVariantHash.contains(f)) {
				tfout.writeln(Strings.concat(elems, Strings.tab));
			}
			elems = tfin.readLineElems(TextFile.tab);

		}

		tfin.close();
		tfout.close();
	}

	public void rewritePlinkResults(String assocfile, String assocout) throws IOException {
		// SNP CHR BP      P
		TextFile tfIn = new TextFile(assocfile, TextFile.R);

		tfIn.readLine();
		String ln = tfIn.readLine();
		TextFile out = new TextFile(assocout, TextFile.W);
		out.writeln("SNP\tCHR\tBP\tP\tzscore");
		while (ln != null) {

			while (ln.contains("  ")) {
				ln = ln.replaceAll("  ", " ");
			}

			while (ln.startsWith(" ")) {
				ln = ln.substring(1, ln.length());
			}

			// # CHR                       SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
			String[] elems = ln.split(" ");
			String snp = elems[1];
			String chr = elems[0];
			String pos = elems[2];
			String p = elems[8];

			try {
				Double pd = Double.parseDouble(p);

				double z = ZScores.pToZ(pd);

				if (!chr.equals("X") && !chr.equals("23")) {
					out.writeln(snp + "\t" + chr + "\t" + pos + "\t" + p + "\t" + z);
				}


			} catch (NumberFormatException e) {

			}

			ln = tfIn.readLine();
		}

		out.close();
		tfIn.close();


	}

	public void preparePlinkDataset(String origmap, String removeFile, String outdir) throws IOException {
		PedAndMapFunctions pedAndMapFunctions = new PedAndMapFunctions();
		VCFFunctions vcfFunctions = new VCFFunctions();
		ProcessBuilder pb = null;

		String referenceVCF = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/MergedWithIC/merged-lowfreqvariantsremoved-sorted-phased-sorted.vcf";
//
		String map = origmap + ".map";
		String refmap = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC.original.map";
//			String chrupdate = "/Data/ImmunoChip/US/raci_us-raciChrNames.map";
//			pedAndMapFunctions.rewriteMapFileChromosomeNames(refmap, map, chrupdate);
//
		String rsnameref = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC-RsIDsFromT1DStudy.map";
		String rsupdate = origmap + "-raciChrNames.map";
//		pedAndMapFunctions.updateMapFileRsIdsUsingMapFile(rsnameref, map, rsupdate);


		String rsupdatededup = origmap + "-raciChrNames-dedup.map";
//		pedAndMapFunctions.deduplicateMAP(rsupdate, rsupdatededup);

		String bedout = origmap + "-raciChrNames-dedup.bed";
//		pedAndMapFunctions.rewriteMapToBed(rsupdatededup, bedout);

		// liftover

		String lifted = origmap + "-lifted.bed";
		String unlifted = origmap + "-unlifted.bed";

		pb = new ProcessBuilder("/Data/Projects/2014-FR-Reseq/ImmunoChip/liftOver", bedout,
				"/Data/Projects/2014-FR-Reseq/ImmunoChip/hg18ToHg19.over.chain.gz", lifted, unlifted);
//		run(pb);

		String hg19map = origmap + "-raciChrNames-hg19.map";

		pedAndMapFunctions.convertPostLiftOverMAP(rsupdatededup, hg19map, lifted);
//		String hg19mapdedup = origmap+"-raciChrNames-hg19-dedup.map";
//		pedAndMapFunctions.identifyDuplicatesInMap(hg19map, hg19mapdedup);

		String hg19mapupd = origmap + "-raciChrNames-hg19-updRS.map";
		String dbsnpvcf = "/Data/Projects/2014-FR-Reseq/GATKResources/dbsnp_138.b37.vcf.gz";
//		pedAndMapFunctions.updateRSNames(dbsnpvcf, hg19map, hg19mapupd);

//
		String variantSelect = origmap + "-raciChrNames-hg19-dedup-selectVariants.txt";
		String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/loci.txt";
//		pedAndMapFunctions.filterMap(hg19mapupd, regionFile, variantSelect);

		Gpio.copyFile(origmap + ".map", origmap + ".hg18map");
		Gpio.copyFile(hg19mapupd, origmap + ".map");

		// select variants from plink


		if (removeFile != null) {

			pb = new ProcessBuilder("/Data/Tools/plink-1.07-mac-intel/plink1.9",
					"--extract", "/Data/ImmunoChip/US/raci_us-raciChrNames-hg19-dedup-selectVariants.txt",
					"--file", origmap,
					"--recode", "--out", origmap + "-filtered",
					"--remove", removeFile);
		} else {
			pb = new ProcessBuilder("/Data/Tools/plink-1.07-mac-intel/plink1.9",
					"--extract", "/Data/ImmunoChip/US/raci_us-raciChrNames-hg19-dedup-selectVariants.txt",
					"--file", origmap,
					"--recode", "--out", origmap + "-filtered");
		}


//		run(pb);

		Gpio.copyFile(origmap + ".hg18map", origmap + ".map");

		String ped = origmap + "-filtered";
		String vcf = origmap + "-filtered.vcf";
		//	vcfFunctions.convertPEDToVCF(ped, vcf);

		String sortedvcf = origmap + "-filtered-sorted.vcf";
		String bashCommand = "cat " + vcf + " | /Data/Tools/vcftools/bin/vcf-sort > " + sortedvcf;
		String bashfilename = origmap + "-filter-sort.sh";
		TextFile tf = new TextFile(bashfilename, TextFile.W);
		tf.writeln("#!/bin/bash\n" + bashCommand);
		tf.close();
		pb = new ProcessBuilder("bash", bashfilename);
		run(pb);

		// next up: beagle compare

		String tmpout = outdir + "beaglecompare/";
		if (Gpio.exists(tmpout)) {
// empty directory
			String[] files = Gpio.getListOfFiles(tmpout);
			for (String f : files) {

				File file = new File(tmpout + f);
				if (file.isFile()) {
					file.delete();
				}

			}

		}
		Gpio.createDir(tmpout);

		String files = "";
		System.out.println("Making genotypes conform to reference");
		for (int i = 1; i < 23; i++) {
			System.out.println("chr " + i);
			System.out.println(sortedvcf);
			System.out.println(referenceVCF);
			pb = new ProcessBuilder("java",
					"-jar",
					"/Data/Tools/beagle/conform-gt.r1174.jar",
					"ref=" + referenceVCF,
					"gt=" + sortedvcf,
					"chrom=" + i,
					"out=" + tmpout + "chr" + i,
					"match=POS");
			run(pb);

			// gunzip and sort the resulting file
			bashfilename = tmpout + "gunzip.sh";
			String sortedvcf2 = tmpout + "chr" + i + "-sorted.vcf";
			bashCommand = "gzcat " + tmpout + "chr" + i + ".vcf | /Data/Tools/vcftools/bin/vcf-sort > " + sortedvcf2;
			tf = new TextFile(bashfilename, TextFile.W);
			tf.writeln("#!/bin/bash\n" + bashCommand);
			tf.close();
			pb = new ProcessBuilder("bash", bashfilename);
			run(pb);

			// sort the filtered vcf
			files += " " + sortedvcf;
		}

		// merge the new vcf files
		// use vcf tools merge
		String merged = tmpout + "merged.vcf";
		String mergedsorted = tmpout + "merged-sorted.vcf";
		bashfilename = tmpout + "merged.sh";
		bashCommand = "/Data/Tools/vcftools/bin/vcf-concat " + files + " > " + merged + "\n" +
				"cat " + merged + " | /Data/Tools/vcftools/bin/vcf-sort > " + mergedsorted;

		tf = new TextFile(bashfilename, TextFile.W);
		tf.writeln("#!/bin/bash\n" + bashCommand);
		tf.close();

		pb = new ProcessBuilder("bash", bashfilename);
		run(pb);

		// impute filtered
//		pb = new ProcessBuilder("java",
//				"-jar",
//				"/Data/Tools/beagle/beagle.r1399.jar",
//				"ref=" + referenceVCF,
//				"gt=" + mergedsorted,
//				"out=" + tmpout + "phasedAndImputed",
//				"nthreads=4");
//		run(pb);
//
//		// impute unfiltered
//		String sortedvcfphased = origmap + "-filtered-sorted-phasedAndImputed";
//		pb = new ProcessBuilder("java",
//				"-jar",
//				"/Data/Tools/beagle/beagle.r1399.jar",
//				"ref=" + referenceVCF,
//				"gt=" + sortedvcf,
//				"out=" + sortedvcfphased,
//				"nthreads=4");
//		run(pb);

	}

	public void sortVCF(String vcf, String sortedvcf, String bashfilename) throws IOException {


		String cat = "cat";
		if (vcf.endsWith(".gz")) {
			cat = "gzcat";
		}
		String bashCommand = cat + " " + vcf + " | /Data/Tools/vcftools/bin/vcf-sort > " + sortedvcf;


		TextFile tf = new TextFile(bashfilename, TextFile.W);
		tf.writeln("#!/bin/bash\n" + bashCommand);
		tf.close();
		ProcessBuilder pb = new ProcessBuilder("bash", bashfilename);
		run(pb);
	}

	public void concatVCF(String files, String merged, String mergedsorted, String bashfilename) throws IOException {
		String bashCommand = "/Data/Tools/vcftools/bin/vcf-concat " + files + " > " + merged + "\n" +
				"cat " + merged + " | /Data/Tools/vcftools/bin/vcf-sort > " + mergedsorted;

		TextFile tf = new TextFile(bashfilename, TextFile.W);
		tf.writeln("#!/bin/bash\n" + bashCommand);
		tf.close();

		ProcessBuilder pb = new ProcessBuilder("bash", bashfilename);
		run(pb);
	}

	public void bgzipAndIndex(String vcf, String bashfilename) throws IOException {
		String bashCommand = "/Data/Tools/tabix/tabix-0.2.6/bgzip " + vcf + "\n" +
				"/Data/Tools/tabix/tabix-0.2.6/tabix -p vcf " + vcf + ".gz";

		if (Gpio.exists( vcf + ".gz")) {
			new File( vcf + ".gz").delete();
		}

		TextFile tf = new TextFile(bashfilename, TextFile.W);
		tf.writeln("#!/bin/bash\n" + bashCommand);
		tf.close();

		ProcessBuilder pb = new ProcessBuilder("bash", bashfilename);
		run(pb);
	}

	public void mergeVCF(String files, String mergedout, String bashfilename) throws IOException {
		String bashCommand = "/Data/Tools/bcftools/bcftools-1.2/bcftools merge -m all -O v -o  " + mergedout + " " + files + "\n";

		TextFile tf = new TextFile(bashfilename, TextFile.W);
		tf.writeln("#!/bin/bash\n" + bashCommand);
		tf.close();

		ProcessBuilder pb = new ProcessBuilder("bash", bashfilename);
		run(pb);

		sortVCF(mergedout, mergedout + "sorted.vcf", bashfilename + "-sort.sh");
	}


	public void run(ProcessBuilder builder) throws IOException {
		Map<String, String> environ = builder.environment();

		final Process process = builder.start();

		java.io.SequenceInputStream is = new SequenceInputStream(process.getInputStream(), process.getErrorStream());

		InputStreamReader isr = new InputStreamReader(is);
		BufferedReader br = new BufferedReader(isr);

		String line;
		while ((line = br.readLine()) != null) {
			System.out.println(line);
		}
	}

	public void replaceSampleNames(String sequencingVCF, String sampletoSampleList, String vcfOut) throws IOException {

//		HashMap<String, String> set =
		TextFile tf = new TextFile(sampletoSampleList, TextFile.R);
		Map<String, String> select = tf.readAsHashMap(0, 1);
		tf.close();

		TextFile tfout = new TextFile(vcfOut, TextFile.W);

		TextFile in = new TextFile(sequencingVCF, TextFile.R);
		String ln = in.readLine();
		boolean[] selectSample = new boolean[0];
		int nrHeaderElems = 9;
		while (ln != null) {


			if (ln.startsWith("##")) {
				tfout.writeln(ln);
			} else if (ln.startsWith("#CHROM")) {
				String[] elems = ln.split("\t");
				selectSample = new boolean[elems.length];
				String header = elems[0];
				for (int i = 1; i < elems.length; i++) {
					if (i < nrHeaderElems) {
						header += "\t" + elems[i];
					} else {
						String newName = select.get(elems[i]);
						if (newName != null) {
							selectSample[i] = true;
							header += "\t" + newName;
						}
					}
				}
				tfout.writeln(header);
			} else {
				String[] elems = ln.split("\t");
				String out = elems[0];
				for (int i = 1; i < elems.length; i++) {
					if (i < nrHeaderElems) {
						out += "\t" + elems[i];
					} else {
						if (selectSample[i]) {
							out += "\t" + elems[i];
						}
					}
				}
				tfout.writeln(out);
			}

			ln = in.readLine();
		}
		in.close();

		tfout.close();
	}

	public void performPCAOnVCFSampleCorrelationMatrix(String sequencingVCFGenotypes) throws IOException, Exception {


		// load samples
		TextFile tf = new TextFile(sequencingVCFGenotypes, TextFile.R);

		String[] header = tf.readLineElems(TextFile.tab);

		ArrayList<String> samples = new ArrayList<String>();
		for (int i = 2; i < header.length; i++) {
			samples.add(header[i]);
		}
		System.out.println(samples.size());

		int nrCols = header.length - 2;
		int nrRows = tf.countLines();

		tf.close();
		tf.open();
		tf.readLine();


		double[][] data = new double[nrRows][nrCols];
		String[] elems = tf.readLineElems(TextFile.tab);
		ArrayList<String> variants = new ArrayList<String>();
		int vctr = 0;
		while (elems != null) {

			variants.add(elems[0]);
			for (int i = 2; i < elems.length; i++) {
				data[vctr][i - 2] = Double.parseDouble(elems[i]);
			}
			vctr++;
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>();
		ds.setMatrix(data);
		ds.setColObjects(samples);
		ds.setRowObjects(variants);

		ds.save(sequencingVCFGenotypes + "-matrix.txt");


	}

	public void removeVariantsFromList(String includeTheseVariantsFilter, String exclusionList, String newVariantList) throws IOException {
		TextFile tf = new TextFile(includeTheseVariantsFilter, TextFile.R);
		ArrayList<String> list = tf.readAsArrayList();
		tf.close();

		TextFile tf2 = new TextFile(exclusionList, TextFile.R);
		Set<String> set = tf2.readAsSet(0, TextFile.tab);
		tf2.close();

		TextFile tf3 = new TextFile(newVariantList, TextFile.W);
		for (String s : list) {
			if (!set.contains(s)) {
				tf3.writeln(s);
			}
		}
		tf3.close();

	}


	public void mergeICPedFiles(String[] immunoChipPEDFiles, String[] immunoChipMAPFiles, String variantSelection, String mergedImmunoChipOutput) throws IOException {

		VCFFunctions vcfFunctions = new VCFFunctions();
		PedAndMapFunctions pedAndMapFunctions = new PedAndMapFunctions();
		// first, we need to define a new order for the variants...
		// read the map files and take the intersect
		int snpcounter = 0;
		HashMap<String, Integer> countPerSNP = new HashMap<String, Integer>();

		int ftctr = 0;
		for (String mapfile : immunoChipMAPFiles) {
			ArrayList<String> featureList = pedAndMapFunctions.readMAPFile(mapfile);
			for (String f : featureList) {
				Integer count = countPerSNP.get(f);
				if (count == null) {
					count = 1;
					if (count == 1 && ftctr > 0) {
						System.out.println(f + " not found in first mapfile");
					}
				} else {
					count++;
				}
				countPerSNP.put(f, count);
			}
			ftctr++;
		}

		Set<String> keys = countPerSNP.keySet();
		ArrayList<String> overlappingFeatures = new ArrayList<String>();
		for (String f : keys) {
			if (countPerSNP.get(f) == immunoChipMAPFiles.length) {
				overlappingFeatures.add(f);
			} else {
//				System.out.println(f+" not found in all datasets");
			}
		}

		System.out.println(overlappingFeatures.size() + " overlapping SNPs");

		// index the new variants
		HashMap<String, Integer> snpIndex = new HashMap<String, Integer>();
		int snpCounter = 0;


		HashSet<String> variantsToSelect = new HashSet<String>();
		TextFile tfVariant = new TextFile(variantSelection, TextFile.R);
		variantsToSelect.addAll(tfVariant.readAsArrayList());
		tfVariant.close();

		ArrayList<String> snpsInTheirFinalOrder = new ArrayList<String>();
		for (String snp : overlappingFeatures) {
			String[] snpelems = snp.split("@");
			String actualVariantName = snpelems[snpelems.length - 1];
			if (variantsToSelect.contains(actualVariantName)) {
				snpIndex.put(snp, snpCounter);
				snpsInTheirFinalOrder.add(snp);
				snpCounter++;
			}
		}

		System.out.println(snpIndex.size() + " final SNPs");

		// read the sample selection


		// parse the ped files
		Pattern space = Pattern.compile(" ");
		// 14364 TB03070470 0 0 1 1
		int elemsInPEDHeader = 6;
		int nrOfElemsRequiredForSNPs = snpsInTheirFinalOrder.size() * 2;
		int sizeOfOutput = elemsInPEDHeader + nrOfElemsRequiredForSNPs;
		TextFile out = new TextFile(mergedImmunoChipOutput + ".ped", TextFile.W);
		for (int f = 0; f < immunoChipPEDFiles.length; f++) {
			ArrayList<String> map = pedAndMapFunctions.readMAPFile(immunoChipMAPFiles[f]);

			// determine which SNPs to load and how to reorder them
			int[] newPosition = new int[map.size() * 2];
			for (int i = 0; i < map.size(); i++) {
				String snp = map.get(i);
				Integer newIndex = snpIndex.get(snp);
				if (newIndex == null) {
					newPosition[i * 2] = -1;
					newPosition[(i * 2) + 1] = -1;
				} else {
					newPosition[i * 2] = newIndex * 2;
					newPosition[(i * 2) + 1] = (newIndex * 2) + 1;
				}
			}


			// parse the ped file
			TextFile tf = new TextFile(immunoChipPEDFiles[f], TextFile.R);
			System.out.println("parsing " + immunoChipPEDFiles[f]);
			String[] elems = tf.readLineElems(Strings.whitespace);
			System.out.println(elems.length + "\t" + sizeOfOutput + "\t" + snpIndex.size());

			while (elems != null) {

				String sample = elems[1];


				String[] output = new String[sizeOfOutput];
				for (int i = 0; i < elemsInPEDHeader; i++) {
					output[i] = elems[i];
				}

				int ctr = 0;
				for (int i = elemsInPEDHeader; i < elems.length; i += 2) {

					int newIndex = newPosition[i - elemsInPEDHeader];
					if (newIndex >= 0) {

//						String snp = map.get((i - elemsInPEDHeader) / 2);
//						System.out.println(i + "\t" + newIndex + "\t" + snp);

						newIndex += elemsInPEDHeader;
						output[newIndex] = elems[i];
						output[(newIndex + 1)] = elems[i + 1];
					}
				}

				out.writeln(Strings.concat(output, space));

				elems = tf.readLineElems(space);
			}

			tf.close();
		}

		out.close();

		// write map file
		TextFile mapout = new TextFile(mergedImmunoChipOutput + ".map", TextFile.W);
		for (int i = 0; i < snpsInTheirFinalOrder.size(); i++) {
			String str = snpsInTheirFinalOrder.get(i);
			String[] strElems1 = str.split("@");
			String[] strElems2 = strElems1[0].split("_");
			mapout.writeln(strElems2[0] + "\t" + strElems1[strElems1.length - 1] + "\t0\t" + strElems2[1]);
		}
		mapout.close();

	}


	public HashMap<String, Integer> index(String[] s) {
		int c = 0;
		HashMap<String, Integer> index = new HashMap<String, Integer>();
		for (String s1 : s) {
			index.put(s1, c);
			c++;
		}
		return index;
	}

	public void determineVCFSummaryStatistics(String vcf, String out) throws IOException {
		VCFFunctions vcfFunctions = new VCFFunctions();


		ArrayList<String> vcfSamples = vcfFunctions.getVCFSamples(vcf);
		System.out.println(vcfSamples.size() + " samples loaded from VCF");

		HashMap<String, Integer> sampleMap = index(vcfSamples.toArray(new String[0]));
// get variants from VCF
		ArrayList<Feature> variantsOnVCF = vcfFunctions.getVariantsFromVCF(vcf);
		System.out.println(variantsOnVCF.size() + " variants on VCF");

		ArrayList<String> autosomal = new ArrayList<String>();
		HashMap<Feature, Integer> variantMap = new HashMap<Feature, Integer>();
		int ctr = 0;
		for (Feature f : variantsOnVCF) {
			if (!f.getChromosome().equals(Chromosome.X)) {
				variantMap.put(f, ctr);
				autosomal.add(f.getChromosome().getName() + ":" + f.getStart());
				ctr++;
			}
		}
		Pair<byte[][], String[]> genotypePair = vcfFunctions.loadVCFGenotypes(vcf, sampleMap, variantMap);
		// byte[sample][variant]

		writeGenotypes(genotypePair, vcfSamples, autosomal, out + "allVCFGenotypes.txt");
		byte[][] genotypes = genotypePair.getLeft();
		String[] alleles = genotypePair.getRight();
		Triple<Pair<double[], String[]>, double[], double[]> params = getGenotypeParams(genotypePair);
		double[] mafs = params.getLeft().getLeft();


		TextFile outf = new TextFile(out, TextFile.W);
		double[] xvals = new double[genotypes.length];
		double[] sampleCR = new double[genotypes.length];
		for (int sample = 0; sample < genotypes.length; sample++) {
			int nrTotal = 0;
			int nrCalled = 0;
			for (int gt = 0; gt < genotypes[sample].length; gt++) {
				//if (mafs[gt] > 0.01) {
				if (genotypes[sample][gt] == -1) {

				} else {
					nrCalled++;
				}
				nrTotal++;
				//}
			}
			double cr = (double) nrCalled / nrTotal;
			sampleCR[sample] = cr;
			outf.writeln(vcfSamples.get(sample) + "\t" + cr + "\t" + nrCalled + "\t" + nrTotal);
		}
		outf.close();

		Arrays.sort(sampleCR);

		Coverage c = new Coverage();
		Chart barchart = c.initBarChart("Call-rate per sample", "Sample", "Call-rate", 1200, 800, true, 50);
		barchart.addSeries("Series1", xvals, sampleCR);
		barchart.getStyleManager().setLegendVisible(false);
		BitmapEncoder.saveBitmapWithDPI(barchart, out + "SampleCallRate", BitmapEncoder.BitmapFormat.PNG, 300);


	}


	/*
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
	 */

	public Triple<Pair<double[], String[]>, double[], double[]> getGenotypeParams(Pair<byte[][], String[]> genotypePair) {
		byte[][] genotypes = genotypePair.getLeft();
		String[] alleles = genotypePair.getRight();

		String[] minorAlleles = new String[alleles.length];
		double[] mafs = new double[alleles.length];
		double[] crs = new double[alleles.length];
		double[] hwes = new double[alleles.length];
		for (int gt = 0; gt < genotypes[0].length; gt++) {
			String[] gtAlleles = alleles[gt].split("/");
			String allele1 = gtAlleles[0];
			String allele2 = gtAlleles[1];
			int nrCalled = 0;
			int nrTotal = 0;
			int nrHomAA = 0;
			int nrHomBB = 0;
			int nrHets = 0;
			for (int sample = 0; sample < genotypes.length; sample++) {
				if (genotypes[sample][gt] == -1) {
					// missing
				} else if (genotypes[sample][gt] == 0) {
					nrHomAA++;
					nrCalled += 2;
				} else if (genotypes[sample][gt] == 1) {
					nrHets++;
					nrCalled += 2;
				} else {
					nrHomBB++;
					nrCalled += 2;
				}
				nrTotal += 2;
			}
			int freqAlleleA = 2 * nrHomAA + nrHets; // genotypeFreq[0] + genotypeFreq[1];
			int freqAlleleB = 2 * nrHomBB + nrHets; // genotypeFreq[2] + genotypeFreq[1];
			double maf = (double) (freqAlleleA) / (nrCalled);
			String minorAllele = allele1;
			if (freqAlleleA > freqAlleleB) {
				maf = 1 - maf;
				minorAllele = allele2;
			}
			minorAlleles[gt] = minorAllele;
			double cr = (double) nrCalled / nrTotal;
			double hwep = HWE.calculateExactHWEPValue(nrHets, nrHomAA, nrHomBB);
			mafs[gt] = maf;
			crs[gt] = cr;
			hwes[gt] = hwep;
			// System.out.println(gt + "\t" + cr + "\t" + maf);
		}
		return new Triple<Pair<double[], String[]>, double[], double[]>(new Pair<double[], String[]>(mafs, minorAlleles), crs, hwes);
	}


	public void compareVCFGenotypesToPedAndMap(String vcf, String plinkfile, String out, boolean onlySharedSamples) throws IOException {

		// get samples from plink data
		VCFFunctions vcfFunctions = new VCFFunctions();
		PedAndMapFunctions pedAndMapFunctions = new PedAndMapFunctions();

		ArrayList<Pair<String, String>> famSamples = pedAndMapFunctions.parseFam(plinkfile + ".fam");
		System.out.println(famSamples.size() + " samples loaded from PED/MAP");
		// get samples from VCF
		ArrayList<String> vcfSamples = vcfFunctions.getVCFSamples(vcf);
		System.out.println(vcfSamples.size() + " samples loaded from VCF");

		// intersect samples
		ArrayList<String> samplesInFam = new ArrayList<String>();
		for (Pair<String, String> p : famSamples) {
			samplesInFam.add(p.getRight());
		}
		HashSet<String> intersectedSamples = intersectSamples(samplesInFam, vcfSamples);

		System.out.println(intersectedSamples.size() + " samples shared");
		writeHash(intersectedSamples, out + "sharedSamples.txt");

		// get variants from VCF

		ArrayList<Feature> variantsOnVCF = vcfFunctions.getVariantsFromVCF(vcf);
		System.out.println(variantsOnVCF.size() + " variants on VCF");

		// get variants from PLINK file that have same position
		ArrayList<Feature> variantsOnMap = pedAndMapFunctions.getVariantsFromMapFile(plinkfile + ".map");
		System.out.println(variantsOnMap.size() + " variants on MAP");

		// intersect variants
		HashSet<Feature> intersectedVariants = intersectVariants(variantsOnMap, variantsOnVCF);
		System.out.println(intersectedVariants.size() + " shared variants");
		writeFeatureHash(intersectedVariants, out + "sharedVariants.txt");

		ArrayList<String> sharedSamples = new ArrayList<String>();
		sharedSamples.addAll(intersectedSamples);
		Collections.sort(sharedSamples);

		HashMap<String, Integer> sampleMapPED = new HashMap<String, Integer>();
		HashMap<String, Integer> sampleMapVCF = new HashMap<String, Integer>();
		int s = 0;

		ArrayList<String> newSampleOrderPed = new ArrayList<String>();
		ArrayList<String> newSampleOrderVCF = new ArrayList<String>();
		for (String sample : sharedSamples) {
			sampleMapPED.put(sample, s);
			sampleMapVCF.put(sample, s);
			newSampleOrderPed.add(sample);
			newSampleOrderVCF.add(sample);
			s++;
		}

		if (!onlySharedSamples) {
			int s2 = sampleMapPED.size();

			for (String sample : samplesInFam) {
				if (!sampleMapPED.containsKey(sample)) {
					sampleMapPED.put(sample, s2);
					newSampleOrderPed.add(sample);
					s2++;
				}
			}

			s2 = sampleMapVCF.size();
			for (String sample : vcfSamples) {
				if (!sampleMapVCF.containsKey(sample)) {
					sampleMapVCF.put(sample, s2);
					newSampleOrderVCF.add(sample);
					s2++;
				}
			}
		}

		HashMap<Feature, Integer> variantMap = new HashMap<Feature, Integer>();

		ArrayList<String> variants = new ArrayList<String>();
		HashMap<String, Feature> strToFeat = new HashMap<String, Feature>();
		for (Feature f : intersectedVariants) {
			variants.add(f.getChromosome().getName() + ":" + f.getStart());
			strToFeat.put(f.getChromosome().getName() + ":" + f.getStart(), f);
		}
		Collections.sort(variants);

		int v = 0;
		for (String var : variants) {
			Feature f = strToFeat.get(var);
			variantMap.put(f, v);
			v++;
		}


		// load genotypes for VCF
		// format: byte[sample][genotype] String[nrVariants];


		Pair<byte[][], String[]> vcfGenotypesPair = vcfFunctions.loadVCFGenotypes(vcf, sampleMapVCF, variantMap);
		writeGenotypes(vcfGenotypesPair, newSampleOrderVCF, variants, out + "vcfGenotypes.txt");


		// load genotypes for PLINK

		Pair<byte[][], String[]> pedGenotypesPair = pedAndMapFunctions.loadPedGenotypes(plinkfile + ".ped", plinkfile + ".map", sampleMapPED, variantMap);
		writeGenotypes(pedGenotypesPair, newSampleOrderPed, variants, out + "pedGenotypes.txt");

		// compare alleles
		String[] allelesVCF = vcfGenotypesPair.getRight();
		String[] allelesPED = pedGenotypesPair.getRight();


//		System.out.println("vcf vcf correlations");
//		correlateSamples(vcfGenotypesPair, vcfGenotypesPair, samples, out + "vcfSampleCorrelations.txt", excludeVariant, flipGenotypes);
//		System.out.println("ped ped correlations");
//		correlateSamples(pedGenotypesPair, pedGenotypesPair, samples, out + "pedSampleCorrelations.txt", excludeVariant, flipGenotypes);


		Triple<Pair<double[], String[]>, double[], double[]> vcfParams = getGenotypeParams(vcfGenotypesPair);
		Triple<Pair<double[], String[]>, double[], double[]> pedParams = getGenotypeParams(pedGenotypesPair);
		double mafThreshold = 0.01;
		double crThreshold = 0.05;
		Pair<boolean[], boolean[]> flipvariants = determineAlleleFlipsAndExcludeVariants(allelesPED, pedParams.getLeft(), allelesVCF, vcfParams.getLeft(), variants, out + "vcfpedAlleleFlips.txt",
				mafThreshold, crThreshold);
		boolean[] excludeVariant = flipvariants.getLeft();
		boolean[] flipGenotypes = flipvariants.getRight();

		// correlate samples
		System.out.println("vcf ped sample correlations");


		Pair<double[][], double[][]> correlationOutputSamples = correlateSamples(vcfGenotypesPair, pedGenotypesPair, newSampleOrderVCF, newSampleOrderPed, out + "vcfpedSampleCorrelations.txt", excludeVariant, flipGenotypes);


		boolean[] sampleSwapped = writeBestCorrelations(correlationOutputSamples, newSampleOrderVCF, newSampleOrderPed, out + "vcfpedBestSampleCorrelations.txt");

		// correlate variants
		System.out.println("vcf ped genotype correlations");


		Pair<double[][], double[][]> correlationOutputVariants = correlateVariants(vcfGenotypesPair, pedGenotypesPair, variants, out + "vcfpedGenotypeCorrelations.txt", excludeVariant, flipGenotypes, sampleSwapped);
		writeBestCorrelations(correlationOutputVariants, variants, variants, out + "vcfpedBestGenotypeCorrelations.txt");

		// write discordantGenotypes
		writeDiscordantGenotypes(vcfGenotypesPair, pedGenotypesPair, variants, out + "discordantGenotypes.txt", excludeVariant, flipGenotypes, sampleSwapped, sharedSamples);


	}

	public void writeDiscordantGenotypes(Pair<byte[][], String[]> vcfGenotypesPair, Pair<byte[][], String[]> pedGenotypesPair, ArrayList<String> variants, String outputFileName, boolean[] excludeVariant, boolean[] flipGenotypes, boolean[] sampleSwapped, ArrayList<String> samples) throws IOException {

		TextFile out = new TextFile(outputFileName, TextFile.W);

		String header = "Variant\tAllelesVCF\tAllelesPED (after flip)\tFlipped\tSample\tGenotypeVCF\tGenotypePED (after flip)\tGenotypeVCFAlleles\tGenotypePEDAlleles (after flip)";
		out.writeln(header);

		byte[][] genotypesVCF = vcfGenotypesPair.getLeft();
		String[] allelesVCF = vcfGenotypesPair.getRight();
		byte[][] genotypesPED = pedGenotypesPair.getLeft();
		String[] allelesPED = pedGenotypesPair.getRight();

		for (int v = 0; v < variants.size(); v++) {
			if (!excludeVariant[v]) {
				String variantOut = variants.get(v);

				String variantAllelesVCF = allelesVCF[v];
				String[] allelesVCFArr = variantAllelesVCF.split("/");
				String variantAllelesPED = allelesPED[v];
				String[] allelesPEDArr = variantAllelesPED.split("/");

				boolean flipped = flipGenotypes[v];
				if (flipped) {
					String all1 = allelesPEDArr[0];
					String all2 = allelesPEDArr[1];
					allelesPEDArr[0] = all2;
					allelesPEDArr[1] = all1;
				}


				variantOut += "\t" + variantAllelesVCF + "\t" + allelesPEDArr[0] + "/" + allelesPEDArr[1] + "\t" + flipped;

				for (int s = 0; s < genotypesPED.length; s++) {
					if (s < sampleSwapped.length && !sampleSwapped[s]) {

						int gt1 = genotypesVCF[s][v];
						int gt2 = genotypesPED[s][v];

						if (flipped) {
							gt2 = Math.abs(gt2 - 2);
							if (gt2 == 3) {
								gt2 = -1;
							}
						}


						if (gt1 != gt2 || variants.get(v).equals("Chr15:79152806")) {
							String sample = samples.get(s);

							String ln = variantOut + "\t" +
									sample + "\t" +
									gt1 + "\t" +
									gt2 + "\t" +
									gtToAlleles(allelesVCFArr, gt1) + "\t" +
									gtToAlleles(allelesPEDArr, gt2);
							out.writeln(ln);
						}

					}
				}
			}
		}

		out.close();

	}

	public String gtToAlleles(String[] alleles, int gt) {
		if (gt == 0) {
			return alleles[0] + "/" + alleles[0];
		} else if (gt == 1) {
			return alleles[0] + "/" + alleles[1];
		} else if (gt == 2) {
			return alleles[1] + "/" + alleles[1];
		} else {
			return "NA";
		}
	}

	public boolean[] writeBestCorrelations(Pair<double[][], double[][]> correlationInput, ArrayList<String> labels1, ArrayList<String> labels2, String out) throws IOException {
		double[][] correlations = correlationInput.getLeft();
		double[][] n = correlationInput.getRight();
		TextFile outf = new TextFile(out + "-rows.txt", TextFile.W);

		outf.writeln("Label\t" +
				"Correlation\t" +
				"n\t" +
				"BestMatch\t" +
				"BestCorrelation\t" +
				"n\t" +
				"Match");

		//  this works only for the diagonal
		boolean[] sampleSwapped = new boolean[labels1.size()];
		for (int i = 0; i < labels1.size(); i++) {
			String ln = "";
			if (i < labels2.size()) {
				ln = labels1.get(i) + "\t"
						+ correlations[i][i]
						+ "\t" + n[i][i];
			} else {
				ln = labels1.get(i) + "\tnull\tnull";
			}
			double max = -1;
			double origMax = -1;
			int maxS = i;
			for (int j = 0; j < labels2.size(); j++) {
				double corr = correlations[i][j];
				if (!Double.isNaN(corr)) {
					double absCorr = Math.abs(corr);
					if (absCorr > max) {
						max = absCorr;
						origMax = corr;
						maxS = j;
					}
				}
			}
			if (max == -1) {
				ln += "\tUnknown" +
						"\tNaN" +
						"\t0" +
						"\tUnknown";
				sampleSwapped[i] = true;
			} else {
				ln += "\t" + labels2.get(maxS) + "\t"
						+ origMax + "\t"
						+ n[i][maxS] + "\t"
						+ (maxS == i);
				sampleSwapped[i] = !(maxS == i);
			}
			outf.writeln(ln);
		}
		outf.close();
		return sampleSwapped;
	}

	public Pair<boolean[], boolean[]> determineAlleleFlipsAndExcludeVariants(String[] allelesPED, Pair<double[], String[]> pedMinorAlleles,
																			 String[] allelesVCF, Pair<double[], String[]> vcfMinorAlleles,
																			 ArrayList<String> variants, String out,
																			 double mafThreshold, double crThreshold) throws IOException {
		boolean[] excludeVariant = new boolean[allelesPED.length];
		boolean[] flipGenotypes = new boolean[allelesPED.length];

		double[] mafsPED = pedMinorAlleles.getLeft();
		String[] minorAllelesPED = pedMinorAlleles.getRight();

		double[] mafsVCF = vcfMinorAlleles.getLeft();
		String[] minorAllelesVCF = vcfMinorAlleles.getRight();

		/*
		for (int variant = 0; variant < excludeVariant.length; variant++) {
			if (mafs[variant] < mafThreshold || cr[variant] < crTheshold)) {
				excludeVariant[variant] = true;
				System.out.println("excluding variant: " + variants.get(variant) + "\t" + mafs[variant] + "\t" + cr[variant]);
			}
		}
		 */
		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln("variant\t" +
				"allelesPED\tminorAllelePED\tmafPED\t" +
				"allelesVCF\tminorAlleleVCF\tmafVCF\t" +
				"deltaMAF\tflip\tReason\tExcludeVariant");


		for (int a = 0; a < allelesPED.length; a++) {
			String allelePED = allelesPED[a];
			String alleleVCF = allelesVCF[a];
			String allele1ped = allelePED.split("/")[0];
			String allele1vcf = alleleVCF.split("/")[0];
			String allele2ped = allelePED.split("/")[1];
			String allele2vcf = alleleVCF.split("/")[1];


			String output = variants.get(a) + "\t" +
					allelesPED[a] + "\t" + minorAllelesPED[a] + "\t" + mafsPED[a] + "\t" +
					allelesVCF[a] + "\t" + minorAllelesVCF[a] + "\t" + mafsVCF[a] + "\t" + Math.abs(mafsVCF[a] - mafsPED[a]);


			if (allele2ped.equals("null") || allele1ped.equals("null") || allele1vcf.equals("null") || allele2vcf.equals("null")) {
				output += "\tfalse\tNullAllele\ttrue";
				excludeVariant[a] = true;
			} else {
				if (allele1ped.equals(BaseAnnot.getComplement(allele2ped)) && allele1vcf.equals(BaseAnnot.getComplement(allele2vcf))) {


					// A/T or G/C SNP
					// try to figure out allele flip on the basis of minor alleles

					String minorAllelePED = minorAllelesPED[a];
					String minorAlleleVCF = minorAllelesVCF[a];
					boolean flipAlleles = false;
					if ((allele1ped.equals(minorAllelePED) && allele1vcf.equals(minorAlleleVCF)) || (allele2ped.equals(minorAllelePED) && allele2vcf.equals(minorAlleleVCF))) {
						// same direction

					} else {
						flipAlleles = true;
					}

					output += "\t" + flipAlleles + "\tComplement\t" + false;
					flipGenotypes[a] = flipAlleles;

					excludeVariant[a] = false;
				} else {
					Boolean flipAlleles = BaseAnnot.flipalleles(allelePED, allele1ped, alleleVCF, allele1vcf);
					if (flipAlleles == null) {
						excludeVariant[a] = true;
						output += "\tfalse\tIncompatibleAlleles\ttrue";
					} else {
						output += "\t" + flipAlleles + "\tOK\tfalse";
						flipGenotypes[a] = flipAlleles;
					}
				}
			}
			outf.writeln(output);
		}
		outf.close();
		return new Pair<boolean[], boolean[]>(excludeVariant, flipGenotypes);
	}


	public Pair<double[][], double[][]> correlateSamples(Pair<byte[][], String[]> genotypesPair1,
														 Pair<byte[][], String[]> genotypesPair2,
														 ArrayList<String> samples1,
														 ArrayList<String> samples2,
														 String s,
														 boolean[] excludeCols,
														 boolean[] flipAlleles) throws IOException {


		byte[][] genotypes1 = genotypesPair1.getLeft();
		byte[][] genotypes2 = genotypesPair2.getLeft();
		double[][] n = new double[samples1.size()][samples2.size()];
		double[][] correlations = new double[samples1.size()][samples2.size()];

		for (int i = 0; i < samples1.size(); i++) {

			for (int j = 0; j < samples2.size(); j++) {
				ArrayList<Double> d1 = new ArrayList<Double>();
				ArrayList<Double> d2 = new ArrayList<Double>();
				for (int g = 0; g < genotypes1[i].length; g++) {
					boolean flip = flipAlleles[g];
					if (genotypes1[i][g] != -1 && genotypes2[j][g] != -1 && !excludeCols[g]) {
						d1.add((double) genotypes1[i][g]);
						double gt = (double) genotypes2[j][g];
						if (flip) {
							gt = Math.abs(gt - 2);
						}
						d2.add(gt);

					}
				}

				double[] g1 = Primitives.toPrimitiveArr(d1.toArray(new Double[0]));
				double[] g2 = Primitives.toPrimitiveArr(d2.toArray(new Double[0]));
				correlations[i][j] = Correlation.correlate(g1, g2);

				n[i][j] = d1.size();

			}
		}

		writeCorrelationMatrix(correlations, samples1, samples2, s);
		writeCorrelationMatrix(n, samples1, samples2, s + "-n.txt");

		return new Pair<double[][], double[][]>(correlations, n);
	}

	public Pair<double[][], double[][]> correlateVariants(Pair<byte[][], String[]> genotypesPair1,
														  Pair<byte[][], String[]> genotypesPair2, ArrayList<String> variants, String s,
														  boolean[] excludeVariants, boolean[] flipAlleles, boolean[] sampleSwapped) throws IOException {

		double[][] correlations = new double[variants.size()][variants.size()];
		double[][] n = new double[variants.size()][variants.size()];
		byte[][] genotypes1 = genotypesPair1.getLeft();
		byte[][] genotypes2 = genotypesPair2.getLeft();

		for (int variantId1 = 0; variantId1 < variants.size(); variantId1++) {

			if (!excludeVariants[variantId1]) {
				for (int variantId2 = 0; variantId2 < variants.size(); variantId2++) {
					if (!excludeVariants[variantId2]) {


						ArrayList<Double> d1 = new ArrayList<Double>();
						ArrayList<Double> d2 = new ArrayList<Double>();

						boolean flip2 = flipAlleles[variantId2];

						// this can only use the diagonal, where samples are supposedly equal
						for (int sample = 0; sample < genotypes1.length && sample < genotypes2.length; sample++) {

							double gt1 = (double) genotypes1[sample][variantId1];
							double gt2 = (double) genotypes2[sample][variantId2];

							if (!sampleSwapped[sample] && gt1 != -1 && gt2 != -1) {
								d1.add(gt1);

								if (flip2) {
									gt2 = Math.abs(gt2 - 2);
								}
								d2.add(gt2);
							}
						}

						double[] g1 = Primitives.toPrimitiveArr(d1.toArray(new Double[0]));
						double[] g2 = Primitives.toPrimitiveArr(d2.toArray(new Double[0]));

						double c = JSci.maths.ArrayMath.correlation(g1, g2);
						if (d1.isEmpty()) {
							if (!excludeVariants[variantId1] && !excludeVariants[variantId2]) {
								System.err.println("ERROR: no genotypes for variants: " + variants.get(variantId1) + " or " + variants.get(variantId2));
								c = Double.NaN;
							} else {
								c = Double.NaN;
							}
						}
						correlations[variantId1][variantId2] = c;
						n[variantId1][variantId2] = d1.size();
						n[variantId2][variantId1] = d1.size();
						correlations[variantId2][variantId1] = c;

						if (variants.get(variantId1).equals("Chr6:167433729") && variants.get(variantId2).equals("Chr6:167433729")) {
							for (int q = 0; q < d1.size(); q++) {
								System.out.println(q + "\t" + d1.get(q) + "\t" + d2.get(q));
							}

							System.out.println(c);
						}


					} else {
						correlations[variantId1][variantId2] = Double.NaN;
						correlations[variantId2][variantId1] = correlations[variantId1][variantId2];
					}
				}
			} else {
				for (int variantId2 = 0; variantId2 < variants.size(); variantId2++) {
					correlations[variantId1][variantId2] = Double.NaN;
					correlations[variantId2][variantId1] = correlations[variantId1][variantId2];
				}
			}

		}


		writeCorrelationMatrix(correlations, variants, variants, s);
		writeCorrelationMatrix(n, variants, variants, s + "-n.txt");

		return new Pair<double[][], double[][]>(correlations, n);
	}

	public void writeCorrelationMatrix(double[][] correlations, ArrayList<String> labels1, ArrayList<String> labels2, String s) throws IOException {
		TextFile out = new TextFile(s, TextFile.W);
		String header = "-";
		for (int i = 0; i < labels2.size(); i++) {
			header += "\t" + labels2.get(i);
		}
		out.writeln(header);
		for (int i = 0; i < labels1.size(); i++) {
			String ln = labels1.get(i);
			for (int j = 0; j < correlations[i].length; j++) {
				ln += "\t" + correlations[i][j];
			}
			out.writeln(ln);
		}
		out.close();
	}


	public void writeGenotypes(Pair<byte[][], String[]> vcfGenotypesPair, ArrayList<String> samples, ArrayList<String> variants, String s) throws IOException {
		TextFile out = new TextFile(s, TextFile.W);

		byte[][] genotypes = vcfGenotypesPair.getLeft();
		String[] alleles = vcfGenotypesPair.getRight();
		String header = "variant\talleles";
		for (int i = 0; i < genotypes.length; i++) {
			header += "\t" + samples.get(i);
		}
		out.writeln(header);
		for (int i = 0; i < genotypes[0].length; i++) {
			String ln = variants.get(i) + "\t" + alleles[i];
			for (int j = 0; j < genotypes.length; j++) {
				ln += "\t" + genotypes[j][i];
			}
			out.writeln(ln);
		}
		out.close();

	}

	public void writeHash(HashSet<String> set, String s) throws IOException {
		TextFile out = new TextFile(s, TextFile.W);
		for (String v : set) {
			out.writeln(v);
		}
		out.close();
	}

	public void writeFeatureHash(HashSet<Feature> set, String s) throws IOException {
		TextFile out = new TextFile(s, TextFile.W);
		for (Feature v : set) {
			out.writeln(v.getChromosome().getName() + "\t" + v.getStart());
		}
		out.close();
	}


	public HashSet<String> intersectSamples(ArrayList<String> s1, ArrayList<String> s2) {
		HashSet<String> intersect = new HashSet<String>();
		HashSet<String> set1 = new HashSet<String>();
		set1.addAll(s1);
		for (String s : s2) {
			if (set1.contains(s)) {
				intersect.add(s);
			}
		}
		return intersect;
	}

	public HashSet<Feature> intersectVariants(ArrayList<Feature> f1, ArrayList<Feature> f2) {
		HashSet<Feature> shared = new HashSet<Feature>();
		for (Feature feat1 : f1) {
			for (Feature feat2 : f2) {
				if (feat1.overlaps(feat2)) {
					shared.add(feat1);
				}
			}
		}
		return shared;

	}


}
