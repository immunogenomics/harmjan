package nl.harmjanwestra.vcfutils;

import com.gs.collections.api.list.MutableList;
import com.gs.collections.impl.multimap.list.FastListMultimap;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.genotypes.GenotypeTools;
import nl.harmjanwestra.utilities.vcf.VCFFunctions;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.*;

/**
 * Created by hwestra on 9/25/15.
 */
public class VCFMerger {
	public static void main(String[] args) {
		try {

			if (args.length < 6) {
				System.out.println("Usage: linuxbool chr vcfsort refvcf testvcf matchedout sepstr");
			} else {


				String linuxStr = args[0];
				boolean linux = Boolean.parseBoolean(linuxStr);
				String chrStr = args[1];
				int chr = Integer.parseInt(chrStr);
				String vcfsort = args[2];
				String refvcf = args[3];
				String testvcf = args[4];
				String matchedout = args[5];
				String separatorStr = args[6];

				String separator = "\t";
				if (separatorStr.equals("tab")) {
					separator = "\t";
				} else if (separatorStr.equals("space")) {
					separator = " ";
				}

				VCFMerger m = new VCFMerger();

				m.mergeAndIntersect(linux, chr, vcfsort, refvcf, testvcf, matchedout, separator);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	public void mergeAndIntersect(boolean linux, int chrint, String vcfsort, String refVCF, String testVCF, String matchedPanelsOut, String separator) throws IOException {
		Chromosome chr = Chromosome.parseChr("" + chrint);

		mergeAndIntersectVCFVariants(
				refVCF,
				testVCF,
				matchedPanelsOut + "ref-matched-" + chr.getName() + ".vcf.gz",
				matchedPanelsOut + "test-matched-" + chr.getName() + ".vcf.gz",
				matchedPanelsOut + "ref-test-merged-" + chr.getName() + ".vcf.gz",
				separator,
				matchedPanelsOut + "mergelog-" + chr.getName() + ".txt",
				true);

		VCFFunctions t = new VCFFunctions();
		t.sortVCF(linux, vcfsort, matchedPanelsOut + "ref-matched-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "ref-matched-sorted-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "ref-sort-" + chr.getName() + ".sh");
		t.sortVCF(linux, vcfsort, matchedPanelsOut + "test-matched-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "test-matched-sorted-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "test-sort-" + chr.getName() + ".sh");
		t.sortVCF(linux, vcfsort, matchedPanelsOut + "ref-test-merged-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "ref-test-merged-sorted-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "test-sort-" + chr.getName() + ".sh");
	}

	private void mergeAndIntersectVCFVariants(String refVCF,
											  String testVCF,
											  String vcf1out,
											  String vcf2out,
											  String vcfmergedout,
											  String separatorInMergedFile,
											  String logoutfile,
											  boolean keepNonOverlapping) throws IOException {

		System.out.println("Merging: ");
		System.out.println("ref: " + refVCF);
		System.out.println("test: " + testVCF);
		System.out.println("out: " + vcfmergedout);


		TextFile tf = new TextFile(refVCF, TextFile.R);
		TextFile mergedOut = new TextFile(vcfmergedout, TextFile.W);
		TextFile vcf1OutTf = new TextFile(vcf1out, TextFile.W);

		String header = "";
		String ln = tf.readLine();
		HashSet<String> headerLines = new HashSet<String>();
		while (ln != null) {
			if (ln.startsWith("##")) {


				if (!headerLines.contains(ln)) {
					if (header.length() == 0) {
						header += ln;
					} else {
						header += "\n" + ln;
					}
					vcf1OutTf.writeln(ln);
					headerLines.add(ln);
				}
			} else if (ln.startsWith("#")) {
				vcf1OutTf.writeln(ln);
			} else {
				break;
			}

			ln = tf.readLine();
		}
		tf.close();

		TextFile tf2 = new TextFile(testVCF, TextFile.R);
		TextFile vcf2OutTf = new TextFile(vcf2out, TextFile.W);
		ln = tf2.readLine();
		while (ln != null) {
			if (ln.startsWith("##")) {

				if (!headerLines.contains(ln)) {
					if (header.length() == 0) {
						header += ln;
					} else {
						header += "\n" + ln;
					}
					vcf2OutTf.writeln(ln);
					headerLines.add(ln);
				}

			} else if (ln.startsWith("#")) {
				vcf2OutTf.writeln(ln);
			} else {
				break;
			}

			ln = tf2.readLine();
		}
		tf2.close();

		VCFFunctions t = new VCFFunctions();

		ArrayList<String> samples1 = t.getVCFSamples(refVCF);
		ArrayList<String> samples2 = t.getVCFSamples(testVCF);

		// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1   SAMPLE2
		String sampleheaderLn = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
		for (String s : samples1) {
			sampleheaderLn += "\t" + s;
		}
		for (String s : samples2) {
			sampleheaderLn += "\t" + s;
		}

		mergedOut.writeln(header + "\n" + sampleheaderLn);


		// reload the variants, get the positions with multiple variants at that position


		VCFGenotypeData iterator = new VCFGenotypeData(refVCF);
		VCFGenotypeData iterator2 = new VCFGenotypeData(testVCF);

		HashSet<Chromosome> chromosomes = new HashSet<Chromosome>();

		HashSet<Feature> allFeatures = new HashSet<Feature>();
		System.out.println("Inventorizing variants in: " + refVCF);
		while (iterator.hasNext()) {
			VCFVariant f = iterator.next();
			Feature feat = new Feature();
			feat.setChromosome(Chromosome.parseChr(f.getChr()));
			feat.setStart(f.getPos());
			feat.setStop(f.getPos());
			chromosomes.add(feat.getChromosome());
			allFeatures.add(feat);
		}

		System.out.println(allFeatures.size() + " features after loading ref vcf");
		iterator.close();

		System.out.println("Inventorizing variants in: " + testVCF);

		while (iterator2.hasNext()) {
			VCFVariant f = iterator2.next();
			Feature feat = new Feature();
			feat.setChromosome(Chromosome.parseChr(f.getChr()));
			feat.setStart(f.getPos());
			feat.setStop(f.getPos());
			chromosomes.add(feat.getChromosome());
			allFeatures.add(feat);
		}
		iterator2.close();
		System.out.println(allFeatures.size() + " features after loading test vcf");

		ArrayList<Feature> allFeatureArr = new ArrayList<Feature>();
		allFeatureArr.addAll(allFeatures);
		Collections.sort(allFeatureArr, new FeatureComparator(false));

		TextFile logOut = new TextFile(logoutfile, TextFile.W);
		logOut.writeln("chr\t" +
				"pos\t" +
				"refVariant#atPos\t" +
				"refVariantId\t" +
				"refVariantRefAllele\t" +
				"refVariantAltAlleles\t" +
				"refVariantMinorAllele\t" +
				"refVariantMAF\t" +
				"restVariant#atPos\t" +
				"testVariantId\t" +
				"testVariantRefAllele\t" +
				"testVariantAltAlleles\t" +
				"testVariantMinorAllele\t" +
				"testVariantMAF\t" +
				"testVariantRefAlleleAM\t" +
				"testVariantAltAllelesAM\t" +
				"testVariantMinorAlleleAM\t" +
				"Reason");

		TextFile uniqueRef = new TextFile(vcf1out + "-uniqueVars.txt", TextFile.W);
		TextFile uniqueTest = new TextFile(vcf2out + "-uniqueVars.txt", TextFile.W);
		for (Chromosome chr : chromosomes) {
			ArrayList<Feature> chrFeatures = new ArrayList<Feature>();
			HashSet<Feature> chrFeatureSet = new HashSet<Feature>();
			for (Feature f : allFeatureArr) {
				if (f.getChromosome().equals(chr)) {
					chrFeatures.add(f);
					chrFeatureSet.add(f);
				}
			}


			// load the genotypes..
			System.out.println("loading variants from: " + refVCF);
			FastListMultimap<Feature, VCFVariant> refVariantMap = t.getVCFVariants(refVCF, chrFeatureSet, true);
			System.out.println("loading variants from: " + testVCF);
			FastListMultimap<Feature, VCFVariant> testVariantMap = t.getVCFVariants(testVCF, chrFeatureSet, true);

			System.out.println(chr.getName() + "\t" + chrFeatures.size() + " features " + refVariantMap.size() + " in ref and " + testVariantMap.size() + " in test");

			for (Feature f : chrFeatures) {


				MutableList<VCFVariant> refVariants = refVariantMap.get(f);

				MutableList<VCFVariant> testVariants = testVariantMap.get(f);


				if (!(refVariantMap.get(f).size() > 0 && testVariantMap.get(f).size() > 0)) {
					// variant unique to one

					if (keepNonOverlapping) {
						if (refVariantMap.containsKey(f)) {
							for (int x = 0; x < refVariants.size(); x++) {
								VCFVariant var1 = refVariants.get(x);
								vcf1OutTf.writeln(var1.toVCFString());
								uniqueRef.writeln(var1.getChr() + "\t" + var1.getPos() + "\t" + var1.getId());
							}

						} else if (testVariantMap.containsKey(f)) {
							for (int x = 0; x < testVariants.size(); x++) {
								VCFVariant var2 = testVariants.get(x);
								vcf2OutTf.writeln(var2.toVCFString());
								uniqueTest.writeln(var2.getChr() + "\t" + var2.getPos() + "\t" + var2.getId());
							}
						}
					}
				} else {
					boolean[] testVariantAlreadyWritten = new boolean[testVariants.size()];


					for (int x = 0; x < refVariants.size(); x++) {
						VCFVariant refVariant = refVariants.get(x);

						boolean refVariantAlreadyMerged = false;
						if (refVariant.getPos() == 2985620) {
							System.out.println(x + " - ref: " + Strings.concat(refVariant.getAlleles(), Strings.comma, 0, refVariant.getAlleles().length));
						}

						String logln = refVariant.getChr()
								+ "\t" + refVariant.getPos()
								+ "\t" + x
								+ "\t" + refVariant.getId()
								+ "\t" + refVariant.getAlleles()[0]
								+ "\t" + Strings.concat(refVariant.getAlleles(), Strings.comma, 1, refVariant.getAlleles().length)
								+ "\t" + refVariant.getMinorAllele()
								+ "\t" + refVariant.getMAF();

						for (int y = 0; y < testVariants.size(); y++) {

							String logoutputln = logln;


							VCFVariant testVariant = testVariants.get(y);
							if (testVariant.getPos() == 2985620) {
								System.out.println(y + " " + testVariantAlreadyWritten[y] + " test: " + Strings.concat(testVariant.getAlleles(), Strings.comma, 0, testVariant.getAlleles().length));
							}
							if (testVariantAlreadyWritten[y] || refVariantAlreadyMerged) {
								logoutputln += "\t" + y
										+ "\t" + testVariant.getId()
										+ "\t" + testVariant.getAlleles()[0]
										+ "\t" + Strings.concat(testVariant.getAlleles(), Strings.comma, 1, testVariant.getAlleles().length)
										+ "\t" + testVariant.getMinorAllele()
										+ "\t" + testVariant.getMAF() + "\t-\t-\tRefVariantAlreadyMerged";
								// don't write again //
							} else {


								logoutputln += "\t" + y
										+ "\t" + testVariant.getId()
										+ "\t" + testVariant.getAlleles()[0]
										+ "\t" + Strings.concat(testVariant.getAlleles(), Strings.comma, 1, testVariant.getAlleles().length)
										+ "\t" + testVariant.getMinorAllele()
										+ "\t" + testVariant.getMAF();


								if (!refVariant.getId().equals(testVariant.getId())) {
									// check whether the name is equal (this may matter for positions with multiple variants).
									logoutputln += "\t-\t-\tDifferentNames";
								} else {

									if (refVariant.getId().equals("rs2240498")) {
										System.out.println("Got it!");
									}

									Pair<String, String> outputpair = mergeVariants(refVariant, testVariant, separatorInMergedFile);
									if (outputpair.getRight() != null) {
										vcf1OutTf.writeln(refVariant.toVCFString());
										vcf2OutTf.writeln(testVariant.toVCFString());
										mergedOut.writeln(outputpair.getRight());
										testVariantAlreadyWritten[y] = true;
										refVariantAlreadyMerged = true;
									}
									logoutputln += outputpair.getLeft();

								}
							}
							logOut.writeln(logoutputln);


						}

						if (refVariant.getPos() == 2985620) {
							for (int y = 0; y < testVariants.size(); y++) {
								System.out.println(testVariantAlreadyWritten[y] + " -- " + refVariantAlreadyMerged);
							}
						}

					}
				}

			}

		}
		uniqueRef.close();
		uniqueTest.close();

		logOut.close();
		vcf1OutTf.close();
		vcf2OutTf.close();
		mergedOut.close();


	}

	private Pair<String, String> mergeVariants(VCFVariant refVariant, VCFVariant testVariant,
											   String separatorInMergedFile
	) {
		int nridenticalalleles = 0;

		String[] refAlleles = refVariant.getAlleles();
		String refMinorAllele = refVariant.getMinorAllele();
		String[] testVariantAlleles = testVariant.getAlleles();
		String testVariantMinorAllele = testVariant.getMinorAllele();

		for (int i = 0; i < refAlleles.length; i++) {
			String allele1 = refAlleles[i];
			for (int j = 0; j < testVariantAlleles.length; j++) {
				if (testVariantAlleles[j].equals(allele1)) {
					nridenticalalleles++;
				}
			}
		}

		GenotypeTools gtools = new GenotypeTools();

		boolean complement = false;
		if (nridenticalalleles == 0) {
			// try complement
			complement = true;
			String[] complementAlleles2 = gtools.convertToComplement(testVariantAlleles);
			testVariantMinorAllele = gtools.getComplement(testVariantMinorAllele);
			nridenticalalleles = 0;

			for (int i = 0; i < refAlleles.length; i++) {
				String allele1 = refAlleles[i];
				for (int j = 0; j < testVariantAlleles.length; j++) {
					if (complementAlleles2[j].equals(allele1)) {
						nridenticalalleles++;
					}
				}
			}
		}

		VCFFunctions t = new VCFFunctions();

		String logoutputln = "";
		boolean flipped = false;
		if (refVariant.getAlleles().length == 2 && testVariant.getAlleles().length == 2) {
			// simple case: both are biallelic..
			// check if the minor alleles are equal. else, skip the variant.
			if (nridenticalalleles == 2) {
//				if (refAlleles[0].equals(BaseAnnot.getComplement(refAlleles[1]))
//						&& !refMinorAllele.equals(testVariantMinorAllele)) {
//					// both variants ar AT or GC snps
//					logoutputln += "\t-\t-\tAT or CG with DiffMinor";
//
//				} else

				if (testVariantMinorAllele.equals(refMinorAllele) || (testVariant.getMAF() > 0.45 && refVariant.getMAF() > 0.45)) {
					// check whether the reference allele is equal
					String[] tmpAlleles = testVariantAlleles;
					if (complement) {
						testVariant.convertAllelesToComplement();
						tmpAlleles = testVariant.getAlleles();
					}

					if (!refAlleles[0].equals(tmpAlleles[0])) {
						testVariant.flipReferenceAlelele();
						flipped = true;
					}

					logoutputln += "\t" + testVariant.getAlleles()[0] + "\t" + Strings.concat(testVariant.getAlleles(), Strings.comma, 1, testVariant.getAlleles().length);


					// merge
					String mergeStr = t.mergeVariants(refVariant, testVariant, separatorInMergedFile);
//					mergedOut.writeln(mergeStr);

					if (complement) {
						logoutputln += "\tOK-Complement";
					} else {
						logoutputln += "\tOK";
					}
					if (flipped) {
						logoutputln += "-flippedAlleles";
					}

					return new Pair<String, String>(logoutputln, mergeStr);

				} else {
					// write to log?
					logoutputln += "\t-\t-\tNotOK-DiffMinor";
					return new Pair<String, String>(logoutputln, null);
				}
			} else {
				// write to log?
				logoutputln += "\t-\t-\tNotOK-IncompatibleAlleles";
				return new Pair<String, String>(logoutputln, null);
			}

		} else if (nridenticalalleles > 1) {

			// recode the genotypes towards the joint set of alleles
			// get a list of all alleles at this locus...
			HashSet<String> uniqueAlleles = new HashSet<String>();
			uniqueAlleles.addAll(Arrays.asList(refAlleles));
			uniqueAlleles.addAll(Arrays.asList(testVariantAlleles));

			HashMap<String, Integer> alleleMap = new HashMap<String, Integer>();
			ArrayList<String> newAlleles = new ArrayList<String>();
			for (int i = 0; i < refAlleles.length; i++) {
				alleleMap.put(refAlleles[i], i);
				newAlleles.add(refAlleles[i]);
			}

			for (int i = 0; i < testVariantAlleles.length; i++) {
				String alleleStr = testVariantAlleles[i];
				if (!alleleMap.containsKey(alleleStr)) {
					alleleMap.put(alleleStr, alleleMap.size());
					newAlleles.add(alleleStr);
				}
			}


			// recode testVariant
			refVariant.recodeAlleles(alleleMap, newAlleles.toArray(new String[0]));
			testVariant.recodeAlleles(alleleMap, newAlleles.toArray(new String[0]));

			logoutputln += "\t" + testVariant.getAlleles()[0] + "\t" + Strings.concat(testVariant.getAlleles(), Strings.comma, 1, testVariant.getAlleles().length);

			// merge
			String mergeStr = t.mergeVariants(refVariant, testVariant, separatorInMergedFile);
			logoutputln += "\tOK-MultiAllelic-AllelesRecoded";
			return new Pair<String, String>(logoutputln, mergeStr);

		} else {
			// variant we can't fix
			logoutputln += "\t-\t-\tNotOK-CantFix";
			return new Pair<String, String>(logoutputln, null);
		}
	}


	public void merge(String vcf1, String vcf2, String out) throws IOException {

		VCFGenotypeData data1 = new VCFGenotypeData(vcf1);
		VCFGenotypeData data2 = new VCFGenotypeData(vcf2);

		ArrayList<String> samples1 = data1.getSamples();
		System.out.println(samples1.size() + " samples in " + vcf1);
		ArrayList<String> samples2 = data2.getSamples();
		System.out.println(samples2.size() + " samples in " + vcf2);

		HashMap<String, Integer> samples1Hash = new HashMap<String, Integer>();
		for (int i = 0; i < samples1.size(); i++) {
			samples1Hash.put(samples1.get(i), i);
		}

		ArrayList<String> sharedSamples = new ArrayList<String>();
		ArrayList<Integer> sharedSamplesIndex = new ArrayList<Integer>();
		boolean[] includeSample1 = new boolean[samples1.size()];
		for (int i = 0; i < samples2.size(); i++) {
			String sample2 = samples2.get(i);
			Integer sample1Index = samples1Hash.get(sample2);
			if (sample1Index != null) {
				sharedSamples.add(sample2);
				sharedSamplesIndex.add(sample1Index);
				includeSample1[sample1Index] = true;
			} else {
				sharedSamplesIndex.add(-1);
			}
		}

		System.out.println(sharedSamples.size() + " samples shared between VCFs");


		if (sharedSamples.size() == 0) {

			HashMap<String, VCFVariant> variantMap = new HashMap<String, VCFVariant>();
			while (data1.hasNext()) {
				VCFVariant var = data1.next();
				variantMap.put(var.toString(), var);
			}

			System.out.println(variantMap.size() + " variants loaded from: " + vcf1);

			ArrayList<String> mergedSamples = new ArrayList<String>();
			mergedSamples.addAll(samples1);
			mergedSamples.addAll(samples2);

			TextFile outf = new TextFile(out, TextFile.W);

			outf.writeln("##fileformat=VCFv4.1");
			String header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
			for (int i = 0; i < samples1.size(); i++) {
				header += "\t" + samples1.get(i);
			}
			for (int i = 0; i < samples2.size(); i++) {
				header += "\t" + samples2.get(i);
			}
			outf.writeln(header);

			HashSet<String> vars2 = new HashSet<String>();
			int shared = 0;
			int sharedWritten = 0;
			int vcf2specific = 0;

			TextFile mergelog = new TextFile(out + "-mergelog.txt.gz", TextFile.W);

			while (data2.hasNext()) {
				VCFVariant var2 = data2.next();
				VCFVariant var1 = variantMap.get(var2.toString());
				StringBuilder builder = new StringBuilder(mergedSamples.size() * 3 + 250);
				builder.append(var2.getChr());
				builder.append("\t").append(var2.getPos());
				builder.append("\t").append(var2.getId());


				if (var1 == null) {
					// only present in ds2
					builder.append("\t").append(var2.getAlleles()[0]);
					builder.append("\t").append(Strings.concat(var2.getAlleles(), Strings.comma, 1, var2.getAlleles().length));
					builder.append("\t").append(".");
					builder.append("\t").append(".");
					builder.append("\t").append(".");
					builder.append("\t").append("GT");

					for (int i = 0; i < samples1.size(); i++) {
						builder.append("\t").append("./.");
					}
					byte[][] alleles = var2.getGenotypeAllelesNew();
					for (int i = 0; i < samples2.size(); i++) {

						if (alleles[0][i] == -1) {
							builder.append("\t").append("./.");
						} else {
							builder.append("\t").append(alleles[0][i])
									.append("/")
									.append(alleles[1][i]);
						}
					}
					outf.writeln(builder.toString());
					vcf2specific++;
				} else {

					String logln = var1.getChr()
							+ "\t" + var1.getPos()
							+ "\t" + var1.getId()
							+ "\t" + var1.getAlleles()[0]
							+ "\t" + Strings.concat(var1.getAlleles(), Strings.comma, 1, var1.getAlleles().length)
							+ "\t" + var1.getMinorAllele()
							+ "\t" + var1.getMAF();
					logln += "\t" + var2.getId()
							+ "\t" + var2.getAlleles()[0]
							+ "\t" + Strings.concat(var2.getAlleles(), Strings.comma, 1, var2.getAlleles().length)
							+ "\t" + var2.getMinorAllele()
							+ "\t" + var2.getMAF();


					// merge
					Pair<String, String> outputpair = mergeVariants(var1, var2, "/");
					String mergeStr = outputpair.getLeft();
					logln += mergeStr;
					mergelog.writeln(logln);
					if (outputpair.getRight() != null) {
						outf.writeln(outputpair.getRight());
						sharedWritten++;
					}
					shared++;
				}
				vars2.add(var2.toString());
			}

			mergelog.close();


			System.out.println(vars2.size() + " variants in: " + vcf2);
			System.out.println(vcf2specific + " specific variants in " + vcf2);
			System.out.println(shared + " shared variants");
			double percWritten = (double) sharedWritten / shared;
			System.out.println(sharedWritten + " shared variants written (" + percWritten + ")");


			int vcf1specific = 0;
			Set<String> keyset = variantMap.keySet();

			for (String s : keyset) {
				if (!vars2.contains(s)) {

					VCFVariant var1 = variantMap.get(s);

					StringBuilder builder = new StringBuilder(mergedSamples.size() * 3 + 250);
					builder.append(var1.getChr());
					builder.append("\t").append(var1.getPos());
					builder.append("\t").append(var1.getId());

					// only present in ds1
					builder.append("\t").append(var1.getAlleles()[0]);
					builder.append("\t").append(Strings.concat(var1.getAlleles(), Strings.comma, 1, var1.getAlleles().length));
					builder.append("\t").append(".");
					builder.append("\t").append(".");
					builder.append("\t").append(".");
					builder.append("\t").append("GT");

					byte[][] alleles = var1.getGenotypeAllelesNew();
					for (int i = 0; i < samples1.size(); i++) {

						if (alleles[0][i] == -1) {
							builder.append("\t").append("./.");
						} else {
							builder.append("\t").append(alleles[0][i])
									.append("/")
									.append(alleles[1][i]);
						}
					}
					for (int i = 0; i < samples2.size(); i++) {
						builder.append("\t").append("./.");
					}
					outf.writeln(builder.toString());
					vcf1specific++;
				}
			}

			System.out.println(vcf1specific + " variants specific to: " + vcf1);
			outf.close();
		} else {
			System.out.println("Not supported yet");
		}

	}
}