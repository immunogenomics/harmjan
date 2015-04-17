package nl.harmjanwestra.ngs.GenotypeFormats;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.VCFVariantType;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.HWE;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 * Created by hwestra on 4/8/15.
 */
public class VCFFunctions {

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
								// System.out.println("allele 1 mising");
							} else {
								nrCalled++;
								allele1b = Byte.parseByte(allele1);

							}
							if (allele2.equals(".")) {
								// missing
								sampleAlleles[samplePos][1] = -1;
								// System.out.println("allele 2 mising");
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

//						if (adCol != -1) {
//							// Allelic depths for the ref and alt alleles in the order listed
//							String ad = sampleElems[adCol];
//							String[] adElems = ad.split(",");
//
//						}
//
//						if (dpCol != -1&& dpCol < sampleElems.length) {
//							// Approximate read depth (reads with MQ=255 or with bad mates are filtered)
//							String dp = sampleElems[dpCol];
//							String[] dpElems = dp.split(",");
//						}
//
						if (gqCol != -1 && gqCol < sampleElems.length) {
							// Genotype Quality
//							String gq = sampleElems[gqCol];
//							String[] gqElems = gq.split(",");
							if (gqCol >= sampleElems.length) {
								System.err.println("ERROR: GQCol>length: " + gqCol + "\t" + Strings.concat(sampleElems, Strings.comma));
							}
						}
//
//						if (plCol != -1) {
//							// Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification
//							String pl = sampleElems[plCol];
//							String[] plElems = pl.split(",");
//						}


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

		// determine call rates per sample

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

	public Pair<byte[][], String[]> loadVCFGenotypes(String vcf,
													 HashMap<String, Integer> sampleMap,
													 HashMap<Feature, Integer> variantMap) throws IOException {

		TextFile tf = new TextFile(vcf, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);

		int[] colToNewSample = null;

		int nrVariantsFound = 0;

		byte[][] genotypes = new byte[sampleMap.size()][variantMap.size()];
		String[] alleles = new String[variantMap.size()];

		int nrSamplesFound = 0;
		while (elems != null) {

			Chromosome chr = Chromosome.NA;
			if (elems[0].startsWith("##")) {
				// superheader

			} else if (elems[0].startsWith("#CHROM")) {
				// header
				// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1   SAMPLE
				// header
				colToNewSample = new int[elems.length];
				for (int i = 9; i < elems.length; i++) {
					// T1D/RA specific stuff
					String[] sampleNameElems = elems[i].split("/");
					String sample = "";
					if (sampleNameElems.length == 1) {

						sample = sampleNameElems[0];

					} else if (sampleNameElems.length == 2) {

						// split further
						String tmpSampleName = sampleNameElems[1];
						String[] tmpSampleNameElems = tmpSampleName.split("_");
						sample = tmpSampleNameElems[0];

					} else {

						System.err.println("Error: unexpected sample name");
					}


					Integer id = sampleMap.get(sample);
					// System.out.println(sample + "\t" + id);
					if (id == null) {
						colToNewSample[i] = -1;
					} else {
						colToNewSample[i] = id;
						nrSamplesFound++;
					}
				}
			} else {

				// line
				String ref = elems[3];
				String alt = elems[4];
				VCFVariantType type = VCFVariantType.parseType(ref, alt);
				boolean biAllelic = type.isBiallelic();
				if (biAllelic) {
					chr = Chromosome.parseChr(elems[0]);

					// determine MAF
					Integer pos = Integer.parseInt(elems[1]);
					Feature f = new Feature();
					f.setChromosome(chr);
					f.setStart(pos);
					f.setStop(pos);
					Integer newVariantId = variantMap.get(f);
					if (newVariantId == null) {
						// don't do anything just yet

					} else {
						nrVariantsFound++;

						// load the genotypes
						String[] alternates = alt.split(",");
						String[] currentalleles = new String[1 + alternates.length];
						currentalleles[0] = ref.intern();
						for (int i = 0; i < alternates.length; i++) {
							currentalleles[i + 1] = alternates[i].intern();
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

						alleles[newVariantId] = ref + "/" + alt;

						for (int e = 9; e < elems.length; e++) {
							int newSampleId = colToNewSample[e];
							if (newSampleId > -1) {
								String sampleColumn = elems[e];
								if (sampleColumn.equals("./.")) {
									genotypes[newSampleId][newVariantId] = -1;
								} else {
									String[] sampleElems = sampleColumn.split(":");

									if (gtCol != -1) {
										String gt = sampleElems[gtCol];
										String[] gtElems = gt.split("/");

										if (gtElems[0].equals(".") || gtElems[1].equals(".")) {
											genotypes[newSampleId][newVariantId] = -1;
										} else {
											Byte g1 = Byte.parseByte(gtElems[0]);
											Byte g2 = Byte.parseByte(gtElems[1]);

											byte dosage = g1;
											dosage += g2;
											genotypes[newSampleId][newVariantId] = dosage;
										}


									}
								}
							}
						}
					}
				}
			}


			elems = tf.readLineElems(TextFile.tab);
		}

		System.out.println("genotypes loaded for " + nrVariantsFound + " variants using VCF for " + nrSamplesFound + " samples");
		tf.close();
		return new Pair<byte[][], String[]>(genotypes, alleles);
	}

	public ArrayList<Feature> getVariantsFromVCF(String vcf) throws IOException {
		TextFile tf = new TextFile(vcf, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		ArrayList<Feature> output = new ArrayList<Feature>();
		while (elems != null) {

			Chromosome chr = Chromosome.NA;
			if (elems[0].startsWith("##")) {
				// superheader

			} else if (elems[0].startsWith("#CHROM")) {
				// header
				// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1   SAMPLE
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
				boolean biAllelic = type.isBiallelic();
				if (biAllelic) {
					// determine MAF
					Integer pos = Integer.parseInt(elems[1]);
					Feature f = new Feature();
					if (!chr.equals(Chromosome.X)) {
						f.setChromosome(chr);
						f.setStart(pos);
						f.setStop(pos);
						output.add(f);
					}
				}

			}


			elems = tf.readLineElems(TextFile.tab);
		}


		tf.close();
		return output;
	}

	public ArrayList<String> getVCFSamples(String vcf) throws IOException {

		ArrayList<String> samples = new ArrayList<String>();
		TextFile tf = new TextFile(vcf, TextFile.R);
		String[] lnElems = tf.readLineElems(TextFile.tab);
		boolean headerFound = false;
		while (lnElems != null && !headerFound) {

			if (lnElems[0].startsWith("##")) {
				// metadata
			} else if (lnElems[0].startsWith("#CHROM")) {

				// header
				for (int i = 9; i < lnElems.length; i++) {
					// T1D/RA specific stuff
					String[] sampleNameElems = lnElems[i].split("/");
					if (sampleNameElems.length == 1) {
						samples.add(sampleNameElems[0]);
					} else if (sampleNameElems.length == 2) {

						// split further
						String tmpSampleName = sampleNameElems[1];
						String[] tmpSampleNameElems = tmpSampleName.split("_");
						samples.add(tmpSampleNameElems[0]);

					} else {

						System.err.println("Error: unexpected sample name");
					}
				}
				headerFound = true;


			}

			lnElems = tf.readLineElems(TextFile.tab);
		}

		tf.close();
		return samples;
	}

	private void replaceVCFHeaderWithImmunoChipIDs(String vcfStart, String seqIdToImmunoChipID, String vcfImmunoChipIDs) throws IOException {

		TextFile tf = new TextFile(seqIdToImmunoChipID, TextFile.R);
		Map<String, String> map = tf.readAsHashMap(0, 1);
		tf.close();


		TextFile tf2 = new TextFile(vcfStart, TextFile.R);
		TextFile tf3 = new TextFile(vcfImmunoChipIDs, TextFile.W);
		String ln = tf2.readLine();
		while (ln != null) {


			if (ln.startsWith("#CHROM")) {

				String[] elems = ln.split("\t");
				for (int i = 9; i < elems.length; i++) {

					String[] sampleElems = elems[i].split("/");
					String sampleName = sampleElems[sampleElems.length - 1];
					String replacement = map.get(sampleName);

					if (replacement != null) {
						System.out.println(elems[i] + " --> " + replacement);
						elems[i] = replacement;
					} else {
						elems[i] = sampleName;
						System.out.println(elems[i] + " --> " + sampleName);
					}
				}
				ln = Strings.concat(elems, Strings.tab);
			}

			if (!ln.startsWith("##GATK")) {
				tf3.writeln(ln);
			}


			ln = tf2.readLine();
		}
		tf3.close();
		tf2.close();

	}


	public void filterLowFrequencyVariants(String sequencingVCF, String outputdir, boolean filterGT) throws IOException {

		TextFile in = new TextFile(sequencingVCF, TextFile.R);
		String ln = in.readLine();
		int nrHeaderElems = 9;

		int minimalReadDepth = 10;
		int minimalGenotypeQual = 30;
		int minimalGenotypePHRED = 10;
		if (!filterGT) {
			minimalReadDepth = 0;
			minimalGenotypeQual = 0;
			minimalGenotypePHRED = Integer.MAX_VALUE;
		}

		double callratethreshold = 0.95;
		double mafthreshold = 0.01;


		int[] allelesCalledPerSample = new int[0];
		ArrayList<String> sampleNames = new ArrayList<String>();

		int[] readDepthDist = new int[500];
		int[] genotypeQualDist = new int[101];
		int[] genotypePHREDQualDist = new int[1000];
		double minAlleleFrequency = 0.01;
		double minObservationsPerAllele = 5d;

		int nrVariants = 0;
		TextFile variantQC = new TextFile(outputdir + "variantqc.txt", TextFile.W);
		String head = "Pos\tRs\tFilterOut\tReason\tBiallelic\tRef\tAlt\tMinorAllele\tMinorAlleleFreq\tNrCalled\tCallRate\tEachAlleleObservedFreq>" + minAlleleFrequency + "\tNumberOfTimesEachAlleleIsObserved\tAlleleFrequencies\tQual\tFilter\tInfo";
		variantQC.writeln(head);

		TextFile filteredVCF = new TextFile(outputdir + "filtered.vcf", TextFile.W);

		while (ln != null) {
			if (ln.startsWith("##")) {
// metadata
				filteredVCF.writeln(ln);
			} else if (ln.startsWith("#CHROM")) {
// header
				filteredVCF.writeln(ln);
				String[] elems = ln.split("\t");
				allelesCalledPerSample = new int[elems.length - nrHeaderElems];
				for (int i = nrHeaderElems; i < elems.length; i++) {
					sampleNames.add(elems[i]);
				}
				minObservationsPerAllele = (int) Math.ceil(minAlleleFrequency * sampleNames.size());
				mafthreshold = minObservationsPerAllele / sampleNames.size();
				System.out.println(minObservationsPerAllele);
			} else {
				nrVariants++;
				String[] elems = ln.split("\t");
				String variant = elems[0] + "_" + elems[1];
				String rs = elems[2];
				// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
				String ref = elems[3];
				String alt = elems[4];
				String[] alternateAlleles = alt.split(",");
				String qual = elems[5];
				String filter = elems[6];
				String info = elems[7];

				String[] format = elems[8].split(":");
				int gtCol = -1; // genotype
				int adCol = -1; // Allelic depths for the ref and alt alleles in the order listed
				int dpCol = -1; // Approximate read depth (reads with MQ=255 or with bad mates are filtered)
				int gqCol = -1; // Genotype Quality
				int plCol = -1; // Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification
				int pidCol = -1; // ?
				int pgtCol = -1; // ?
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
					} else if (format[c].equals("PGT")) {
						pgtCol = c;
					} else if (format[c].equals("PID")) {
						pidCol = c;
					}

					// GT:AD:DP:GQ:PGT:PID:PL

				}

				int nrCalled = 0;

				// count nr possible alleles
				int[] nrAllelesObserved = new int[alternateAlleles.length + 1];

				// count number of alleles with certain read depth
				for (int i = nrHeaderElems; i < elems.length; i++) {
					String sampleColumn = elems[i];
					if (sampleColumn.equals("./.")) {
						// not called
					} else {
						String[] sampleElems = sampleColumn.split(":");
						if (gtCol != -1) {
							String gt = sampleElems[gtCol];
							String[] gtElems = gt.split("/");

							int depth = 0;
							int gq = 0;


							if (gtElems[0].equals(".") || gtElems[1].equals(".")) {
								// not called
							} else {
								try {
									if (gqCol != -1) {

										gq = Integer.parseInt(sampleElems[gqCol]);


									}
								} catch (NumberFormatException e) {
								}

								if (gq > genotypeQualDist.length - 1) {
									genotypeQualDist[genotypeQualDist.length - 1]++;
								} else {
									genotypeQualDist[gq]++;
								}

								try {
									if (dpCol != -1) {
										if (sampleElems.length < dpCol + 1) {
											System.out.println(sampleColumn);
											System.exit(-1);
										}
										depth = Integer.parseInt(sampleElems[dpCol]);
									}
								} catch (NumberFormatException e) {
								}

								if (depth > readDepthDist.length - 1) {
									readDepthDist[readDepthDist.length - 1]++;
								} else {
									readDepthDist[depth]++;
								}

								int gt1 = 0;
								int gt2 = 0;
								try {
									gt1 = Integer.parseInt(gtElems[0]);
									gt2 = Integer.parseInt(gtElems[1]);
								} catch (NumberFormatException e) {
								}

								// called
								// check coverage
								if (depth < minimalReadDepth) {
									// not called
								} else {
									// called
									// check the genotype quality
									if (gq < minimalGenotypeQual) {
										// not called
									} else {
										// called
										int gtqual0 = 0;
										String[] plElems = sampleElems[plCol].split(",");
										try {

											int genotype = gt1 + gt2;
											// 0,99,1238 AA AB BB
											if (gt1 != gt2) {
												// heterozygote
											}
											gtqual0 = Integer.parseInt(plElems[genotype]);

										} catch (NumberFormatException e) {
										}
										if (gtqual0 > genotypePHREDQualDist.length - 1) {
											genotypePHREDQualDist[genotypePHREDQualDist.length - 1]++;
										} else {
											genotypePHREDQualDist[gtqual0]++;
										}
										if (gtqual0 > minimalGenotypePHRED) {
											// not called
										} else {
											// called.. this is the final genotype
											allelesCalledPerSample[i - nrHeaderElems]++;
											nrCalled++;
											nrAllelesObserved[gt1]++;
											nrAllelesObserved[gt2]++;
											if (variant.equals("X_153330275")) {
												System.out.println(gt1 + "\t" + gt2);
											}
										}
									}
								}
							}
						}
					}
				}


				double cr = (double) nrCalled / (elems.length - nrHeaderElems);

				String reason = "";
				boolean filterout = false;
				if (cr < callratethreshold) {
					filterout = true;
					reason += "lowCallRate";
				}

				boolean allAllelesObservedFrequently = true;
				int totalAllelesObs = nrCalled * 2;
				String[] allelefreqOut = new String[nrAllelesObserved.length];
				String[] allelefreqOut2 = new String[nrAllelesObserved.length];

				int nrAllelesAboveFreq = 0;
				int nrAllelesThatHaveAlleleFrequency = 0;
				double minAlleleFreq = 1;
				String minorAllele = "";
				for (int i = 0; i < nrAllelesObserved.length; i++) {
					double alleleFreq = (double) nrAllelesObserved[i] / totalAllelesObs;
					allelefreqOut[i] = "" + alleleFreq;
					allelefreqOut2[i] = "" + nrAllelesObserved[i];
					if (nrAllelesObserved[i] > 0) {
						nrAllelesThatHaveAlleleFrequency++;
						if (alleleFreq < minAlleleFreq) {
							if (i == 0) {
								minorAllele = ref;
							} else {
								minorAllele = alternateAlleles[i - 1];
							}
							minAlleleFreq = alleleFreq;
						}
					}
					if (alleleFreq > minAlleleFrequency) {
						nrAllelesAboveFreq++;
					}
				}
				boolean isBiAllelic = false;
				if (nrAllelesThatHaveAlleleFrequency == 1) {
					filterout = true; // monomorphic variants
					reason += ";monomorphic";
				} else {
					if (nrAllelesThatHaveAlleleFrequency == 2) {
						isBiAllelic = true;
						if (minAlleleFreq < mafthreshold) {
							filterout = true;
							reason += ";mafbelowthreshold";
						}
					}
				}
//
// else if (nrAllelesThatHaveAlleleFrequency == 2) {
//					if (nrAllelesAboveFreq == 2) {
//						allAllelesObservedFrequently = true; // biallelic variants
//						if (minAlleleFreq < mafthreshold) {
//							filterout = true;
//							reason += ";biallelic-belowmafthreshold";
//						}
//					}
//					isBiAllelic = true;
//
//				} else if (nrAllelesObserved.length > 2 && nrAllelesAboveFreq >= 2) {
//					if (nrAllelesAboveFreq == 2) {
//						// this is actually a biallelic
//						isBiAllelic = true;
//						if (minAlleleFreq < mafthreshold) {
//							filterout = true;
//							reason += ";multi-biallelic-belowmafthreshold";
//						}
//					} else {
//						if (nrAllelesAboveFreq > 2) {
//							allAllelesObservedFrequently = true; // multi-allelic variants that
//						} else {
//							filterout = true;
//							reason += ";multi-biallelic-belowobservedallelesthreshold";
//						}
//					}
//
//				}

				String variantOut = variant
						+ "\t" + rs
						+ "\t" + filterout
						+ "\t" + reason
						+ "\t" + isBiAllelic
						+ "\t" + ref
						+ "\t" + alt
						+ "\t" + minorAllele
						+ "\t" + minAlleleFreq
						+ "\t" + nrCalled
						+ "\t" + cr
						+ "\t" + allAllelesObservedFrequently
						+ "\t" + Strings.concat(allelefreqOut2, Strings.semicolon)
						+ "\t" + Strings.concat(allelefreqOut, Strings.semicolon);

				variantQC.writeln(variantOut + "\t" + qual + "\t" + filter + "\t" + info);
				if (!filterout) {
					filteredVCF.writeln(ln);
				}
			}

			ln = in.readLine();
		}
		variantQC.close();
		in.close();
		filteredVCF.close();

		TextFile outf1 = new TextFile(outputdir + "sampleCallRate.txt", TextFile.W);
		for (int i = 0; i < allelesCalledPerSample.length; i++) {
			double cr = (double) allelesCalledPerSample[i] / nrVariants;
			outf1.writeln(sampleNames.get(i) + "\t" + allelesCalledPerSample[i] + "\t" + cr);
		}
		outf1.close();

		TextFile outf3 = new TextFile(outputdir + "gq.txt", TextFile.W);
		for (int i = 0; i < genotypeQualDist.length; i++) {
			outf3.writeln(i + "\t" + genotypeQualDist[i]);
		}
		outf3.close();

		TextFile outf4 = new TextFile(outputdir + "gqphred.txt", TextFile.W);
		for (int i = 0; i < genotypePHREDQualDist.length; i++) {
			outf4.writeln(i + "\t" + genotypePHREDQualDist[i]);
		}
		outf4.close();

		TextFile outf5 = new TextFile(outputdir + "readdepth.txt", TextFile.W);
		for (int i = 0; i < readDepthDist.length; i++) {
			outf5.writeln(i + "\t" + readDepthDist[i]);
		}
		outf5.close();


	}

	public void filterVCF(String vcfIn, String vcfVariantsToFilter, String vcfOut) throws IOException {


		TextFile in = new TextFile(vcfVariantsToFilter, TextFile.R);
		String ln = in.readLine();
		HashSet<String> variantsToFilter = new HashSet<String>();
		while (ln != null) {
			if (ln.startsWith("##")) {
// metadata

			} else if (ln.startsWith("#CHROM")) {
// header

			} else {
				String[] elems = ln.split("\t");
				String variant = elems[0] + "_" + elems[1];
				variantsToFilter.add(variant);
			}

			ln = in.readLine();
		}
		in.close();


		TextFile in2 = new TextFile(vcfVariantsToFilter, TextFile.R);
		ln = in2.readLine();
		TextFile filteredVCF = new TextFile(vcfOut, TextFile.W);
		while (ln != null) {
			if (ln.startsWith("##")) {
// metadata
				filteredVCF.writeln(ln);
			} else if (ln.startsWith("#CHROM")) {
// header
				filteredVCF.writeln(ln);
			} else {
				String[] elems = ln.split("\t");
				String variant = elems[0] + "_" + elems[1];
				if (!variantsToFilter.contains(variant)) {
					filteredVCF.writeln(ln);
				}
			}

			ln = in2.readLine();
		}
		in2.close();
		filteredVCF.close();
	}

	public void mergeWithPed(String ped, String vcfIn, String vcfOut) throws IOException {

		PedAndMapFunctions pedAndMapFunctions = new PedAndMapFunctions();
		TextFile tf = new TextFile(ped + ".map", TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);
		ArrayList<String> rsNames = new ArrayList<String>();
		ArrayList<String> pos = new ArrayList<String>();
		HashMap<Feature, Integer> variantMap = new HashMap<Feature, Integer>();
		int ctr = 0;
		while (elems != null) {
			rsNames.add(elems[1]);
			pos.add(elems[0] + "_" + elems[3]);
			Feature f = new Feature();
			f.setChromosome(Chromosome.parseChr(elems[0]));
			f.setStart(Integer.parseInt(elems[3]));
			f.setStop(Integer.parseInt(elems[3]));
			if (variantMap.containsKey(f)) {
				System.out.println(elems[1]);
			}
			variantMap.put(f, ctr);
			ctr++;
			elems = tf.readLineElems(TextFile.tab);
		}
		System.out.println(variantMap.size() + "\t" + rsNames.size());

		Pair<byte[][], String[]> pedGenotypes = null;

		TextFile vcftf = new TextFile(vcfIn, TextFile.R);
		int nrHeaderElems = 9;
		String ln = vcftf.readLine();
		TextFile vcfoutf = new TextFile(vcfOut, TextFile.W);
		while (ln != null) {
			if (ln.startsWith("##")) {
				vcfoutf.writeln(ln);
			} else if (ln.startsWith("#CHROM")) {
				String[] sampleElems = ln.split("\t");
				int ctr2 = 0;
				HashMap<String, Integer> sampleMap = new HashMap<String, Integer>();
				for (int i = nrHeaderElems; i < sampleElems.length; i++) {
					sampleMap.put(sampleElems[i], ctr2);
					ctr2++;
				}
				pedGenotypes = pedAndMapFunctions.loadPedGenotypes(ped + ".ped", ped + ".map", sampleMap, variantMap);
				vcfoutf.writeln(ln);
			} else {
				vcfoutf.writeln(ln);
			}
			ln = vcftf.readLine();
		}


		// now write the ped genotypes
		if (pedGenotypes != null) {
			byte[][] genotypes = pedGenotypes.getLeft(); // [samples][variants]
			String[] alleles = pedGenotypes.getRight();

			// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
			for (int variant = 0; variant < alleles.length; variant++) {
				String posStr = pos.get(variant);
				String[] posStrElems = posStr.split("_");

				String variantAlleles = alleles[variant];
				String[] variantAllelesElems = variantAlleles.split("/");

				if (!variantAllelesElems[0].equals("null") && !variantAllelesElems[1].equals("null")) {
					String outln = posStrElems[0]
							+ "\t" + posStrElems[1]
							+ "\t" + rsNames.get(variant)
							+ "\t" + variantAllelesElems[0]
							+ "\t" + variantAllelesElems[1]
							+ "\t99"
							+ "\t."
							+ "\t."
							+ "\tGT";

					int allelesA = 0;
					int nrCalled = 0;

					for (int i = 0; i < genotypes.length; i++) {

						byte gt = genotypes[i][variant];
						if (gt == -1) {
							outln += "\t./.";
						} else {
							nrCalled++;
							if (gt == 0) {
								allelesA += 2;
								outln += "\t0/0";
							} else if (gt == 1) {
								allelesA++;
								outln += "\t0/1";
							} else {
								outln += "\t1/1";
							}
						}

					}

					double af = (double) allelesA / (nrCalled * 2);
					if (af > 0.5) {
						af = 1 - af;
					}
					if (af > (5d / genotypes.length)) {
						double cr = (double) nrCalled / genotypes.length;
						if (cr >= 0.95) {
							vcfoutf.writeln(outln);
						}
					}


				}


			}


		}


		vcftf.close();
		vcfoutf.close();


	}

	public void convertPEDToVCF(String ped, String vcf) throws IOException {

		PedAndMapFunctions pedAndMapFunctions = new PedAndMapFunctions();

		TextFile tf = new TextFile(ped + ".map", TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);
		ArrayList<String> rsNames = new ArrayList<String>();
		ArrayList<String> pos = new ArrayList<String>();
		HashMap<Feature, Integer> variantMap = new HashMap<Feature, Integer>();
		int ctr = 0;
		while (elems != null) {
			rsNames.add(elems[1]);
			pos.add(elems[0] + "_" + elems[3]);
			Feature f = new Feature();
			f.setChromosome(Chromosome.parseChr(elems[0]));
			f.setStart(Integer.parseInt(elems[3]));
			f.setStop(Integer.parseInt(elems[3]));
			if (variantMap.containsKey(f)) {
				System.out.println(elems[1]);
			}
			variantMap.put(f, ctr);
			ctr++;
			elems = tf.readLineElems(TextFile.tab);
		}
		System.out.println(variantMap.size() + "\t" + rsNames.size());


		// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
		TextFile vcfout = new TextFile(vcf, TextFile.W);
		vcfout.writeln("##fileformat=VCFv4.1");
		vcfout.writeln("##source=HJ");
		String vcfheader = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";


		// inventorize samples the quick and dirty way.

		TextFile tfped = new TextFile(ped + ".ped", TextFile.R);

		String[] pedelems = tfped.readLineElems(Strings.whitespace);
		ArrayList<String> samples = new ArrayList<String>();
		HashMap<String, Integer> sampleMap = new HashMap<String, Integer>();
		int samplectr = 0;
		while (pedelems != null) {
			samples.add(pedelems[1]);
			vcfheader += "\t" + pedelems[1];
			sampleMap.put(pedelems[1], samplectr);
			samplectr++;
			pedelems = tfped.readLineElems(Strings.whitespace);
		}

		tfped.close();

		vcfout.writeln(vcfheader);
		Pair<byte[][], String[]> pedGenotypes = pedAndMapFunctions.loadPedGenotypes(ped + ".ped", ped + ".map", sampleMap, variantMap);

		byte[][] genotypes = pedGenotypes.getLeft(); // [samples][variants]
		String[] alleles = pedGenotypes.getRight();

		for (int variant = 0; variant < alleles.length; variant++) {
			String posStr = pos.get(variant);
			String[] posStrElems = posStr.split("_");

			String variantAlleles = alleles[variant];
			String[] variantAllelesElems = variantAlleles.split("/");

			if (!variantAllelesElems[0].equals("null") && !variantAllelesElems[1].equals("null")) {
				String outln = posStrElems[0]
						+ "\t" + posStrElems[1]
						+ "\t" + rsNames.get(variant)
						+ "\t" + variantAllelesElems[0]
						+ "\t" + variantAllelesElems[1]
						+ "\t."
						+ "\tPASS"
						+ "\t."
						+ "\tGT";

				int allelesA = 0;
				int nrCalled = 0;

				for (int i = 0; i < genotypes.length; i++) {

					byte gt = genotypes[i][variant];
					if (gt == -1) {
						outln += "\t./.";
					} else {
						nrCalled++;
						if (gt == 0) {
							allelesA += 2;
							outln += "\t0/0";
						} else if (gt == 1) {
							allelesA++;
							outln += "\t0/1";
						} else {
							outln += "\t1/1";
						}
					}

				}

				double af = (double) allelesA / (nrCalled * 2);
				if (af > 0.5) {
					af = 1 - af;
				}
				if (af > (5d / genotypes.length)) {
					double cr = (double) nrCalled / genotypes.length;
					if (cr >= 0.95) {
						vcfout.writeln(outln);
					}
				}
			}
		}
		vcfout.close();
	}

	public void filterBeagleReference(String vcfIn, String referenceOut, String regionFile) throws IOException {

		PedAndMapFunctions fun = new PedAndMapFunctions();
		ArrayList<Feature> set = fun.readRegionFile(regionFile);

		TextFile vcftf = new TextFile(vcfIn, TextFile.R);
		int nrHeaderElems = 9;
		String ln = vcftf.readLine();
		TextFile vcfoutf = new TextFile(referenceOut, TextFile.W);
		while (ln != null) {
			if (ln.startsWith("##")) {
				vcfoutf.writeln(ln);
			} else if (ln.startsWith("#CHROM")) {
				vcfoutf.writeln(ln);
			} else {

				String[] elems = ln.split("\t");
				Feature f = new Feature();
				f.setChromosome(Chromosome.parseChr(elems[0]));
				f.setStart(Integer.parseInt(elems[1]));
				f.setStop(Integer.parseInt(elems[1]));
				boolean overlap = false;
				for (Feature region : set) {
					if (f.overlaps(region)) {
						overlap = true;
						break;
					}
				}

				if (overlap) {
					vcfoutf.writeln(ln);
				}
			}
			ln = vcftf.readLine();
		}
		vcfoutf.close();
		vcftf.close();
	}
}
