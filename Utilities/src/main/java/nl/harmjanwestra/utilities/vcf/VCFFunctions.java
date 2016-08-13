package nl.harmjanwestra.utilities.vcf;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import com.gs.collections.impl.multimap.list.FastListMultimap;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.Gender;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.genotypes.GenotypeTools;
import nl.harmjanwestra.utilities.individuals.Individual;
import nl.harmjanwestra.utilities.plink.PedAndMapFunctions;
import nl.harmjanwestra.utilities.plink.PlinkFamFile;
import nl.harmjanwestra.utilities.sets.StringSets;
import nl.harmjanwestra.utilities.shell.ProcessRunner;
import nl.harmjanwestra.utilities.vcf.filter.GenotypeQualityFilter;
import nl.harmjanwestra.utilities.vcf.filter.ReadDepthFilter;
import nl.harmjanwestra.utilities.vcf.filter.VCFGenotypeFilter;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.HWE;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.*;

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
				int dpCol = -1; // Approximate readAsTrack depth (reads with MQ=255 or with bad mates are filtered)
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
							if (gtElems.length == 1) {
								gtElems = gt.split("\\|");
							}
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
//							// Approximate readAsTrack depth (reads with MQ=255 or with bad mates are filtered)
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
								System.out.println("New sample detected: " + sample + " in path: " + vcffile);
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
			System.out.println("Total samples: " + sampleToId.size() + " after path: " + vcffile);
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
							System.err.println("Error in path - missing sample:" + vcffile);
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

	public void rewriteVariantsAtSamePositionAsMultiAllelic(String vcfIn, String vcfOut) throws IOException {

		TextFile tf = new TextFile(vcfIn, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);

		HashSet<Feature> visitedPositions = new HashSet<Feature>();
		HashSet<Feature> nonUniqueFeatures = new HashSet<Feature>();
		HashMap<Feature, String> posToName = new HashMap<Feature, String>();
		HashMap<Feature, Pair<String, HashSet<String>>> posToAlleles = new HashMap<Feature, Pair<String, HashSet<String>>>();

		while (elems != null) {

			Chromosome chr = Chromosome.NA;
			if (elems[0].startsWith("##")) {
				// superheader

			} else if (elems[0].startsWith("#CHROM")) {

			} else {

				String ref = new String(elems[3]).intern();
				String alt = elems[4];
				String[] alleles = alt.split(",");

				ArrayList<String> allelesList = new ArrayList<String>();
				for (String al : alleles) {
					allelesList.add(new String(al).intern());
				}

				String name = elems[2];
				Integer pos = Integer.parseInt(elems[1]);
				chr = Chromosome.parseChr(elems[0]);

				Feature f = new Feature();
				f.setChromosome(chr);
				f.setStart(pos);
				f.setStop(pos);

				if (visitedPositions.contains(f)) {
					nonUniqueFeatures.add(f);
					String name2 = posToName.get(f);
					if (!name.equals(name2)) {
						System.out.println("different rs id for same pos: " + name2);
					}

					// get alleles
					Pair<String, HashSet<String>> alleles2 = posToAlleles.get(f);
					String ref2 = alleles2.getLeft();
					if (!ref.equals(ref2)) {
						System.out.println("reference alleles different at same base pos: " + ref + "\t" + ref2 + "\t" + elems[0] + ":" + pos);
					}

					for (String s : alleles) {
						alleles2.getRight().add(new String(s).intern());
					}


				} else {
					HashSet<String> alleleHash = new HashSet<String>();
					Pair<String, HashSet<String>> allelePair = new Pair<String, HashSet<String>>(ref, alleleHash);

					for (String s : alleles) {
						alleleHash.add(new String(s).intern());
					}

					posToAlleles.put(f, allelePair);
				}


				visitedPositions.add(f);
				posToName.put(f, name);


			}

			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();


		System.out.println(nonUniqueFeatures.size() + " multi-allelic variants");


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
						int dpCol = -1; // Approximate readAsTrack depth (reads with MQ=255 or with bad mates are filtered)
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
										if (gtElems.length == 1) {
											gtElems = gt.split("\\|");
										}

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

	public ArrayList<Feature> getVariantsFromVCF(String vcf, boolean includeVariantName, boolean biAllelicOnly) throws IOException {
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
				if (chr.equals(Chromosome.NA)) {
					System.out.println("Could not parse chr: " + elems[0] + " for variant: " + elems[2]);
				}
				String[] alternates = alt.split(",");
				String[] alleles = new String[1 + alternates.length];
				alleles[0] = ref.intern();
				for (int i = 0; i < alternates.length; i++) {
					alleles[i + 1] = alternates[i].intern();
				}

				VCFVariantType type = VCFVariantType.parseType(ref, alt);
				boolean biAllelic = type.isBiallelic();

				// determine MAF
				Integer pos = Integer.parseInt(elems[1]);
				Feature f = new Feature();
				if (!biAllelicOnly || (biAllelic && biAllelicOnly)) {
//					if (!chr.equals(Chromosome.X)) {
					f.setChromosome(chr);

					if (includeVariantName) {
						f.setName(elems[2]);
					}

					f.setStart(pos);
					f.setStop(pos);
					output.add(f);
//					}
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


	public void filterLowFrequencyVariants(String sequencingVCF,
										   String outputdir,
										   String famfile, boolean filterGT,
										   int minimalReadDepth,
										   int minimalGenotypeQual,
										   double callratethreshold,
										   int minObservationsPerAllele) throws IOException {


		System.out.println("Filtering for low freq variants: " + sequencingVCF);
		TextFile in = new TextFile(sequencingVCF, TextFile.R);
		String ln = in.readLine();
		int nrHeaderElems = 9;

		double mafthreshold = 1;
		int minimalGenotypePHRED = 10;

		if (!filterGT) {
			minimalReadDepth = 0;
			minimalGenotypeQual = 0;
			minimalGenotypePHRED = Integer.MAX_VALUE;
		}

		ArrayList<String> sampleNames = new ArrayList<String>();

		int[] readDepthDist = new int[500];
		int[] genotypeQualDist = new int[101];
		int[] genotypePHREDQualDist = new int[1000];


		int[] allelesCalledPerSample = null;
		int nrVariants = 0;
		TextFile variantQC = new TextFile(outputdir + "variantqc.txt", TextFile.W);
		String head = "Chr\t" +
				"Pos\t" +
				"Rs\t" +
				"FilterOut\t" +
				"Reason\t" +
				"Biallelic\t" +
				"Multiallelic\t" +
				"Ref\t" +
				"Alt\t" +
				"MinorAllele\t" +
				"MinorAlleleFreq\t" +
				"NrCalled\t" +
				"CallRate\t" +
				"NumberOfTimesEachAlleleIsObserved\t" +
				"AlleleFrequencies\t" +
				"Qual\t" +
				"Filter\t" +
				"Info";
		variantQC.writeln(head);


		TextFile filteredVCF = new TextFile(outputdir + "filtered.vcf.gz", TextFile.W);

		HashMap<String, Integer> sampleMap = new HashMap<String, Integer>();
		Gender[] sampleIsFemale = null;
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
					sampleMap.put(elems[i], i - nrHeaderElems);
				}

				if (famfile != null) {
					PlinkFamFile pf = new PlinkFamFile(famfile);
					sampleIsFemale = new Gender[sampleMap.size()];
					for (Individual ind : pf.getSamples()) {
						Integer id = sampleMap.get(ind.getName());
						if (id != null) {
							sampleIsFemale[id] = ind.getGender();
						}
					}
				}

				mafthreshold = (double) minObservationsPerAllele / sampleNames.size();
				System.out.println("maf threshold: " + mafthreshold + " min obs: " + minObservationsPerAllele);
			} else {
				nrVariants++;

				ArrayList<VCFGenotypeFilter> filters = new ArrayList<>();
				if (minimalGenotypeQual > 0) {
					filters.add(new GenotypeQualityFilter(minimalGenotypeQual));
				}
				if (minimalReadDepth > 0) {
					filters.add(new ReadDepthFilter(minimalReadDepth));
				}

				VCFVariant variant = new VCFVariant(ln, filters, true);

				short[] depths = variant.getApproximateDepth();

				for (int depth : depths) {
					if (depth > readDepthDist.length - 1) {
						readDepthDist[readDepthDist.length - 1]++;
					} else {
						readDepthDist[depth]++;
					}
				}

				short[] gqs = variant.getGenotypeQuals();
				for (int gq : gqs) {
					if (gq > genotypeQualDist.length - 1) {
						genotypeQualDist[genotypeQualDist.length - 1]++;
					} else {
						genotypeQualDist[gq]++;
					}
				}


//				if (gtqual0 > genotypePHREDQualDist.length - 1) {
//					genotypePHREDQualDist[genotypePHREDQualDist.length - 1]++;
//				} else {
//					genotypePHREDQualDist[gtqual0]++;
//				}

				String reason = "";
				boolean filterout = false;
				String alt = "";
				Chromosome chr = Chromosome.parseChr(variant.getChr());
				if (sampleIsFemale != null && (chr.equals(Chromosome.X) || chr.equals(Chromosome.Y))) {
					variant.recalculateMAFAndCallRate(sampleIsFemale, null);
				}

				if (variant.getCallrate() < callratethreshold) {
					filterout = true;
					reason += "lowCallRate";
				}
				if (variant.isMonomorphic()) {
					filterout = true; // monomorphic variants
					reason += ";monomorphic";
					alt = "null";
				} else if (variant.isBiallelic()) {
					double maf = variant.getMAF();
					if (maf < mafthreshold) {
						filterout = true;
						reason += ";mafbelowthreshold";
					}

					alt = Strings.concat(variant.getAlleles(), Strings.comma, 1, variant.getAlleles().length);

				} else {

					double maf = variant.getMAF();
					if (maf < mafthreshold) {
						filterout = true;
						reason += ";mafbelowthreshold";
					}
					alt = Strings.concat(variant.getAlleles(), Strings.comma, 1, variant.getAlleles().length);

				}

				double maf = variant.getMAF();
				double[] alleleFrequencies = variant.getAlleleFrequencies();


				DoubleMatrix2D genotypes = variant.getGenotypeAllelesAsMatrix2D();
				int[] variantAllelesObserved = variant.getNrAllelesObserved();
				int nrCalled = 0;
				for (int i = 0; i < genotypes.rows(); i++) {
					if (genotypes.getQuick(i, 0) != -1) {
						nrCalled++;
						allelesCalledPerSample[i]++;
					}
					if (genotypes.getQuick(i, 1) != -1) {
						nrCalled++;
						allelesCalledPerSample[i]++;
					}
				}

				String variantOut = variant.getChr()
						+ "\t" + variant.getPos()
						+ "\t" + variant.getId()
						+ "\t" + filterout
						+ "\t" + reason
						+ "\t" + variant.isBiallelic()
						+ "\t" + variant.isMultiallelic()
						+ "\t" + variant.getAlleles()[0]
						+ "\t" + alt
						+ "\t" + variant.getMinorAllele()
						+ "\t" + maf
						+ "\t" + nrCalled
						+ "\t" + variant.getCallrate()
						+ "\t" + Strings.concat(variantAllelesObserved, Strings.semicolon)
						+ "\t" + Strings.concat(alleleFrequencies, Strings.semicolon)
						+ "\t" + variant.getQual()
						+ "\t" + variant.getFilter()
						+ "\t" + variant.getInfo();


				variantQC.writeln(variantOut);
				if (!filterout && filterGT) {
					// rewrite alt allele if there are only two alleles observed, while VCF defines > 2 alt alleles
					if (variant.isBiallelic() && variantAllelesObserved.length > 2) {
						ArrayList<String> alleles = new ArrayList<String>();

						HashMap<String, Integer> alleleMap = new HashMap<String, Integer>();
						int ctr = 0;
						for (int i = 0; i < variantAllelesObserved.length; i++) {
							if (variantAllelesObserved[i] > 0) {
								alleles.add(variant.getAlleles()[i]);
								alleleMap.put(variant.getAlleles()[i], ctr);
								ctr++;
							}
						}

						// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
//						String[] lnElems = ln.split("\t");
//						lnElems[3] = alleles.get(0);
//						lnElems[4] = alleles.get(1);
//						String[] format = lnElems[8].split(":");
//						int gtCol = -1; // genotype
//
//						for (int c = 0; c < format.length; c++) {
//							if (format[c].equals("GT")) {
//								gtCol = c;
//							}
//						}

//						for (int i = nrHeaderElems; i < lnElems.length; i++) {
//							String[] sampleElems = lnElems[i].split(":");
//							if (gtCol != -1) {
//								String gt = sampleElems[gtCol];
//								String[] gtElems = gt.split("/");
//								if (gtElems.length == 1) { // phased genotypes
//									gtElems = gt.split("\\|");
//								}
//
//								if (gtElems[0].equals(".") || gtElems[1].equals(".")) {
//									lnElems[i] = "./.";
//								} else {
//									byte gt1 = 0;
//									byte gt2 = 0;
//									try {
//										gt1 = Byte.parseByte(gtElems[0]);
//										gt2 = Byte.parseByte(gtElems[1]);
//
//										lnElems[i] = gtToGt[gt1] + "/" + gtToGt[gt2];
//
//									} catch (NumberFormatException e) {
//										lnElems[i] = "./.";
//									}
//								}
//							}
//						}

						String[] newAlleles = alleles.toArray(new String[0]);
						variant.recodeAlleles(alleleMap, newAlleles);
						// filteredVCF.writeln(Strings.concat(lnElems, Strings.tab));
						filteredVCF.writeln(variant.toVCFString());

					} else {
						filteredVCF.writeln(variant.toVCFString());
					}
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

	public void filterMendelianErrors(String vcfIn, String famIn, String vcfOut, String outdir) throws IOException {
		HashMap<String, Pair<String, String>> trios = new HashMap<String, Pair<String, String>>();
		TextFile fam = new TextFile(famIn, TextFile.R);
		String[] elems = fam.readLineElems(TextFile.tab);
		System.out.println("Filtering mendel errors:");
		System.out.println("in: " + vcfIn);
		System.out.println("out: " + vcfOut);
		System.out.println("fam: " + famIn);

		while (elems != null) {

			String sample = elems[1];
			String p1 = elems[2];
			String p2 = elems[3];
			if (!p1.equals("0") && !p2.equals("0")) {
				if (sample.equals("40203302")) {
					System.out.println(sample + "\t" + p1 + "\t" + p2);
				}
				Pair p = new Pair<String, String>(p1, p2);
				trios.put(sample, p);
			}
			elems = fam.readLineElems(TextFile.tab);
		}
		fam.close();

		System.out.println(trios.size() + " family relationships/trios");

		String[] samples = new String[]{};
		boolean[] hasParents = new boolean[]{};
		int[][] sampleToParentId = new int[][]{};

		TextFile vcf = new TextFile(vcfIn, TextFile.R);

		TextFile out = new TextFile(vcfOut, TextFile.W);
		TextFile log = new TextFile(outdir + "mendelianerrorlog.txt", TextFile.W);
		TextFile logmendel = new TextFile(outdir + "mendelianerrors.txt", TextFile.W);


		String ln = vcf.readLine();
		while (ln != null) {
			if (ln.startsWith("#")) {

				if (ln.startsWith("#CHROM")) {
					HashMap<String, Integer> sampleIndex = new HashMap<String, Integer>();
					int sampleCtr = 0;
					String[] vcfElemes = ln.split("\t");
					samples = new String[vcfElemes.length - 9];
					hasParents = new boolean[vcfElemes.length - 9];
					sampleToParentId = new int[vcfElemes.length - 9][2];
					int nrSamplesWithParents = 0;
					System.out.println(vcfElemes.length + " columns on the VCF");
					System.out.println(hasParents.length + " spaces for parents.");
					for (int i = 9; i < vcfElemes.length; i++) {
						String sample = vcfElemes[i];
						sampleIndex.put(sample, sampleCtr);
						samples[sampleCtr] = sample;
						Pair<String, String> parents = trios.get(sample);
						if (parents != null) {
							hasParents[sampleCtr] = true;
							nrSamplesWithParents++;
						}
						sampleCtr++;
					}
					System.out.println(nrSamplesWithParents + " samples with parents in dataset");

					for (int i = 0; i < samples.length; i++) {
						Pair<String, String> parents = trios.get(samples[i]);
						if (parents != null) {
							hasParents[i] = true;
							String p1 = parents.getLeft();
							String p2 = parents.getRight();

							Integer i1 = sampleIndex.get(p1);
							Integer i2 = sampleIndex.get(p2);
							if (i1 == null) {
								i1 = -1;
							}
							if (i2 == null) {
								i2 = -1;
							}
							sampleToParentId[i][0] = i1;
							sampleToParentId[i][1] = i2;
						} else {
							sampleToParentId[i][0] = -1;
							sampleToParentId[i][1] = -1;
						}
					}
				}
				out.writeln(ln);
			} else {
				VCFVariant var = new VCFVariant(ln);
				if (!var.isBiallelic()) {
					// multi allelic variants complicate things..... let's just output these
					out.writeln(var.toVCFString());
				} else {
					DoubleMatrix2D alleles = var.getGenotypeAllelesAsMatrix2D();
					String[] allelesStr = var.getAlleles();
					int nrErrors = 0;
					int nrTested = 0;
					for (int s = 0; s < samples.length; s++) {
						if (hasParents[s]) {

							byte s1allele = (byte) alleles.getQuick(s, 0); //alleles[0][s];
							byte s2allele = (byte) alleles.getQuick(s, 1); //alleles[1][s];
							if (s1allele != -1) {
								int p1 = sampleToParentId[s][0];
								int p2 = sampleToParentId[s][1];


								HashSet<Integer> allowedAlleles = new HashSet<Integer>();
								boolean p1hasgenotype = false;
								if (p1 != -1) {
									byte allele1p1 = (byte) alleles.getQuick(p1, 0); //alleles[0][p1];
									if (allele1p1 == -1) {
										p1hasgenotype = false;
									} else {
										p1hasgenotype = true;
									}
								}

								boolean p2hasgenotype = false;
								if (p2 != -1) {
									byte allele1p2 = (byte) alleles.getQuick(p2, 0); // alleles[0][p2];
									if (allele1p2 == -1) {
										p2hasgenotype = false;
									} else {
										p2hasgenotype = true;
									}
								}


								if (p1hasgenotype || p2hasgenotype) {
									int ctr = 0;
									nrTested++;
									if (p1hasgenotype && p2hasgenotype) {
										for (int a = 0; a < 2; a++) {
											for (int b = 0; b < 2; b++) {

												allowedAlleles.add((byte) alleles.getQuick(p1, a) + (byte) alleles.getQuick(p2, b)); // 0+0 0+1 1+0 1+1
												ctr++;
											}

										}
									} else if (p1hasgenotype) { // other parent missing?

										int sum = (byte) alleles.getQuick(p1, 0) + (byte) alleles.getQuick(p1, 1);
										if (sum == 0 || sum == 2) {
											allowedAlleles.add(sum); // allow homozygous parent allele
											allowedAlleles.add(1); // allow hets if homozygous

										} else {
											// all combinations possible
											allowedAlleles.add(0); // allow homozygous parent allele
											allowedAlleles.add(1); // allow hets if homozygous
											allowedAlleles.add(2); // allow homozygous other allele
										}

									} else if (p2hasgenotype) {
										int sum = (byte) alleles.getQuick(p2, 0) + (byte) alleles.getQuick(p2, 1);

										if (sum == 0 || sum == 2) {
											allowedAlleles.add(sum); // allow homozygous parent allele
											allowedAlleles.add(1); // allow hets if homozygous
										} else {
											// all combinations possible
											allowedAlleles.add(0); // allow homozygous parent allele
											allowedAlleles.add(1); // allow hets if homozygous
											allowedAlleles.add(2); // allow homozygous other allele
										}
									}

									int sampleSum = s1allele + s2allele;


									Integer[] allowedCodes = allowedAlleles.toArray(new Integer[allowedAlleles.size()]);
									String allowedStr = Strings.concat(allowedCodes, Strings.comma);

									String p1Str = "Unknown";
									if (p1 != -1) {
										p1Str = samples[p1];
									}
									String p2Str = "Unknown";
									if (p2 != -1) {
										p2Str = samples[p2];
									}
									String p1GenotypeStr = "./.";
									int sum1 = -1;
									if (p1hasgenotype) {
										sum1 = (byte) alleles.getQuick(p1, 0) + (byte) alleles.getQuick(p1, 1);
										p1GenotypeStr = allelesStr[(byte) alleles.getQuick(p1, 0)] + "/" + allelesStr[(byte) alleles.getQuick(p1, 1)];
									}
									String p2GenotypeStr = "./.";
									if (p2hasgenotype) {
										p2GenotypeStr = allelesStr[(byte) alleles.getQuick(p2, 0)] + "/" + allelesStr[(byte) alleles.getQuick(p2, 1)];
									}
									String sampleStr = allelesStr[s1allele] + "/" + allelesStr[s2allele];


									if (!allowedAlleles.contains(sampleSum)) {
// errorrrrr
										var.getGenotypeAlleles()[0][s] = -1;
										var.getGenotypeAlleles()[1][s] = -1;
										nrErrors++;
										logmendel.writeln(var.getChr() + "\t" + var.getPos() + "\t" + var.getId() + "\t" + Strings.concat(allelesStr, Strings.comma)
												+ "\t" + samples[s] + "\t" + sampleStr + "\t" + s1allele + "\t" + s2allele
												+ "\t" + p1Str + "\t" + p1GenotypeStr + "\t"
												+ p2Str + "\t" + p2GenotypeStr + "\t" + false);
									} else {
//										logmendel.writeln(var.getChr() + "\t" + var.getPos() + "\t" + var.getId() + "\t" + Strings.concat(allelesStr, Strings.comma)
//												+ "\t" + samples[s] + "\t" + sampleStr + "\t" + s1allele + "\t" + s2allele
//												+ "\t" + p1Str + "\t" + p1GenotypeStr + "\t"
//												+ p2Str + "\t" + p2GenotypeStr + "\t" + true);
									}


								}
							}
						}
					}

					double perc = ((double) nrErrors / nrTested);


					boolean includedAfterFilter = false;

					Chromosome chr = Chromosome.parseChr(var.getChr());
					if (chr.equals(Chromosome.X)) {
						out.writeln(var.toVCFString());
						includedAfterFilter = true;
					} else {
						var.recalculateMAFAndCallRate();
						if (perc < 0.1 && var.getMAF() > 0.005 && var.getCallrate() > 0.90) {
							out.writeln(var.toVCFString());
							includedAfterFilter = true;
						}
					}


					log.writeln(var.getChr() + "\t" + var.getPos() + "\t" + var.getId() + "\t" + nrErrors + "\t" + nrTested + "\t" + perc + "\t" + includedAfterFilter);

				}


			}
			ln = vcf.readLine();

		}

		logmendel.close();
		log.close();
		out.close();
		vcf.close();
	}

	public void filterVCFVariants(String vcfIn, String vcfVariantsToFilter, String vcfOut) throws IOException {


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


		TextFile in2 = new TextFile(vcfIn, TextFile.R);
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
			// f.setName(elems[1]);
			if (variantMap.containsKey(f)) {
				System.out.println("Duplicate SNP!: " + elems[1] + " replacing chromosome with chr NA");
				f.setChromosome(Chromosome.NA);
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
		convertPEDToVCF(ped, vcf, true);
	}

	public void convertPEDToVCF(String ped, String vcf, boolean applyfilter) throws IOException {

		System.out.println("converting " + ped + " to " + vcf);
		PedAndMapFunctions pedAndMapFunctions = new PedAndMapFunctions();


		ArrayList<String> rsNames = new ArrayList<String>();
		ArrayList<String> pos = new ArrayList<String>();
		HashMap<Feature, Integer> variantMap = new HashMap<Feature, Integer>();

		HashSet<Feature> nonUniqueFeatures = new HashSet<Feature>();
		HashSet<Feature> allFeatures = new HashSet<Feature>();
		TextFile tf = new TextFile(ped + ".map", TextFile.R);

		// get a list of duplicated features?
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			Feature f = new Feature();
			f.setChromosome(Chromosome.parseChr(elems[0]));
			f.setStart(Integer.parseInt(elems[3]));
			f.setStop(Integer.parseInt(elems[3]));

			if (allFeatures.contains(f)) {
				nonUniqueFeatures.add(f);
			}
			allFeatures.add(f);

			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();
		tf.open();
		elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		// index the unique features.
		while (elems != null) {

			Feature f = new Feature();
			f.setChromosome(Chromosome.parseChr(elems[0]));
			f.setStart(Integer.parseInt(elems[3]));
			f.setStop(Integer.parseInt(elems[3]));

			if (!nonUniqueFeatures.contains(f)) {
				rsNames.add(elems[1]);
				pos.add(elems[0] + "_" + elems[3]);
				variantMap.put(f, ctr);
				ctr++;
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		System.out.println(variantMap.size() + " variant positions \t" + rsNames.size() + " rs names");


		// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
		TextFile vcfout = new TextFile(vcf, TextFile.W);
		vcfout.writeln("##fileformat=VCFv4.1");
		vcfout.writeln("##source=HJ");
		vcfout.writeln("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
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

				if (applyfilter) {


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
				} else {
					vcfout.writeln(outln);
				}
			}
		}
		vcfout.close();
	}


//	public void compareAndCorrectVCFVariants(String vcf1,
//											 String vcf1out,
//											 String vcf2,
//											 String vcf2out,
//											 String log,
//											 boolean keepNonIntersectingVariants,
//											 boolean keepMonomorphic) throws IOException {
//
//		ArrayList<Feature> vcf1variants = getVariantsFromVCF(vcf1);
//		ArrayList<Feature> vcf2variants = getVariantsFromVCF(vcf2);
//		GenotypeTools tools = new GenotypeTools();
//		HashSet<Feature> intersectVariants = tools.intersectVariants(vcf1variants, vcf2variants);
//
//
//		System.out.println("checking " + vcf1 + " vs " + vcf2);
//		System.out.println(intersectVariants.size() + " variants intersect");
//
//		ArrayList<Feature> variantsToExclude = new ArrayList<Feature>();
//
//		TextFile logfile = new TextFile(log, TextFile.W);
//		logfile.writeln(intersectVariants.size() + " variants shared");
//
//		logfile.writeln("chr\tpos\trs\talleles2\tminor2\tmaf2\tallleles-ref\tminor-ref\tmaf-ref\treason\tinfo2");
//		if (intersectVariants.size() > 0) {
//			// get genotypes for intersected variants..
//			FastListMultimap<Feature, VCFVariant> variantsvcf1 = getVCFVariants(vcf1, intersectVariants, false);
//
//			TextFile outvcf = new TextFile(vcf2out, TextFile.W);
//
//			TextFile vcfIn = new TextFile(vcf2, TextFile.R);
//
//			String ln = vcfIn.readLine();
//			while (ln != null) {
//
//				if (ln.startsWith("#")) {
//					outvcf.writeln(ln);
//				} else {
//					VCFVariant variant2 = new VCFVariant(ln);
//					Feature f = new Feature();
//					f.setChromosome(Chromosome.parseChr(variant2.getChr()));
//					f.setStart(variant2.getPos());
//					f.setStop(variant2.getPos());
//
//
//					if (intersectVariants.contains(f)) {
//						String[] alleles2 = variant2.getAlleles();
//						String minorAllele2 = variant2.getMinorAllele();
//
//						ListIterable<VCFVariant> refVariants = variantsvcf1.get(f); // get all variants at this position..
//						// compare
//
//						boolean variantHasMatch = false;
//						int nrMatches = 0;
//						for (VCFVariant refVariant : refVariants) {
//							// check the alleles
//							String[] refAlleles = refVariant.getAlleles();
//
//							String refMinorAllele = refVariant.getMinorAllele();
//
//							String logstr = variant2.getChr()
//									+ "\t" + variant2.getPos()
//									+ "\t" + variant2.getId()
//
//									+ "\t" + Strings.concat(alleles2, Strings.comma)
//									+ "\t" + minorAllele2
//									+ "\t" + variant2.getMAF()
//									+ "\t" + Strings.concat(refAlleles, Strings.comma)
//									+ "\t" + refMinorAllele
//									+ "\t" + refVariant.getMAF()
//									+ "\t" + variant2.getInfo();
//
//
//							int nridenticalalleles = 0;
//
//							for (int i = 0; i < refAlleles.length; i++) {
//								String allele1 = refAlleles[i];
//								for (int j = 0; j < alleles2.length; j++) {
//									if (alleles2[j].equals(allele1)) {
//										nridenticalalleles++;
//									}
//								}
//							}
//
//							if (refVariant.isBiallelic() && variant2.isBiallelic() ||
//									(keepMonomorphic && (refVariant.isMultiallelic() || variant2.isMonomorphic()))) {
//								String allele11 = refAlleles[0];
//								String allele12 = refAlleles[1];
//								String allele21 = alleles2[0];
//								String allele22 = alleles2[1];
//
//								// this is true if variants have identical alleles
//								// and is always true for A/T and C/G variants
//								if (nridenticalalleles == 2) {
//
//									// check A/T G/C variants
//									if (allele11.equals(BaseAnnot.getComplement(allele12)) && !refMinorAllele.equals(minorAllele2)) {
//										// both variants ar AT or GC snps
//										logstr += "\tAT or CG with DiffMinor";
//										// System.out.println("AT/CG: minor alleles different: " + allele11 + "/" + allele12 + "-" + refMinorAllele + "\t" + allele21 + "/" + allele22 + "-" + minorAllele2);
//									} else {
//										// alleles are identical. check the direction of the alleles
//										if (refMinorAllele.equals(minorAllele2)) {
//											if (allele11.equals(allele21)) {
//												// reference alleles are identical
//
//												variantHasMatch = true;
//												nrMatches++;
//												logstr += "\tOK";
//											} else {
//												// flip the genotypes
//												ln = flipReferenceAllele(ln);
//												variantHasMatch = true;
//												nrMatches++;
//												logstr += "\tFLIPOK";
//											}
//										} else {
//											// System.out.println("minor alleles different: " + allele11 + "/" + allele12 + "-" + refMinorAllele + "\t" + allele21 + "/" + allele22 + "-" + minorAllele2);
//											logstr += "\tDiffMinor";
//										}
//									}
//								} else {
//									// try complement
//									String[] complementAlleles2 = convertToComplement(alleles2);
//									String complementMinorAllele2 = getComplement(minorAllele2);
//									allele21 = complementAlleles2[0];
//									allele22 = complementAlleles2[1];
//									nridenticalalleles = 0;
//
//									for (int i = 0; i < refAlleles.length; i++) {
//										String allele1 = refAlleles[i];
//										for (int j = 0; j < alleles2.length; j++) {
//											if (complementAlleles2[j].equals(allele1)) {
//												nridenticalalleles++;
//											}
//										}
//									}
//
//									if (nridenticalalleles == 2) {
//										if (refMinorAllele.equals(complementMinorAllele2)) {
//											if (allele11.equals(allele21)) {
//												// reference alleles are identical
//												// rewrite allele codes to reflect complement
//												// #CHROM  POS     ID REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1   SAMPLE2
//												String[] elems = ln.split("\t");
//												elems[3] = allele21;
//												elems[4] = allele22;
//												ln = Strings.concat(elems, Strings.tab);
//												variantHasMatch = true;
//												nrMatches++;
//												logstr += "\tOK-Complement";
//											} else {
//												// flip the genotypes
//												ln = flipReferenceAllele(ln);
//												String[] elems = ln.split("\t");
//												elems[4] = allele21;
//												elems[3] = allele22;
//												ln = Strings.concat(elems, Strings.tab);
//												variantHasMatch = true;
//												nrMatches++;
//												logstr += "\tFLIPOK-Complement";
//											}
//										} else {
//											// System.out.println("minor alleles different: " + allele11 + "/" + allele12 + "-" + refMinorAllele + "\t" + allele21 + "/" + allele22 + "-" + minorAllele2);
//											logstr += "\tDiffMinor-Complement";
//										}
//									} else {
//										//System.out.println("incompatible alleles: " + allele11 + "/" + allele12 + "-" + refMinorAllele + "\t" + allele21 + "/" + allele22 + "-" + minorAllele2);
//										logstr += "\tIncompatibleAlleles";
//									}
//								}
//
//							} else if (refVariant.isMonomorphic() || variant2.isMonomorphic()) {
//								// exclude variant from both datasets
//
////							System.out.println("monomorphic: " + Strings.concat(refAlleles, Strings.comma) + "-" + refMinorAllele
////									+ Strings.concat(alleles2, Strings.comma) + "-" + minorAllele2);
//								logstr += "\tMonomorph";
//							} else if (refVariant.isMultiallelic() && variant2.isMultiallelic()) {
//								// check whether the alleles are the same and their allele frequencies compare as well.
//
//								if (nridenticalalleles == refVariant.getAlleles().length && nridenticalalleles == variant2.getAlleles().length) {
//									logstr += "\tOK-MultiAllelic";
//								} else {
////								System.out.println("multi-allelic: " + Strings.concat(refAlleles, Strings.comma) + "-" + refMinorAllele
////										+ Strings.concat(alleles2, Strings.comma) + "-" + minorAllele2);
//									logstr += "\tIncompatibleAlleles-MultiAllelic";
//								}
//							} else {
//
////							System.out.println("inconclusive: " + Strings.concat(refAlleles, Strings.comma) + "-" + refMinorAllele
////									+ Strings.concat(alleles2, Strings.comma) + "-" + minorAllele2);
//								logstr += "\tInconclusive";
//							}
//
//
//							logfile.writeln(logstr + "\t" + variantHasMatch);
//						}
//
//						if (!variantHasMatch) {
//							variantsToExclude.add(f);
//						} else {
//							outvcf.writeln(ln);
//						}
//
//						if (refVariants.size() > 1) {
//							if (nrMatches == 0) {
//								System.err.println(variant2.getId() + " has " + refVariants.size() + " counterparts in ref. " + nrMatches + " of which match alleles");
//							}
////								System.out.println();
//						}
//
//
//					} else if (keepNonIntersectingVariants) {
//						outvcf.writeln(ln);
//
//					}
//				}
//
//				ln = vcfIn.readLine();
//			}
//
//			vcfIn.close();
//			outvcf.close();
//
//		}
//		logfile.close();
//
//		System.out.println("About to exclude: " + variantsToExclude.size() + " variants from reference dataset...");
//
//		TextFile vcfoutref = new TextFile(vcf1out, TextFile.W);
//		TextFile vcfinref = new TextFile(vcf1, TextFile.R);
//		String ln = vcfinref.readLine();
//		while (ln != null) {
//
//			if (ln.startsWith("#")) {
//				vcfoutref.writeln(ln);
//			} else {
//				VCFVariant variant = new VCFVariant(ln);
//				Feature f = new Feature();
//				f.setChromosome(Chromosome.parseChr(variant.getChr()));
//				f.setStart(variant.getPos());
//				f.setStop(variant.getPos());
//				if (!variantsToExclude.contains(f)) {
//					if (keepNonIntersectingVariants) {
//						vcfoutref.writeln(ln);
//					} else if (intersectVariants.contains(f)) {
//						vcfoutref.writeln(ln);
//					}
//				}
//			}
//			ln = vcfinref.readLine();
//		}
//		vcfoutref.close();
//		vcfinref.close();
//
//	}


	public String flipReferenceAllele(String ln) {

		// #CHROM  POS     ID REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1   SAMPLE2

		String[] elems = ln.split("\t");
		String ref = elems[3];
		String alt = elems[4];

		elems[3] = alt;
		elems[4] = ref;

		String[] format = elems[8].split(":");
		int gtCol = -1; // genotype
		int adCol = -1; // Allelic depths for the ref and alt alleles in the order listed
		int dpCol = -1; // Approximate readAsTrack depth (reads with MQ=255 or with bad mates are filtered)
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


		for (int i = 9; i < elems.length; i++) {
			String sample = elems[i];
			String[] sampleElems = sample.split(":");
			String gt = sampleElems[gtCol];
			String[] gtElems = gt.split("/");
			String sep = "/";
			if (gtElems.length == 1) {
				gtElems = gt.split("\\|");
				sep = "|";
			}

			if (gtElems[0].equals("0")) {
				gtElems[0] = "1";
			} else {
				gtElems[0] = "0";
			}
			if (gtElems[1].equals("0")) {
				gtElems[1] = "1";
			} else {
				gtElems[1] = "0";
			}

			sampleElems[gtCol] = gtElems[0] + sep + gtElems[1];
			elems[i] = Strings.concat(sampleElems, Strings.colon);
		}

		return Strings.concat(elems, Strings.tab);

	}

	public FastListMultimap<Feature, VCFVariant> getVCFVariants(String vcf, HashSet<Feature> map, boolean loadGenotypes) throws IOException {
		TextFile vcftf = new TextFile(vcf, TextFile.R);

		String ln = vcftf.readLine();

		FastListMultimap<Feature, VCFVariant> output = new FastListMultimap<Feature, VCFVariant>();

		while (ln != null) {
			if (!ln.startsWith("##") && !ln.startsWith("#CHROM")) {
				VCFVariant variant = null;
				if (loadGenotypes) {
					variant = new VCFVariant(ln, VCFVariant.PARSE.GENOTYPES);
				} else {
					variant = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
				}


				Feature f = new Feature();
				f.setChromosome(Chromosome.parseChr(variant.getChr()));
				f.setStart(variant.getPos());
				f.setStop(variant.getPos());
				if (map == null || map.contains(f)) {
					output.put(f, variant);
				}

			}
			ln = vcftf.readLine();
		}
		vcftf.close();
		return output;
	}

	public void determineVCFSummaryStatistics(String vcf, String out) throws IOException {
		VCFFunctions vcfFunctions = new VCFFunctions();

		GenotypeTools t = new GenotypeTools();

		ArrayList<String> vcfSamples = vcfFunctions.getVCFSamples(vcf);
		System.out.println(vcfSamples.size() + " samples loaded from VCF");

		HashMap<String, Integer> sampleMap = StringSets.index(vcfSamples.toArray(new String[0]));
// get variants from VCF
		ArrayList<Feature> variantsOnVCF = vcfFunctions.getVariantsFromVCF(vcf, false, true);
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


		t.writeGenotypes(genotypePair, vcfSamples, autosomal, out + "allVCFGenotypes.txt");
		byte[][] genotypes = genotypePair.getLeft();
		String[] alleles = genotypePair.getRight();
		Triple<Pair<double[], String[]>, double[], double[]> params = t.getGenotypeParams(genotypePair);
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

//		Coverage c = new Coverage();
//		Chart barchart = c.initBarChart("Call-rate per sample", "Sample", "Call-rate", 1200, 800, true, 50);
//		barchart.addSeries("Series1", xvals, sampleCR);
//		barchart.getStyleManager().setLegendVisible(false);
//		BitmapEncoder.saveBitmapWithDPI(barchart, out + "SampleCallRate", BitmapEncoder.BitmapFormat.PNG, 300);


	}

	public void replaceHeader(String vcfin, String ref, String vcfout) throws IOException {

		TextFile tf = new TextFile(ref, TextFile.R);
		String header = "";
		String ln = tf.readLine();
		while (ln != null) {
			if (ln.startsWith("##")) {
				if (header.length() == 0) {
					header += ln;
				} else {
					header += "\n" + ln;

				}


			}
			ln = tf.readLine();
		}
		tf.close();

		TextFile tfin = new TextFile(vcfin, TextFile.R);
		TextFile tfout = new TextFile(vcfout, TextFile.W);
		ln = tfin.readLine();
		while (ln != null) {
			if (ln.startsWith("##")) {

			} else if (ln.startsWith("#")) {

				tfout.writeln(header);
				tfout.writeln(ln);
			} else {
				tfout.writeln(ln);
			}
			ln = tfin.readLine();
		}
		tfin.close();
		tfout.close();
	}


	public String mergeVariants(VCFVariant var1, VCFVariant var2, String separator) {
		StringBuilder output = new StringBuilder();
		output.append(var1.getChr());
		output.append("\t");
		output.append(var1.getPos());
		output.append("\t");
		output.append(var1.getId());
		output.append("\t");
		output.append(var1.getAlleles()[0]);
		output.append("\t");
		output.append(Strings.concat(var1.getAlleles(), Strings.comma, 1, var1.getAlleles().length));
		output.append("\t.\t.\t.\tGT");

		DoubleMatrix2D genotypeAlleles1 = var1.getGenotypeAllelesAsMatrix2D();
		DoubleMatrix2D genotypeAlleles2 = var2.getGenotypeAllelesAsMatrix2D();
		if (separator == null) {
			separator = "/";
		}

		for (int i = 0; i < genotypeAlleles1.rows(); i++) {

			String allele1 = "" + (byte) genotypeAlleles1.get(i, 0);
			String allele2 = "" + (byte) genotypeAlleles1.get(i, 1);

			if (genotypeAlleles1.get(i, 0) == -1) {
				allele1 = ".";
			}
			if (genotypeAlleles1.get(i, 1) == -1) {
				allele2 = ".";
			}

			output.append("\t").append(allele1).append(separator).append(allele2);
		}

		for (int i = 0; i < genotypeAlleles2.rows(); i++) {
			String allele1 = "" + genotypeAlleles2.getQuick(i, 0);
			String allele2 = "" + genotypeAlleles2.getQuick(i, 1);

			if (genotypeAlleles2.getQuick(i, 0) == -1) {
				allele1 = ".";
			}
			if (genotypeAlleles2.getQuick(i, 1) == -1) {
				allele2 = ".";
			}

			output.append("\t").append(allele1).append(separator).append(allele2);
		}

		return output.toString();
	}

	public void splitMultipleAllelicVariants(String vcfIn, String vcfOut) throws IOException {
		System.out.println("Splitting alleles in: " + vcfIn);
		System.out.println("Writing to: " + vcfOut);
		TextFile tfIn = new TextFile(vcfIn, TextFile.R);
		TextFile tfOut = new TextFile(vcfOut, TextFile.W);

		String ln = tfIn.readLine();
		while (ln != null) {
			if (ln.startsWith("#")) {
				tfOut.writeln(ln);
			} else {
				VCFVariant var = new VCFVariant(ln);
				String[] alleles = var.getAlleles();
				if (alleles.length <= 2) {
					tfOut.writeln(ln);
				} else {

					// nrVariants = alleles.length-1
					DoubleMatrix2D genotypeAlleles = var.getGenotypeAllelesAsMatrix2D();


					int nrLinesWritten = 0;
					for (byte i = 1; i < alleles.length; i++) {

						String ref = alleles[0];
						String alt = alleles[i];
						String lnout = var.getChr()
								+ "\t" + var.getPos()
								+ "\t" + var.getId()
								+ "\t" + ref
								+ "\t" + alt
								+ "\t.\t.\t.\tGT";


						for (int ind = 0; ind < genotypeAlleles.rows(); ind++) {
							byte gt1 = (byte) genotypeAlleles.getQuick(ind, 0);
							byte gt2 = (byte) genotypeAlleles.getQuick(ind, 1);


							if (gt1 == i) {
								gt1 = 1;
							} else {
								gt1 = 0;
							}

							if (gt2 == i) {
								gt2 = 1;
							} else {
								gt2 = 0;
							}

							lnout += "\t" + gt1 + "/" + gt2;
						}
						tfOut.writeln(lnout);
						nrLinesWritten++;
					}

					System.out.println(alleles.length + " alleles detected for variant: " + var.getId() + "/" + var.getChr() + ":" + var.getPos() + " \tnrlines: " + nrLinesWritten);
				}
			}
			ln = tfIn.readLine();
		}
		tfIn.close();
		tfOut.close();
	}

	public void splitPerChromosome(String vcfIn, String outdirprefix) throws IOException {
		System.out.println("Splitting: " + vcfIn + " to " + outdirprefix);

		HashMap<Chromosome, TextFile> textFiles = new HashMap<Chromosome, TextFile>();

		VCFGenotypeData d = new VCFGenotypeData(vcfIn);
		int ctr = 0;
		while (d.hasNext()) {
			VCFVariant v = d.next();
			Feature f = v.asFeature();
			ctr++;
			if (f.getChromosome().equals(Chromosome.NA)) {
				System.out.println("unknown chromosome found..?");
			} else {
				if (!textFiles.containsKey(f.getChromosome())) {
					textFiles.put(f.getChromosome(), new TextFile(outdirprefix + "-" + f.getChromosome().getName() + ".vcf.gz", TextFile.W));
				}
			}

		}
		d.close();

		Set<Chromosome> chromosomes = textFiles.keySet();

		for (Chromosome chr : chromosomes) {
			System.out.println(chr.getName() + " found");
		}

		System.out.println(ctr + " over " + textFiles.size() + " chromosomes");
		TextFile tf = new TextFile(vcfIn, TextFile.R);

		String ln = tf.readLine();
		while (ln != null) {

			if (ln.startsWith("#")) {
				for (Chromosome chr : chromosomes) {
					textFiles.get(chr).writeln(ln);
				}
			} else {
				VCFVariant var = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
				Chromosome chr = Chromosome.parseChr(var.getChr());
				TextFile tf2 = textFiles.get(chr);
				if (tf2 != null) {
					tf2.writeln(ln);
				} else {
					System.err.println("ERRORRRRR: chr: " + chr.getName() + " not found? for variant: " + var.getChr() + "-" + var.getPos() + "/" + var.getId());
				}

			}
			ln = tf.readLine();
		}

		tf.close();


		for (Chromosome chr : chromosomes) {
			textFiles.get(chr).close();
		}

	}

	public void filterSamples(String vcfIn, String vcfOut, String sampleList) throws IOException {

		TextFile sampleFile = new TextFile(sampleList, TextFile.R);
		HashSet<String> allowedSamples = (HashSet<String>) sampleFile.readAsSet(0, TextFile.tab);
		sampleFile.close();

		TextFile out = new TextFile(vcfOut, TextFile.W);
		TextFile tfIn = new TextFile(vcfIn, TextFile.R);
		String ln = tfIn.readLine();
		boolean[] includeCol = new boolean[0];
		while (ln != null) {
			if (ln.startsWith("##")) {
				out.writeln(ln);
			} else if (ln.startsWith("#CHROM")) {
				// populate index
				String[] elems = ln.split("\t");
				String outHeader = Strings.concat(elems, Strings.tab, 0, 9);
				includeCol = new boolean[elems.length];
				for (int i = 9; i < elems.length; i++) {
					String sample = elems[i];
					if (allowedSamples.contains(sample)) {
						outHeader += "\t" + sample;
						includeCol[i] = true;
					}
				}
				out.writeln(outHeader);
			} else {
				String[] elems = ln.split("\t");
				StringBuilder builder = new StringBuilder();
				builder.append(Strings.concat(elems, Strings.tab, 0, 9));
				for (int i = 9; i < elems.length; i++) {
					if (includeCol[i]) {
						builder.append("\t").append(elems[i]);
					}
				}
				out.writeln(builder.toString());
			}
			ln = tfIn.readLine();
		}

		out.close();
		tfIn.close();


	}


	public void isolateGT(String vcfIn, String vcfOut) throws IOException {
		TextFile tf = new TextFile(vcfIn, TextFile.R);
		TextFile tfout = new TextFile(vcfOut, TextFile.W);
		String[] elems = tf.readLineElems(TextFile.tab);

		while (elems != null) {
			if (elems[0].startsWith("#")) {

				tfout.writeln(Strings.concat(elems, Strings.tab));
			} else {

				boolean include = true;

				// // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1   SAMPLE2
				String[] infoElems = elems[7].split(";");
				for (int i = 0; i < infoElems.length; i++) {
					String info = infoElems[i];
					if (info.startsWith("AR2")) {
						String[] infoStr = info.split("=");
						if (Double.parseDouble(infoStr[1]) < 0.8) {
							include = false;
						}
						break;
					}
				}

				if (include) {
					String[] formatElems = elems[8].split(":");
					int gtCol = -1;
					for (int f = 0; f < formatElems.length; f++) {
						if (formatElems[f].equals("GT")) {
							gtCol = f;
							break;
						}
					}

					StringBuilder output = new StringBuilder();
					output.append(Strings.concat(elems, Strings.tab, 0, 8)).append("\t").append("GT");
					for (int e = 9; e < elems.length; e++) {
						String[] sampleElems = elems[e].split(":");
						output.append("\t").append(sampleElems[gtCol]);
					}
					tfout.writeln(output.toString());
				}

			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		tfout.close();

	}


	public void filterVariantsWithWeirdEncoding(String s) throws IOException {
		System.out.println("Looking for encoding issues in path: " + s);
		TextFile tmp = new TextFile(s + "-tmp.vcf.gz", TextFile.W);
		TextFile tf = new TextFile(s, TextFile.R);

		String ln = tf.readLine();
		while (ln != null) {

			if (ln.startsWith("#")) {
				tmp.writeln(ln);
			} else {

				try {
					VCFVariant variants = new VCFVariant(ln);
					tmp.writeln(ln);
				} catch (ArrayIndexOutOfBoundsException e) {
					String[] bleh = ln.split("\t");
					if (bleh.length > 3) {
						System.err.println(bleh[0] + "" + bleh[1] + "" + bleh[2] + " has issues with encoding. in path: " + s);
					}
				}
			}
			ln = tf.readLine();
		}

		tf.close();
		tmp.close();
		Gpio.moveFile(s + "-tmp.vcf.gz", s);
	}

	public void sortVCF(boolean linux, String vcfsort, String vcf, String sortedvcf, String bashfilename) throws IOException {


		String cat = "cat";
		if (vcf.endsWith(".gz")) {
			cat = "gzcat";
			if (linux) {
				cat = "zcat";
			}
		}

		boolean gzipoutput = false;
		if (sortedvcf.endsWith(".gz")) {
			gzipoutput = true;
			sortedvcf = sortedvcf.replaceAll(".gz", "");
		}
		String bashCommand = cat + " " + vcf + " | " + vcfsort + " > " + sortedvcf;

		if (gzipoutput) {
			if (Gpio.exists(sortedvcf + ".gz")) {
				File f = new File(sortedvcf + ".gz");
				f.delete();
			}
			bashCommand += "\n";
			bashCommand += "gzip " + sortedvcf + "\n";

		}

		TextFile tf = new TextFile(bashfilename, TextFile.W);
		tf.writeln("#!/bin/bash\n" + bashCommand);
		tf.close();
		ProcessBuilder pb = new ProcessBuilder("bash", bashfilename);
		ProcessRunner.run(pb);
	}

	public void concatVCF(String vcfConcat, String vcfSort, String files, String merged, String mergedsorted, String bashfilename) throws IOException {


		boolean gzipoutput = false;
		if (mergedsorted.endsWith(".gz")) {
			gzipoutput = true;
			mergedsorted = mergedsorted.replaceAll(".gz", "");
		}
		String bashCommand = vcfConcat + " " + files + " > " + merged + "\n" +
				"cat " + merged + " | " + vcfSort + " > " + mergedsorted + "\n" +
				"rm " + merged + "\n";

		if (gzipoutput) {
			if (Gpio.exists(mergedsorted + ".gz")) {
				File f = new File(mergedsorted + ".gz");
				f.delete();
			}
			bashCommand += "gzip " + mergedsorted + "\n";
		}

		TextFile tf = new TextFile(bashfilename, TextFile.W);
		tf.writeln("#!/bin/bash\n" + bashCommand);
		tf.close();

		ProcessBuilder pb = new ProcessBuilder("bash", bashfilename);
		ProcessRunner.run(pb);
	}

	public void bgzipAndIndex(String vcf, String bashfilename) throws IOException {
		String bashCommand = "/Data/Tools/tabix/tabix-0.2.6/bgzip " + vcf + "\n" +
				"/Data/Tools/tabix/tabix-0.2.6/tabix -p vcf " + vcf + ".gz";

		if (Gpio.exists(vcf + ".gz")) {
			new File(vcf + ".gz").delete();
		}

		TextFile tf = new TextFile(bashfilename, TextFile.W);
		tf.writeln("#!/bin/bash\n" + bashCommand);
		tf.close();

		ProcessBuilder pb = new ProcessBuilder("bash", bashfilename);
		ProcessRunner.run(pb);
	}

	public void mergeVCF(boolean linux, String vcfsort, String files, String mergedout, String bashfilename) throws IOException {
		String bashCommand = "/Data/Tools/bcftools/bcftools-1.2/bcftools merge -m all -O v -o  " + mergedout + " " + files + "\n";

		TextFile tf = new TextFile(bashfilename, TextFile.W);
		tf.writeln("#!/bin/bash\n" + bashCommand);
		tf.close();

		ProcessBuilder pb = new ProcessBuilder("bash", bashfilename);
		ProcessRunner.run(pb);

		sortVCF(linux, vcfsort, mergedout, mergedout + "sorted.vcf", bashfilename + "-sort.sh");
	}

	public void replaceMissingCodes(String vcfIn, String vcfOut) throws IOException {
		TextFile tf = new TextFile(vcfIn, TextFile.R);
		TextFile out = new TextFile(vcfOut, TextFile.W);

		String ln = tf.readLine();
		while (ln != null) {

			if (ln.startsWith("#")) {
				out.writeln(ln);
			} else {
				ln = ln.replaceAll("-1/-1", "./.");
				out.writeln(ln);
			}
			ln = tf.readLine();
		}

		tf.close();
		out.close();

	}

	public void removedups(String vcfIn, String vcfOut) throws IOException {
		System.out.println("Dedup");
		System.out.println("in: " + vcfIn);
		System.out.println("out: " + vcfOut);
		TextFile tf = new TextFile(vcfIn, TextFile.R);
		String ln = tf.readLine();
		int linenumber = 0;

		ArrayList<Double> callrates = new ArrayList<>();
		HashMap<Integer, ArrayList<Integer>> posToIndex = new HashMap<>();
		int nrvars = 0;
		while (ln != null) {
			if (ln.startsWith("#")) {
			} else {
				VCFVariant var = new VCFVariant(ln, true);

				ArrayList<Integer> linenumbers = posToIndex.get(var.getPos());
				if (linenumbers == null) {
					linenumbers = new ArrayList<>();
				}
				linenumbers.add(linenumber);
				callrates.add(var.getCallrate());
				posToIndex.put(var.getPos(), linenumbers);
				nrvars++;
			}
			linenumber++;
			ln = tf.readLine();
		}

		tf.close();
		tf.open();
		TextFile out = new TextFile(vcfOut, TextFile.W);
		linenumber = 0;
		ln = tf.readLine();
		HashSet<Integer> positionsWritten = new HashSet<>();
		int nrwritten = 0;
		while (ln != null) {
			if (ln.startsWith("#")) {
				out.writeln(ln);
			} else {
				VCFVariant var = new VCFVariant(ln, true);
				ArrayList<Integer> linenumbers = posToIndex.get(var.getPos());

				if (linenumbers.size() <= 1) {
					out.writeln(ln);
					nrwritten++;
				} else {

					int position = var.getPos();
					if (!positionsWritten.contains(position)) {
						int maxcrln = linenumber;
						double maxcallrate = -1;
						for (int lnum : linenumbers) {
							double cr = callrates.get(lnum);
							if (cr > maxcallrate) {
								maxcallrate = cr;
								maxcrln = lnum;
							}
						}

						if (maxcrln == linenumber) {
							out.writeln(ln);
							positionsWritten.add(position);
							nrwritten++;
						}

					}
				}

			}

			linenumber++;
			ln = tf.readLine();
		}
		out.close();
		System.out.println(nrwritten + " / " + nrvars + " vars written");
	}
}
