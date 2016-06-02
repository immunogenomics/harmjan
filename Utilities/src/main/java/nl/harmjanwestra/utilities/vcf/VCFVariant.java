package nl.harmjanwestra.utilities.vcf;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.matrix.ByteMatrix2D;
import nl.harmjanwestra.utilities.matrix.ShortMatrix2D;
import nl.harmjanwestra.utilities.vcf.filter.VCFGenotypeFilter;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.HWE;
import umcg.genetica.text.Strings;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * Created by hwestra on 4/21/15.
 */
public class VCFVariant {

	private static final int nrHeaderElems = 9;
	private static final Pattern slash = Pattern.compile("/");
	private static final Pattern pipe = Pattern.compile("\\|");
	private static final Pattern nullGenotype = Pattern.compile("\\./\\.");
	private final HashMap<String, String> info = new HashMap<String, String>();
	private ArrayList<VCFGenotypeFilter> filters;
	private ByteMatrix2D genotypeAlleles; // format [alleles][individuals] (save some memory by making only two individual-sized arrays)
	private DoubleMatrix2D imputedDosages; // format [alleles][individuals] (save some memory by making only two individual-sized arrays)
	private int[] nrAllelesObserved;
	private double[] genotypeDosages;
//	private DoubleMatrix2D genotypeProbabilies;

	private short[] genotypeQuals;
	private short[] approximateDepth;
	private boolean monomorphic;
	private double callrate;
	private boolean multiallelic;
	private double hwep;
	private boolean biallelic = false;
	private double[] alleleFrequencies;
	private String minorAllele;
	private String[] alleles = null;
	private String chr = null;
	private int pos = -1;
	private String id = null;
	private int qual = -1;
	private String filter = null;
	private double MAF;
	private String separator = "/";


	// GT:AB:AD:DP:GQ:PL
	private int gtCol = -1; // genotype
	private int abCol = -1; // allelic balance
	private int adCol = -1; // Allelic depths for the ref and alt alleles in the order listed
	private int dpCol = -1; // Approximate readAsTrack depth (reads with MQ=255 or with bad mates are filtered)
	private int gqCol = -1; // Genotype Quality
	private int plCol = -1; // Normalized, Phred-scaled likelihoods for genotypes
	private int pidCol = -1; // ?
	private int pgtCol = -1; // ?
	private int dsCol = -1;
	private int gpCol = -1;

	public DenseDoubleMatrix2D getGenotypeProbabilies() {
		return genotypeProbabilies;
	}

	private boolean ignoregender;
	private ShortMatrix2D allelicDepth;
	private int totalCalledAlleles;
	private DenseDoubleMatrix2D genotypeProbabilies;

	public VCFVariant(String ln) {
		this(ln, PARSE.ALL);
	}

	public VCFVariant(String ln, boolean ignoregender) {
		this.ignoregender = ignoregender;

		parse(ln, PARSE.ALL);
	}

	public VCFVariant(String ln, ArrayList<VCFGenotypeFilter> filters, boolean ignoregender) {
		this.ignoregender = ignoregender;
		this.filters = filters;

		parse(ln, PARSE.ALL);
	}

	public VCFVariant(String ln, PARSE p) {

		parse(ln, p);
	}

	public VCFVariant(String chr, Integer pos, String id, String alleleStr, String info, ByteMatrix2D alleles, DoubleMatrix2D dosages) {
		this.chr = new String(chr).intern();
		this.pos = pos;
		this.id = id;
		String[] allelesElems = alleleStr.split(",");
		this.alleles = allelesElems;
		parseInfoString(info);
		this.genotypeAlleles = alleles;
		this.imputedDosages = dosages;

		recalculateMAFAndCallRate();

	}

	public short[][] getAllelicDepth() {
		return allelicDepth.toArray();
	}

	private void parse(String ln, PARSE p) {

		// count number of alleles with certain readAsTrack depth
		parseHeader(ln);

		if (p.equals(PARSE.ALL) || p.equals(PARSE.GENOTYPES)) {
			parseGenotypes(ln, p);
		}

	}

	public Double getImputationQualityScore() {
		// BEAGLE VCF qual score
		String output = info.get("AR2");
		if (output == null) {
			// PBWT / Impute2
			output = info.get("INFO");
		}
		if (output != null) {
			return Double.parseDouble(output);
		}
		return null;
	}

	public void parseGenotypes(String ln, PARSE p) {

		if (ln != null) {

			String[] tokenArr = Strings.tab.split(ln);

			if (tokenArr.length > 9) { // allow VCFs without any actual genotypes
				// parse actual genotypes.
				int nrTokens = tokenArr.length;
				int nrSamples = nrTokens - nrHeaderElems;

				genotypeAlleles = new ByteMatrix2D(2, nrSamples);
				approximateDepth = new short[nrSamples];
				genotypeQuals = new short[nrSamples];

				DenseDoubleMatrix2D genotypeProbabilies = null;

				for (int t = 9; t < nrTokens; t++) {
					String token = tokenArr[t];
					int indPos = t - nrHeaderElems;
					String sampleColumn = token;
					if (nullGenotype.equals(sampleColumn)) {
						// not called
						genotypeAlleles.setQuick(0, indPos, (byte) -1);
						genotypeAlleles.setQuick(1, indPos, (byte) -1);

					} else {
						//String[] sampleElems = Strings.colon.split(sampleColumn);

						String[] sampleTokens = Strings.colon.split(sampleColumn);

						for (int s = 0; s < sampleTokens.length; s++) {
							String sampleToken = sampleTokens[s];

							if (s == gtCol) {
								String gt = sampleToken;
								if (nullGenotype.equals(gt)) {
									genotypeAlleles.setQuick(0, indPos, (byte) -1);
									genotypeAlleles.setQuick(1, indPos, (byte) -1);

								} else {
									String[] gtElems = slash.split(gt);
									separator = "/";
									if (gtElems.length == 1) { // phased genotypes ??
										gtElems = pipe.split(gt);
										separator = "|";
									}

									byte gt1 = 0;
									byte gt2 = 0;

									if (gtElems[0].equals(".")) {
										genotypeAlleles.setQuick(0, indPos, (byte) -1);
										genotypeAlleles.setQuick(1, indPos, (byte) -1);

									} else {

										if (gtElems.length < 2) {
											System.out.println(ln);
											System.out.println(Strings.concat(gtElems, Strings.tab));
											System.exit(-1);
										}

										try {
											gt1 = Byte.parseByte(gtElems[0]);
											gt2 = Byte.parseByte(gtElems[1]);

											genotypeAlleles.setQuick(0, indPos, gt1);
											genotypeAlleles.setQuick(1, indPos, gt2);


										} catch (NumberFormatException e) {
											System.out.println("Cannot parse genotype string: " + token + " nr elems: " + gtElems.length);
											genotypeAlleles.setQuick(0, indPos, (byte) -1);
											genotypeAlleles.setQuick(1, indPos, (byte) -1);
										}
									}
								}
							} else if (s == dsCol) {
								// dosages
								if (p.equals(PARSE.ALL)) {
									try {
										if (genotypeDosages == null) {
											genotypeDosages = new double[nrSamples];
										}
										genotypeDosages[indPos] = Double.parseDouble(sampleToken);
									} catch (NumberFormatException e) {

									}
								}
							} else if (s == gpCol) {
								// genotype probs
								if (p.equals(PARSE.ALL)) {
									String[] gpElems = Strings.comma.split(sampleToken);

									try {

										if (genotypeProbabilies == null) {
											genotypeProbabilies = new DenseDoubleMatrix2D(gpElems.length, nrSamples);
										}

										for (int g = 0; g < gpElems.length; g++) {
											genotypeProbabilies.setQuick(g, indPos, Double.parseDouble(gpElems[g]));

										}

									} catch (NumberFormatException e) {

									}
								}
							} else if (s == adCol) {
								// depth of sequencing per allele
								String[] adElems = Strings.comma.split(sampleToken);
								try {
									if (allelicDepth == null) {
										allelicDepth = new ShortMatrix2D(alleles.length, nrSamples);
									}

//									try {
									for (int g = 0; g < adElems.length; g++) {
										allelicDepth.setQuick(g, indPos, Short.parseShort(adElems[g]));
									}
//									} catch (ArrayIndexOutOfBoundsException e) {
//										e.printStackTrace();
//
//
////										System.out.println("Error with: " + ln);
//
//										System.out.println();
//										System.out.println(adElems.length + "\t" + alleles.length + "\t" + nrSamples + "\t" + Strings.concat(alleles, Strings.comma));
//										System.exit(-1);
//									}

								} catch (NumberFormatException e) {

								}
							} else if (s == gqCol) {
								// genotype quals
								short gq = 0;

								try {
									if (gqCol != -1) {
										gq = Short.parseShort(sampleToken);
									}
								} catch (NumberFormatException e) {
								}
								genotypeQuals[indPos] = gq;

							} else if (s == dpCol) {
								// approximate depth of sequencing
								short depth = 0;

								try {
									if (dpCol != -1) {
										depth = Short.parseShort(sampleToken);
									}
								} catch (NumberFormatException e) {
								}
								approximateDepth[indPos] = depth;
							}
						}
					}
				}

				if (genotypeProbabilies != null) {
					int nrAlleles = alleles.length;
					imputedDosages = new DenseDoubleMatrix2D(genotypeProbabilies.columns(), nrAlleles - 1);
					for (int i = 0; i < genotypeProbabilies.columns(); i++) {
						int alctr = 0;
						for (int a1 = 0; a1 < nrAlleles; a1++) {
							for (int a2 = a1; a2 < nrAlleles; a2++) {
								double dosageval = genotypeProbabilies.getQuick(alctr, i);
								if (a1 > 0) {
									imputedDosages.setQuick(i, a1 - 1, imputedDosages.getQuick(i, a1 - 1) + dosageval);
								}
								if (a2 > 0) {
									imputedDosages.setQuick(i, a2 - 1, imputedDosages.getQuick(i, a2 - 1) + dosageval);
								}
								alctr++;
							}
						}
					}
					this.genotypeProbabilies = genotypeProbabilies;
				}

				if (filters != null) {
					for (VCFGenotypeFilter filter : filters) {
						filter.filter(this);
					}
				}
				recalculateMAFAndCallRate();
			}
		}


	}


	private void parseHeader(String ln) {
		// parse line header
		String ref = "";

		if (ln != null) {
			int strlen = ln.length();
			int substrlen = 1500; // this should capture most annotation //
			if (strlen < substrlen) {
				substrlen = strlen;
			}
			String lnheader = ln.substring(0, substrlen);
			String[] tokenArr = Strings.tab.split(lnheader);

			for (int t = 0; t < 9; t++) {
				if (t < tokenArr.length) {
					String token = tokenArr[t];

					switch (t) {
						case 0:
							this.chr = new String(token).intern();
							break;
						case 1:
							pos = Integer.parseInt(token);
							break;
						case 2:
							id = new String(token);
							break;
						case 3:
							ref = token;
							break;
						case 4:
							String alt = token;
							String[] alternateAlleles = alt.split(",");
							alleles = new String[1 + alternateAlleles.length];
							alleles[0] = new String(ref).intern();
							for (int i = 0; i < alternateAlleles.length; i++) {
								alleles[1 + i] = new String(alternateAlleles[i]).intern();
							}
//							System.out.println(id + "\t" + Strings.concat(alleles, Strings.comma));
							nrAllelesObserved = new int[alternateAlleles.length + 1];
							break;
						case 5:
							String qualStr = token;
							try {
								qual = Integer.parseInt(qualStr);
							} catch (NumberFormatException e) {

							}
							break;
						case 6:
							filter = new String(token).intern();
							break;
						case 7:
							String infoStr = token;
							parseInfoString(infoStr);
							break;
						case 8:
							String[] format = Strings.colon.split(token);

							for (int c = 0; c < format.length; c++) {
								if (format[c].equals("GT")) {
									gtCol = c;
								} else if (format[c].equals("AD")) {
									adCol = c;
								} else if (format[c].equals("AB")) {
									abCol = c;
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
								} else if (format[c].equals("DS")) {
									dsCol = c;
								} else if (format[c].equals("GP")) {
									gpCol = c;
								}
								// GT:AD:DP:GQ:PGT:PID:PL
							}

							if (gtCol == -1) {
								System.out.println("No GT COL: " + token);
								System.exit(-1);
							}
							break;
					}
				}

			}
		}

	}

	private void parseInfoString(String infoStr) {
		String[] infoElems = Strings.semicolon.split(infoStr);

		if (!infoStr.equals(".")) {
			for (int e = 0; e < infoElems.length; e++) {
				String[] infoElemElems = Strings.equalssign.split(infoElems[e]);
				String id = new String(infoElemElems[0]).intern();
				if (infoElemElems.length >= 2) {
					String val = new String(infoElemElems[1]).intern();
					info.put(id, val);
				} else {
					info.put(id, null);
				}
			}
		}
	}

	public boolean isMonomorphic() {
		return monomorphic;
	}

	public double getCallrate() {
		return callrate;
	}

	public boolean isMultiallelic() {
		return multiallelic;
	}

	public double getHwep() {
		return hwep;
	}

	public boolean isBiallelic() {
		return biallelic;
	}

	public double[] getAlleleFrequencies() {
		return alleleFrequencies;
	}

	public int getTotalAlleleCount() {
		return totalCalledAlleles;
	}

	public String[] getAlleles() {
		return alleles;
	}

	public String getMinorAllele() {
		return minorAllele;
	}

	public String getChr() {
		return chr;
	}

	public int getPos() {
		return pos;
	}

	public String getId() {
		return id;
	}

	public int getQual() {
		return qual;
	}

	public String getFilter() {
		return filter;
	}

	public short[] getGenotypeQuals() {
		return genotypeQuals;
	}

	public short[] getApproximateDepth() {
		return approximateDepth;
	}

	public double getMAF() {
		return MAF;
	}

	public int[] getNrAllelesObserved() {
		return nrAllelesObserved;
	}

	public byte[][] getGenotypeAlleles() {
		return genotypeAlleles.toArray();
	}

	public void flipReferenceAlelele() {

		String allele0 = alleles[0];
		String allele1 = alleles[1];
		alleles[0] = allele1;
		alleles[1] = allele0;

		if (genotypeAlleles != null) {

			// only works for biallelic variants!
			for (int i = 0; i < genotypeAlleles.columns(); i++) {
				for (int j = 0; j < genotypeAlleles.rows(); j++) {
					if (genotypeAlleles.getQuick(j, i) == -1) {
						genotypeAlleles.setQuick(j, i, (byte) -1);
					} else {
						genotypeAlleles.setQuick(j, i, (byte) Math.abs(genotypeAlleles.getQuick(j, i) - 1));
					}
				}
			}
		}
	}

	public double[][] getImputedDosages() {
		return getImputedDosagesAsMatrix2D().toArray();
	}

	public Boolean alleleFlip(VCFVariant var2) {

		if (isBiallelic() && var2.isBiallelic()) {
			return BaseAnnot.flipalleles(allelesAsString(), minorAllele, var2.allelesAsString(), var2.getMinorAllele());
		} else {
			System.out.println("WARNING: multi allelic allele flip not implemented at this point!");
			return false;
		}


	}


	//???
	public double[][] getDosages() {
		int nrAlleles = getAlleles().length;
		double[][] output = new double[genotypeAlleles.columns()][nrAlleles - 1]; // allow for multiple alleles
		for (int i = 0; i < output.length; i++) {
			for (int j = 0; j < genotypeAlleles.rows(); j++) {
				int a = (int) genotypeAlleles.getQuick(j, i);
				if (a == -1) {
					// missing genotype
					output[i][0] = Double.NaN;
				} else if (a > 0) {
					if (a > 0) {
						a -= 1;
					}
					output[i][a]++;
				}

			}
		}

		return output;
	}

	public String allelesAsString() {
		return Strings.concat(alleles, Pattern.compile("/"));
	}

	public void recodeAlleles(HashMap<String, Integer> alleleMap, String[] newAlleles) {


		int[] alleleRecode = new int[alleles.length];
		boolean allelesremoved = false;
		int allelesremain = 0;
		ArrayList<String> remainingalleles = new ArrayList<String>(3);
		for (int i = 0; i < alleles.length; i++) {
			Integer alleleCd = alleleMap.get(alleles[i]);
			;
			if (alleleCd != null) {
				alleleRecode[i] = alleleCd;
				allelesremain++;
				remainingalleles.add(alleles[i]);
			} else {
				System.out.println("Removing allele: " + alleles[i] + " from variant: " + chr + ":" + pos);
				alleleRecode[i] = -1;
				allelesremoved = true;
			}
		}
		if (allelesremoved) {
			System.out.println(allelesremain + " alleles remain for variant " + chr + ":" + pos + " prev: " + Strings.concat(alleles, Strings.comma) + "\tnew: " + Strings.concat(remainingalleles, Strings.comma));
		}


		if (genotypeAlleles != null) {
			for (int i = 0; i < genotypeAlleles.columns(); i++) {
				for (int j = 0; j < genotypeAlleles.rows(); j++) {
					if (genotypeAlleles.getQuick(j, i) == -1) {
						genotypeAlleles.setQuick(j, i, (byte) -1);
					} else {
						if (alleleRecode[genotypeAlleles.getQuick(j, i)] == -1) {
							System.err.println("Allele " + alleleRecode[genotypeAlleles.getQuick(j, i)] + " removed!");
							System.exit(-1);
						}
						genotypeAlleles.setQuick(j, i, (byte) alleleRecode[genotypeAlleles.getQuick(j, i)]);
					}
				}
			}
		}

		alleles = newAlleles;
	}

	public void convertAllelesToComplement() {
		alleles = convertToComplement(alleles);
	}

	private String[] convertToComplement(String[] alleles2) {
		String[] complement = new String[alleles2.length];
		for (int i = 0; i < complement.length; i++) {
			String allele = alleles2[i];
			complement[i] = getComplement(allele);

		}
		return complement;
	}

	private String getComplement(String allele) {
		String out = "";
		for (int j = 0; j < allele.length(); j++) {
			char c = allele.charAt(j);
			if (c == 'A') {
				out += "T";
			} else if (c == 'T') {
				out += "A";
			} else if (c == 'G') {
				out += "C";
			} else if (c == 'C') {
				out += "G";
			} else {
				out += "N";
			}
		}
		return out;
	}

	public String toVCFString() {
		StringBuilder builder = new StringBuilder(100000);
		builder.append(chr);
		builder.append("\t");
		builder.append(pos);
		builder.append("\t");
		builder.append(id);
		builder.append("\t");
		builder.append(alleles[0]);
		builder.append("\t");
		builder.append(Strings.concat(alleles, Strings.comma, 1, alleles.length));
		builder.append("\t.\t.\t").append(getInfoString()).append("\tGT");

		if (hasImputationDosages()) {
			builder.append(":DS");
		}

		for (int i = 0; i < genotypeAlleles.columns(); i++) {
			byte a1 = genotypeAlleles.getQuick(0, i);
			byte a2 = genotypeAlleles.getQuick(1, i);
			String al1 = "" + a1;
			String al2 = "" + a2;
			if (a1 == -1) {
				al1 = ".";
			}
			if (a2 == -1) {
				al2 = ".";
			}
			builder.append("\t");
			builder.append(al1).append(separator).append(al2);
			if (hasImputationDosages()) {

				// samples x alleles
				builder.append(":");
				for (int a = 0; a < imputedDosages.columns(); a++) {
					if (a == 0) {
						builder.append(imputedDosages.getQuick(i, a));
					} else {
						builder.append(",").append(imputedDosages.getQuick(i, a));
					}
				}

			}

		}
		return builder.toString();
	}

	@Override
	public String toString() {
		return this.chr + "_" + this.pos + "_" + this.id;
//		return this.chr + "-" + this.pos;
	}

	public double[] getGenotypeDosages() {
		return genotypeDosages;
	}

//	public double[][] getGenotypeProbabilies() {
//		return genotypeProbabilies.toArray();
//	}

	public HashMap<String, String> getInfo() {
		return info;
	}

	public boolean hasImputationDosages() {
		return imputedDosages != null;
	}

	public void recalculateMAFAndCallRate() {
		recalculateMAFAndCallRate(null);
	}

	public void recalculateMAFAndCallRate(Boolean[] individualIsFemale) {

		int nrCalled = 0;
		nrAllelesObserved = new int[nrAllelesObserved.length];
		int[] nrAllelesObservedLocal = new int[nrAllelesObserved.length];
		int nrIndividuals = genotypeAlleles.columns();
		Chromosome chromosome = Chromosome.parseChr(chr);
		if (individualIsFemale != null) {
			int nrFemales = 0;
			int nrMales = 0;
			for (int i = 0; i < nrIndividuals; i++) {
				if (individualIsFemale[i] != null) {
					if (individualIsFemale[i]) {
						nrFemales++;
					} else {
						nrMales++;
					}
				}
			}
			if (chromosome.equals(Chromosome.X)) {
				nrIndividuals = nrFemales;
			} else if (chromosome.equals(Chromosome.Y)) {
				nrIndividuals = nrMales;
			} else {
				nrIndividuals = nrIndividuals;
			}
		}


		for (int i = 0; i < genotypeAlleles.columns(); i++) {
			if (chromosome.equals(Chromosome.X)) {
				if (ignoregender || (individualIsFemale != null && individualIsFemale[i] != null && individualIsFemale[i])) {
					if (genotypeAlleles.getQuick(0, i) != -1) {
						nrCalled++;
						nrAllelesObservedLocal[genotypeAlleles.getQuick(0, i)]++;
						nrAllelesObservedLocal[genotypeAlleles.getQuick(1, i)]++;
					}
				} else if (individualIsFemale == null) {
//					System.err.println("ERROR: cannot calculate chr X MAF if gender information unavailable.");
					throw new IllegalArgumentException("ERROR: cannot calculate chr X MAF if gender information unavailable.");
				}
			} else if (chromosome.equals(Chromosome.Y)) {
				if (ignoregender || (individualIsFemale != null && individualIsFemale[i] != null && !individualIsFemale[i])) {
					if (genotypeAlleles.getQuick(0, i) != -1) {
						nrCalled++;
						nrAllelesObservedLocal[genotypeAlleles.getQuick(0, i)]++;
						nrAllelesObservedLocal[genotypeAlleles.getQuick(1, i)]++;
					}
				} else if (individualIsFemale == null) {
					System.err.println("ERROR: cannot calculate chr Y MAF if gender information unavailable.");
					throw new IllegalArgumentException("ERROR: cannot calculate chr X MAF if gender information unavailable.");
				}
			} else {
				if (genotypeAlleles.getQuick(0, i) != -1) {
					nrCalled++;
					nrAllelesObservedLocal[genotypeAlleles.getQuick(0, i)]++;
					nrAllelesObservedLocal[genotypeAlleles.getQuick(1, i)]++;

				}
			}

			if (genotypeAlleles.getQuick(0, i) != -1) {
				nrAllelesObserved[genotypeAlleles.getQuick(0, i)]++;
				nrAllelesObserved[genotypeAlleles.getQuick(1, i)]++;
			}
		}

		callrate = (double) nrCalled / nrIndividuals;

		int totalAllelesObs = nrCalled * 2;

		int nrAllelesThatHaveAlleleFrequency = 0;
		double minAlleleFreq = 2;
		alleleFrequencies = new double[nrAllelesObserved.length];
		minorAllele = null;
		totalCalledAlleles = totalAllelesObs;

		for (int i = 0; i < nrAllelesObserved.length; i++) {
			double alleleFreq = (double) nrAllelesObservedLocal[i] / totalAllelesObs;
			alleleFrequencies[i] = alleleFreq;

			if (nrAllelesObservedLocal[i] > 0) {
				nrAllelesThatHaveAlleleFrequency++;
				if (alleleFreq < minAlleleFreq) {
					if (i == 0) {
						minorAllele = alleles[0];
					} else {
						minorAllele = alleles[i];
					}
					minAlleleFreq = alleleFreq;
				}
			}
		}

		MAF = minAlleleFreq;
		if (MAF == 1) { // flip alleles if monomorphic
			MAF = 0;
			if (minorAllele.equals(alleles[0])) {
				minorAllele = alleles[1];
			} else {
				minorAllele = alleles[0];
			}
		}

		calculateHWEP();

		if (nrAllelesThatHaveAlleleFrequency == 2) {
			biallelic = true;

			// TODO: calculate HWE P
			hwep = 0;
		} else if (nrAllelesThatHaveAlleleFrequency > 2) {
			multiallelic = true;
		} else {
			monomorphic = true;
		}
	}

	public void calculateHWEP() {


		int nrAlleles = alleles.length;
		if (nrAlleles == 2 && 1 == 2) {

			int hets = 0;
			int homs1 = 0;
			int homs2 = 0;

			for (int i = 0; i < genotypeAlleles.columns(); i++) {
				byte a1 = genotypeAlleles.getQuick(0, i);
				if (a1 != -1) {
					byte a2 = genotypeAlleles.getQuick(1, i);
					if (a1 == a2) {
						if (a1 == 0) {
							homs1++;
						} else {
							homs2++;
						}
					} else {
						hets++;
					}
				}
			}

			hwep = HWE.calculateExactHWEPValue(hets, homs1, homs2);
		} else {
			int nrCombinations = (nrAlleles * (nrAlleles + 1)) / 2;
			int[] obs = new int[nrCombinations];
			int[] nrHomozygous = new int[nrAlleles];
			int nrCalled = 0;
			int[][] index = new int[nrAlleles][nrAlleles];

			int ctr = 0;
			for (int i = 0; i < nrAlleles; i++) {
				for (int j = i; j < nrAlleles; j++) {
					index[i][j] = ctr;
					index[j][i] = ctr;
					ctr++;
				}
			}

			// count the homozygous
			for (int i = 0; i < genotypeAlleles.columns(); i++) {
				byte a1 = genotypeAlleles.getQuick(0, i);
				if (a1 != -1) {
					byte a2 = genotypeAlleles.getQuick(1, i);
					if (a1 == a2) {
						nrHomozygous[a1]++;
					}

					int id = index[a1][a2];
					obs[id]++;
					nrCalled++;
				}
			}

			double[] freqs = new double[nrAlleles];
			for (int i = 0; i < nrAlleles; i++) {
				for (int j = 0; j < nrAlleles; j++) {
					int id = index[i][j];
					if (i == j) {
						freqs[i] += (2 * obs[id]);
					} else {
						freqs[i] += obs[id];
					}
				}
			}
			for (int i = 0; i < nrAlleles; i++) {
				freqs[i] /= (nrCalled * 2);
			}

			ctr = 0;
			double chisq = 0;
			for (int i = 0; i < nrAlleles; i++) {
				for (int j = i; j < nrAlleles; j++) {
					double expectedFreq;
					if (i == j) {
						expectedFreq = (freqs[i] * freqs[i]) * nrCalled; // homozygote freq
					} else {
						expectedFreq = (2 * (freqs[i] * freqs[j])) * nrCalled; // heterozygote freq
					}
					double observedFreq = obs[ctr];
					double oe = (observedFreq - expectedFreq);
					if (oe != 0 && expectedFreq != 0) {
						chisq += ((oe * oe) / expectedFreq);
					}
					ctr++;
				}
			}

			int df = (nrCombinations - nrAlleles);
			hwep = umcg.genetica.math.stats.ChiSquare.getP(df, chisq);
		}

	}

	public String getInfoString() {

		if (info == null || info.isEmpty()) {
			return ".";
		} else {
			Set<String> keys = info.keySet();
			String infoStr = "";
			for (String key : keys) {
				Object infoObj = info.get(key);
				if (infoObj == null) {
					if (infoStr.length() == 0) {
						infoStr += key;
					} else {
						infoStr += ";" + key;
					}
				} else {
					if (infoStr.length() == 0) {
						infoStr += key + "=" + infoObj.toString();
					} else {
						infoStr += ";" + key + "=" + infoObj.toString();
					}
				}
			}
			return infoStr;
		}
	}

	public Feature asFeature() {
		Feature output = new Feature(Chromosome.parseChr(chr), pos, pos);
		output.setName(id);
		return output;

	}


	public int getGTCol() {
		return gtCol;
	}

	public String getSeparator() {
		return separator;
	}

	public boolean isIndel() {
		String[] alleles = getAlleles();
		boolean isIndel = false;
		for (String s : alleles) {
			if (s.length() > 1) {
				return true;
			}
		}
		return false;
	}

	public Chromosome getChrObj() {
		if (chr == null) {
			return Chromosome.NA;
		} else {
			return Chromosome.parseChr(chr);
		}

	}

	public byte[] getGenotypesAsByteVector() {
		byte[][] alleles = getGenotypeAlleles();
		byte[] output = new byte[alleles[0].length];
		for (int i = 0; i < alleles[0].length; i++) {
			if (alleles[0][i] == -1) {
				output[i] = -1;
			} else {
				output[i] = (byte) (alleles[0][i] + alleles[1][i]);
			}
		}
		return output;
	}

	public String getMinorAlleleFromInfoField() {

		String alleleCounts = info.get("AC");
		String totalCounts = info.get("AN");

		if (alleleCounts != null && totalCounts != null) {
			String[] elems = alleleCounts.split(",");
			double[] freqs = new double[elems.length];
			double totalAlleles = Double.parseDouble(totalCounts);
			double refFreq = 1;
			for (int i = 0; i < elems.length; i++) {
				double ct = Double.parseDouble(elems[i]);
				freqs[i] = ct / totalAlleles;
				refFreq -= freqs[i];
			}

			int minorAlleleNr = 0;
			double minorAlleleFreq = refFreq;
			for (int i = 0; i < elems.length; i++) {
				if (freqs[i] < minorAlleleFreq) {
					minorAlleleNr = i + 1;
					minorAlleleFreq = freqs[i];
				}
			}

			MAF = minorAlleleFreq;
			minorAllele = alleles[minorAlleleNr];
			return minorAllele;
		}

		return null;

	}

	public ByteMatrix2D getGenotypeAllelesAsMatrix2D() {
		return genotypeAlleles;
	}

	public DoubleMatrix2D getImputedDosagesAsMatrix2D() {
		// parse genotype probs
		int nrAlleles = getAlleles().length;
		DoubleMatrix2D dosages = null;
		if (imputedDosages == null) {
			dosages = new DenseDoubleMatrix2D(genotypeAlleles.columns(), nrAlleles - 1);
			for (int i = 0; i < genotypeAlleles.columns(); i++) {
				for (int j = 0; j < genotypeAlleles.rows(); j++) {
					dosages.setQuick(i, 0, genotypeAlleles.getQuick(j, i));
				}
			}
			return dosages;
		} else {
			return imputedDosages;
		}
	}

	public enum PARSE {
		HEADER,
		GENOTYPES,
		ALL
	}
}
