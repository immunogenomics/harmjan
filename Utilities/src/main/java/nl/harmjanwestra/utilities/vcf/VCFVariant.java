package nl.harmjanwestra.utilities.vcf;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.enums.Gender;
import nl.harmjanwestra.utilities.features.Feature;
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

enum PHASE {
	UNPHASED("/", Pattern.compile("/")),
	PHASED("|", Pattern.compile("\\|"));

	private final String ph;
	private final Pattern pattern;

	PHASE(String ph, Pattern pattern) {
		this.ph = ph;
		this.pattern = pattern;
	}

	public static PHASE getPhase(String character) {
		if (character.contains("|")) {
			return PHASE.PHASED;
		} else {
			return PHASE.UNPHASED;
		}
	}

	public String getSep() {
		return ph;
	}

	public Pattern getPattern() {
		return pattern;
	}

	public String toString() {
		return ph;
	}
}

public class VCFVariant {

	private static final int nrHeaderElems = 9;

	private static final Pattern nullGenotype = Pattern.compile("\\./\\.");
	private final HashMap<String, String> info = new HashMap<String, String>();
	private boolean[] samplesToInclude;
	private ArrayList<VCFGenotypeFilter> filters;


	private DoubleMatrix2D genotypeAlleles; // format [individuals][alleles] (save some memory by making only two individual-sized arrays)
	private DoubleMatrix2D genotypeProbabilies;
	private DoubleMatrix2D dosages; // this can hold the imputed dosages, or if they are not set, the dosages from the genotypes
	private ShortMatrix2D allelicDepth;
	private int[] nrAllelesObserved;
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
	// private String separator = new String("/").intern();

	private PHASE phase;
	private boolean ignoregender;
	private int totalCalledAlleles;


	private int[] nrAllelesObservedCases;
	private int[] nrAllelesObservedControls;
	private double[] alleleFrequenciesCases;
	private double[] alleleFrequenciesControls;
	private double hwepCases;
	private double hwepControls;
	int constructor = 0;

	private SampleAnnotation sampleAnnotation;
	private double MAFControls;

	public VCFVariant(String ln) {
		constructor = 1;
		parse(ln, PARSE.ALL);
	}

	public VCFVariant(String ln, boolean ignoregender) {
		constructor = 2;
		this.ignoregender = ignoregender;
		parse(ln, PARSE.ALL);

	}

	public VCFVariant(String ln, ArrayList<VCFGenotypeFilter> filters, boolean ignoregender, SampleAnnotation annotation) {
		constructor = 3;
		this.ignoregender = ignoregender;
		this.sampleAnnotation = annotation;
		this.filters = filters;
		parse(ln, PARSE.ALL);

	}

	public VCFVariant(String ln, ArrayList<VCFGenotypeFilter> filters, boolean ignoregender) {
		constructor = 4;
		this.ignoregender = ignoregender;
		this.filters = filters;
		parse(ln, PARSE.ALL);
	}

	public VCFVariant(String ln, PARSE p) {
		constructor = 5;
		parse(ln, p);
	}

	public VCFVariant(String ln, PARSE p, boolean[] samplesToInclude) {
		constructor = 6;
		this.samplesToInclude = samplesToInclude;
		parse(ln, p);
	}

	public VCFVariant(String ln, PARSE p, boolean[] samplesToInclude, SampleAnnotation sampleAnnotation) {
		constructor = 7;
		this.samplesToInclude = samplesToInclude;
		this.sampleAnnotation = sampleAnnotation;
		parse(ln, p);
	}

	public VCFVariant(String ln, PARSE p, SampleAnnotation sampleAnnotation) {
		constructor = 8;
		this.sampleAnnotation = sampleAnnotation;
		parse(ln, p);
	}

	public VCFVariant(String chr, Integer pos, String id, String alleleStr, String info, DoubleMatrix2D alleles, DoubleMatrix2D dosages, SampleAnnotation annotation) {
		constructor = 9;
		this.chr = new String(chr).intern();
		this.pos = pos;
		this.id = id;
		String[] allelesElems = alleleStr.split(",");
		this.alleles = allelesElems;
		parseInfoString(info);
		this.genotypeAlleles = alleles;
		this.dosages = dosages;
		this.sampleAnnotation = annotation;
		recalculateMAFAndCallRate();

	}

	public DoubleMatrix2D getGenotypeProbabilies() {
		return genotypeProbabilies;
	}


	public short[][] getAllelicDepth() {
		return allelicDepth.toArray();
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

	public void parse(String ln, PARSE p) {

		// GT:AB:AD:DP:GQ:PL
		int gtCol = -1; // genotype
		int abCol = -1; // allelic balance
		int adCol = -1; // Allelic depths for the ref and alt alleles in the order listed
		int dpCol = -1; // Approximate readAsTrack depth (reads with MQ=255 or with bad mates are filtered)
		int gqCol = -1; // Genotype Quality
		int plCol = -1; // Normalized, Phred-scaled likelihoods for genotypes
		int pidCol = -1; // ?
		int pgtCol = -1; // ?
		int dsCol = -1;
		int gpCol = -1;

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

		if (p.equals(PARSE.ALL) || p.equals(PARSE.GENOTYPES)) {
			if (ln != null) {

				String[] tokenArr = Strings.tab.split(ln);

				if (tokenArr.length > 9) { // allow VCFs without any actual genotypes
					// parse actual genotypes.
					int nrTokens = tokenArr.length;
					int nrSamples = nrTokens - nrHeaderElems;
					if (samplesToInclude != null) {
						nrSamples = 0;
						for (int i = 0; i < samplesToInclude.length; i++) {
							if (samplesToInclude[i]) {
								nrSamples++;
							}
						}
					}


					genotypeAlleles = DoubleFactory2D.dense.make(nrSamples, 2, -1);// new DenseDoubleMatrix2D(nrSamples, 2);


					genotypeProbabilies = null;

					int includedSampleCtr = 0;
					for (int t = 9; t < nrTokens; t++) {
						String token = tokenArr[t];
						int indPos = t - nrHeaderElems;
						if (samplesToInclude == null || samplesToInclude[indPos]) {
							String sampleColumn = token;
							if (nullGenotype.equals(sampleColumn)) {
								// not called
								genotypeAlleles.setQuick(includedSampleCtr, 0, -1);
								genotypeAlleles.setQuick(includedSampleCtr, 1, -1);

							} else {
								//String[] sampleElems = Strings.colon.split(sampleColumn);

								String[] sampleTokens = Strings.colon.split(sampleColumn);
								for (int s = 0; s < sampleTokens.length; s++) {
									String sampleToken = sampleTokens[s];

									if (s == gtCol) {
										String gt = sampleToken;
										if (!nullGenotype.equals(gt)) {

											PHASE phase = PHASE.getPhase(gt);
											String[] gtElems = phase.getPattern().split(gt);
											this.phase = phase;
//											if (gtElems.length == 1) { // phased genotypes ??
//												gtElems = pipe.split(gt);
//												phase = PHASE.PHASED;
//											}

											if (!gtElems[0].equals(".")) {

												if (gtElems.length < 2) {
													System.out.println(ln);
													System.out.println();
													System.out.println("gtelems: " + Strings.concat(gtElems, Strings.tab));
													System.out.println("length of GT Elems < 2");
													System.out.println("Actual length: " + gtElems.length);
													System.out.println("Number of total sample elems: " + sampleTokens.length);
													System.out.println("Number of total tokens: " + nrTokens);
													System.exit(-1);
												}

												try {
													byte gt1 = 0;
													byte gt2 = 0;
													gt1 = Byte.parseByte(gtElems[0]);
													gt2 = Byte.parseByte(gtElems[1]);

													genotypeAlleles.setQuick(includedSampleCtr, 0, gt1);
													genotypeAlleles.setQuick(includedSampleCtr, 1, gt2);


												} catch (NumberFormatException e) {
													System.out.println("Cannot parse genotype string: " + token + " nr elems: " + gtElems.length);

												}
											}
										}
									} else if (s == dsCol) {
										// dosages
										if (p.equals(PARSE.ALL)) {
											try {
												if (dosages == null) {
													dosages = new DenseDoubleMatrix2D(nrSamples, alleles.length - 1);
												}

												String[] dsElems = Strings.comma.split(sampleToken);
												for (int q = 0; q < dsElems.length; q++) {
													dosages.set(includedSampleCtr, q, Double.parseDouble(dsElems[q]));
												}

											} catch (NumberFormatException e) {

											}
										}
									} else if (s == gpCol) {
										// genotype probs
										if (p.equals(PARSE.ALL)) {
											String[] gpElems = Strings.comma.split(sampleToken);

											try {
												if (genotypeProbabilies == null) {
													genotypeProbabilies = new DenseDoubleMatrix2D(nrSamples, gpElems.length);
												}

												for (int g = 0; g < gpElems.length; g++) {
													genotypeProbabilies.setQuick(includedSampleCtr, g, Double.parseDouble(gpElems[g]));
												}
											} catch (NumberFormatException e) {

											}
										}
									} else if (s == adCol) {
										// depth of sequencing per allele
										String[] adElems = Strings.comma.split(sampleToken);
										try {
											if (allelicDepth == null) {
												allelicDepth = new ShortMatrix2D(nrSamples, alleles.length);
											}

//									try {
											for (int g = 0; g < adElems.length; g++) {
												allelicDepth.setQuick(includedSampleCtr, g, Short.parseShort(adElems[g]));
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
										if (genotypeQuals == null) {
											genotypeQuals = new short[nrSamples];
										}
										try {
											if (gqCol != -1) {
												gq = Short.parseShort(sampleToken);
											}
										} catch (NumberFormatException e) {
										}
										genotypeQuals[includedSampleCtr] = gq;

									} else if (s == dpCol) {
										// approximate depth of sequencing
										short depth = 0;

										if (approximateDepth == null) {
											approximateDepth = new short[nrSamples];

										}

										try {
											if (dpCol != -1) {
												depth = Short.parseShort(sampleToken);
											}
										} catch (NumberFormatException e) {
										}
										approximateDepth[includedSampleCtr] = depth;
									}
								}
							}
							includedSampleCtr++;
						}

					}

					if (genotypeProbabilies != null && dosages == null) {
						int nrAlleles = alleles.length;
						dosages = new DenseDoubleMatrix2D(genotypeProbabilies.rows(), nrAlleles - 1);
						for (int i = 0; i < genotypeProbabilies.rows(); i++) {
							int alctr = 0;
							for (int a1 = 0; a1 < nrAlleles; a1++) {
								for (int a2 = a1; a2 < nrAlleles; a2++) {
									double dosageval = genotypeProbabilies.getQuick(i, alctr);
									if (a1 > 0) {
										dosages.setQuick(i, a1 - 1, dosages.getQuick(i, a1 - 1) + dosageval);
									}
									if (a2 > 0) {
										dosages.setQuick(i, a2 - 1, dosages.getQuick(i, a2 - 1) + dosageval);
									}
									alctr++;
								}
							}
						}
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

	public double[][] getGenotypeAlleles() {
		return genotypeAlleles.toArray();
	}

	public void flipReferenceAlelele() {

		String allele0 = alleles[0];
		String allele1 = alleles[1];
		alleles[0] = allele1;
		alleles[1] = allele0;

		if (genotypeAlleles != null) {

			// only works for biallelic variants!
			for (int i = 0; i < genotypeAlleles.rows(); i++) {
				for (int j = 0; j < genotypeAlleles.columns(); j++) {
					if (genotypeAlleles.getQuick(i, j) == -1) {
						genotypeAlleles.setQuick(i, j, (byte) -1);
					} else {
						genotypeAlleles.setQuick(i, j, (byte) Math.abs(genotypeAlleles.getQuick(i, j) - 1));
					}
				}
			}
		}
	}


	public Boolean alleleFlip(VCFVariant var2) {

		if (isBiallelic() && var2.isBiallelic()) {
			return BaseAnnot.flipalleles(allelesAsString(), minorAllele, var2.allelesAsString(), var2.getMinorAllele());
		} else {
			System.out.println("WARNING: multi allelic allele flip not implemented at this point!");
			return false;
		}
	}

	public DoubleMatrix2D getGenotypeDosagesAsMatrix2D() {
		DoubleMatrix2D gtdosage = new DenseDoubleMatrix2D(genotypeAlleles.rows(), alleles.length - 1);

		for (int i = 0; i < genotypeAlleles.rows(); i++) {
			for (int j = 0; j < genotypeAlleles.columns(); j++) {
				int allele = (int) genotypeAlleles.getQuick(i, j);
				if (allele == -1) {
					for (int q = 0; q < gtdosage.columns(); q++) {
						gtdosage.setQuick(i, q, -1);
					}
				} else if (allele > 0) {
					allele -= 1;
					double ct = gtdosage.getQuick(i, allele);
					ct++;
					gtdosage.setQuick(i, allele, ct);
				}
			}
		}

		return gtdosage;
	}

	public double[][] getDosage() {
		return getDosagesAsMatrix2D().toArray();
	}

	public double[][] getGenotypeDosage() {
		return getGenotypeDosagesAsMatrix2D().toArray();
	}

	public DoubleMatrix2D getDosagesAsMatrix2D() {
		// parse genotype probs
		if (this.dosages == null) {
			return getGenotypeDosagesAsMatrix2D();
		} else {
			return this.dosages;
		}
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
			for (int i = 0; i < genotypeAlleles.rows(); i++) {
				for (int j = 0; j < genotypeAlleles.columns(); j++) {
					if (genotypeAlleles.getQuick(i, j) == -1) {
						genotypeAlleles.setQuick(i, j, (byte) -1);
					} else {
						int gt = (int) genotypeAlleles.getQuick(i, j);
						if (alleleRecode[gt] == -1) {
							System.err.println("Allele " + alleleRecode[gt] + " removed!");
							System.exit(-1);
						}
						genotypeAlleles.setQuick(i, j, (byte) alleleRecode[gt]);
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
		return toVCFString(true);
	}

	public String toVCFString(boolean includeHeader) {
		StringBuilder builder = new StringBuilder(100000);
		if (includeHeader) {
			builder.append(toVCFHeader());
		}


		for (int i = 0; i < genotypeAlleles.rows(); i++) {
			byte a1 = (byte) genotypeAlleles.getQuick(i, 0);
			byte a2 = (byte) genotypeAlleles.getQuick(i, 1);
			String al1 = "" + a1;
			String al2 = "" + a2;
			if (a1 == -1) {
				al1 = ".";
			}
			if (a2 == -1) {
				al2 = ".";
			}
			if (builder.length() > 0) {
				builder.append("\t");
			}
			builder.append(al1).append(phase.toString()).append(al2);
			if (genotypeProbabilies != null) {

				// samples x alleles
				builder.append(":");
				for (int a = 0; a < genotypeProbabilies.columns(); a++) {
					if (a == 0) {
						builder.append(genotypeProbabilies.getQuick(i, a));
					} else {
						builder.append(",").append(genotypeProbabilies.getQuick(i, a));
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


//	public double[][] getGenotypeProbabilies() {
//		return genotypeProbabilies.toArray();
//	}

	public HashMap<String, String> getInfo() {
		return info;
	}

	public boolean hasImputationDosages() {
		return dosages != null;
	}

	public int[] getNrAllelesObservedCases() {
		return nrAllelesObservedCases;
	}

	public int[] getNrAllelesObservedControls() {
		return nrAllelesObservedControls;
	}

	public double[] getAlleleFrequenciesCases() {
		return alleleFrequenciesCases;
	}

	public double[] getAlleleFrequenciesControls() {
		return alleleFrequenciesControls;
	}

	public void recalculateMAFAndCallRate() {

		DiseaseStatus[][] sampleDiseaseStatus = null;
		Gender[] individualGender = null;
		if (sampleAnnotation != null) {

			if (sampleAnnotation.getIndividualGender() != null && sampleAnnotation.getIndividualGender().length != genotypeAlleles.rows()) {
				throw new IllegalArgumentException("Sample gender status length does not match number of loaded genotypes. " + sampleAnnotation.getIndividualGender().length + " gender status vs " + genotypeAlleles.rows() + " genotypes");
			}
			if (sampleAnnotation.getSampleDiseaseStatus() != null && sampleAnnotation.getSampleDiseaseStatus().length != genotypeAlleles.rows()) {
				throw new IllegalArgumentException("Sample disease status length does not match number of loaded genotypes. " + sampleAnnotation.getSampleDiseaseStatus().length + " disease status vs " + genotypeAlleles.rows() + " genotypes");
			}

			sampleDiseaseStatus = sampleAnnotation.getSampleDiseaseStatus();
			individualGender = sampleAnnotation.getIndividualGender();
		} else {
//			System.out.println("sampleAnnotation is null");
//			Exception e = new Exception("");
//			System.out.println("Constructor: " + constructor);
//			e.printStackTrace();
//			System.exit(-1);
		}

		int nrCalled = 0;
		int nrCalledCases = 0;
		int nrCalledControls = 0;
		nrAllelesObserved = new int[alleles.length];

		if (sampleDiseaseStatus != null) {
			nrAllelesObservedCases = new int[alleles.length];
			nrAllelesObservedControls = new int[alleles.length];
		}

		int nrIndividuals = genotypeAlleles.rows();
		Chromosome chromosome = getChrObj();
		if (individualGender != null) {
			int nrFemales = 0;
			int nrMales = 0;
			for (int i = 0; i < nrIndividuals; i++) {
				Gender gender = individualGender[i];
				if (gender != null && gender == Gender.MALE) {
					nrMales++;
				} else if (gender != null && gender == Gender.FEMALE) {
					nrFemales++;
				}
			}
			if (chromosome.equals(Chromosome.X)) {
				nrIndividuals = nrFemales;
			} else if (chromosome.equals(Chromosome.Y)) {
				nrIndividuals = nrMales;
			}
		}


		for (int i = 0; i < genotypeAlleles.rows(); i++) {
			Gender gender = null;
			DiseaseStatus diseaseStatus = null;
			if (individualGender != null) {
				gender = individualGender[i];
			}
			if (sampleDiseaseStatus != null) {
				diseaseStatus = sampleDiseaseStatus[i][0];
			}

			int gt1 = (int) genotypeAlleles.getQuick(i, 0);
			if (gt1 != -1) {
				int gt2 = (int) genotypeAlleles.getQuick(i, 1);
				if (chromosome.isAutosome()) {
					nrCalled++;
					nrAllelesObserved[gt1]++;
					nrAllelesObserved[gt2]++;
					if (diseaseStatus != null) {
						if (diseaseStatus == DiseaseStatus.CASE) {
							nrAllelesObservedCases[gt1]++;
							nrAllelesObservedCases[gt2]++;
							nrCalledCases++;
						} else if (diseaseStatus == DiseaseStatus.CONTROL) {
							nrAllelesObservedControls[gt1]++;
							nrAllelesObservedControls[gt2]++;
							nrCalledControls++;
						}
					}
				} else {
					if (chromosome.equals(Chromosome.X)) {
						if (ignoregender || (gender != null && gender == Gender.FEMALE)) {
							nrCalled++;
							nrAllelesObserved[gt1]++;
							nrAllelesObserved[gt2]++;
							if (diseaseStatus != null) {
								if (diseaseStatus == DiseaseStatus.CASE) {
									nrAllelesObservedCases[gt1]++;
									nrAllelesObservedCases[gt2]++;
									nrCalledCases++;
								} else if (diseaseStatus == DiseaseStatus.CONTROL) {
									nrAllelesObservedControls[gt1]++;
									nrAllelesObservedControls[gt2]++;
									nrCalledControls++;
								}
							}

						} else if (individualGender == null) {
//						System.err.println("ERROR: cannot calculateWithinDataset chr X MAF if gender information unavailable.");
//						throw new IllegalArgumentException("ERROR: cannot calculateWithinDataset chr X MAF if gender information unavailable.");
						}
					} else if (chromosome.equals(Chromosome.Y)) {
						if (ignoregender || (gender != null && gender == Gender.MALE)) {
							nrCalled++;
							nrAllelesObserved[gt1]++;
							nrAllelesObserved[gt2]++;
							if (diseaseStatus != null) {
								if (diseaseStatus == DiseaseStatus.CASE) {
									nrAllelesObservedCases[gt1]++;
									nrAllelesObservedCases[gt2]++;
									nrCalledCases++;
								} else if (diseaseStatus == DiseaseStatus.CONTROL) {
									nrAllelesObservedControls[gt1]++;
									nrAllelesObservedControls[gt2]++;
									nrCalledControls++;
								}
							}
						} else if (individualGender == null) {
//						System.err.println("ERROR: cannot calculateWithinDataset chr Y MAF if gender information unavailable.");
						}
					}
				}


			}
		}

		alleleFrequencies = new double[nrAllelesObserved.length];
		if (sampleDiseaseStatus != null) {
			alleleFrequenciesCases = new double[nrAllelesObserved.length];
			alleleFrequenciesControls = new double[nrAllelesObserved.length];
		}

		if (nrCalled == 0) {
			callrate = 0;
			MAF = 0;
			MAFControls = 0;
			minorAllele = null;
			totalCalledAlleles = 0;

			hwep = 0;
			hwepCases = 0;
			hwepControls = 0;
			for (int i = 0; i < nrAllelesObserved.length; i++) {
				alleleFrequencies[i] = 0;
				if (sampleDiseaseStatus != null) {
					alleleFrequenciesCases[i] = 0;
					alleleFrequenciesControls[i] = 0;
				}
			}

			monomorphic = true;
		} else {
			callrate = (double) nrCalled / nrIndividuals;
			int totalAllelesObs = nrCalled * 2;
			int nrAllelesThatHaveAlleleFrequency = 0;
			double minAlleleFreq = 2;
			minorAllele = null;
			totalCalledAlleles = totalAllelesObs;
			for (int i = 0; i < nrAllelesObserved.length; i++) {
				double alleleFreq = (double) nrAllelesObserved[i] / totalAllelesObs;
				alleleFrequencies[i] = alleleFreq;

				if (sampleDiseaseStatus != null) {
					alleleFrequenciesCases[i] = (double) nrAllelesObservedCases[i] / (nrCalledCases * 2);
					alleleFrequenciesControls[i] = (double) nrAllelesObservedControls[i] / (nrCalledControls * 2);
				}


				if (nrAllelesObserved[i] > 0) {
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
			if (nrAllelesThatHaveAlleleFrequency == 2) {
				biallelic = true;
			} else if (nrAllelesThatHaveAlleleFrequency > 2) {
				multiallelic = true;
			} else {
				monomorphic = true;
			}
			calculateHWEP();
		}


	}

	public double getHwepCases() {
		return hwepCases;
	}

	public double getHwepControls() {
		return hwepControls;
	}

	public void calculateHWEP() {


		DiseaseStatus[][] sampleDiseaseStatus = null;
		if (sampleAnnotation != null) {
			sampleDiseaseStatus = sampleAnnotation.getSampleDiseaseStatus();
			if (sampleAnnotation.getSampleDiseaseStatus() != null && sampleAnnotation.getSampleDiseaseStatus().length != genotypeAlleles.rows()) {
				throw new IllegalArgumentException("Sample disease status length does not match number of loaded genotypes. " + sampleAnnotation.getSampleDiseaseStatus().length + " disease status vs " + genotypeAlleles.rows() + " genotypes");
			}
		}

		if (sampleDiseaseStatus != null && sampleDiseaseStatus.length != genotypeAlleles.rows()) {
			throw new IllegalArgumentException("Sample disease status length does not match number of loaded genotypes. " + sampleDiseaseStatus.length + " disease status vs " + genotypeAlleles.rows() + " genotypes");
		}

		int nrAlleles = alleles.length;
		if (nrAlleles == 2) {

			int hets = 0;
			int homs1 = 0;
			int homs2 = 0;

			int hetsCases = 0;
			int homs1Cases = 0;
			int homs2Cases = 0;

			int hetsControls = 0;
			int homs1Controls = 0;
			int homs2Controls = 0;

			for (int i = 0; i < genotypeAlleles.rows(); i++) {
				int a1 = (int) genotypeAlleles.getQuick(i, 0);
				DiseaseStatus diseaseStatus = null;
				if (sampleDiseaseStatus != null) {
					diseaseStatus = sampleDiseaseStatus[i][0];
				}
				if (a1 != -1) {
					int a2 = (int) genotypeAlleles.getQuick(i, 1);
					if (a1 == a2) {
						if (a1 == 0) {
							homs1++;
							if (diseaseStatus != null && diseaseStatus == DiseaseStatus.CASE) {
								homs1Cases++;
							} else if (diseaseStatus != null && diseaseStatus == DiseaseStatus.CONTROL) {
								homs1Controls++;
							}
						} else {
							homs2++;
							if (diseaseStatus != null && diseaseStatus == DiseaseStatus.CASE) {
								homs2Cases++;
							} else if (diseaseStatus != null && diseaseStatus == DiseaseStatus.CONTROL) {
								homs2Controls++;
							}
						}
					} else {
						hets++;
						if (diseaseStatus != null && diseaseStatus == DiseaseStatus.CASE) {
							hetsCases++;
						} else if (diseaseStatus != null && diseaseStatus == DiseaseStatus.CONTROL) {
							hetsControls++;
						}
					}
				}
			}

			hwep = HWE.calculateExactHWEPValue(hets, homs1, homs2);
			if (sampleDiseaseStatus != null) {
				hwepCases = HWE.calculateExactHWEPValue(hetsCases, homs1Cases, homs2Cases);
				hwepControls = HWE.calculateExactHWEPValue(hetsControls, homs1Controls, homs2Controls);
			}
		} else {
			// use a chi-square test for multi-allelic variants (less precise with lower allele frequencies)
			int nrCombinations = (nrAlleles * (nrAlleles + 1)) / 2;
			int nrCalled = 0;
			int nrCalledCases = 0;
			int nrCalledControls = 0;

			int[] obs = new int[nrCombinations];
			int[] nrHomozygous = new int[nrAlleles];
			double[] freqs = new double[nrAlleles];

			double[] freqsCases = null;
			double[] freqsControls = null;

			int[] obsCases = null;
			int[] nrHomozygousCases = null;
			int[] obsControls = null;
			int[] nrHomozygousControls = null;

			if (sampleDiseaseStatus != null) {
				obsCases = new int[nrCombinations];
				nrHomozygousCases = new int[nrAlleles];
				obsControls = new int[nrCombinations];
				nrHomozygousControls = new int[nrAlleles];
				freqsCases = new double[nrAlleles];
				freqsControls = new double[nrAlleles];
			}

			// index the allele combinations
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
			for (int i = 0; i < genotypeAlleles.rows(); i++) {
				DiseaseStatus diseaseStatus = null;
				if (sampleDiseaseStatus != null) {
					diseaseStatus = sampleDiseaseStatus[i][0];
				}
				int a1 = (int) genotypeAlleles.getQuick(i, 0);
				if (a1 != -1) {
					int a2 = (int) genotypeAlleles.getQuick(i, 1);
					if (a1 == a2) {
						if (diseaseStatus != null && diseaseStatus == DiseaseStatus.CASE) {
							nrHomozygousCases[a1]++;
						} else if (diseaseStatus != null && diseaseStatus == DiseaseStatus.CONTROL) {
							nrHomozygousControls[a1]++;
						}
						nrHomozygous[a1]++;
					}

					int id = index[a1][a2];
					obs[id]++;
					if (diseaseStatus != null && diseaseStatus == DiseaseStatus.CASE) {
						obsCases[id]++;
						nrCalledCases++;
					} else if (diseaseStatus != null && diseaseStatus == DiseaseStatus.CONTROL) {
						obsControls[id]++;
						nrCalledControls++;
					}
					nrCalled++;
				}
			}


			for (int i = 0; i < nrAlleles; i++) {
				for (int j = 0; j < nrAlleles; j++) {
					int id = index[i][j];
					if (i == j) {
						freqs[i] += (2 * obs[id]);
						if (freqsCases != null) {
							freqsCases[i] += (2 * obsCases[id]);
							freqsControls[i] += (2 * obsControls[id]);
						}
					} else {
						freqs[i] += obs[id];
						if (freqsCases != null) {
							freqsCases[i] += (obsCases[id]);
							freqsControls[i] += (obsControls[id]);
						}
					}
				}
			}
			for (int i = 0; i < nrAlleles; i++) {
				freqs[i] /= (nrCalled * 2);
				if (freqsCases != null) {
					freqsCases[i] /= (nrCalledCases * 2);
					freqsControls[i] /= (nrCalledControls * 2);
				}
			}

			ctr = 0;
			double chisq = 0;
			double chisqCases = 0;
			double chisqControls = 0;
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

					if (freqsCases != null) {
						double expectedFreqCases;
						double expectedFreqControls;
						if (i == j) {
							expectedFreqCases = (freqsCases[i] * freqsCases[i]) * nrCalledCases; // homozygote freq
							expectedFreqControls = (freqsControls[i] * freqsControls[i]) * nrCalledControls; // homozygote freq
						} else {
							expectedFreqCases = (2 * (freqsCases[i] * freqsCases[j])) * nrCalledCases; // heterozygote freq
							expectedFreqControls = (2 * (freqsControls[i] * freqsControls[j])) * nrCalledControls; // heterozygote freq
						}
						double observedFreqCases = obsCases[ctr];
						double observedFreqControls = obsControls[ctr];
						double oeCases = (observedFreqCases - expectedFreqCases);
						double oeControls = (observedFreqControls - expectedFreqControls);
						if (oeCases != 0 && expectedFreqCases != 0) {
							chisqCases += ((oeCases * oeCases) / expectedFreqCases);
						}

						if (oeControls != 0 && expectedFreqControls != 0) {
							chisqControls += ((oeControls * oeControls) / expectedFreqControls);
						}

					}
					ctr++;
				}
			}

			int df = (nrCombinations - nrAlleles);
			hwep = umcg.genetica.math.stats.ChiSquare.getP(df, chisq);
			if (freqsCases != null) {
				hwepCases = umcg.genetica.math.stats.ChiSquare.getP(df, chisqCases);
				hwepControls = umcg.genetica.math.stats.ChiSquare.getP(df, chisqControls);
			}
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

//	public int getGTCol() {
//		return gtCol;
//	}

	public String getSeparator() {
		return phase.toString();
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

//	public byte[] getGenotypesAsByteVector() {
//		byte[][] alleles = getGenotypeAlleles();
//		byte[] output = new byte[alleles[0].length];
//		for (int i = 0; i < alleles[0].length; i++) {
//			if (alleles[0][i] == -1) {
//				output[i] = -1;
//			} else {
//				output[i] = (byte) (alleles[0][i] + alleles[1][i]);
//			}
//		}
//		return output;
//	}

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

	public DoubleMatrix2D getGenotypeAllelesAsMatrix2D() {
		return genotypeAlleles;
	}

	public byte[] getGenotypesAsByteVector() {
		if (alleles.length == 2) {
			DoubleMatrix2D d = getGenotypeDosagesAsMatrix2D();
			byte[] arr = new byte[d.rows()];
			for (int i = 0; i < arr.length; i++) {
				arr[i] = (byte) d.get(i, 0);
			}
			return arr;
		} else {
			return null;
		}

	}

	public int getNrAlleles() {
		return alleles.length;
	}

	public int getNrSamples() {
		return genotypeAlleles.rows();
	}

	public VCFVariant variantFromAllele(int allele) {
		if (allele == 0) {
			throw new IllegalArgumentException("allele should be bigger than 0");
		}
		if (allele > alleles.length) {
			throw new IllegalArgumentException("allele not present");
		}



		/*
			private DoubleMatrix2D genotypeAlleles; // format [individuals][alleles] (save some memory by making only two individual-sized arrays)
	private DoubleMatrix2D genotypeProbabilies;
	private DoubleMatrix2D dosages; // this can hold the imputed dosages, or if they are not set, the dosages from the genotypes
	private ShortMatrix2D allelicDepth;
	private int[] nrAllelesObserved;
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
	private boolean ignoregender;
	private int totalCalledAlleles;
		 */


		DenseDoubleAlgebra dda = new DenseDoubleAlgebra();


		DoubleMatrix2D odosages = dda.subMatrix(dosages, 0, genotypeProbabilies.rows() - 1, allele - 1, allele - 1);

		DoubleMatrix2D ogenotypeAlleles = new DenseDoubleMatrix2D(genotypeAlleles.rows(), genotypeAlleles.columns());
		int alleleindex = allele - 1;
		for (int r = 0; r < genotypeAlleles.rows(); r++) {
			for (int c = 0; c < genotypeAlleles.columns(); c++) {
				if (genotypeAlleles.getQuick(r, c) == alleleindex) {
					ogenotypeAlleles.set(r, c, 1);
				}
			}
		}
		String[] oalleles = new String[]{alleles[0], alleles[allele]};

		VCFVariant output = new VCFVariant(chr,
				pos,
				id + "_" + alleles[allele],
				Strings.concat(oalleles, Strings.comma),
				getInfoString(),
				ogenotypeAlleles,
				odosages,
				sampleAnnotation);

		return output;
	}

	public boolean isPhased() {
		return phase.equals(PHASE.PHASED);
	}

	public boolean isAutosomal() {
		return getChrObj().isAutosome();
	}

	public String toVCFHeader() {
		StringBuilder builder = new StringBuilder();
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

		if (genotypeProbabilies != null) {
			builder.append(":GP");
		}

		return builder.toString();

	}

	public double getMAFControls() {
		if (alleleFrequenciesControls != null) {
			double min = 1;
			for (double d : alleleFrequenciesControls) {
				if (d < min) {
					min = d;
				}
			}
			return min;
		} else {
			return getMAF();
		}
	}

	public double getMAFCases() {
		if (alleleFrequenciesCases != null) {
			double min = 1;
			for (double d : alleleFrequenciesCases) {
				if (d < min) {
					min = d;
				}
			}
			return min;
		} else {
			return getMAF();
		}
	}

	public enum PARSE {
		HEADER,
		GENOTYPES,
		ALL
	}
}
