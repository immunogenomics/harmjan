package nl.harmjanwestra.ngs.GenotypeFormats;

import umcg.genetica.text.Strings;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;

/**
 * Created by hwestra on 4/21/15.
 */
public class VCFVariant {

	private static final int nrHeaderElems = 9;
	private int[] nrAllelesObserved;
	private final byte[][] genotypeAlleles;
	private final HashMap<String, Double> info = new HashMap<String, Double>();
	private final double[] genotypeDosages;
	private final double[][] genotypeProbs;
	private int[] genotypeQuals;
	private int[] allelicDepths;
	private boolean monomorphic;
	private double callrate;
	private boolean multiallelic;
	private double hwep;
	private boolean biallelic = false;
	private double[] allelefrequencies;
	private String minorAllele;

	private String[] alleles = null;
	private String chr = null;
	private int pos = -1;
	private String id = null;
	private int qual = -1;
	private String filter = null;
	private double MAF;
	private String separator = "/";

	private static final Pattern slash = Pattern.compile("/");
	private static final Pattern pipe = Pattern.compile("\\|");


	public VCFVariant(String ln, boolean skipLoadingGenotypes, boolean skipsplittinggenotypes) {
		this(ln, 0, 0, skipLoadingGenotypes, skipsplittinggenotypes);
	}

	public VCFVariant(String ln, boolean skipLoadingGenotypes) {
		this(ln, 0, 0, skipLoadingGenotypes, false);
	}

	public VCFVariant(String ln) {
		this(ln, 0, 0, false, false);
	}


	public VCFVariant(String ln, int minimalReadDepth, int minimalGenotypeQual, boolean skipLoadingGenotypes, boolean skipsplittinggenotypes) {

		String[] elems = ln.split("\t");

		if (elems.length < 9) {
			System.err.println("ERROR in vcf line: " + ln);
			System.exit(-1);
		}
		// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
		this.chr = new String(elems[0]).intern();
		pos = Integer.parseInt(elems[1]);
		id = new String(elems[2]);

		String ref = elems[3];
		String alt = elems[4];
		String[] alternateAlleles = alt.split(",");
		alleles = new String[1 + alternateAlleles.length];
		alleles[0] = new String(ref).intern();
		for (int i = 0; i < alternateAlleles.length; i++) {
			alleles[1 + i] = new String(alternateAlleles[i]).intern();
		}

		String qualStr = elems[5];
		try {
			qual = Integer.parseInt(qualStr);
		} catch (NumberFormatException e) {

		}

		filter = new String(elems[6]).intern();

		String infoStr = elems[7];
		String[] infoElems = Strings.semicolon.split(infoStr);

		if (!infoStr.equals(".")) {
			for (int e = 0; e < infoElems.length; e++) {
				String[] infoElemElems = Strings.equalssign.split(infoElems[e]);
				String id = new String(infoElemElems[0]).intern();
				if (infoElemElems.length > 1) {
					try {
						Double val = Double.parseDouble(infoElemElems[1]);
						info.put(id, val);
					} catch (NumberFormatException ex) {

					}
				} else {
					if (infoElems[e].equals("DB")) {
						info.put(id, 1d);
					} else {
						System.out.println("info: " + infoElems[e] + " not splitable");
					}

				}
			}
		}


		String[] format = Strings.colon.split(elems[8]);
		int gtCol = -1; // genotype

		int adCol = -1; // Allelic depths for the ref and alt alleles in the order listed
		int dpCol = -1; // Approximate read depth (reads with MQ=255 or with bad mates are filtered)
		int gqCol = -1; // Genotype Quality
		int plCol = -1; // Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification
		int pidCol = -1; // ?
		int pgtCol = -1; // ?
		int dsCol = -1;
		int gpCol = -1;
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
			} else if (format[c].equals("DS")) {
				dsCol = c;
			} else if (format[c].equals("GP")) {
				gpCol = c;
			}
			// GT:AD:DP:GQ:PGT:PID:PL
		}

		if (gtCol == -1) {
			System.out.println("No GT COL: " + elems[8]);
			System.exit(-1);
		}

		int nrCalled = 0;

		// count nr possible alleles
		nrAllelesObserved = new int[alternateAlleles.length + 1];

		byte[][] tmpgenotypes = new byte[elems.length - nrHeaderElems][2];

		allelicDepths = new int[elems.length - nrHeaderElems];
		genotypeQuals = new int[elems.length - nrHeaderElems];
		// count number of alleles with certain read depth

		if (dsCol != -1) {
			genotypeDosages = new double[elems.length - nrHeaderElems];
		} else {
			genotypeDosages = null;
		}

		if (gpCol != -1) {
			if (alleles.length == 3) {
				genotypeProbs = new double[elems.length - nrHeaderElems][6]; // assume diploid for now
			} else {
				genotypeProbs = new double[elems.length - nrHeaderElems][3]; // assume diploid for now
			}
		} else {
			genotypeProbs = null;
		}

		for (int i = nrHeaderElems; i < elems.length; i++) {
			int indPos = i - nrHeaderElems;
			String sampleColumn = elems[i];
			if (sampleColumn.equals("./.")) {
				// not called
				tmpgenotypes[indPos][0] = -1;
				tmpgenotypes[indPos][1] = -1;
			} else {
				String[] sampleElems = Strings.colon.split(sampleColumn);

				if (dsCol != -1) {
					try {
						genotypeDosages[indPos] = Double.parseDouble(sampleElems[dsCol]);
					} catch (NumberFormatException e) {

					}
				}
				if (gpCol != -1) {
					String[] gpElems = Strings.comma.split(sampleElems[gpCol]);

					try {

						if (gpElems.length > genotypeProbs[indPos].length) {
							// System.err.println("More genotype probs than expected: " + sampleElems[gpCol] + " for variant " + chr + ":" + pos);

						} else {

							for (int g = 0; g < gpElems.length; g++) {
								genotypeProbs[indPos][g] = Double.parseDouble(gpElems[g]);
							}
						}

					} catch (NumberFormatException e) {

					}
				}
				if (gtCol != -1) {
					String gt = sampleElems[gtCol];
					String[] gtElems = slash.split(gt);
					separator = "/";
					if (gtElems.length == 1) { // phased genotypes
						gtElems = pipe.split(gt);
						separator = "|";
					}

					int depth = 0;
					int gq = 0;

					if (Strings.dot.equals(gtElems[0]) || Strings.dot.equals(gtElems[1])) {
						// not called
						tmpgenotypes[indPos][0] = -1;
						tmpgenotypes[indPos][1] = -1;
					} else {
						try {
							if (gqCol != -1) {
								gq = Integer.parseInt(sampleElems[gqCol]);
							}
						} catch (NumberFormatException e) {
						}
						genotypeQuals[indPos] = gq;

						try {
							if (dpCol != -1) {
								depth = Integer.parseInt(sampleElems[dpCol]);
							}
						} catch (NumberFormatException e) {
						}
						allelicDepths[indPos] = depth;

						byte gt1 = 0;
						byte gt2 = 0;
						try {
							gt1 = Byte.parseByte(gtElems[0]);
							gt2 = Byte.parseByte(gtElems[1]);

							if (dpCol == -1) {
								// depth variable unavailable
								nrCalled++;

								if (gt1 >= nrAllelesObserved.length || gt2 >= nrAllelesObserved.length) {
									System.err.println("Found more alleles than expected: " + gt1 + "/" + gt2 + " " + Strings.concat(alleles, Strings.forwardslash) + " for variant " + chr + ":" + pos + ":" + id);
									System.exit(-1);
								}
								nrAllelesObserved[gt1]++;
								nrAllelesObserved[gt2]++;

								tmpgenotypes[indPos][0] = gt1;
								tmpgenotypes[indPos][1] = gt2;
							} else {
								// called
								// check coverage
								if (depth < minimalReadDepth) {
									// not called
									tmpgenotypes[indPos][0] = -1;
									tmpgenotypes[indPos][1] = -1;
								} else {
									// called
									// check the genotype quality
									if (gq < minimalGenotypeQual) {
										// not called
										tmpgenotypes[indPos][0] = -1;
										tmpgenotypes[indPos][1] = -1;
									} else {
										// called
//								int gtqual0 = 0;
//								String[] plElems = sampleElems[plCol].split(",");
//								try {
//
//									int genotype = gt1 + gt2;
//									// 0,99,1238 AA AB BB
//									if (gt1 != gt2) {
//										// heterozygote
//									}
//									gtqual0 = Integer.parseInt(plElems[genotype]);
//
//								} catch (NumberFormatException e) {
//								}
										// called.. this is the final genotype
										nrCalled++;
										nrAllelesObserved[gt1]++;
										nrAllelesObserved[gt2]++;

										tmpgenotypes[indPos][0] = gt1;
										tmpgenotypes[indPos][1] = gt2;
									}
								}
							}

						} catch (NumberFormatException e) {
							System.out.println("Cannot parse genotype string: " + sampleElems[gtCol] + " nr elems: " + gtElems.length);
							tmpgenotypes[indPos][0] = -1;
							tmpgenotypes[indPos][1] = -1;
						}

					}
				}
			}
		}

		if (!skipLoadingGenotypes) {
			genotypeAlleles = tmpgenotypes;
		} else {
			genotypeAlleles = null;
		}

		callrate = (double) nrCalled / (elems.length - nrHeaderElems);

		int totalAllelesObs = nrCalled * 2;

		int nrAllelesThatHaveAlleleFrequency = 0;
		double minAlleleFreq = 2;
		allelefrequencies = new double[nrAllelesObserved.length];
		minorAllele = null;

		for (int i = 0; i < nrAllelesObserved.length; i++) {
			double alleleFreq = (double) nrAllelesObserved[i] / totalAllelesObs;
			allelefrequencies[i] = alleleFreq;

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
		}

		MAF = minAlleleFreq;
		if (MAF == 1) {
			MAF = 0;
			if (minorAllele.equals(ref)) {
				minorAllele = alt;
			} else {
				minorAllele = ref;
			}
		}

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

	public double[] getAllelefrequencies() {
		return allelefrequencies;
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

	public int[] getGenotypeQuals() {
		return genotypeQuals;
	}

	public int[] getAllelicDepths() {
		return allelicDepths;
	}

	public double getMAF() {
		return MAF;
	}

	public int[] getNrAllelesObserved() {
		return nrAllelesObserved;
	}

	public byte[][] getGenotypeAlleles() {
		return genotypeAlleles;
	}

	public void flipReferenceAlelele() {

		String allele0 = alleles[0];
		String allele1 = alleles[1];
		alleles[0] = allele1;
		alleles[1] = allele0;

		if (genotypeAlleles != null) {

			// only works for biallelic variants!
			for (int i = 0; i < genotypeAlleles.length; i++) {
				for (int j = 0; j < genotypeAlleles[i].length; j++) {
					if (genotypeAlleles[i][j] == -1) {
						genotypeAlleles[i][j] = -1;
					} else {
						genotypeAlleles[i][j] = (byte) Math.abs(genotypeAlleles[i][j] - 1);
					}
				}
			}


		}
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
			for (int i = 0; i < genotypeAlleles.length; i++) {
				for (int j = 0; j < genotypeAlleles[i].length; j++) {
					if (genotypeAlleles[i][j] == -1) {
						genotypeAlleles[i][j] = -1;
					} else {
						if (alleleRecode[genotypeAlleles[i][j]] == -1) {
							System.err.println("Allele " + alleleRecode[genotypeAlleles[i][j]] + " removed!");
							System.exit(-1);
						}
						genotypeAlleles[i][j] = (byte) alleleRecode[genotypeAlleles[i][j]];
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
		builder.append("\t.\t.\t.\tGT");


		for (int i = 0; i < genotypeAlleles.length; i++) {
			String al1 = "" + genotypeAlleles[i][0];
			String al2 = "" + genotypeAlleles[i][1];
			if (genotypeAlleles[i][0] == -1) {
				al1 = ".";
			}
			if (genotypeAlleles[i][1] == -1) {
				al2 = ".";
			}
			builder.append("\t");
			builder.append(al1).append(separator).append(al2);

		}
		return builder.toString();
	}

	public double[] getGenotypeDosages() {
		return genotypeDosages;
	}

	public double[][] getGenotypeProbs() {
		return genotypeProbs;
	}

	public HashMap<String, Double> getInfo() {
		return info;
	}

	public void recalculateMAFAndCallRate() {


		int nrCalled = 0;
		nrAllelesObserved = new int[nrAllelesObserved.length];
		for (int i = 0; i < genotypeAlleles.length; i++) {
			if (genotypeAlleles[i][0] != -1) {
				nrCalled++;
				nrAllelesObserved[genotypeAlleles[i][0]]++;
				nrAllelesObserved[genotypeAlleles[i][1]]++;
			}
		}

		callrate = (double) nrCalled / (genotypeAlleles.length);

		int totalAllelesObs = nrCalled * 2;

		int nrAllelesThatHaveAlleleFrequency = 0;
		double minAlleleFreq = 2;
		allelefrequencies = new double[nrAllelesObserved.length];
		minorAllele = null;

		for (int i = 0; i < nrAllelesObserved.length; i++) {
			double alleleFreq = (double) nrAllelesObserved[i] / totalAllelesObs;
			allelefrequencies[i] = alleleFreq;

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

			// TODO: calculate HWE P
			hwep = 0;
		} else if (nrAllelesThatHaveAlleleFrequency > 2) {
			multiallelic = true;
		} else {
			monomorphic = true;
		}
	}
}
