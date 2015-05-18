package nl.harmjanwestra.ngs.GenotypeFormats;

import umcg.genetica.text.Strings;

import java.util.HashMap;

/**
 * Created by hwestra on 4/21/15.
 */
public class VCFVariant {

	private static final int nrHeaderElems = 9;
	private final int[] nrAllelesObserved;
	private final byte[][] genotypeAlleles;
	private final String info;
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

	public VCFVariant(String ln) {
		this(ln, 0, 0, false);
	}

	public VCFVariant(String ln, int minimalReadDepth, int minimalGenotypeQual, boolean skipLoadingGenotypes) {

		String[] elems = ln.split("\t");

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

		// TODO: parse info elements
		info = elems[7];
		String[] infoElems = info.split(";");

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
		for (int i = nrHeaderElems; i < elems.length; i++) {
			int indPos = i - nrHeaderElems;
			String sampleColumn = elems[i];
			if (sampleColumn.equals("./.")) {
				// not called
				tmpgenotypes[indPos][0] = -1;
				tmpgenotypes[indPos][1] = -1;
			} else {
				String[] sampleElems = sampleColumn.split(":");
				if (gtCol != -1) {
					String gt = sampleElems[gtCol];
					String[] gtElems = gt.split("/");
					separator = "/";
					if (gtElems.length == 1) { // phased genotypes
						gtElems = gt.split("\\|");
						separator = "|";
					}

					int depth = 0;
					int gq = 0;

					if (gtElems[0].equals(".") || gtElems[1].equals(".")) {
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
							System.out.println("Cannot parse genotype string: " + sampleElems[gtCol]);
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

	public String getInfo() {
		return info;
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
		for (int i = 0; i < alleles.length; i++) {
			alleleRecode[i] = alleleMap.get(alleles[i]);
		}


		if (genotypeAlleles != null) {
			for (int i = 0; i < genotypeAlleles.length; i++) {
				for (int j = 0; j < genotypeAlleles[i].length; j++) {
					if (genotypeAlleles[i][j] == -1) {
						genotypeAlleles[i][j] = -1;
					} else {
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
		String output = chr + "\t" + pos + "\t" + id + "\t" + alleles[0] + "\t" + Strings.concat(alleles, Strings.comma, 1, alleles.length) + "\t.\t.\t.\tGT";
		for (int i = 0; i < genotypeAlleles.length; i++) {
			String al1 = "" + genotypeAlleles[i][0];
			String al2 = "" + genotypeAlleles[i][1];
			if (genotypeAlleles[i][0] == -1) {
				al1 = ".";
			}
			if (genotypeAlleles[i][1] == -1) {
				al2 = ".";
			}
			output += "\t" + al1 + separator + al2;
		}
		return output;
	}
}
