package nl.harmjanwestra.gwas;

import cern.colt.matrix.tbit.BitVector;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.text.Strings;

import java.util.ArrayList;

/**
 * Created by hwestra on 11/3/16.
 */
public class HaplotypeGroup {

	private final DoubleMatrix2D haplotypeMatrixCollapsed;
	private final ArrayList<BitVector> group0;
	private final ArrayList<BitVector> group1;
	private final ArrayList<VCFVariant> variants;
	private final DiseaseStatus[][] diseaseStatus;
	double freqCases;
	double freqControls;
	int minorGroup = 0;

	public HaplotypeGroup(DoubleMatrix2D haplotypeMatrixCollapsed,
						  ArrayList<BitVector> group0,
						  ArrayList<BitVector> group1,
						  ArrayList<VCFVariant> variants,
						  DiseaseStatus[][] diseaseStatuses) {
		this.haplotypeMatrixCollapsed = haplotypeMatrixCollapsed;
		this.group0 = group0;
		this.group1 = group1;
		this.variants = variants;
		this.diseaseStatus = diseaseStatuses;
		calcFreq();
	}

	public double getFreqCases() {
		return freqCases;
	}

	public double getFreqControls() {
		return freqControls;
	}

	private void calcFreq() {
		int chrCases = 0;
		int chrControls = 0;
		double allelesCases = 0;
		double allelesControls = 0;
//		TextFile outf = null;
//		try {
//			outf = new TextFile("/Data/tmp/tnfaip3/hap.txt", TextFile.W);
			for (int i = 0; i < haplotypeMatrixCollapsed.size(); i++) {
				double hap = haplotypeMatrixCollapsed.get(i, 0);
				if (hap == 2) {
					if (diseaseStatus[i][0].equals(DiseaseStatus.CASE)) {
						chrCases += 2;
						allelesCases += 2;
					} else {
						chrControls += 2;
						allelesControls += 2;
					}
				} else if (hap == 1) {
					if (diseaseStatus[i][0].equals(DiseaseStatus.CASE)) {
						chrCases += 2;
						allelesCases++;
					} else {
						chrControls += 2;
						allelesControls++;
					}
				} else {
					if (diseaseStatus[i][0].equals(DiseaseStatus.CASE)) {
						chrCases += 2;
					} else {
						chrControls += 2;
					}
				}
//				outf.writeln(i + "\t" + hap);
			}
//			outf.close();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		freqCases = allelesCases / chrCases;

		freqControls = allelesControls / chrControls;
		if (freqControls < 0.5) {
			minorGroup = 1;
		}
		//System.exit(-1);
	}


	public String[] getAlleles() {
		return new String[]{getHaplotypeDesc(group0), getHaplotypeDesc(group1)};
	}

	public String getHaplotypeDesc(ArrayList<BitVector> group) {
		String groupdesc = "";

		for (int a = 0; a < group.size(); a++) {
			BitVector bitVector = group.get(a);
			String hapStr = "";
			for (int b = 0; b < bitVector.size(); b++) {
				String allele = "";
				if (bitVector.get(b)) {
					allele = variants.get(b).getAlleles()[1];
				} else {
					allele = variants.get(b).getAlleles()[0];
				}
				if (hapStr.length() == 0) {
					hapStr += allele;
				} else {
					hapStr += "|" + allele;
				}
			}

			if (groupdesc.length() == 0) {
				groupdesc += hapStr;
			} else {
				groupdesc += "+" + hapStr;
			}
		}

		return groupdesc;
	}

	public DoubleMatrix2D getHaplotypeMatrixCollapsed() {
		return haplotypeMatrixCollapsed;
	}


	public SNPFeature asFeature() {
		Chromosome chr = variants.get(0).getChrObj();
		int start = variants.get(0).getPos();
		int stop = variants.get(variants.size() - 1).getPos();
		SNPFeature snp = new SNPFeature(chr, start, stop);
		snp.setAFCases(freqCases);
		snp.setAFControls(freqControls);
		if (freqControls > 0.5) {
			snp.setMaf(1 - freqControls);
		} else {
			snp.setMaf(freqControls);
		}

		snp.setName(getVariantsAsString());
		snp.setAlleles(getAlleles());
		if (minorGroup == 0) {
			snp.setMinorAllele(getHaplotypeDesc(group0));
		} else {
			snp.setMinorAllele(getHaplotypeDesc(group1));
		}

		return snp;
	}

	public String getVariantsAsString() {
		String variantStr = "";
		for (VCFVariant variant : variants) {
			if (variantStr.length() == 0) {
				variantStr += variant.getId();
			} else {
				variantStr += "," + variant.getId();
			}
		}

		return variantStr;
	}

	@Override
	public String toString() {

		String str = getVariantsAsString() + "\t" + Strings.concat(getAlleles(), Strings.forwardslash) + "\t" + freqControls;
		return str;
	}

	public ArrayList<VCFVariant> getVariants() {
		return variants;
	}

	public ArrayList<BitVector> getGroup0() {
		return group0;
	}

	public ArrayList<BitVector> getGroup1() {
		return group1;
	}
}
