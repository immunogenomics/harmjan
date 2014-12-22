package nl.harmjanwestra.utilities.vcf;

/**
 * Created by hwestra on 12/18/14.
 */
public enum VCFVariantType {

	Unknown(false),
	UnknownBiallelic(true),
	MultiAllelicSNP(false),
	MultiAllelicDeletion(false),
	MultiAllelicInsert(false),
	MultiAllelicIndel(false),
	Deletion(true),
	InDel(true),
	BiallelicSNP(true),
	Insertion(true);

	boolean b = false;

	private VCFVariantType(boolean b) {
		this.b = b;
	}

	public static VCFVariantType parseType(String ref, String alt) {

		String[] altelems = alt.split(",");
		boolean altSingleBp = true;
		if (altelems.length > 1) {

			for (String s : altelems) {
				if (s.length() > 1) {
					altSingleBp = false;
				}
			}

			if (altSingleBp && ref.length() == 1) {
				return MultiAllelicSNP;
			} else if (altSingleBp && ref.length() > 1) {
				return MultiAllelicDeletion;
			} else if (!altSingleBp && ref.length() == 1) {
				return MultiAllelicInsert;
			} else if (!altSingleBp && ref.length() > 1) {
				return MultiAllelicIndel;
			} else {
				return Unknown;
			}

		} else {
			if (ref.length() > 1 && alt.length() == 1) {
				return Deletion;
			} else if (ref.length() > 1 && alt.length() > 1) {
				return InDel;
			} else if (ref.length() == 1 && alt.length() == 1) {
				return BiallelicSNP;
			} else if (ref.length() == 1 && alt.length() > 1) {
				return Insertion;
			} else {
				return UnknownBiallelic;
			}
		}


	}

	public boolean isBiallelic() {
		return b;
	}
}

