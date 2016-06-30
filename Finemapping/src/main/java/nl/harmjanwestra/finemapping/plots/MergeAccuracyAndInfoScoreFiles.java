package nl.harmjanwestra.finemapping.plots;

import nl.harmjanwestra.finemapping.VariantCounter;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 6/30/16.
 */
public class MergeAccuracyAndInfoScoreFiles {


	public static void main(String[] args) {

		String seqpanelvcf = "/Data/tmp/2016-05-28/seqpanelfiltered-maf0005-cr0950-rd10-gq30-runNamesFixed-RASampleNamesFixed-badSamplesRemoved-mixupsFixed.vcf.gz";
		String variantsOnIC = "/Data/tmp/2016-05-20/T1D-recode-stats.vcf.gz";

		String[] labels = new String[]{
				"EUR", "COSMO", "HRC"
		};
		String[] accfiles = new String[]{
				"/Data/tmp/2016-06-29-quals/Acc/T1D-EUR.txt",
				"/Data/tmp/2016-06-29-quals/Acc/T1D-COSMO.txt",
				"/Data/tmp/2016-06-29-quals/Acc/T1D-HRC-w100kb.txt"
		};
		String[] infofiles = new String[]{
				"/Data/tmp/2016-06-29-quals/INFO/T1D-Beagle1kg-regionfiltered-EUR-ImpQualsReplaced-stats.vcf.gz",
				"/Data/tmp/2016-06-29-quals/INFO/T1D-Beagle1kg-regionfiltered-COSMO-ImpQualsReplaced-stats.vcf.gz",
				"/Data/tmp/2016-06-29-quals/INFO/T1D-HRC-w100kb.vcf.gz"
		};

		String outfile = "/Data/tmp/2016-06-29-quals/T1D-INFOAndAccMerged.txt";

		try {
			MergeAccuracyAndInfoScoreFiles r = new MergeAccuracyAndInfoScoreFiles();
			r.merge(seqpanelvcf, variantsOnIC, accfiles, infofiles, labels, outfile);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void merge(String seqpanelvcf, String variantsOnIC, String[] accfiles, String[] infofiles, String[] labels, String outfile) throws IOException {

		// get a list of maf > 0.005 variants that on the sequencingpanel
		TextFile tf = new TextFile(seqpanelvcf, TextFile.R);
		ArrayList<VCFVariant> seqpanel = new ArrayList<>();
		String ln = tf.readLine();

		double mafthreshold = 0.01;
		double upperthreshold = 1;
//		double infothreshold = 0.5;
		boolean includeICVariants = false;

		boolean includeId = true;
		boolean includeIndels = true;

		VariantCounter counter = new VariantCounter();
		HashSet<String> variantsOnICHash = counter.loadVariantHash(variantsOnIC, includeId);
		System.out.println(variantsOnICHash.size() + " total on IC");
		ArrayList<VCFVariant> variantsNotOnImmunoChip = new ArrayList<>();
		ArrayList<VCFVariant> variantsOnImmunoChip = new ArrayList<>();
		while (ln != null) {
			if (!ln.startsWith("#")) {
				String[] elems = ln.split("\t");
				Chromosome chr = Chromosome.parseChr(elems[0]);
				if (chr.isAutosome()) {
					VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.ALL);
					if (variant.getMAF() > mafthreshold && variant.getMAF() < upperthreshold) {
						seqpanel.add(variant);
						boolean varOnIc = counter.isVariantOnIC(elems, variantsOnICHash, includeId);
						boolean indel = counter.isIndel(elems);

						if (includeIndels || (!includeIndels && !indel)) {
							if (!varOnIc) {
								variantsNotOnImmunoChip.add(variant);
							} else {
								variantsOnImmunoChip.add(variant);
							}
						}
					}
				}
			}
			ln = tf.readLine();
		}
		tf.close();

		System.out.println(seqpanel.size() + " variants in VCF");
		System.out.println(variantsNotOnImmunoChip.size() + " not on IC");
		System.out.println(variantsOnImmunoChip.size() + " on IC");
		System.out.println("MAF> " + mafthreshold);
//		System.out.println("INFO> " + infothreshold);


		HashSet<String> sequencedVariantsHash = null;
		ArrayList<VCFVariant> allvars = null;
		if (includeICVariants) {
			allvars = new ArrayList<>();
			allvars.addAll(variantsNotOnImmunoChip);
			allvars.addAll(variantsOnImmunoChip);
			sequencedVariantsHash = counter.hashit(allvars, includeId);
		} else {
			allvars = new ArrayList<>();
			allvars.addAll(variantsNotOnImmunoChip);
			sequencedVariantsHash = counter.hashit(variantsNotOnImmunoChip, includeId);
		}

		System.out.println(allvars.size() + " variants to evaluate");

		HashMap<String, Integer> variantToId = new HashMap<String, Integer>();
		for (int v = 0; v < allvars.size(); v++) {
			ln = allvars.get(v).toString();
			variantToId.put(ln, v);
		}

		double[][][] scores = new double[infofiles.length][sequencedVariantsHash.size()][2]; // [referencepanels][variants][info/acc]

		for (int q = 0; q < accfiles.length; q++) {
			TextFile tf2 = new TextFile(accfiles[q], TextFile.R);
			String[] elems = tf2.readLineElems(TextFile.tab);
			while (elems != null) {
				String variantName = elems[0];
				Integer varId = variantToId.get(variantName);
				if (varId != null) {
					scores[q][varId][0] = Double.parseDouble(elems[12]);
				}
				elems = tf2.readLineElems(TextFile.tab);
			}
			tf2.close();
		}

		for (int q = 0; q < infofiles.length; q++) {
			TextFile tf2 = new TextFile(infofiles[q], TextFile.R);
			System.out.println("parsing: " + accfiles[q]);
			ln = tf2.readLine();
			while (ln != null) {
				if (!ln.startsWith("#")) {
					VCFVariant var = new VCFVariant(ln);
					String varname = var.toString();
					Integer varId = variantToId.get(varname);
					if (varId != null) {
						scores[q][varId][1] = var.getImputationQualityScore();
					}
				}
				ln = tf2.readLine();
			}
			tf.close();
		}

		TextFile out = new TextFile(outfile, TextFile.W);

		String header = "variant\tmafInSeqPanel";
		for (int l = 0; l < labels.length; l++) {
			header += "\t" + labels[l] + "-Accuracy\t" + labels[l] + "-INFO";
		}
		out.writeln(header);

		for (int v = 0; v < allvars.size(); v++) {
			ln = allvars.get(v).toString() + "\t" + allvars.get(v).getMAF();
			for (int l = 0; l < labels.length; l++) {
				ln += "\t" + scores[l][v][0] + "\t" + scores[l][v][1];
			}
			out.writeln(ln);
		}
		out.close();
	}

}
