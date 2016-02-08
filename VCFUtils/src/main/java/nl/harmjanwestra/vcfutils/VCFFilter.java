package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by Harm-Jan on 02/04/16.
 */
public class VCFFilter {


	public void filter(String in, String out, double mafthreshold, int readdepth, int gqual, double callratethreshold, boolean onlyAutosomes) throws IOException {

		TextFile tf1 = new TextFile(in, TextFile.R);
		TextFile tf2 = new TextFile(out, TextFile.W);
		String ln = tf1.readLine();


		System.out.println("Filtering: " + in);
		System.out.println("Out: " + out);
		System.out.println("mafthreshold > " + mafthreshold);
		System.out.println("readdepth > " + readdepth);
		System.out.println("gqual > " + gqual);
		System.out.println("callrate > " + callratethreshold);

		TextFile filterlog = new TextFile(out + "_log.txt", TextFile.W);

		int read = 0;
		int filtered = 0;
		int lowcr = 0;
		int lowmaf = 0;
		int lowcrnf = 0;
		int lowmafnf = 0;

		int distlen = 100;
		int[] crf = new int[distlen];
		int[] crfnf = new int[distlen];
		int[] mff = new int[distlen];
		int[] mffnf = new int[distlen];

		int[] gqf = new int[distlen];
		int[] rdf = new int[distlen];

		while (ln != null) {
			if (ln.startsWith("#")) {
				tf2.writeln(ln);
			} else {
				VCFVariant var = new VCFVariant(ln, readdepth, gqual, true);
				VCFVariant varnofilter = new VCFVariant(ln, 0, 0, true);

				if ((onlyAutosomes && !var.getChr().equals(Chromosome.X) && !var.getChr().equals(Chromosome.Y)) || !onlyAutosomes) {
					double m = var.getMAF() * 2;
					double c = var.getCallrate();
					int bin = (int) Math.floor(m * distlen);
					if (bin >= distlen) {
						bin = distlen - 1;
					}
					mff[bin]++;
					bin = (int) Math.floor(c * distlen);
					if (bin >= distlen) {
						bin = distlen - 1;
					}
					crf[bin]++;

					m = varnofilter.getMAF() * 2;
					c = varnofilter.getCallrate();
					bin = (int) Math.floor(m * distlen);
					if (bin >= distlen) {
						bin = distlen - 1;
					}
					mffnf[bin]++;
					bin = (int) Math.floor(c * distlen);
					if (bin >= distlen) {
						bin = distlen - 1;
					}
					crfnf[bin]++;

					double meanDepth = 0;
					double meanQual = 0;

					short[] rd = var.getApproximateDepth();
					for (int i = 0; i < rd.length; i++) {
						meanDepth += rd[i];
						if (rd[i] >= distlen) {
							rdf[distlen - 1]++;
						} else {
							rdf[rd[i]]++;
						}
					}

					short[] gq = var.getGenotypeQuals();
					for (int i = 0; i < gq.length; i++) {
						meanQual += gq[i];

						if (gq[i] >= crf.length) {
							gqf[gqf.length - 1]++;
						} else {
							gqf[gq[i]]++;
						}
					}

					meanDepth /= rd.length;
					meanQual /= gq.length;

					if (var.getCallrate() >= callratethreshold && var.getMAF() >= mafthreshold) {
						tf2.writeln(var.toVCFString());
					} else {


						filterlog.writeln(var.toString()
								+ "\t" + varnofilter.getMAF()
								+ "\t" + var.getMAF()
								+ "\t" + (var.getMAF() < mafthreshold)
								+ "\t" + varnofilter.getCallrate()
								+ "\t" + var.getCallrate()
								+ "\t" + (var.getCallrate() < callratethreshold)
								+ "\t" + meanDepth
								+ "\t" + meanQual);

						if (var.getCallrate() < callratethreshold) {
							lowcr++;
						}
						if (var.getMAF() < mafthreshold) {
							lowmaf++;
						}

						if (varnofilter.getCallrate() < callratethreshold) {
							lowcrnf++;
						}
						if (varnofilter.getMAF() < mafthreshold) {
							lowmafnf++;
						}


						filtered++;
					}
					read++;
				}
			}
			ln = tf1.readLine();
		}
		System.out.println(filtered + " variants filtered out of: " + read + " (" + ((double) filtered / read) + "). " + (read - filtered) + " remain.");
		System.out.println(lowcr + " low call rate " + ((double) lowcr / filtered));
		System.out.println(lowmaf + " low mafthreshold " + ((double) lowmaf / filtered));
		System.out.println(lowcrnf + " low call rate wo filtering " + ((double) lowcrnf / filtered));
		System.out.println(lowmafnf + " low mafthreshold  wo filtering " + ((double) lowmafnf / filtered));
		tf1.close();
		tf2.close();

		System.out.println();
		for (int i = 0; i < crf.length; i++) {
			System.out.println(i + "\t" + crf[i] + "\t" + mff[i] + "\t" + crfnf[i] + "\t" + mffnf[i] + "\t" + gqf[i] + "\t" + rdf[i]);
		}
		filterlog.close();
	}

}
