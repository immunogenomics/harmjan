package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.vcf.filter.genotypefilters.AllelicDepthFilter;
import nl.harmjanwestra.utilities.vcf.filter.genotypefilters.GenotypeQualityFilter;
import nl.harmjanwestra.utilities.vcf.filter.genotypefilters.VCFGenotypeFilter;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by Harm-Jan on 02/04/16.
 */
public class VCFFilter {


	public void filter(String in, String out, double mafthreshold, double callratethreshold, Integer readdepth, Integer gqual, Double allelicBalance, boolean onlyAutosomes) throws IOException {

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

		ArrayList<VCFGenotypeFilter> filters = new ArrayList<>();
		if (gqual != null) {
			filters.add(new GenotypeQualityFilter(gqual));
			System.out.println("Adding genotyping filter: " + gqual);
		}

		if (allelicBalance != null) {
			if (readdepth == null) {
				filters.add(new AllelicDepthFilter(allelicBalance, 0));
				System.out.println("Adding allele balance filter: ab" + allelicBalance + "\tread depth: " + 0);
			} else {
				filters.add(new AllelicDepthFilter(allelicBalance, readdepth));
				System.out.println("Adding allele balance filter: ab" + allelicBalance + "\tread depth: " + readdepth);
			}
		} else if (readdepth != null) {
			System.out.println("Adding read depth filter: " + readdepth);
		}


		while (ln != null) {
			if (ln.startsWith("##")) {
				tf2.writeln(ln);
			} else if (ln.startsWith("#CHROM")) {
				tf2.writeln("##VCFGenotypeFilter=maf>" + mafthreshold + ",cr>" + callratethreshold + ",dp>" + readdepth + ",gq>" + gqual + ",ab>" + allelicBalance);
				tf2.writeln(ln);
			} else {
				VCFVariant varFiltered = new VCFVariant(ln, filters, true);
				VCFVariant varnofilter = new VCFVariant(ln, true);

				if ((onlyAutosomes && !varFiltered.getChr().equals(Chromosome.X) && !varFiltered.getChr().equals(Chromosome.Y)) || !onlyAutosomes) {
					double m = varFiltered.getMAF() * 2;
					double c = varFiltered.getCallrate();
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

					short[] rd = varFiltered.getApproximateDepth();
					if(rd!=null) {
						for (int i = 0; i < rd.length; i++) {
							meanDepth += rd[i];
							if (rd[i] >= distlen) {
								rdf[distlen - 1]++;
							} else {
								rdf[rd[i]]++;
							}
						}
						meanDepth /= rd.length;
					}

					short[] gq = varFiltered.getGenotypeQuals();
					if(gq!=null) {
						for (int i = 0; i < gq.length; i++) {
							meanQual += gq[i];

							if (gq[i] >= crf.length) {
								gqf[gqf.length - 1]++;
							} else {
								gqf[gq[i]]++;
							}
						}
						meanQual /= gq.length;
					}



					if (varFiltered.getCallrate() >= callratethreshold && varFiltered.getMAF() >= mafthreshold) {
						String[] elems = ln.split("\t");

						String infoblock = elems[8];
						String[] format = Strings.colon.split(infoblock);

						int gtcol = -1;
						for (int q = 0; q < format.length; q++) {
							if (format[q].equals("GT")) {
								gtcol = q;
							}
						}

						String sep = varFiltered.getSeparator();
						double[][] alleles = varFiltered.getGenotypeAlleles();
						// write back, keep annotations..
						for (int i = 9; i < elems.length; i++) {
							int indid = i - 9;
							String[] subElems = elems[i].split(":");
							if (alleles[0][indid] == -1) {
								subElems[gtcol] = "." + sep + ".";
							} else {
								subElems[gtcol] = alleles[0][indid] + sep + alleles[1][indid];
							}
							elems[i] = Strings.concat(subElems, Strings.colon);
						}
						tf2.writeln(Strings.concat(elems, Strings.tab));
					} else {
						filterlog.writeln(varFiltered.toString()
								+ "\t" + varnofilter.getMAF()
								+ "\t" + varFiltered.getMAF()
								+ "\t" + (varFiltered.getMAF() < mafthreshold)
								+ "\t" + varnofilter.getCallrate()
								+ "\t" + varFiltered.getCallrate()
								+ "\t" + (varFiltered.getCallrate() < callratethreshold)
								+ "\t" + meanDepth
								+ "\t" + meanQual);

						if (varFiltered.getCallrate() < callratethreshold) {
							lowcr++;
						}
						if (varFiltered.getMAF() < mafthreshold) {
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

					if (read % 1000 == 0) {
						System.out.print(filtered + " filtered out / " + read + " total filtered\r");

					}
				}
			}
			ln = tf1.readLine();
		}
		System.out.println();
		System.out.println(filtered + " variants filtered out of: " + read + " (" + ((double) filtered / read) + "). " + (read - filtered) + " remain.");
		System.out.println(lowcr + " low call rate " + ((double) lowcr / filtered));
		System.out.println(lowmaf + " low mafthreshold " + ((double) lowmaf / filtered));
		System.out.println(lowcrnf + " low call rate wo filtering " + ((double) lowcrnf / filtered));
		System.out.println(lowmafnf + " low mafthreshold  wo filtering " + ((double) lowmafnf / filtered));
		tf1.close();
		tf2.close();

		System.out.println();
		System.out.println("Bin\tCallRate\tMaf\tCallRateNonFiltered\tMafNonFiltered\tGQ\tReadDepth");
		for (int i = 0; i < crf.length; i++) {
			System.out.println(i + "\t" + crf[i] + "\t" + mff[i] + "\t" + crfnf[i] + "\t" + mffnf[i] + "\t" + gqf[i] + "\t" + rdf[i]);
		}
		filterlog.close();
	}

	/*

		// Maybe this filtering should be moved to a separate class or
		// something?
		if (minimalReadDepth > 0 && dpCol != -1) {

		}
		if (minimalGenotypeQual > 0 && gqCol != -1) {

		}
	 */

	public void filterNonACTGVariants(String vcfin, String vcfout) throws IOException {
		System.out.println("in: " + vcfin);
		System.out.println("out: " + vcfout);
		TextFile out = new TextFile(vcfout, TextFile.W);
		TextFile tf = new TextFile(vcfin, TextFile.R);

		int written = 0;
		int total = 0;
		String ln = tf.readLine();
		while (ln != null) {
			if (ln.startsWith("#")) {
				out.writeln(ln);
			} else {
				VCFVariant var = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
				String[] alleles = var.getAlleles();
				boolean allelesok = true;
				for (String allele : alleles) {
					if (allele.toUpperCase().contains("A") ||
							allele.toUpperCase().contains("C") ||
							allele.toUpperCase().contains("T") ||
							allele.toUpperCase().contains("G")) {
					} else {
						allelesok = false;
					}
				}
				if (allelesok) {
					out.writeln(ln);
					written++;
				} else {
					System.out.println("Weird alleles: " + var.getId() + "\t" + Strings.concat(alleles, Strings.forwardslash));
				}
				total++;
			}
			ln = tf.readLine();
		}

		System.out.println("weird alleles: " + written + "/" + total + " written..");
		tf.close();
		out.close();
	}
}
