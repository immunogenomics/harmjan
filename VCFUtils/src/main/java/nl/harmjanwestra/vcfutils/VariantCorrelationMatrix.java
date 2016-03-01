package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by Harm-Jan on 02/27/16.
 */
public class VariantCorrelationMatrix {

	public void correlate(String vcf, String bedregion, String out) throws IOException {


		Feature region = new Feature();
		String[] regionelems = bedregion.split(":");
		Chromosome chr = Chromosome.parseChr(regionelems[0]);
		String[] poselems = regionelems[1].split("-");
		int start = Integer.parseInt(poselems[0]);
		int stop = Integer.parseInt(poselems[1]);

		region.setChromosome(chr);
		region.setStart(start);
		region.setStop(stop);

		System.out.println("parsing: " + vcf + " for region: " + region.toString());

		ArrayList<VCFVariant> variants = new ArrayList<VCFVariant>();
		int ctr = 0;
		TextFile vcfin = new TextFile(vcf, TextFile.R);
		String ln = vcfin.readLine();
		while (ln != null) {
			if (!ln.startsWith("#")) {
				VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
				if (region.overlaps(variant.asFeature())) {
					variant = new VCFVariant(ln, VCFVariant.PARSE.GENOTYPES);
					if (variant.getMAF() > 0.005 && variant.getAlleles().length == 2) {
						variants.add(variant);
					}
				}
				ctr++;
				if (ctr % 1000 == 0) {
					System.out.println(ctr + " variants parsed. " + variants.size() + " within region: " + region.toString());
				}
			}

			ln = vcfin.readLine();
		}
		vcfin.close();

		System.out.println(variants.size() + " variants in region " + region.toString());

		double[][] matrix = new double[variants.size()][variants.size()];

		for (int i = 0; i < variants.size(); i++) {
			double[] genotypes1 = convertToDouble(variants.get(i));
			for (int j = i + 1; j < variants.size(); j++) {
				double[] genotypes2 = convertToDouble(variants.get(j));
				matrix[i][j] = Correlation.correlate(genotypes1, genotypes2);
				matrix[j][i] = matrix[i][j];
			}
			matrix[i][i] = 1;
			if (i % 10 == 0) {
				System.out.println(i + " out of " + variants.size() + " variants correlated.");
			}
		}

		savematrix(matrix, out);
		savelist(variants, out);

	}

	private void savelist(ArrayList<VCFVariant> variants, String out) throws IOException {
		TextFile tfout = new TextFile(out + "-variantlist.txt", TextFile.W);
		for (int i = 0; i < variants.size(); i++) {
			tfout.writeln(variants.get(i).toString() + "-" + variants.get(i).getId());
		}
		tfout.close();
	}

	private void savematrix(double[][] matrix, String out) throws IOException {
		TextFile tfout = new TextFile(out, TextFile.W);
		for (int row = 0; row < matrix.length; row++) {
			tfout.writeln(Strings.concat(matrix[row], Strings.tab));
		}
		tfout.close();
	}

	public double[] convertToDouble(VCFVariant vcfVariant) {
		byte[][] alleles = vcfVariant.getGenotypeAlleles();
		double[] output = new double[alleles[0].length];
		for (int i = 0; i < alleles[0].length; i++) {
			if (alleles[0][i] == -1) {
				output[i] = -1;
			} else {
				output[i] = (alleles[0][i] + alleles[1][i]);
			}
		}

		return output;
	}
}
