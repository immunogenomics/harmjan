package nl.harmjanwestra.vcfutils;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.LDPanel;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.vcf.VCFVariantComparator;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Created by Harm-Jan on 02/27/16.
 */
public class VariantLDMatrix {

	public void correlate(String vcf, String bedregion, String out, boolean printheaders, boolean dprime, boolean useCorrelation) throws IOException {
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

		// sort variants by position
		Collections.sort(variants, new VCFVariantComparator());

		System.out.println(variants.size() + " variants in region " + region.toString());

		System.out.println("Sorting variants");
		Collections.sort(variants, new VCFVariantComparator());

		double[][] matrix = new double[variants.size()][variants.size()];
		DetermineLD ldcalc = new DetermineLD();
		for (int i = 0; i < variants.size(); i++) {
			matrix[i][i] = 1d;
			double[] genotypes1 = convertToDouble(variants.get(i));
			for (int j = i + 1; j < variants.size(); j++) {
				double[] genotypes2 = convertToDouble(variants.get(j));

				if (useCorrelation) {
					matrix[i][j] = Correlation.correlate(genotypes1, genotypes2);
				} else {
					Pair<Double, Double> ldInfo = ldcalc.getLD(variants.get(i), variants.get(j));

					if (ldInfo != null) {
						if (dprime) {
							matrix[i][j] = ldInfo.getLeft();
						} else {
							matrix[i][j] = ldInfo.getRight();
						}
					} else {
						matrix[i][j] = 0;
					}
				}
				matrix[j][i] = matrix[i][j];
			}
			matrix[i][i] = 1;
			if (i % 10 == 0) {
				System.out.println(i + " out of " + variants.size() + " variants correlated.");
			}
		}

		savematrix(matrix, variants, printheaders, out);
		savelist(variants, out);

	}

	public void plot(String correlationMatrix, String output) throws IOException {

		Grid grid = new Grid(1000, 1000, 1, 1, 10, 10);


		LDPanel panel = new LDPanel(0, 0);
		DoubleMatrixDataset<String, String> set = new DoubleMatrixDataset<String, String>();
		set.loadDoubleData(correlationMatrix);
		panel.setData(set);
		grid.addPanel(panel);

		try {
			grid.draw(output);
		} catch (DocumentException e) {
			e.printStackTrace();
		}


	}

	private void savelist(ArrayList<VCFVariant> variants, String out) throws IOException {
		TextFile tfout = new TextFile(out + "-variantlist.txt", TextFile.W);
		for (int i = 0; i < variants.size(); i++) {
			tfout.writeln(variants.get(i).toString());
		}
		tfout.close();
	}

	private void savematrix(double[][] matrix, ArrayList<VCFVariant> variants, boolean printheaders, String out) throws IOException {
		TextFile tfout = new TextFile(out, TextFile.W);
		if (printheaders) {
			String header = "-";
			for (int i = 0; i < variants.size(); i++) {
				header += "\t" + variants.get(i).toString();
			}
			tfout.writeln(header);
		}
		for (int row = 0; row < matrix.length; row++) {
			String header = "";
			if (printheaders) {
				header = variants.get(row).toString();
			}
			tfout.writeln(header + "\t" + Strings.concat(matrix[row], Strings.tab));
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
