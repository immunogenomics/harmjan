package nl.harmjanwestra.finemappingtools.gwas;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.finemappingtools.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.vcf.VCFVariantComparator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

/**
 * Created by hwestra on 7/7/17.
 */
public class GUESS extends LRTest {

	public GUESS(LRTestOptions options) throws IOException {
		super(options);
		convert();
	}

	private void convert() throws IOException {
		System.out.println("Testing haplotypes");
		options.setSplitMultiAllelic(true);

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(options.getBedfile());
		ArrayList<VCFVariant> variants = readVariants(options.getOutputdir() + "variantlog.txt", regions);
		Collections.sort(variants, new VCFVariantComparator());

		if (variants.isEmpty()) {
			System.out.println("Sorry no variants found.");
			System.exit(-1);
		}

		for (Feature region : regions) {
			// make output directory
			String dirGuessScript = options.getOutputdir() + "/" + region.toString() + "/";
			String dirOut = options.getOutputdir() + "/" + region.toString() + "/input/";
			String dirGuessOut = options.getOutputdir() + "/" + region.toString() + "/output/";
			Gpio.createDir(dirGuessScript);
			Gpio.createDir(dirOut);
			Gpio.createDir(dirGuessOut);


			// filter variants
			ArrayList<VCFVariant> regionVariants = getRegionVariants(variants, region);

			TextFile xOut = new TextFile(dirOut + "x.txt", TextFile.W);
			// write number of observations
			int nrSamples = sampleAnnotation.getSampleName().length;
			xOut.writeln("" + nrSamples);

			// write number of predictors
			DoubleMatrix2D covariates = sampleAnnotation.getCovariates();
			DiseaseStatus[][] diseaseStatuses = sampleAnnotation.getSampleDiseaseStatus();

			// prune for correlated covariates
			covariates = prune(covariates);

			int totalPredictors = covariates.columns() + regionVariants.size();
			xOut.writeln("" + totalPredictors);

			// one row per sample
			String tab = "\t";
			for (int i = 0; i < nrSamples; i++) {
				// append covariates first
				StringBuilder sampleOut = new StringBuilder();
				for (int j = 0; j < covariates.columns(); j++) {
					if (j == 0) {
						sampleOut.append(covariates.getQuick(i, j));
					} else {
						sampleOut.append(tab).append(covariates.getQuick(i, j));
					}
				}

				// now append genotypes
				for (int j = 0; j < variants.size(); j++) {
					double dosage = variants.get(j).getDosagesAsMatrix2D().getQuick(i, 0);
					if (covariates.columns() == 0) {
						sampleOut.append(dosage);
					} else {
						sampleOut.append(tab).append(dosage);
					}
				}
				sampleOut.append("\n");
				xOut.append(sampleOut);
			}
			xOut.close();

			// write list of variants
			TextFile varOut = new TextFile(dirOut + "variants.txt", TextFile.W);
			for (int i = 0; i < regionVariants.size(); i++) {
				VCFVariant variant = regionVariants.get(i);
				varOut.writeln(variant.getChr().toString()
						+ "\t" + variant.getPos()
						+ "\t" + variant.getId()
						+ "\t" + Strings.concat(variant.getAlleles(), Strings.comma)
						+ "\t" + variant.getMinorAllele());
			}
			varOut.close();

			// write list of phenotypes
			TextFile yOut = new TextFile(dirOut + "y.txt", TextFile.W);
			yOut.writeln("" + nrSamples);
			yOut.writeln("" + 1);

			for (int i = 0; i < diseaseStatuses.length; i++) {
				yOut.writeln("" + diseaseStatuses[i][0].getNumber());
			}
			yOut.close();

			// write GUESS command line
			String guessshell = "GUESS -X " + dirOut + "x.txt -Y " + dirOut + "y.txt -out_full " + dirGuessOut;
			TextFile scriptout = new TextFile(dirGuessScript + "bashguess.sh", TextFile.W);
			scriptout.writeln(guessshell);
			scriptout.close();

		}
	}

	private DoubleMatrix2D prune(DoubleMatrix2D covariates) {

		// run backwards
		HashSet<Integer> colsToRemove = new HashSet<>();
		for (int i = covariates.columns() - 1; i > -1; i--) {
			double[] x = covariates.viewColumn(i).toArray();

			// check if this covariate is correlated with anything
			for (int j = i - 1; j > -1; j--) {
				double[] y = covariates.viewColumn(j).toArray();
				// correlate
				double r = JSci.maths.ArrayMath.correlation(x, y);
				if (Math.abs(r) >= 1d) {
					// remove i
					colsToRemove.add(i);
					break;
				}
			}
		}

		if (!colsToRemove.isEmpty()) {

			int[] rowIndexes = new int[covariates.rows()];
			for (int i = 0; i < rowIndexes.length; i++) {
				rowIndexes[i] = i;
			}

			int[] colIndexes = new int[covariates.columns() - colsToRemove.size()];
			int ctr = 0;
			for (int i = 0; i < covariates.columns(); i++) {
				if (!colsToRemove.contains(i)) {
					colIndexes[ctr] = i;
					ctr++;
				}
			}

			// return a submatrix with the selected columns removed
			return covariates.viewSelection(rowIndexes, colIndexes);

		} else {
			return covariates;
		}
	}

	private ArrayList<VCFVariant> getRegionVariants(ArrayList<VCFVariant> variants, Feature region) {
		ArrayList<VCFVariant> output = new ArrayList<VCFVariant>();
		for (VCFVariant v : variants) {
			if (!v.isMultiallelic() && v.asFeature().overlaps(region)) {
				output.add(v);
			}
		}
		return output;
	}

}
