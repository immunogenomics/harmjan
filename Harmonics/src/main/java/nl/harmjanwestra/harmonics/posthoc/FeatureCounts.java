package nl.harmjanwestra.harmonics.posthoc;

import com.itextpdf.text.DocumentException;

import eqtlmappingpipeline.normalization.Normalizer;
import nl.harmjanwestra.utilities.coverage.CoverageMeasures;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.BoxPlotPanel;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

/**
 * Created by hwestra on 9/14/15.
 */
public class FeatureCounts {

	public static void main(String[] args) {
		try {
			FeatureCounts f = new FeatureCounts();
			String sampleFileName = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/FileToSampleName.txt";
			String outdir = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/FeatureCounts/top500ExpressedGenes/output/";
			Gpio.createDir(outdir);
			f.mergeAndProcessFeaturesCountSummaryFiles(sampleFileName, 0, 11, outdir + "summary.txt");
			String covariates = outdir + "summary.txt_cov.txt";
			f.mergeAndProcessFeatureCountsFiles(sampleFileName, 0, 10, 12, covariates, outdir);


//			String outfilename = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/FeatureCounts/tpm-qqnorm-l2t.pdf";
//			try {
//				f.customBoxPlot(outfilename,
//						"/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/FeatureCounts/topPeaksNorm/topPeaks_tpm.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.txt.gz");
//
//				outfilename = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/FeatureCounts/tpm-qqnorm-l2t-Z.pdf";
//				f.customBoxPlot(outfilename,
//						"/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/FeatureCounts/topPeaksNorm/topPeaks_tpm.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz");
//			} catch (DocumentException e) {
//				e.printStackTrace();
//			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void customBoxPlot(String outfilename, String file) throws IOException, DocumentException {
		DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(file);
		ds.transposeDataset();
		String[] q = ds.rowObjects.toArray(new String[0]);

		boxPlotTPM(ds.rawData, outfilename, q);

	}

	public double[][] correlate(String matrix, boolean includeZeros, boolean nonparametric, boolean correlateRows) throws IOException {
		umcg.genetica.math.matrix.DoubleMatrixDataset<String, String> ds = new umcg.genetica.math.matrix.DoubleMatrixDataset<String, String>(matrix);

		if (!correlateRows) {
			ds.transposeDataset();
		}

		double[][] correlationMatrix = new double[ds.nrRows][ds.nrRows];
		for (int r1 = 0; r1 < ds.nrRows; r1++) {
			correlationMatrix[r1][r1] = 1;
			double[] xtmp = ds.rawData[r1];

			for (int r2 = r1 + 1; r2 < ds.nrRows; r2++) {
				double[] ytmp = ds.rawData[r2];
				double[] y = null;
				double[] x = null;
				if (includeZeros) {
					y = ytmp;
					x = xtmp;
				} else {
					Pair<double[], double[]> pair = removeZeroes(xtmp, ytmp);
					x = pair.getLeft();
					y = pair.getRight();
				}

				double r = 0;
				if (nonparametric) {
					SpearmansCorrelation corr = new SpearmansCorrelation();
					r = corr.correlation(x, y);
				} else {
					r = JSci.maths.ArrayMath.correlation(x, y);
				}
				correlationMatrix[r1][r2] = r;
				correlationMatrix[r2][r1] = r;
			}
		}
		if (!correlateRows) {
			ds.transposeDataset();
		}
		return correlationMatrix;
	}

	private Pair<double[], double[]> removeZeroes(double[] xtmp, double[] ytmp) {

		int nrzero = 0;
		for (int i = 0; i < xtmp.length; i++) {
			if (xtmp[i] == 0 && ytmp[i] == 0) {
				nrzero++;
			}
		}

		int ctr = 0;
		double[] y = new double[xtmp.length - nrzero];
		double[] x = new double[xtmp.length - nrzero];
		for (int i = 0; i < xtmp.length; i++) {
			if (xtmp[i] == 0 && ytmp[i] == 0) {

			} else {
				x[ctr] = xtmp[i];
				y[ctr] = ytmp[i];
				ctr++;
			}
		}

		return new Pair<double[], double[]>(x, y);
	}

	private void mergeAndProcessFeaturesCountSummaryFiles(String sampleFilename,
														  int samplecol,
														  int filecol,
														  String out) throws IOException {

		PeakMerge pm = new PeakMerge();
		Pair<String[], String[]> samplecombo = pm.loadSamplePairs(sampleFilename, samplecol, filecol);
		String[] samplefiles = samplecombo.getRight();
		String[] sampleNames = samplecombo.getLeft();


//		Status / medpop / srlab / cd4timelinepilot / atac / bwa - aligned / 20141202 - 0 - hrs - sorted - dedup.bam
//		Assigned 0
//		Unassigned_Ambiguity 0
//		Unassigned_MultiMapping 0
//		Unassigned_NoFeatures 13388577
//		Unassigned_Unmapped 3478752
//		Unassigned_MappingQuality 7996091
//		Unassigned_FragmentLength 4219645
//		Unassigned_Chimera 0
//		Unassigned_Secondary 8885624
//		Unassigned_Nonjunction 0
//		Unassigned_Duplicate 25238000

		TextFile outtf = new TextFile(out, TextFile.W);
		TextFile outtf2 = new TextFile(out + "_cov.txt", TextFile.W);

		outtf2.writeln("Sample\tAssigned");
		outtf.writeln("Sample" +
				"\tAssigned" +
				"\tUnassigned_Ambiguity" +
				"\tUnassigned_MultiMapping" +
				"\tUnassigned_NoFeatures" +
				"\tUnassigned_Unmapped" +
				"\tUnassigned_MappingQuality" +
				"\tUnassigned_FragmentLength" +
				"\tUnassigned_Chimera" +
				"\tUnassigned_Secondary" +
				"\tUnassigned_Nonjunction" +
				"\tUnassigned_Duplicate");

		for (int i = 0; i < samplefiles.length; i++) {
			String f = samplefiles[i];

			TextFile sf = new TextFile(f, TextFile.R);
			sf.readLine();
			String[] elems = sf.readLineElems(TextFile.tab);
			String ln = sampleNames[i];

			String[] lnElems = new String[11];

			while (elems != null) {
				if (elems[0].startsWith("Assigned")) {
					lnElems[0] = elems[1];
				} else if (elems[0].startsWith("Unassigned_Ambiguity")) {
					lnElems[1] = elems[1];
				} else if (elems[0].startsWith("Unassigned_MultiMapping")) {
					lnElems[2] = elems[1];
				} else if (elems[0].startsWith("Unassigned_NoFeatures")) {
					lnElems[3] = elems[1];
				} else if (elems[0].startsWith("Unassigned_Unmapped")) {
					lnElems[4] = elems[1];
				} else if (elems[0].startsWith("Unassigned_MappingQuality")) {
					lnElems[5] = elems[1];
				} else if (elems[0].startsWith("Unassigned_FragmentLength")) {
					lnElems[6] = elems[1];
				} else if (elems[0].startsWith("Unassigned_Chimera")) {
					lnElems[7] = elems[1];
				} else if (elems[0].startsWith("Unassigned_Secondary")) {
					lnElems[8] = elems[1];
				} else if (elems[0].startsWith("Unassigned_Nonjunction")) {
					lnElems[9] = elems[1];
				} else if (elems[0].startsWith("Unassigned_Duplicate")) {
					lnElems[10] = elems[1];
				}
				elems = sf.readLineElems(TextFile.tab);
			}
			outtf2.writeln(ln + "\t" + lnElems[0]);
			outtf.writeln(ln + "\t" + Strings.concat(lnElems, Strings.tab));
		}

		outtf2.close();

		outtf.close();


	}

	private void mergeAndProcessFeatureCountsFiles(String sampleFileName,
												   int samplecol,
												   int filecol,
												   int fraglencol,
												   String covariates,
												   String outdir) throws IOException {
		PeakMerge pm = new PeakMerge();
		Pair<String[], String[]> samplecombo = pm.loadSamplePairs(sampleFileName, samplecol, filecol);
		String[] samplefiles = samplecombo.getRight();
		String[] sampleNames = samplecombo.getLeft();

		// Geneid	Chr	Start	End	Strand	Length	/medpop/srlab/cd4timelinepilot/atac/bwa-aligned/20141202-0-hrs-sorted-dedup.bam

		// assume the same data is in each file
		TextFile tf = new TextFile(samplefiles[0], TextFile.R);
		int peakCtr = 0;
		tf.readLine();
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		ArrayList<String> peakNames = new ArrayList<String>();
		while (elems != null) {

			String name = elems[1] + ":" + elems[2] + "-" + elems[3];
			peakNames.add(name);

			elems = tf.readLineElems(TextFile.tab);
		}

		System.out.println(peakNames.size() + " peaks found");


		tf.close();

		double[][] dataMatrix = new double[peakNames.size()][sampleNames.length];
		for (int d = 0; d < samplefiles.length; d++) {
			tf = new TextFile(samplefiles[d], TextFile.R);
			peakCtr = 0;
			tf.readLine();
			tf.readLine();
			elems = tf.readLineElems(TextFile.tab);

			while (elems != null) {
				dataMatrix[peakCtr][d] = Double.parseDouble(elems[6]);
				peakCtr++;
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
		}

		umcg.genetica.math.matrix.DoubleMatrixDataset<String, String> ddm = new umcg.genetica.math.matrix.DoubleMatrixDataset<String, String>();
		List<String> sampleList = Arrays.asList(sampleNames);
		System.out.println("saving to: " + outdir + "raw.txt");
		try {
			ddm.setRawData(dataMatrix);
			ddm.colObjects = sampleList;
			ddm.rowObjects = peakNames;
			ddm.recalculateHashMaps();
			ddm.save(outdir + "raw.txt");
		} catch (Exception e) {
			e.printStackTrace();
		}

		// TPM that thing
		// load fragment lengths
		double[] fraglen = new double[sampleNames.length];
		TextFile sfile = new TextFile(sampleFileName, TextFile.R);
		sfile.readLine();

		int fctr = 0;
		HashSet<String> sampleNamesHash = new HashSet<String>();
		sampleNamesHash.addAll(sampleList);

		String[] selems = sfile.readLineElems(TextFile.tab);
		while (selems != null) {
			String sample = selems[samplecol];
			String file = selems[filecol];
			if (sampleNamesHash.contains(sample)) {

				double d = Double.parseDouble(selems[fraglencol]);
				fraglen[fctr] = d;
				fctr++;
			}
			selems = sfile.readLineElems(TextFile.tab);
		}

		CoverageMeasures cov = new CoverageMeasures();


		Pair<double[][], ArrayList<String>> output = cov.convertToTPM(dataMatrix, peakNames, fraglen);
		ddm = new umcg.genetica.math.matrix.DoubleMatrixDataset<String, String>();

		System.out.println("saving to: " + outdir + "tpm.txt");
		try {
			ddm.setRawData(output.getLeft());
			ddm.colObjects = sampleList;
			ddm.rowObjects = output.getRight();
			ddm.recalculateHashMaps();
			ddm.save(outdir + "tpm.txt");
		} catch (Exception e) {
			e.printStackTrace();
		}

		Normalizer p = new Normalizer();


		boolean orthogonalizecovariates = false;
		String normoutdir = outdir + "/normalized/";
		Gpio.createDir(outdir);
		String probeIncludeList = null;
		String sampleIncludeList = null;
		int nrPCAsOverSamplesToRemove = 10;
		int nrIntermediatePCAsOverSamplesToRemoveToOutput = 5;

		boolean runQQNorm = true;
		boolean runLog2Transform = true;
		boolean runCenterScale = true;
		boolean runMTransform = false;
		boolean runPCA = true;
		boolean adjustCovariates = true;
		if (covariates == null) {
			adjustCovariates = false;
		}

		boolean forceMissingValues = true;
		boolean forceReplacementOfMissingValues = false;
		boolean forceReplacementOfMissingValues2 = false;
		boolean treatZerosAsNulls = true;
		boolean forceNormalDistribution = false;

		p.normalize(outdir + "tpm.txt",
				probeIncludeList,
				sampleIncludeList,
				nrPCAsOverSamplesToRemove,
				nrIntermediatePCAsOverSamplesToRemoveToOutput,
				covariates,
				orthogonalizecovariates,
				normoutdir,
				runQQNorm,
				runLog2Transform,
				runMTransform,
				runCenterScale,
				runPCA,
				adjustCovariates,
				forceMissingValues,
				forceReplacementOfMissingValues,
				forceReplacementOfMissingValues2,
				treatZerosAsNulls,
				forceNormalDistribution);


		ArrayList<String> newRowNames = output.getRight();
		double[][] tpm = output.getLeft();

		ddm.transposeDataset();
		try {
			boxPlotTPM(ddm.rawData, outdir + "tpm-boxplots.pdf", sampleNames);
		} catch (DocumentException e) {
			e.printStackTrace();
		}


//		// log transform
//
//		Log2Transform.log2transform(tpm);
//
//		ddm = new DoubleMatrixDataset<String, String>();
//
//		System.out.println("saving to: " + out);
//		try {
//			ddm.setColObjects(sampleList);
//			ddm.setRowObjects(output.getRight());
//			ddm.setMatrix(tpm);
//			ddm.save(out + "_tpm_logtransform.txt");
//		} catch (Exception e) {
//			e.printStackTrace();
//		}


	}

	private void boxPlotTPM(double[][] data, String outfilename, String[] sampleNames) throws IOException, DocumentException {
		Grid grid = new Grid(1200, 600, 1, 1, 100, 100);
		BoxPlotPanel boxplotpanel = new BoxPlotPanel(1, 1);
		Range r = new Range(0, 0, 0, 300);
		boxplotpanel.setDrawDataPoints(false);
		boxplotpanel.setLabels(sampleNames);
		boxplotpanel.setOutputIQRS(outfilename + "-iqrs.txt");
		boxplotpanel.setRange(r);
		boxplotpanel.setData(data);
		grid.addPanel(boxplotpanel);
		grid.draw(outfilename);
	}

}
