package nl.harmjanwestra;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.ProbeAnnotation;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.RankArray;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Created by hwestra on 2/9/16.
 */
public class BroRunner {


	public static void main(String[] args) {

		BroRunner b = new BroRunner();
		try {
			boolean removezeros = true;
			String homfile = "/Data/Homozygosity/EGCUT/hg19/roh/st01_03_unpruned.hom.indiv";
			String linkFile = "/Sync/Dropbox/ROH_Expression/Linker/Prote_Vcode_to_Pcode,txt.txt";
			String expfile = "/Volumes/Datasets/Datasets/GeneticalGenomicsDatasets/EGCUT/2014-10-28-NormalizerOutput/ExpressionData.txt.gz.PCAOverSamplesEigenvectorsTransposed.txt.gz";//"/Data/Homozygosity/EGCUT/ExpressionData.QuantileNormalized.Log2Transformed.txt.gz";
			String probeannot = "/Data/GGD/2013-07-18-HT12v3.txt";
			String out1 = "/Data/Homozygosity/EGCUT/hg19/roh/st01_03_unpruned.hom.indiv_EigV_woZeros.txt";

			b.runHomInDiv(homfile, removezeros, linkFile, expfile, probeannot, out1);
//			homfile = "/Data/Homozygosity/GRNG/hg19/roh/st01_03_unpruned.hom.indiv";
//			linkFile = null;// "/Sync/Dropbox/ROH_Expression/Linker/Prote_Vcode_to_Pcode,txt.txt";
//			expfile = "/Data/Homozygosity/GRNG/ExpressionData.txt.QuantileNormalized.Log2Transformed.txt.gz";
//			probeannot = "/Data/GGD/2013-07-18-HT12v3.txt";
//			String out2 = "/Data/Homozygosity/GRNG/hg19/roh/st01_03_unpruned.hom.indiv_QQL2T_woZeros.txt";
//
//
//			b.runHomInDiv(homfile, removezeros, linkFile, expfile, probeannot, out2);
////
//
//			String file1 = out1;//"/Data/Homozygosity/EGCUT/hg19/roh/st01_03_unpruned.hom.indiv_expCorr40PCs_woZeros.txt";
//			int n1 = 856;
//			String file2 = out2; // "/Data/Homozygosity/GRNG/hg19/roh/st01_03_unpruned.hom.indiv_expCorr40PCs_woZeros.txt";
//			int n2 = 1240;
//			String outf = "/Data/Homozygosity/st01_03_unpruned.hom.indiv_QQL2T_woZeros_meta.txt";
//			b.meta(file1, n1, file2, n2, outf);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	public void runHomInDiv(String homindivFile, boolean removezeros, String linkFile, String expFile, String probeAnnot, String out) throws IOException {

		HashMap<String, String> sampleLink = null;
		if (linkFile != null) {
			sampleLink = new HashMap<String, String>();
			TextFile tf = new TextFile(linkFile, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				sampleLink.put(elems[1], elems[0]);
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
		}

		// load homInd
		HomIndivFile homfile = new HomIndivFile(homindivFile, removezeros, sampleLink);

		// expression
		DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(expFile);

		// prune datasets
		String[] homSamples = homfile.getSamples();
		HashMap<String, Integer> newIndex = new HashMap<String, Integer>();

		int sharedCtr = 0;
		for (int i = 0; i < homSamples.length; i++) {
			if (ds.hashCols.containsKey(homSamples[i])) {
				newIndex.put(homSamples[i], sharedCtr);
				sharedCtr++;
			}
		}

		if (sharedCtr < 10) {
			System.err.println("Error: shared samples " + sharedCtr + " < 10");
		} else {

			ProbeAnnotation pb = new ProbeAnnotation(probeAnnot);


			System.out.println(sharedCtr + " samples shared");
			homfile.reOrder(newIndex);
			ds = pruneExpSamples(ds, newIndex);

			double[][] matrix = ds.getRawData();
			List<String> probes = ds.rowObjects;
			double[] homAvg = homfile.getKb();
			double[] homAvgKb = homfile.getAvgkb();
			RankArray rda = new RankArray();
			rda.rank(homAvg, true);
			rda.rank(homAvgKb, true);
			Correlation corr = new Correlation();
			TextFile outfile = new TextFile(out, TextFile.W);
			outfile.writeln("Probe\tChr\tLocation\tGene\tn\tKbR\tAvgKbR\tKbZ\tAvgKbZ\tKbP\tAvgKbP");
			Correlation.correlationToZScore(homAvgKb.length);
			for (int i = 0; i < matrix.length; i++) {
				double[] exp = matrix[i];
				rda.rank(exp, true);
				double r1 = Correlation.correlate(exp, homAvg);
				double r2 = Correlation.correlate(exp, homAvgKb);

				double z1 = Correlation.convertCorrelationToZScore(exp.length, r1);
				double z2 = Correlation.convertCorrelationToZScore(exp.length, r2);

				double p1 = ZScores.zToP(z1);
				double p2 = ZScores.zToP(z2);

				String probe = ds.rowObjects.get(i);
				Integer pid = pb.getProbeToProbeId().get(probe);
				int loc = 0;
				String gene = "";
				short chr = 0;

				if (pid != null) {
					chr = pb.getChr()[pid];
					loc = pb.getChrStart()[pid];
					gene = pb.getProbeAnnotation()[pid];
				}

				outfile.writeln(probe
						+ "\t" + chr
						+ "\t" + loc
						+ "\t" + gene
						+ "\t" + exp.length
						+ "\t" + r1
						+ "\t" + r2
						+ "\t" + z1
						+ "\t" + z2
						+ "\t" + p1
						+ "\t" + p2);


				if (i % 1000 == 0) {
					System.out.println(i + "\t" + ds.rowObjects.size());
				}
			}
			outfile.close();

		}


	}


	public void meta(String in1, int n1, String in2, int n2, String outf) throws IOException {


		ArrayList<String> probes1 = new ArrayList<>();
		ArrayList<Double> vals1 = new ArrayList<>();

		TextFile tf = new TextFile(in1, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, Integer> probeId = new HashMap<>();
		int ctr = 0;
		while (elems != null) {
			String probe = elems[0];
			String val = elems[8];
			probes1.add(probe);

			vals1.add(Double.parseDouble(val));
			probeId.put(probe, ctr);
			elems = tf.readLineElems(TextFile.tab);
			ctr++;
		}
		tf.close();

		TextFile out = new TextFile(outf, TextFile.W);

		tf = new TextFile(in2, TextFile.R);
		tf.readLine();
		elems = tf.readLineElems(TextFile.tab);

		int[] sampleSizes = new int[]{n1, n2};
		out.writeln("probe\tgene\tvals\twz\tp");
		while (elems != null) {
			String probe = elems[0];
			String gene = elems[3];
			String val = elems[8];

			Integer index = probeId.get(probe);
			double val1 = vals1.get(index);
			double[] z = new double[]{val1, Double.parseDouble(val)};
			double wz = ZScores.getWeightedZ(z, sampleSizes);
			double p = ZScores.zToP(wz);

			out.writeln(probe
					+ "\t" + gene
					+ "\t" + Strings.concat(z, Strings.semicolon)
					+ "\t" + wz
					+ "\t" + p);

			elems = tf.readLineElems(TextFile.tab);

		}
		tf.close();

		out.close();


	}

	private DoubleMatrixDataset<String, String> pruneExpSamples(DoubleMatrixDataset<String, String> ds, HashMap<String, Integer> newIndex) {

		double[][] org = ds.getRawData();
		double[][] out = new double[org.length][newIndex.size()];

		List<String> samples = ds.colObjects;
		List<String> newSampls = new ArrayList<>();
		for (int c = 0; c < samples.size(); c++) {
			String sample = samples.get(c);
			Integer idx = newIndex.get(sample);
			if (idx != null) {
				newSampls.add(sample);

				for (int r = 0; r < org.length; r++) {
					out[r][idx] = org[r][c];
				}
			}
		}

		DoubleMatrixDataset<String, String> output = new DoubleMatrixDataset<>();
		output.rawData = out;
		output.rowObjects = ds.rowObjects;
		output.colObjects = newSampls;
		output.recalculateHashMaps();
		return output;
	}

}
