package nl.harmjanwestra.gwas;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.FisherExactTest;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by Harm-Jan on 06/01/16.
 */
public class MafHWEPTest {

//	public static void main(String[] args) {
//		String vcf = "D:\\tmp\\2016-06-01\\RA-controls-chr19.vcf.gz";
//		String outfilename = "D:\\tmp\\2016-06-01\\hwepout.txt";
//		MafHWEPTest t = new MafHWEPTest();
//		try {
//			t.run(vcf, outfilename);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//	}
//
//	public void run(String vcffilename, String outfilename) throws IOException {
//		VCFGenotypeData data = new VCFGenotypeData(vcffilename);
//		ArrayList<String> samples = data.getSamples();
//		data.close();
//		TextFile vcf = new TextFile(vcffilename, TextFile.R);
//
//
//
//		DiseaseStatus[] finalDiseaseStatus = new DiseaseStatus[samples.size()];
//		boolean[] includeGenotype = new boolean[samples.size()];
//		for (int i = 0; i < finalDiseaseStatus.length; i++) {
//			finalDiseaseStatus[i] = DiseaseStatus.CONTROL;
//			includeGenotype[i] = true;
//		}
//
//		String ln = vcf.readLine();
//		LRTestTask tasktmp = new LRTestTask();
//		TextFile out = new TextFile(outfilename, TextFile.W);
//		out.writeln("id\tmissing\tmaf\thwep");
//		int ctr = 0;
//		while (ln != null) {
//			if (!ln.startsWith("#") && ln.contains("rs200458265")) {
//				VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.ALL);
//				if (variant.getId().equals("rs200458265")) {
//					System.out.println("Got it!");
//					variant.calculateHWEP();
//				}
//				Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>> unfilteredGenotypeData = tasktmp.filterAndRecodeGenotypes(
//						includeGenotype,
//						variant.getGenotypeAllelesAsMatrix2D(),
//						finalDiseaseStatus,
//						variant.getAlleles().length,
//						samples.size());
//				Triple<Integer, Double, Double> stats = unfilteredGenotypeData.getRight();
//
//				out.writeln(variant.getId() + "\t" + stats.getLeft() + "\t" + stats.getMiddle() + "\t" + stats.getRight() + "\t" + variant.getHwep());
//				ctr++;
//				if (ctr % 1000 == 0) {
//					System.out.println(ctr);
//				}
//			}
//			ln = vcf.readLine();
//		}
//		out.close();
//
//		vcf.close();
//	}
}
