package nl.harmjanwestra.finemapping.genotypesandimputation;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by Harm-Jan on 07/19/16.
 */
public class MetaAnalysisCovarWriter {

	public static void main(String[] args) {
		MetaAnalysisCovarWriter c = new MetaAnalysisCovarWriter();
//		String remainlistRA = "D:\\tmp\\2016-07-19\\sim\\listofsharedsamples-RA-rewritten-remaining.txt";
//		String remainlistT1D = "D:\\tmp\\2016-07-19\\sim\\listofsharedsamples-T1D-rewritten-remaining.txt";
//		String remainlistMergeRA = "D:\\tmp\\2016-07-19\\sim\\listofsharedsamples-RA-rewritten-remaining-dedup.txt";
//		String remainlistMergeT1D = "D:\\tmp\\2016-07-19\\sim\\listofsharedsamples-T1D-rewritten-remaining-dedup.txt";
//
//		try {
//			c.dedup(remainlistRA, remainlistT1D, remainlistMergeRA, remainlistMergeT1D);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

		String remainlistMergeRA = "D:\\tmp\\2016-08-03-sim\\out\\listofsharedsamples-RA-rewritten-remaining-dedup.txt";
		String remainlistMergeT1D = "D:\\tmp\\2016-08-03-sim\\out\\listofsharedsamples-T1D-rewritten-remaining-dedup.txt";

		String famin = "D:\\tmp\\2016-08-03-sim\\phenotypes\\covarmerged.txtmergedfam.fam";
		String famout = "D:\\tmp\\2016-08-03-sim\\phenotypes\\meta-RA.fam";
		String covin = "D:\\tmp\\2016-08-03-sim\\phenotypes\\covarmerged.txtmergedCovariates.txt";
		String covout = "D:\\tmp\\2016-08-03-sim\\phenotypes\\meta-RA-covar.txt";
		String diseasein = "D:\\tmp\\2016-08-03-sim\\phenotypes\\covarmerged.txtmergeddisease.txt";
		String diseaseout = "D:\\tmp\\2016-08-03-sim\\phenotypes\\meta-RA-disease.txt";
		try {
			c.filterfam(remainlistMergeRA, famin, famout);
			c.filtercov(remainlistMergeRA, covin, covout, true);
			c.filtercov(remainlistMergeRA, diseasein, diseaseout, false);
		} catch (IOException e) {
			e.printStackTrace();
		}

		System.out.println();
		System.out.println("T1D");


		famin = "D:\\tmp\\2016-08-03-sim\\phenotypes\\T1D-recode-maf0005-ICRegions-samplenamefix-pseudo.vcf.gz-filtered-merged.fam";
		famout = "D:\\tmp\\2016-08-03-sim\\phenotypes\\meta-T1D.fam";
		covin = "D:\\tmp\\2016-08-03-sim\\phenotypes\\2016-03-11-T1D-covarmerged.txtmergedCovariates-withPseudos.txt";
		covout = "D:\\tmp\\2016-08-03-sim\\phenotypes\\meta-T1D-covar.txt";
		diseasein = "D:\\tmp\\2016-08-03-sim\\phenotypes\\2016-03-11-T1D-diseaseStatusWithPseudos.txt";
		diseaseout = "D:\\tmp\\2016-08-03-sim\\phenotypes\\meta-T1D-disease.txt";
		try {
			c.filterfam(remainlistMergeT1D, famin, famout);
			c.filtercov(remainlistMergeT1D, covin, covout, true);
			c.filtercov(remainlistMergeT1D, diseasein, diseaseout, false);
		} catch (IOException e) {
			e.printStackTrace();
		}


	}

	private void dedup(String remainlistRA, String remainlistT1D, String remainlistFilteredRA, String remainListFilteredT1D) throws IOException {
		TextFile tf = new TextFile(remainlistRA, TextFile.R);
		ArrayList<String> sampleList = tf.readAsArrayList();
		tf.close();

		TextFile tf2 = new TextFile(remainlistT1D, TextFile.R);
		ArrayList<String> sampleList2 = tf2.readAsArrayList();
		tf2.close();


		HashSet<String> duplicates = new HashSet<>();
		HashSet<String> visited = new HashSet<String>();
		visited.addAll(sampleList);
		for (String sample : sampleList2) {
			if (visited.contains(sample)) {
				duplicates.add(sample);
			}
		}

		System.out.println(duplicates.size() + " duplicates...");

		TextFile out = new TextFile(remainlistFilteredRA, TextFile.W);
		for (String sample : sampleList) {
			if (!duplicates.contains(sample)) {
				out.writeln(sample);
			}
		}

		TextFile out2 = new TextFile(remainListFilteredT1D, TextFile.W);

		for (String sample : sampleList2) {
			if (!duplicates.contains(sample)) {
				out2.writeln(sample);
			}
		}

		for (String sample : duplicates) {
			if (Math.random() > 0.5) {
				out.writeln(sample);
			} else {
				out2.writeln(sample);
			}
		}

		out.close();
		out2.close();
	}


	public void filterfam(String samplenamelist, String famin, String famout) throws IOException {
		TextFile tf = new TextFile(samplenamelist, TextFile.R);
		ArrayList<String> sampleList = tf.readAsArrayList();
		tf.close();

		HashSet<String> set = new HashSet<String>();
		set.addAll(sampleList);

		TextFile fin = new TextFile(famin, TextFile.R);
		TextFile fout = new TextFile(famout, TextFile.W);

		int remain = 0;
		String[] elems = fin.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length >= 2) {
				String sample = elems[1];
				if (set.contains(sample)) {
					fout.writeln(Strings.concat(elems, Strings.tab));
					remain++;
				}
			}
			elems = fin.readLineElems(TextFile.tab);
		}
		fin.close();
		fout.close();
		System.out.println(remain + " samples remain");
	}

	public void filtercov(String samplenamelist, String covin, String covout, boolean writeheader) throws IOException {
		TextFile tf = new TextFile(samplenamelist, TextFile.R);
		ArrayList<String> sampleList = tf.readAsArrayList();
		tf.close();

		HashSet<String> set = new HashSet<String>();
		set.addAll(sampleList);

		TextFile fin = new TextFile(covin, TextFile.R);
		TextFile fout = new TextFile(covout, TextFile.W);
		if (writeheader) {
			fout.writeln(fin.readLine());
		}
		int remain = 0;
		String[] elems = fin.readLineElems(TextFile.tab);
		while (elems != null) {
			String sample = elems[0];
			if (set.contains(sample)) {
				fout.writeln(Strings.concat(elems, Strings.tab));
				remain++;
			}
			elems = fin.readLineElems(TextFile.tab);
		}
		fin.close();
		fout.close();
		System.out.println(remain + " samples remain");
	}

}
