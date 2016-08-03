package nl.harmjanwestra.finemapping;

import nl.harmjanwestra.utilities.individuals.Individual;
import nl.harmjanwestra.utilities.plink.PlinkFamFile;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

/**
 * Created by Harm-Jan on 07/19/16.
 */
public class SamplePicker {

	public static void main(String[] args) {

		String list = "D:\\tmp\\2016-08-03-sim\\sharedsamples-geneticsim0.2.txt";
		String out1 = "D:\\tmp\\2016-08-03-sim\\out\\listofsharedsamples-RA.txt";
		String out2 = "D:\\tmp\\2016-08-03-sim\\out\\listofsharedsamples-T1D.txt";
		String famt1d = "D:\\tmp\\2016-08-03-sim\\phenotypes\\T1D-recode-maf0005-ICRegions-samplenamefix-pseudo.vcf.gz-filtered-merged.fam";
		String famra = "D:\\tmp\\2016-08-03-sim\\phenotypes\\covarmerged.txtmergedfam.fam";
		SamplePicker s = new SamplePicker();
		try {
			s.picksamples(list, out1, out2);
			list = "D:\\tmp\\2016-08-03-sim\\RASamplex.txt";
			String rewriteout = "D:\\tmp\\2016-08-03-sim\\out\\listofsharedsamples-RA-rewritten.txt";
			s.rewrite(out1, list, rewriteout);
			String remainingout = "D:\\tmp\\2016-08-03-sim\\out\\listofsharedsamples-RA-rewritten-remaining.txt";
			s.writeListOfRemainingSamples(list, rewriteout, famra, remainingout);
			System.out.println("Filtering T1D samples");
			list = "D:\\tmp\\2016-08-03-sim\\T1DSamplex.txt";
			rewriteout = "D:\\tmp\\2016-08-03-sim\\out\\listofsharedsamples-T1D-rewritten.txt";
			s.rewrite(out2, list, rewriteout);
			remainingout = "D:\\tmp\\2016-08-03-sim\\out\\listofsharedsamples-T1D-rewritten-remaining.txt";
			s.writeListOfRemainingSamples(list, rewriteout, famt1d, remainingout);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void writeListOfRemainingSamples(String list1, String excludelist, String famfile, String out) throws IOException {
		TextFile tf = new TextFile(list1, TextFile.R);
		ArrayList<String> sampleList = tf.readAsArrayList();
		tf.close();

		TextFile tf2 = new TextFile(excludelist, TextFile.R);
		ArrayList<String> sampleList2 = tf2.readAsArrayList();
		tf2.close();

		HashSet<String> map = new HashSet<String>();
		for (String sample : sampleList2) {
			map.add(sample);
		}

		PlinkFamFile f = new PlinkFamFile(famfile);
		ArrayList<Individual> faminds = f.getSamples();
		HashSet<String> map2 = new HashSet<String>();
		for (Individual i : faminds) {
			if (!map.contains(i.getName())) {
				map2.add(i.getName());
			}
		}


		int ctr = 0;
		TextFile outf = new TextFile(out, TextFile.W);
		for (String sample : sampleList) {
			if (map2.contains(sample)) {
				outf.writeln(sample);
				ctr++;
			}
		}
		outf.close();
		System.out.println(ctr + " remaining.");

	}

	public void rewrite(String list, String samples, String out) throws IOException {
		TextFile tf = new TextFile(list, TextFile.R);
		ArrayList<String> sampleList = tf.readAsArrayList();
		tf.close();

		TextFile tf2 = new TextFile(samples, TextFile.R);
		ArrayList<String> sampleList2 = tf2.readAsArrayList();
		tf2.close();

		HashMap<String, String> map = new HashMap<String, String>();
		for (String sample : sampleList2) {
			map.put(sample + "_" + sample, sample);
		}

		TextFile outf = new TextFile(out, TextFile.W);
		for (String sample : sampleList) {
			String other = map.get(sample);
			if (other == null) {
				String[] elems = sample.split("_");
				if (elems.length > 2) {
					sample = Strings.concat(elems, Pattern.compile("_"), 1, elems.length);
					sample = sample + "_" + sample;
				} else {
					sample = elems[1] + "_" + elems[1];
				}
				other = map.get(sample);
				if (other == null) {
					System.out.println("could not find sample: " + sample);
				} else {
					outf.writeln(other);
				}
			} else {
				outf.writeln(other);
			}
		}
		outf.close();
	}

	public void picksamples(String in, String out1, String out2) throws IOException {
		TextFile tf = new TextFile(in, TextFile.R);

		ArrayList<Pair<String, String>> pairs = new ArrayList<Pair<String, String>>();
		String[] elems = tf.readLineElems(TextFile.tab);

		while (elems != null) {
			if (elems.length == 2) {
				pairs.add(new Pair<String, String>(elems[0], elems[1]));
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();


		HashSet<String> visitedSamples = new HashSet<String>();
		ArrayList<String> duplicates = new ArrayList<>();
		for (Pair<String, String> p : pairs) {
			if (!visitedSamples.contains(p.getLeft())) {
				visitedSamples.add(p.getLeft());
			} else {
				duplicates.add(p.getLeft());
			}

			if (!visitedSamples.contains(p.getRight())) {
				visitedSamples.add(p.getRight());
			} else {
				duplicates.add(p.getRight());
			}
		}


		// select which sample to exclude
		TextFile outtf1 = new TextFile(out1, TextFile.W);
		TextFile outtf2 = new TextFile(out2, TextFile.W);

		for (int i = 0; i < pairs.size(); i++) {

			Pair<String, String> p = pairs.get(i);
			double d = Math.random();
			if (d > 0.5) {
				outtf1.writeln(p.getLeft());
			} else {
				outtf2.writeln(p.getRight());
			}

		}

		outtf1.close();
		outtf2.close();

	}
}
