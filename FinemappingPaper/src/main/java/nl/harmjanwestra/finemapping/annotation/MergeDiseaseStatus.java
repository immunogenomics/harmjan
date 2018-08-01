package nl.harmjanwestra.finemapping.annotation;

import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 7/18/17.
 */
public class MergeDiseaseStatus {
	
	public static void main(String[] args) {
		try {
//			String in1 = "/Data/tmp/2017-07-18/2016-03-11-T1D-diseaseStatusWithPseudos.txt";
//			String in2 = "/Data/tmp/2017-07-18/covarmerged.txtmergeddisease.txt";
//			String out = "/Data/tmp/2017-07-18/diseaseStatusMulti.txt";
//
//			run(in1, in2, out);
			
			String in1 = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\multinomial\\diseaseStatus\\listofsharedsamples-RA-rewritten-remaining-dedup.txt";
			String in2 = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\multinomial\\diseaseStatus\\listofsharedsamples-T1D-rewritten-remaining-dedup.txt";
			String disease = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\multinomial\\diseaseStatus\\diseaseStatusMulti.txt";
			String out = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\multinomial\\diseaseStatus\\controlsCohorts.txt";
//			runSelectControls(in1, in2, disease, out);
			
			String cov = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\multinomial\\meta-pca.txt";
			out = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\multinomial\\meta-pcaWRALabel.txt";
			
			addCohortId(in1, in2, cov, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static void addCohortId(String in1, String in2, String cov, String out) throws IOException {
		TextFile tf1 = new TextFile(in1, TextFile.R);
		ArrayList<String> samples1 = tf1.readAsArrayList();
		tf1.close();
		
		TextFile tf2 = new TextFile(in2, TextFile.R);
		ArrayList<String> samples2 = tf2.readAsArrayList();
		tf1.close();
		
		HashSet<String> set1 = new HashSet<>();
		set1.addAll(samples1);
		
		TextFile output = new TextFile(out, TextFile.W);
		TextFile in = new TextFile(cov, TextFile.R);
		String header = in.readLine() + "\tIsRA";
		output.writeln(header);
		String ln = in.readLine();
		while (ln != null) {
			String[] elems = ln.split("\t");
			String sample = elems[0];
			int isRA = 0;
			if (set1.contains(sample)) {
				isRA = 1;
			}
			output.writeln(ln + "\t" + isRA);
			ln = in.readLine();
		}
		in.close();
		output.close();
	}
	
	private static void runSelectControls(String in1, String in2, String disease, String out) throws IOException {
		
		HashSet<String> controls = new HashSet<String>();
		TextFile tf = new TextFile(disease, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, DiseaseStatus> diseaseStatus1 = new HashMap<String, DiseaseStatus>();
		while (elems != null) {
			String sample = elems[0];
			DiseaseStatus s1 = DiseaseStatus.parseStatus(elems[1]);
			DiseaseStatus s2 = DiseaseStatus.parseStatus(elems[2]);
			if (s1.equals(DiseaseStatus.CONTROL) && s2.equals(DiseaseStatus.CONTROL)) {
				controls.add(sample);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile tf1 = new TextFile(in1, TextFile.R);
		ArrayList<String> samples1 = tf1.readAsArrayList();
		tf1.close();
		
		TextFile tf2 = new TextFile(in2, TextFile.R);
		ArrayList<String> samples2 = tf2.readAsArrayList();
		tf1.close();
		
		TextFile outfd = new TextFile(out, TextFile.W);
		int fromFile1 = 0;
		for (String s : samples1) {
			if (controls.contains(s)) {
				outfd.writeln(s + "\tControl");
				fromFile1++;
			}
		}
		int fromFile2 = 0;
		for (String s : samples2) {
			if (controls.contains(s)) {
				outfd.writeln(s + "\tCase");
				fromFile2++;
			}
		}
		outfd.close();
		
		System.out.println("From file1: " + fromFile1);
		System.out.println("From file2: " + fromFile2);
		
	}
	
	public static void run(String in1, String in2, String out) throws IOException {
		
		HashMap<String, DiseaseStatus> status1 = read(in1);
		HashMap<String, DiseaseStatus> status2 = read(in2);
		HashSet<String> sharedSamples = new HashSet<String>();
		HashSet<String> allSamples = new HashSet<String>();
		allSamples.addAll(status1.keySet());
		allSamples.addAll(status2.keySet());
		System.out.println(allSamples.size() + " unique samples");
		
		for (String sample : status1.keySet()) {
			if (status2.containsKey(sample)) {
				DiseaseStatus s1 = status1.get(sample);
				DiseaseStatus s2 = status2.get(sample);
				if (!s1.equals(s2)) {
					System.out.println("WARNING: different disease status for sample: " + sample + "\t" + status1 + "\t" + status2);
				}
				sharedSamples.add(sample);
				
			}
			
		}
		System.out.println(sharedSamples.size() + " shared samples");
		
		TextFile outf = new TextFile(out, TextFile.W);
		for (String sample : allSamples) {
			
			
			DiseaseStatus s1 = status1.get(sample);
			if (s1 == null) {
				s1 = DiseaseStatus.CONTROL;
			}
			DiseaseStatus s2 = status2.get(sample);
			if (s2 == null) {
				s2 = DiseaseStatus.CONTROL;
			}
			
			String ln = sample + "\t" + s1.toString() + "\t" + s2.toString();
			outf.writeln(ln);
		}
		
		outf.close();
		
	}
	
	public static HashMap<String, DiseaseStatus> read(String in) throws IOException {
		TextFile tf = new TextFile(in, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, DiseaseStatus> diseaseStatus1 = new HashMap<String, DiseaseStatus>();
		while (elems != null) {
			String sample = elems[0];
			DiseaseStatus s = DiseaseStatus.parseStatus(elems[1]);
			diseaseStatus1.put(sample, s);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return diseaseStatus1;
	}
}
