package nl.harmjanwestra.goscializer;

import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 3/28/16.
 */
public class Main {

	public static void main(String[] args) {

		try {

			String snplist = "/Users/hwestra/Downloads/SNPs20160328.txt";
			String dbsnp = "/Data/Projects/2014-FR-Reseq/GATKResources/dbsnp_138.b37.vcf.gz";
			String out = "/Users/hwestra/Downloads/SNPs20160328-snpmap.txt";

			Main m = new Main();
//			m.getAnnotationForSNPList(snplist, dbsnp, out);

			m.makeJobs();

		} catch (IOException e) {
			e.printStackTrace();

		}

	}

	public void getAnnotationForSNPList(String snplist, String dbsnpVCF, String out) throws IOException {

		ArrayList<String> rsIds = new ArrayList<String>();
		TextFile tf = new TextFile(snplist, TextFile.R);
		String ln = tf.readLine();
		while (ln != null) {
			rsIds.add(ln);
			ln = tf.readLine();
		}
		tf.close();
		System.out.println(rsIds.size() + " rsids loaded....");

		HashSet<String> rsIdHash = new HashSet<String>();
		rsIdHash.addAll(rsIds);

		HashMap<String, String> rsIdToAnnotation = new HashMap<String, String>();

		VCFGenotypeData data = new VCFGenotypeData(dbsnpVCF);
		int ctr = 0;
		while (data.hasNext()) {
			VCFVariant variant = data.next();
			if (rsIdHash.contains(variant.getId())) {
				rsIdToAnnotation.put(variant.getId(), variant.getId() + "\t" + variant.getChr().toLowerCase() + "\t" + variant.getPos());
			}

			if (ctr % 10000 == 0) {
				System.out.println(ctr + " parsed");
			}
			ctr++;
		}

		TextFile tfout = new TextFile(out, TextFile.W);
		tfout.writeln("SNP\tChrom\tBP");
		for (String s : rsIds) {
			String output = rsIdToAnnotation.get(s);
			if (output != null) {
				tfout.writeln(output);
			}
		}
		tfout.close();

	}


	public void makeJobs() throws IOException {

		// ./goshifter.py [--snpmap FILE --proxymap FILE] --annotation FILE --permute INT --ld DIR --out FILE [--rsquared NUM --window NUM --min-shift NUM --max-shift NUM --ld-extend NUM --no-ld]


		TextFile tf = new TextFile("/Data/tmp/2016-03-28/h3k4me3.txt", TextFile.R);
		ArrayList<String> annotations = new ArrayList<String>();
		String ln = tf.readLine();
		while (ln != null) {
			annotations.add(ln);
			ln = tf.readLine();
		}
		tf.close();

		TextFile out = new TextFile("/Data/tmp/2016-03-28/h3k4me3-jobs.sh", TextFile.W);
		for (String s : annotations) {
			String ds = "encode";
			if (s.contains("roadmap")) {
				ds = "roadmap";
			}
			String[] annotation = s.split("/");
			String annotname = annotation[annotation.length - 1];
			annotname = ds + "_" + annotname;
			annotname = annotname.replaceAll(".xls.gz", "");
			annotname = annotname.replaceAll(".bed.gz", "");

			System.out.println(annotname);


			String job = "/medpop/srlab/hwestra/tools/python2.7/bin/python /medpop/srlab/hwestra/goshifter/goshifter.py ";
			job += "--snpmap /medpop/srlab/hwestra/chikashi/2016-03-28-SNPS.txt ";
			job += "--permute 10000 ";
			job += "--ld /medpop/srlab/external-data/1000genomes/GoShifterLd/ ";
			job += "--annotation " + s + " ";
			job += "--out /medpop/srlab/hwestra/chikashi/output/" + annotname + " ";
			job += " > /medpop/srlab/hwestra/chikashi/output/" + annotname + ".log ";
			out.writeln("echo " + s);
			out.writeln(job);

		}
		out.close();


	}

}
