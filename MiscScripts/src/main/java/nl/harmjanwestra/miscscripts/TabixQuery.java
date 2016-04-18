package nl.harmjanwestra.miscscripts;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 4/18/16.
 */
public class TabixQuery {


	public static void main(String[] args) {
		try {

			TabixQuery q = new TabixQuery();
			q.convertToProxyFinderFormat();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	public void makeTabixJobs() throws IOException {
		String snpfile = "/Data/tmp/2016-03-25/RA-snpmap.txt";
		String scriptout = "/Data/tmp/2016-03-25/RA-tabixquery.sh";

		TextFile in = new TextFile(snpfile, TextFile.R);
		in.readLine();
		TextFile out = new TextFile(scriptout, TextFile.W);

		String ln = in.readLine();
		while (ln != null) {
			String[] elems = ln.split("\t");
			if (elems.length >= 3) {
				String snp = elems[0];
				String chr = elems[1];
				Integer pos = Integer.parseInt(elems[2]);
				int sta = pos - 500000;
				int sto = pos + 500000;
				double r2min = 0.8;
				String tabixfile = "/Code/Python/GoShifter/ld/" + chr + ".EUR.tsv.gz";
				out.writeln("tabix " + tabixfile + " " + chr + ":" + sta + "-" + sto + " | awk '$6 >= " + r2min + " {{print $0}}' | grep -w " + snp);
			}
			ln = in.readLine();
		}


		in.close();
		out.close();
	}

	public void convertToProxyFinderFormat() throws IOException {


		String insnps = "/Data/tmp/2016-03-25/RA-snpmap.txt";
		String in = "/Data/tmp/2016-03-25/oldProxies.txt";
		String out = "/Data/tmp/2016-03-25/oldProxies-broshifter.txt";

		HashSet<String> querySNPs = new HashSet<String>();
		TextFile insnpstf = new TextFile(insnps, TextFile.R);

		insnpstf.readLine();
		String[] insnpselems = insnpstf.readLineElems(TextFile.tab);

		HashMap<String, String[]> snpInput = new HashMap<String, String[]>();
		while (insnpselems != null) {
			querySNPs.add(insnpselems[0]);
			snpInput.put(insnpselems[0], insnpselems);
			insnpselems = insnpstf.readLineElems(TextFile.tab);
		}
		insnpstf.close();


		TextFile inf = new TextFile(in, TextFile.R);
		TextFile outf = new TextFile(out, TextFile.W);

		String[] elems = inf.readLineElems(TextFile.tab);
		HashSet<String> seenSNPs = new HashSet<String>();
		HashSet<String> outputlns = new HashSet<String>();
		while (elems != null) {

			int pos1 = Integer.parseInt(elems[1]);
			int pos2 = Integer.parseInt(elems[3]);

			int diff = Math.abs(pos1 - pos2);

			String snp1 = elems[2];
			String snp2 = elems[4];

			// chr1	114303808	rs6679677	114303808	rs6679677	1	1
//		Chr1	114303808	rs6679677	Chr1	114303808	rs6679677	0	1.0

			if (querySNPs.contains(snp1)) {
				seenSNPs.add(snp1);
				String outln = elems[0]
						+ "\t" + elems[1]
						+ "\t" + elems[2]
						+ "\t" + elems[0]
						+ "\t" + elems[3]
						+ "\t" + elems[4]
						+ "\t" + diff
						+ "\t" + elems[6];
				if (!outputlns.contains(outln)) {
					outf.writeln(outln);
					outputlns.add(outln);
				} else {
					System.out.println("Already printed line");
				}

			} else {
				if (querySNPs.contains(snp2)) {
					// snp is on second position.. flip around..
					System.out.println("flip da snp " + snp2);
					String outln = elems[0]
							+ "\t" + elems[3]
							+ "\t" + elems[4]
							+ "\t" + elems[0]
							+ "\t" + elems[1]
							+ "\t" + elems[2]
							+ "\t" + diff
							+ "\t" + elems[6];
					if (!outputlns.contains(outln)) {
						outf.writeln(outln);
						outputlns.add(outln);
					} else {
						System.out.println("Already printed swapped line");
					}
				}
			}
			elems = inf.readLineElems(TextFile.tab);
		}
		inf.close();

		// check whether any snps were missed..
		for (String snp : querySNPs) {
			if (!seenSNPs.contains(snp)) {
				String[] snpelems = snpInput.get(snp);
				System.out.println(snp + " not found in proxylist.");
				String outln = snpelems[1]
						+ "\t" + snpelems[2]
						+ "\t" + snpelems[0]
						+ "\t" + snpelems[1]
						+ "\t" + snpelems[2]
						+ "\t" + snpelems[0]
						+ "\t" + 0
						+ "\t" + 1;
				if (!outputlns.contains(outln)) {
					outf.writeln(outln);
					outputlns.add(outln);
				} else {
					System.out.println("Already printed missing line");
				}
			}
		}


		outf.close();
	}


}
