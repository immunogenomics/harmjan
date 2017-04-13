package nl.harmjanwestra.goshifter;

import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

/**
 * Created by hwestra on 2/25/15.
 */
public class GoShifterJobCreator {

	public void run(String[] args) {


		if (args.length < 8) {
			System.out.println("Usage: goshifter chiplist snpmap nrperm lddir out jobout query");
		} else {
			String goshifter = null;
			String chipList = null;
			String snpmap = null;
			String nrPerm = null;
			String ldDir = null;
			String out = null;
			String jobout = null;
			String query = null;

			for (int i = 0; i < args.length; i++) {
				if (args[i].equals("--goshifter")) {
					if (i + 1 < args.length) {
						goshifter = args[i + 1];
					}
				} else if (args[i].equals("--chiplist")) {
					if (i + 1 < args.length) {
						chipList = args[i + 1];
					}
				} else if (args[i].equals("--snpmap")) {
					if (i + 1 < args.length) {
						snpmap = args[i + 1];
					}
				} else if (args[i].equals("--nrperm")) {
					if (i + 1 < args.length) {
						nrPerm = args[i + 1];
					}
				} else if (args[i].equals("--lddir")) {
					if (i + 1 < args.length) {
						ldDir = args[i + 1];
					}
				} else if (args[i].equals("--out")) {
					if (i + 1 < args.length) {
						out = args[i + 1];
					}
				} else if (args[i].equals("--jobout")) {
					if (i + 1 < args.length) {
						jobout = args[i + 1];
					}
				} else if (args[i].equals("--query")) {
					if (i + 1 < args.length) {
						query = args[i + 1];
					}
				}
			}


			boolean ok = true;
			ok &= checkFile("--goshifter", goshifter);
			ok &= checkFile("--chiplist", chipList);
			ok &= checkFile("--snpmap", snpmap);
			if (nrPerm == null) {
				System.out.println("Could not find setting: --nrperm");
				ok &= false;
			}
			if (out == null) {
				System.out.println("Could not find setting: --out");
				ok &= false;
			}

			ok &= checkFile("--lddir", ldDir);
			ok &= checkFile("--jobout", jobout);
			ok &= checkFile("--query", query);

			if (!ok) {
				System.err.println("Some files were not found");
				printUsage();
			} else {
				try {
					if (!Gpio.exists(out)) {
						out = Gpio.formatAsDirectory(out);
						Gpio.createDir(out);
					}

					if (!Gpio.exists(jobout)) {
						jobout = Gpio.formatAsDirectory(jobout);
						Gpio.createDir(jobout);
					}

					createJobs(goshifter, chipList, snpmap, nrPerm, ldDir, out, jobout, query);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		System.exit(0);
	}

	private boolean checkFile(String name, String f) {
		if (f == null) {
			System.out.println("Could not find path: " + f + " for switch " + name);
			return false;
		}
		boolean b = Gpio.exists(f);
		if (!b) {
			System.out.println("could not find path: " + f + " for switch " + name);
		}
		return b;
	}

	private void printUsage() {
		System.out.println("Usage: --jobs --goshifter /path/to/goshifter.py --chiplist /path/to/listofAnnotationFiles.txt --snpmap /path/to/snpmap.txt --nrperm int --lddir /path/to/ldfiles/ --out /path/to/outdir/ --jobout /path/to/joboutput/ --query /path/to/listOfQuerySNPs.txt");

	}

	private void createJobs(String goshifter, String chipList, String snpmap, String nrPerm, String ldDir, String out, String jobout, String query) throws IOException {
		TextFile q = new TextFile(query, TextFile.R);
		HashSet<String> uniqueVariants = new HashSet<String>();
		uniqueVariants.addAll(q.readAsArrayList());
		q.close();

		HashSet<String> variantsFound = new HashSet<String>();


		TextFile snpmapout = new TextFile(snpmap + "_filtered.txt", TextFile.W);
		TextFile q2 = new TextFile(snpmap, TextFile.R);
		snpmapout.writeln(q2.readLine());
		String[] elems = q2.readLineElems(TextFile.tab);
		int lnsWritten = 0;
		while (elems != null) {
			String snp = elems[0];
			if (uniqueVariants.contains(snp)) {
				snpmapout.writeln(Strings.concat(elems, Strings.tab));
				variantsFound.add(snp);
				lnsWritten++;
			}
			elems = q2.readLineElems(TextFile.tab);
		}
		q2.close();
		snpmapout.close();

		System.out.println(lnsWritten + "\t" + variantsFound.size());
		for (String s : variantsFound) {
			System.out.println(s);
		}

		String snpmapoutname = snpmapout.getFileName();

		TextFile tf = new TextFile(chipList, TextFile.R);
		String[] listOfBedFiles = tf.readAsArray();
		tf.close();

		for (String file : listOfBedFiles) {

			String[] fname = file.split("/");
			String fnames = fname[fname.length - 1];

			fnames = fnames.replaceAll(".gz", "");
			fnames = fnames.replaceAll(".bed", "");


			String run = "nice -n 15 python " + goshifter + " --snpmap " + snpmapoutname + " --annotation " + file + " --permute " + nrPerm + " --ld " + ldDir + " --out " + out + fnames + " &> " + out + fnames + ".log";

			TextFile outf = new TextFile(jobout + fnames + ".sh", TextFile.W);
			outf.writeln(run);
			outf.close();

		}
	}
}
