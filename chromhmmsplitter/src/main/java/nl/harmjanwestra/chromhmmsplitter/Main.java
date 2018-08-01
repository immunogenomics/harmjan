package nl.harmjanwestra.chromhmmsplitter;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by Harm-Jan on 02/27/16.
 */
public class Main {

	public static void main(String[] args) {

//		if (args.length < 2) {
//			System.out.println("usage: dirin dirout");
//		} else {
//			Main m = new Main();
//			try {
//
//				m.split(args[0], args[1]);
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//		}

		Main m = new Main();
		try {
			m.annotate("/Data/Enhancers/ChromHMM/celltypenames.txt",
					"/Data/Enhancers/ChromHMM/allChromHMM.txt",
					"/Data/Enhancers/ChromHMM/allChromHMM-desc.txt");
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void annotate(String samplenamefile, String bedfilefile, String out) throws IOException {

		HashMap<String, String> ecodes = new HashMap<String, String>();
		TextFile tf = new TextFile(samplenamefile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 1) {
				String ecode = elems[0];
				String desc = elems[1];
				ecodes.put(ecode, desc);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile bedin = new TextFile(bedfilefile, TextFile.R);
		TextFile outf = new TextFile(out, TextFile.W);

		String ln = bedin.readLine();
		while (ln != null) {

			String f = Gpio.getFileName(ln);


			// E129_25_imputed12marks_mnemonics-9_TxReg.bed.gz
			String[] lnelems = f.split("_");
			String tissue = ecodes.get(lnelems[0]);
			lnelems = f.split("-")[1].split("_");
			String type = lnelems[1].substring(0, lnelems[1].length() - 7);

			if(!type.equals("Quies")) {
				String desc = type + " / " + tissue;
				outf.writeln(desc + "\t" + ln);
			}
			ln = bedin.readLine();
		}

		outf.close();
		bedin.close();

	}

	public void split(String indir, String outdir) throws IOException {

		String[] files = Gpio.getListOfFiles(indir);
		System.out.println(files.length + " files found");
		int ctr = 0;
		for (String file : files) {
			if (file.endsWith(".bed.gz")) {
				File fileobj = new File(file);
				String name = fileobj.getName();
				name = name.replaceAll(".bed.gz", "");
				System.out.print("Processssing " + ctr + "/" + files.length + "\t" + name + " .. ");

				ArrayList<ArrayList<Feature>> data = new ArrayList<ArrayList<Feature>>(26);
				for (int i = 0; i < 26; i++) {
					data.add(new ArrayList<Feature>(1000));
				}

				String[] mnemonics = new String[26];

				TextFile tf = new TextFile(file, TextFile.R);
				String[] elems = tf.readLineElems(Strings.whitespace);
				while (elems != null) {
					Feature f = new Feature();
					f.setChromosome(Chromosome.parseChr(elems[0]));
					int start = Integer.parseInt(elems[1]);
					int stop = Integer.parseInt(elems[2]);
					f.setStart(start);
					f.setStop(stop);

					String mnemonic = elems[3];
					String[] mnemonicelems = mnemonic.split("_");

					int nr = Integer.parseInt(mnemonicelems[0]);
					ArrayList<Feature> storage = data.get(nr);
					storage.add(f);
					mnemonics[nr] = mnemonic;

					elems = tf.readLineElems(Strings.whitespace);
				}
				tf.close();


				int nr = 0;
				for (int i = 0; i < data.size(); i++) {
					ArrayList<Feature> storage = data.get(i);

					if (storage.size() > 0) {
						String mnemonic = mnemonics[i];
						while (mnemonic.contains("/")) {
							mnemonic = mnemonic.replaceAll("/", "");
						}
						while (mnemonic.contains("'")) {
							mnemonic = mnemonic.replaceAll("'", "");
						}
						TextFile tf2 = new TextFile(outdir + name + "-" + mnemonic + ".bed.gz", TextFile.W);
						for (Feature f : storage) {
							tf2.writeln(f.getChromosome().toString() + "\t" + f.getStart() + "\t" + f.getStop());
						}
						tf2.close();
						nr++;
					}

				}
				System.out.println(nr + " files split :)");
				ctr++;
			}

		}


	}
}
