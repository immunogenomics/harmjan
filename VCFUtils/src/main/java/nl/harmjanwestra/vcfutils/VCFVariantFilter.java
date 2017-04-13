package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashSet;

/**
 * Created by Harm-Jan on 07/24/16.
 */
public class VCFVariantFilter {

	enum FILTER {
		positions,
		rsids
	}

	public void filter(String vcf, String out, String listfile, FILTER filter, boolean exclude) throws IOException {
		boolean run = true;
		if (vcf == null) {
			System.err.println("Input file not set");
			run = false;
		}
		if (out == null) {
			System.err.println("Output file not set");
			run = false;
		}
		if (listfile == null) {
			System.err.println("Variant list file not set");
			run = false;
		}

		if (!run) {
			System.exit(-1);
		}
		// expect a bed file?
		TextFile tf = new TextFile(listfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashSet<Feature> variants = new HashSet<Feature>();
		HashSet<String> rsids = new HashSet<String>();
		System.out.println("Loading variant ids from " + listfile);
		while (elems != null) {
			if (elems.length >= 3) {
				Feature f = new Feature();
				Chromosome chrobj = Chromosome.parseChr(elems[0]);
				int pos = Integer.parseInt(elems[1]);

				f.setChromosome(chrobj);
				f.setStart(pos);
				f.setStop(pos);

				if (elems.length >= 4) {
					String name = elems[3];
					f.setName(name);
					rsids.add(name);
				}
				variants.add(f);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println(rsids.size() + " rsids loaded");
		System.out.println(variants.size() + " variants loaded");

		TextFile vcfin = new TextFile(vcf, TextFile.R);
		TextFile vcfout = new TextFile(out, TextFile.W);

		String ln = vcfin.readLine();
		int filtered = 0;
		int written = 0;
		int lnctr = 0;
		while (ln != null) {
			if (ln.startsWith("#")) {
				vcfout.writeln(ln);
				written++;
			} else {
				elems = ln.substring(0, 200).split("\t");
				Feature f = new Feature();
				Chromosome chrobj = Chromosome.parseChr(elems[0]);
				int pos = Integer.parseInt(elems[1]);

				f.setChromosome(chrobj);
				f.setStart(pos);
				f.setStop(pos);

				boolean debug = false;
				if (elems.length >= 3) {
					String name = elems[2];
					f.setName(name);
					if (name.equals("rs3184504")) {
						System.out.println("ping!");
						debug = true;
					}
				}

				boolean present;
				if (filter == FILTER.rsids && f.getName() != null) {
					present = rsids.contains(f.getName());
				} else {
					if (debug) {
						System.out.println("Testing for presence in list...");
					}
					present = isPresent(f, variants);
					if (debug) {
						System.out.println("presence: " + present);
					}
				}

				if ((exclude && !present) || (!exclude && present)) {
					vcfout.writeln(ln);
					written++;
				} else {
					filtered++;
				}
			}
			lnctr++;
			if (lnctr % 1000 == 0) {
				System.out.print(lnctr + " lines,\t" + filtered + " filtered,\t" + written + " written.\r");
			}
			ln = vcfin.readLine();
		}
		vcfin.close();
		vcfout.close();
		System.out.println();
		System.out.println("Done");
	}

	private boolean isPresent(Feature f, HashSet<Feature> variants) {
		for (Feature f2 : variants) {
			if (f.getChromosome().equals(f2.getChromosome()) && f.getStart() == f2.getStart()) {
				return true;
			}
		}
		return false;
	}

}
