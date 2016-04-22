package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 4/22/16.
 */
public class VCFBedFilter {


	public void filterVCFForBedRegions(String vcfIn, String referenceOut, String regionFile, String chr) throws IOException {

		System.out.println("BED filter");
		System.out.println("in: " + vcfIn);
		System.out.println("out: " + referenceOut);
		System.out.println("bed: " + regionFile);

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> set = reader.readAsList(regionFile);
		if (chr != null) {
			ArrayList<Feature> tmp = new ArrayList<>();
			Chromosome chrom = Chromosome.parseChr(chr);
			for (Feature f : set) {
				if (f.getChromosome().equals(chrom)) {
					tmp.add(f);
				}
			}
			set = tmp;
			System.out.println(set.size() + " regions overlapping chr: " + chr);
		}


		TextFile vcfoutf = new TextFile(referenceOut, TextFile.W);
		TextFile vcfoutf2 = new TextFile(referenceOut + "-filterdOutVariants.txt.gz", TextFile.W);
		TextFile vcftf = new TextFile(vcfIn, TextFile.R);
		int nrHeaderElems = 9;
		String ln = vcftf.readLine();

		int parsed = 0;
		int saved = 0;


		while (ln != null) {
			if (ln.startsWith("##")) {
				vcfoutf.writeln(ln);
			} else if (ln.startsWith("#CHROM")) {
				vcfoutf.writeln(ln);
			} else {

				String[] elems = ln.split("\t");
				Feature f = new Feature();
				f.setChromosome(Chromosome.parseChr(elems[0]));
				f.setStart(Integer.parseInt(elems[1]));
				f.setStop(Integer.parseInt(elems[1]));
				boolean overlap = false;
				for (Feature region : set) {
					if (f.overlaps(region)) {
						overlap = true;
						break;
					}
				}

				if (overlap) {
					vcfoutf.writeln(ln);
					saved++;
				} else {
					vcfoutf2.writeln(Strings.concat(elems, Strings.tab, 0, 9));
				}
				parsed++;

				if (parsed % 50000 == 0) {
					System.out.println(saved + "/" + parsed + " written");
				}
			}

			ln = vcftf.readLine();
		}

		vcfoutf.close();
		vcfoutf2.close();
		vcftf.close();

		System.out.println(saved + "/" + parsed + " written");
	}

}
