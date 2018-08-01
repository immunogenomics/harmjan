package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by Harm-Jan on 04/20/16.
 */
public class GWASLogComparison {

	public static void main(String[] args) {

		String loglist = "D:\\tmp\\2016-04-20\\T1DEurLogs.txt";
		String gwaslist = "D:\\tmp\\2016-04-20\\T1DEurGWAS.txt";
		String bedregionfile = "D:\\tmp\\2016-04-20\\";

	}

	public void run(String loglist, String gwaslist, String bedRegionFile) throws IOException {

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> bedRegions = reader.readAsList(bedRegionFile);

		TextFile listtf = new TextFile(loglist, TextFile.R);
		String[] loglistarr = listtf.readAsArray();
		listtf.close();
		listtf = new TextFile(gwaslist, TextFile.R);
		String[] gwaslistarr = listtf.readAsArray();
		listtf.close();


		HashMap<String, String> variantsLog = new HashMap<String, String>();
		HashSet<String> variantsThatShouldBeTested = new HashSet<String>();
		for (String file : loglistarr) {
			TextFile tf = new TextFile(file, TextFile.R);

			tf.readLine();
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				// SNP     Chr     Pos     ImputationQual  MAF     OverlapOK       MAFOk   ImpQualOK

				// rs77716349      1       197752167       0.0     null    true    null    false
				Boolean overlapok = false;
				Boolean mafOK = Boolean.parseBoolean(elems[6]);
				Boolean impQualOK = Boolean.parseBoolean(elems[7]);

				Feature f = new Feature();
				f.setChromosome(Chromosome.parseChr(elems[1]));
				Integer pos = Integer.parseInt(elems[2]);
				f.setStart(pos);
				f.setStop(pos);

				for (Feature r : bedRegions) {
					if (r.overlaps(f)) {
						overlapok = true;
						break;
					}
				}

				String snp = (elems[0] + "_" + elems[1] + "_" + elems[2]).toLowerCase();
				variantsLog.put(snp, Strings.concat(elems, Strings.tab));
				if (overlapok) {
					if (impQualOK) {
						if (mafOK) {
							variantsThatShouldBeTested.add(snp);
						}
					}
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
		}

		for (String file : loglistarr) {
			TextFile tf = new TextFile(file, TextFile.R);
			tf.readLine();
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {

				String snp = (elems[2] + "_" + elems[0] + "_" + elems[1]).toLowerCase();
				if (!variantsThatShouldBeTested.contains(snp)) {
					System.out.println(variantsLog.get(snp));
				}

				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
		}
	}
}
