package nl.harmjanwestra.vcfutils;

import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.vcf.VCFTabix;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Created by hwestra on 5/12/17.
 */
public class TabixRegionFilter {

	public void run(String in, String regionfile, String out) {
		try {

			BedFileReader bfr = new BedFileReader();
			ArrayList<Feature> regions = bfr.readAsList(regionfile);

			Collections.sort(regions, new FeatureComparator());

			VCFTabix t = new VCFTabix(in);
			System.out.println(regions.size() + " regions to filter..");

			TextFile outf = new TextFile(out, TextFile.W);
			TextFile tfin = new TextFile(in, TextFile.R);
			String ln = tfin.readLine();
			while (ln.startsWith("#")) {
				outf.writeln(ln);
				ln = tfin.readLine();
			}
			tfin.close();


			int totallns = 0;
			for (int r = 0; r < regions.size(); r++) {
				Feature region = regions.get(r);
				TabixReader.Iterator iterator = t.query(region);
				ln = iterator.next();
				int lnctr = 0;
				while (ln != null) {
					outf.writeln(ln);
					ln = iterator.next();
					lnctr++;
					totallns++;
					if (lnctr % 1000 == 0) {
						System.out.print(lnctr + " lines written for region: " + region.toString() + ". Total lines written: " + totallns + "\r");
					}
				}
			}
			outf.close();
			System.out.println("Done.. " + totallns + " written.");

		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
