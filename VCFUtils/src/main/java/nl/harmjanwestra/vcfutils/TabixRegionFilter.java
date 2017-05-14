package nl.harmjanwestra.vcfutils;

import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.vcf.VCFTabix;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 5/12/17.
 */
public class TabixRegionFilter {

	public void run(String in, String regionfile, String out) {
		try {
			BedFileReader bfr = new BedFileReader();
			ArrayList<Feature> regions = bfr.readAsList(regionfile);
			VCFTabix t = new VCFTabix(in);

			TextFile outf = new TextFile(out, TextFile.W);
			TextFile tfin = new TextFile(in, TextFile.R);
			String ln = tfin.readLine();
			while (ln.startsWith("#")) {
				outf.writeln(ln);
				ln = tfin.readLine();
			}
			tfin.close();


			for (Feature region : regions) {
				TabixReader.Iterator iterator = t.query(region);
				ln = iterator.next();
				while (ln != null) {
					outf.writeln(ln);
					ln = iterator.next();
				}
			}
			outf.close();

		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
