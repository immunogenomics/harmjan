package nl.harmjanwestra.gwas;

import nl.harmjanwestra.gwas.CLI.BedAssocFilterOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 2/24/16.
 */
public class BedAssocFilter {

	private final BedAssocFilterOptions options;

	public BedAssocFilter(BedAssocFilterOptions options) throws IOException {
		this.options = options;
		this.run();
	}

	public void run() throws IOException {

		// load bed regions to test
		BedFileReader bf = new BedFileReader();
		ArrayList<Feature> regions = bf.readAsList(options.regionFile);

		TextFile out = new TextFile(options.outfile, TextFile.R);

		double threshold = -Math.log10(options.threshold);

		for (int i = 0; i < regions.size(); i++) {
			Feature region = regions.get(i);
			AssociationFile assocFile = new AssociationFile();
			ArrayList<AssociationResult> results = assocFile.read(options.assocFile, region);

			boolean written = false;
			for (int j = 0; j < results.size() && !written; j++) {
				AssociationResult result = results.get(j);
				if (result.getLog10Pval() > threshold) {
					out.writeln(region.toBedString());
					written = true;
				}
			}

		}
		out.close();


	}
}
