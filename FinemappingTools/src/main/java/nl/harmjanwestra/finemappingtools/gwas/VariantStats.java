package nl.harmjanwestra.finemappingtools.gwas;

import nl.harmjanwestra.finemappingtools.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.vcf.VCFVariantComparator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Created by hwestra on 5/5/17.
 */
public class VariantStats extends LRTest {


	public VariantStats(LRTestOptions options) throws IOException {
		super(options);
		run();
	}

	public void run() throws IOException {
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(options.getBedfile());
		ArrayList<VCFVariant> variants = readVariants(options.getOutputdir() + "variantlog.txt", regions);
		Collections.sort(variants, new VCFVariantComparator());

		TextFile out = new TextFile(options.getOutputdir() + "variantstats.txt", TextFile.W);

		String header = "Chr\tPos\tId\tCRCases\tCRControls\tAFCases\tAFControls\tHWE-PCases\tHWE-PControls";
		out.writeln(header);
		for (VCFVariant variant : variants) {
			String ln = variant.getChr() + "\t" + variant.getPos() + "\t" + variant.getId() + "\t" + variant.toString();
			ln += "\t" + variant.getCallrateCases() + "\t" + variant.getCallrateControls() + "\t" + variant.getDiffMissingnessP();
//			ln += "\t" + variant.getNumberCasesCalled() + "\t" + variant.getControlsCalled();
			ln += "\t" + variant.getAlleleFrequenciesCases()[0] + "\t" + variant.getAlleleFrequenciesControls()[0];
			ln += "\t" + variant.getHwepCases() + "\t" + variant.getHwepControls();
			out.writeln(ln);
		}
		out.close();

	}
}
