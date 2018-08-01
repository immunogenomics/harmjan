package nl.harmjanwestra.finemappingtools.gwas;

import nl.harmjanwestra.finemappingtools.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.vcf.filter.variantfilters.*;

import java.io.IOException;
import java.util.ArrayList;

public class VariantFilterList extends LRTest {
	public VariantFilterList(LRTestOptions options) throws IOException {
		super(options);
		run();
	}
	
	private void run() throws IOException {
		
		String bedfile = options.getBedfile();
		ArrayList<Feature> regions = null;
		if (bedfile != null) {
			BedFileReader r = new BedFileReader();
			regions = r.readAsList(bedfile);
		}
		
		VCFVariantFilters filter = new VCFVariantFilters();
		
		// there is a file that limits the snps to include
		if (options.getSnpLimitFile() != null) {
			filter.addFilter(new VCFVariantSetFilter(options.getSnpLimitFile()));
		}
		if (options.getMissingNessP() != null) {
			filter.addFilter(new VCFVariantMissingnessPFilter(options.getMissingNessP()));
		}
		filter.addFilter(new VCFVariantCallRateFilter(options.getCallrateThreshold()));
		filter.addFilter(new VCFVariantImpQualFilter(options.getImputationqualitythreshold(), true));
		filter.addFilter(new VCFVariantMAFFilter(options.getMafthresholdD(), VCFVariantMAFFilter.MODE.OVERALL));
		filter.addFilter(new VCFVariantHWEPFilter(options.getHWEPThreshold(), VCFVariantHWEPFilter.MODE.CONTROLS));
		filter.addFilter(new VCFVariantRegionFilter(regions));
		
		TextFile out = new TextFile(options.getOutputdir() + "-faillist.txt.gz", TextFile.W);
		TextFile out3 = new TextFile(options.getOutputdir() + "-successlist.txt.gz", TextFile.W);
		TextFile out2 = new TextFile(options.getOutputdir() + "-success.txt", TextFile.W);
		TextFile in = new TextFile(options.getVcf(), TextFile.R);
		String ln = in.readLine();
		int ctr = 0;
		int fail = 0;
		int passed = 0;
		System.out.println("Parsing: " + options.getVcf());
		
		while (ln != null) {
			if (!ln.startsWith("#")) {
				String outStr = p(ln, filter, genotypeSamplesWithCovariatesAndDiseaseStatus);
				if (outStr != null) {
					out.writeln(outStr);
					fail++;
				} else {
					VCFVariant v = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
					out3.writeln(v.asSNPFeature().toString());
					passed++;
				}
			}
			ctr++;
			if (ctr % 10 == 0) {
				System.out.print("\r" + ctr + " lines parsed, " + fail + " failed, " + passed + " passed");
			}
			ln = in.readLine();
		}
		out2.writeln(ctr + " lines parsed, " + fail + " failed, " + passed + " passed");
		System.out.println("Done");
		out2.close();
		out.close();
		in.close();
		out3.close();
	}
	
	Feature tyk2 = new Feature(Chromosome.NINETEEN, 10396336, 10628468);
	
	private String p(String s, VCFVariantFilters filters, boolean[] samplesToInclude) {
		
		if (s != null && !s.startsWith("#")) {
			boolean parseln = false;
			if (!filters.hasRegionOrVariantSetFilter()) {
				parseln = true;
			} else {
				String substr = s.substring(0, 200);
				VCFVariant v = new VCFVariant(substr, VCFVariant.PARSE.HEADER);
				parseln = filters.passesRegionOrVariantFilter(v);
				if (!parseln) {
					return "" + v.asSNPFeature().toString();
				}
			}
			
			if (parseln) {
				VCFVariant v = new VCFVariant(s, VCFVariant.PARSE.ALL, samplesToInclude, sampleAnnotation);
				if (v.asFeature().overlaps(tyk2)) {
					if (filters != null) {
						VCFVariantFilters filters2 = new VCFVariantFilters();
						ArrayList<VCFVariantFilter> setFilters = filters.getFilters();
						boolean maffilterset = false;
						for (VCFVariantFilter f : setFilters) {
							if (f instanceof VCFVariantMAFFilter) {
								maffilterset = true;
							} else {
								filters2.add(f);
							}
						}
						if (maffilterset) {
							filters2.addFilter(new VCFVariantMAFFilter(0.005));
						}
						if (!filters2.passesFilters(v)) {
							return "" + v.asSNPFeature().toString();
						} else {
							return null;
						}
					}
				} else if (filters != null) {
					if (!filters.passesFilters(v)) {
						return "" + v.asSNPFeature().toString();
					} else {
						return null;
					}
				}
			}
		}
		return null;
	}
	
	
}
