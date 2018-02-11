package nl.harmjanwestra.finemapping.rebuttal;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class CompareHRCvs1KG {
	
	
	public static void main(String[] args) {
		String disk = "d:";
		String[] hrcfiles = new String[]{
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\hrcimputedassoc\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\hrcimputedassoc\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
		};
		String[] kgfiles = new String[]{
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
			
		};
		String[] out = new String[]{
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\hrcimputedassoc\\comp\\RA.txt",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\hrcimputedassoc\\comp\\T1D.txt",
		};
		String bedregions = disk + "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		
		String tabix = "";
		String tabixSamples = "";
		
		
		CompareHRCvs1KG k = new CompareHRCvs1KG();
		try {
			k.run(hrcfiles, kgfiles, out, bedregions, tabix, tabixSamples);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private void run(String[] hrcfiles, String[] kgfiles, String[] out, String bedfile, String tabix, String tabixSamples) throws IOException {
		
		BedFileReader f = new BedFileReader();
		AssociationFile areader = new AssociationFile();
		ArrayList<Feature> regions = f.readAsList(bedfile);
		// load the 1kg results
		for (int i = 0; i < hrcfiles.length; i++) {
			
			ArrayList<AssociationResult> hrcassoc = areader.read(hrcfiles[i]);
			ArrayList<AssociationResult> kgassoc = areader.read(kgfiles[i]);
			HashSet<String> allsnps = new HashSet<String>();
			for (int q = 0; q < regions.size(); q++) {
				Feature region = regions.get(q);
				
				for (AssociationResult r : kgassoc) {
					if (r.getSnp().overlaps(region)) {
						String snpstr = r.getSnp().getChromosome().toString() + "_" + r.getSnp().getStart();
						allsnps.add(snpstr);
					}
				}
			}
			System.out.println(allsnps.size() + "\t" + kgfiles[i]);
			
			TextFile tfout = new TextFile(out[i], TextFile.W);
			tfout.writeln(areader.getHeader());
			// compare assocs
			for (int q = 0; q < regions.size(); q++) {
				Feature region = regions.get(q);
				HashSet<String> snps = new HashSet<String>();
				for (AssociationResult r : kgassoc) {
					if (r.getSnp().overlaps(region)) {
						String snpstr = r.getSnp().getChromosome().toString() + "_" + r.getSnp().getStart();
						snps.add(snpstr);
					}
				}
				
				for (AssociationResult r : hrcassoc) {
					
					if (r.getSnp().overlaps(region)) {
						if (r.getSnp().getName().equals("rs563485051")) {
							System.out.println("Found it:");
						}
						String snpstr = r.getSnp().getChromosome().toString() + "_" + r.getSnp().getStart();
						if (!snps.contains(snpstr)) {
							String outln = r.toString();
							tfout.writeln(outln);
						}
					}
				}
				
			}
			tfout.close();
		}
	}
	
}
