package nl.harmjanwestra.finemapping.rebuttal2;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;

import java.io.IOException;
import java.util.ArrayList;

public class CountVariantsInRegions {
	
	
	public static void main(String[] args) {
		
		
		String[] assocfiles = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				"d:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Okada\\RA.OKADA.gz",
				"D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
				"d:/Sync/SyncThing/Postdoc/2016-03-RAT1D-Finemapping/Data/ImmunoBase/hg19_gwas_ic_t1d_onengut_meta_4_19_1.tab.gz",
				"D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\ImmunoBase\\hg19_gwas_ic_t1d_onengut_cc_4_19_1.tab.gz"
		};
		String regionfile = "d:/Sync/SyncThing/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		
		CountVariantsInRegions r = new CountVariantsInRegions();
		for (int a = 0; a < assocfiles.length; a++) {
			try {
				r.countAssocs(assocfiles[a], regionfile);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	public void countAssocs(String assocFile, String regionfile) throws IOException {
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(regionfile);
		AssociationFile f = new AssociationFile();
		
		ArrayList<AssociationResult> assocs = f.readRegions(assocFile, regions);
		System.out.println(assocs.size());
		
	}
	
}
