package nl.harmjanwestra.finemapping.comparisons;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 3/25/17.
 */
public class CompareVariantLists {

	public static void main(String[] args) {
		String f1 = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-16-SummaryStats/normal/RA-assoc0.3-COSMO-merged.txt.gz";
		String f2 = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/RA-assoc0.3-COSMO-merged.txt.gz";
		String out = "";
		CompareVariantLists l = new CompareVariantLists();
		try {
			l.run(f1,f2,out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String af1, String af2, String out) throws IOException {
		AssociationFile f = new AssociationFile();
		ArrayList<AssociationResult> r1 = f.read(af1);
		ArrayList<AssociationResult> r2 = f.read(af2);

		HashSet<String> vars2 = new HashSet<String>();
		for (AssociationResult r : r2) {
			vars2.add(r.getSnp().toString());
		}

		int nrShared = 0;
		int nrNotFound = 0;
		for (AssociationResult r : r1) {
			String var = r.getSnp().toString();
			if (vars2.contains(var)) {
				nrShared++;
			} else {
				nrNotFound++;
				System.out.println(var + " not found");
			}
		}
		System.out.println("shared:\t" + nrShared);
		System.out.println("not found:\t" + nrNotFound);


	}
}
