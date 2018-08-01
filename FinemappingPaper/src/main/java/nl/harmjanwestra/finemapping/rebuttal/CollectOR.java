package nl.harmjanwestra.finemapping.rebuttal;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

public class CollectOR {
	
	public static void main(String[] args) {
		String file = "C:\\Sync\\Dropbox\\FineMap\\2018-01-Rebuttal\\tables\\listofsnpswithposterior0.2.txt";
		String order = "C:\\Sync\\Dropbox\\FineMap\\2018-01-Rebuttal\\tables\\order.txt";
		String[] datasets = new String[]{
				"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\META-assoc0.3-COSMO-merged-posterior.txt.gz",
				"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
		};
		
		CollectOR c = new CollectOR();
		try {
			c.filter(file, datasets, order,null);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void filter(String file, String[] datasets, String order, String output) throws IOException {
		
		TextFile t = new TextFile(order, TextFile.R);
		ArrayList<String> snps = t.readAsArrayList();
		
		
		
		AssociationResult[][] results = new AssociationResult[snps.size()][datasets.length];
		
		AssociationFile f = new AssociationFile();
		for (int d = 0; d < datasets.length; d++) {
			ArrayList<AssociationResult> r = f.read(datasets[d]);
			for (int a = 0; a < r.size(); a++) {
				for (int s = 0; s < snps.size(); s++) {
					if (snps.get(s).equals(r.get(a).getSnp().getName())) {
						results[s][d] = r.get(a);
					}
				}
			}
		}
		
		for (int s = 0; s < snps.size(); s++) {
			
			String lnout = "";
			for (int d = 0; d < datasets.length; d++) {
				if (results[s][d] != null) {
					String alleles = Strings.concat(results[s][d].getSnp().getAlleles(), Strings.comma);
					double or = results[s][d].getORs()[0][0];
					lnout += "\t" + alleles + "\t" + or;
				} else {
					lnout += "\t-\t-";
				}
			}
			
			System.out.println(snps.get(s).toString() + lnout);
		}
	}
	
}
