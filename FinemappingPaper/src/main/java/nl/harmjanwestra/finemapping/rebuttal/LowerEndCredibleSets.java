package nl.harmjanwestra.finemapping.rebuttal;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

public class LowerEndCredibleSets {
	public static void main(String[] args) {
		
		String[] assoc = new String[]{
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\META-assoc0.3-COSMO-merged-posterior.txt.gz"
		};
		String[] regions = new String[]{
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-significantregions-75e7.bed",
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-significantregions-75e7.bed",
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\META-assoc0.3-COSMO-merged-significantregions-75e7.bed"
		};
		
		String[] output = new String[]{
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\lowerendoutputcrediblesets\\RA.txt",
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\lowerendoutputcrediblesets\\T1D.txt",
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\lowerendoutputcrediblesets\\META.txt"
		};
		LowerEndCredibleSets c = new LowerEndCredibleSets();
		try {
			c.run(assoc, regions, output);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	
	public void run(String[] assoc, String[] regionfiles, String[] output) throws IOException {
		for (int i = 0; i < assoc.length; i++) {
			int[] bins = new int[100];
			AssociationFile f = new AssociationFile();
			ArrayList<AssociationResult> assocs = f.read(assoc[i]);
			
			BedFileReader reader = new BedFileReader();
			ArrayList<Feature> regions = reader.readAsList(regionfiles[i]);
			ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
			int nrSig = 0;
			int nrCs = 0;
			for (int r = 0; r < regions.size(); r++) {
				ArrayList<AssociationResult> regionassocs = getAssocs(regions.get(r), assocs);
				if (regionassocs != null) {
					nrSig++;
					ArrayList<AssociationResult> credibleset = abp.createCredibleSet(regionassocs, 0.95);
					if (credibleset.size() <= 10) {
						nrCs++;
						// bin that mofo
						for (AssociationResult cr : credibleset) {
							double posterior = cr.getPosterior();
							if (posterior > 0.2) {
								posterior = 0.2;
							}
							double perc = posterior / 0.2;
							int bin = (int) Math.floor(bins.length * perc);
							if (bin < 0) {
								bin = 0;
							}
							if (bin > bins.length - 1) {
								bin = bins.length - 1;
							}
							bins[bin]++;
						}
					}
				}
			}
			
			TextFile tf = new TextFile(output[i], TextFile.W);
			tf.writeln("Posterior\tN");
			for (int j = 0; j < bins.length; j++) {
				double perc = (0.2 / bins.length) * j;
				tf.writeln(perc + "\t" + bins[j]);
			}
			tf.close();
			System.out.println(assoc[i] + "\tnrcs<10: " + nrCs + "\tsig: " + nrSig);
		}
	}
	
	private ArrayList<AssociationResult> getAssocs(Feature feature, ArrayList<AssociationResult> assocs) {
		ArrayList<AssociationResult> output = new ArrayList<>();
		boolean sig = false;
		for (AssociationResult r : assocs) {
			if (r.getSnp().overlaps(feature)) {
				if (r.getPval() < 7.5E-7) {
					sig = true;
				}
				output.add(r);
			}
		}
		if (sig) {
			return output;
		} else {
			return null;
		}
	}
}
