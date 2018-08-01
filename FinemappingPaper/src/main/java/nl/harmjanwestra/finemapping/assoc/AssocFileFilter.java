package nl.harmjanwestra.finemapping.assoc;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 10/13/16.
 */
public class AssocFileFilter {

	public static void main(String[] args) {
		try {
			AssocFileFilter f = new AssocFileFilter();
			double maf = 0.01;
			double hwep = 1E-4;
			double info = 0.3;

			String in = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/Normal/META-assoc0.3-COSMO-merged.txt.gz";
			String out = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/META-assoc0.3-COSMO-merged.txt.gz";

			f.filter(in, out, hwep, maf, info);

			in = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/Normal/RA-assoc0.3-COSMO-merged.txt.gz";
			out = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/RA-assoc0.3-COSMO-merged.txt.gz";
			f.filter(in, out, hwep, maf, info);

			in = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/Normal/T1D-assoc0.3-COSMO-merged.txt.gz";
			out = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/T1D-assoc0.3-COSMO-merged.txt.gz";
			f.filter(in, out, hwep, maf, info);

		} catch (IOException e) {
			e.printStackTrace();
		}

	}


	public void filter(String in, String out, double hwep, double maf, double info) throws IOException {
		TextFile outf = new TextFile(out, TextFile.W);
		AssociationFile f = new AssociationFile();
		outf.writeln(f.getHeader());
		ArrayList<AssociationResult> assoc = f.read(in);
		System.out.println(assoc.size());
		int written = 0;
		for (AssociationResult r : assoc) {
			SNPFeature snp = r.getSnp();
			if (snp.getMaf() > maf && snp.getHwep() > hwep && snp.getImputationQualityScore() > info) {
				outf.writeln(r.toString());
				written++;
			}
		}

		System.out.println(written + " totally written.");

		outf.close();

	}

}
