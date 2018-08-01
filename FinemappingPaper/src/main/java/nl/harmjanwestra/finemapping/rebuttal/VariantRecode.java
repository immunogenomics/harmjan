package nl.harmjanwestra.finemapping.rebuttal;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class VariantRecode {
	
	public static void main(String[] args) {
		String[] assoc = new String[]{
				"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
				"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\META-assoc0.3-COSMO-merged-posterior.txt.gz"
		};
		
		String listin = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\2018-01-16-ListOfVariants.txt";
		String listout = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\2018-01-16-ListOfVariants-associds.txt";
		VariantRecode r = new VariantRecode();
		try {
			r.recodeRsToAssocId(listin, assoc, listout);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void recodeRsToAssocId(String in, String[] in2, String out) throws IOException {
		TextFile tf = new TextFile(in, TextFile.R);
		String ln = tf.readLine();
		ArrayList<String> list = new ArrayList<String>();
		while (ln != null) {
			list.add(ln);
			ln = tf.readLine();
		}
		tf.close();
		
		HashMap<String, String> namemap = new HashMap<String, String>();
		
		AssociationFile f = new AssociationFile();
		for (String in2f : in2) {
			ArrayList<AssociationResult> results = f.read(in2f);
			for (AssociationResult r : results) {
				namemap.put(r.getSnp().getName(), r.getSnp().toString());
			}
		}
		
		TextFile outf = new TextFile(out, TextFile.W);
		for (String s : list) {
			outf.writeln(s + "\t" + namemap.get(s));
		}
		outf.close();
	}
}
