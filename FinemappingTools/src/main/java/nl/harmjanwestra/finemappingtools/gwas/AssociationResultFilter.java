package nl.harmjanwestra.finemappingtools.gwas;

import nl.harmjanwestra.finemappingtools.gwas.CLI.AssociationFilterOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class AssociationResultFilter {
	
	
	public AssociationResultFilter(AssociationFilterOptions options) throws IOException {
		this.filter(options.getInput(), options.getSnploc(), options.getOutputprefix(), options.isPairwise());
	}
	
	public void filter(String in, String snpsToRemove, String out, boolean pairwise) throws IOException {
		
		HashSet<String> combinedIds = new HashSet<String>();
		TextFile tf = new TextFile(snpsToRemove, TextFile.R);
		
		String ln = tf.readLine();
		while (ln != null) {
			combinedIds.add(ln.trim().toLowerCase());
			ln = tf.readLine();
		}
		tf.close();
		System.out.println(combinedIds.size() + " ids to filter read from: " + snpsToRemove);
		
		AssociationFile f = new AssociationFile();
		
		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln(f.getHeader());
		
		TextFile inf = new TextFile(in, TextFile.R);
		inf.readLine(); // skip header
		ln = inf.readLine();
		int nrlines = 0;
		int nrtotal = 0;
		while (ln != null) {
			String[] elems = ln.split("\t");
			String combinedId = elems[3];
			if (!combinedIds.contains(combinedId.trim().toLowerCase())) {
				outf.writeln(ln);
				nrlines++;
			}
			nrtotal++;
			ln = inf.readLine();
		}
		inf.close();
		
		outf.close();
		
		System.out.println(nrlines + " kept \t out of " + nrtotal);
		
	}
}
