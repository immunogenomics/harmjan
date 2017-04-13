package nl.harmjanwestra.finemapping.genotypesandimputation;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by hwestra on 2/21/17.
 */
public class Table1 {

	public static void main(String[] args) {
		String[] files = new String[]{
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/NormalHWEP1e4/RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/NormalHWEP1e4/T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/NormalHWEP1e4/META-assoc0.3-COSMO-merged-posterior.txt.gz"
		};

		try {
			String query = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/NormalHWEP1e4/table1variants.txt";
			Table1 t = new Table1();
			String outf = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/NormalHWEP1e4/table1variants-out.txt";
			t.run(files, query, new SELECT[]{SELECT.PVAL, SELECT.POSTERIOR}, outf);


			String[] filestyk2 = new String[]{
					"/Data/Projects/2016-Finemapping/genotypes/RA/out-assoc/RA-assoc0.3-COSMO-chr19-iter1-gwas-1.txt",
					"/Data/Projects/2016-Finemapping/genotypes/T1D/out-TYK2/T1D-assoc0.3-COSMO-chr19-gwas-1.txt",
					"/Data/Projects/2016-Finemapping/genotypes/META/out-assoc/RA-assoc0.3-COSMO-chr19-iter1-gwas-0.txt",
			};
			outf = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/NormalHWEP1e4/table1variants-out-tyk2sec.txt";
			t.run(filestyk2, query, new SELECT[]{SELECT.PVAL, SELECT.POSTERIOR}, outf);

			String[] filestyk2post = new String[]{
					"/Data/Projects/2016-Finemapping/genotypes/RA/out-assoc/RA-assoc0.3-COSMO-chr19-gwas-0-posterior.txt",
					"/Data/Projects/2016-Finemapping/genotypes/T1D/out-TYK2/T1D-assoc0.3-COSMO-chr19-gwas-0-posterior.txt",
					"/Data/Projects/2016-Finemapping/genotypes/META/out-assoc/RA-assoc0.3-COSMO-chr19-iter1-gwas-0-posterior.txt",
			};
			outf = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/NormalHWEP1e4/table1variants-out-tyk2post.txt";
			t.run(filestyk2post, query, new SELECT[]{SELECT.PVAL, SELECT.POSTERIOR}, outf);


			String[] secondaryFiles = new String[]{
					"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/Conditional/RA-assoc0.3-COSMO-gwas-1-merged.txt.gz",
					"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/Conditional/T1D-assoc0.3-COSMO-gwas-1-merged.txt.gz",
					"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/Conditional/META-assoc0.3-COSMO-gwas-1-merged.txt.gz",
			};
			outf = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/NormalHWEP1e4/table1variants-out-second.txt";
			t.run(secondaryFiles, query, new SELECT[]{SELECT.PVAL, SELECT.POSTERIOR}, outf);

			String[] tertiaryFiles = new String[]{
					"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/Conditional/RA-assoc0.3-COSMO-gwas-2-merged.txt.gz",
					"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/Conditional/T1D-assoc0.3-COSMO-gwas-2-merged.txt.gz",
					"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/Conditional/META-assoc0.3-COSMO-gwas-2-merged.txt.gz",
			};
			outf = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/NormalHWEP1e4/table1variants-out-tert.txt";
			t.run(tertiaryFiles, query, new SELECT[]{SELECT.PVAL, SELECT.POSTERIOR}, outf);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	enum SELECT {
		PVAL,
		OR,
		POSTERIOR,
		BETA
	}

	public void run(String[] files, String queryfile, SELECT[] queue, String outf) throws IOException {
		AssociationFile f = new AssociationFile();


		ArrayList<String> query = new ArrayList<String>();
		TextFile querytf = new TextFile(queryfile, TextFile.R);
		String ln = querytf.readLine();
		int ctr = 0;
		HashMap<String, AssociationResult[]> resultmap = new HashMap<>();

		while (ln != null) {
			String q = ln.trim();
			if (q.length() >= 0) {
				query.add(q);
				resultmap.put(q, new AssociationResult[files.length]);
			}
			ln = querytf.readLine();
		}
		querytf.close();

		for (int i = 0; i < files.length; i++) {
			ArrayList<AssociationResult> results = f.read(files[i]);
			for (AssociationResult r : results) {
				String id = r.getSnp().getName();
//				System.out.println(id);
//				System.exit(-1);
				if (resultmap.containsKey(id)) {
					AssociationResult[] r2 = resultmap.get(id);
					r2[i] = r;
				}
			}
		}

		String header = "id";
		for (SELECT s : queue) {
			for (int q = 0; q < files.length; q++) {
				if (s == SELECT.OR) {
					header += "\tOR";
				} else if (s == SELECT.PVAL) {
					header += "\tPVAL";
				} else if (s == SELECT.POSTERIOR) {
					header += "\tPOSTERIOR";
				} else if (s == SELECT.BETA) {
					header += "\tBETA (SE)";
				}
			}
		}
		TextFile outtf = new TextFile(outf, TextFile.W);
		outtf.writeln(header);
		for (int v = 0; v < query.size(); v++) {
			String id = query.get(v);
			AssociationResult[] r2 = resultmap.get(id);
			String out = id;

			for (SELECT s : queue) {
				for (int q = 0; q < r2.length; q++) {
					if (r2[q] == null) {
						out += "\t-";
						System.out.println("Could not find " + id);
					} else {
						if (s == SELECT.OR) {
							out += "\t" + r2[q].getORs()[0];
						} else if (s == SELECT.PVAL) {
							out += "\t" + r2[q].getPval();
						} else if (s == SELECT.POSTERIOR) {
							out += "\t" + r2[q].getPosterior();
						} else if (s == SELECT.BETA) {
							out += "\t" + r2[q].getBeta() + " (" + r2[q].getSe() + ")";
						}
//						System.out.println("Found " + id);
					}
				}
			}
			outtf.writeln(out);
		}
		outtf.close();


	}
}
