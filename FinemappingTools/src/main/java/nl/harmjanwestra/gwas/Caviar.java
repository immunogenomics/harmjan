package nl.harmjanwestra.gwas;

import nl.harmjanwestra.gwas.CLI.CaviarOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by Harm-Jan on 02/28/16.
 */
public class Caviar {

	public Caviar(CaviarOptions options) throws IOException {
		if (options.convert) {
			this.filterassoc(options.assoc, options.variantlist, options.out);
		} else if (options.filter) {
			this.filterCorrelationMatrix(options.assoc, options.matrix, options.out);
		}
	}

	public void filterassoc(String assocfile, String listfile, String out) throws IOException {

		HashMap<String, Integer> variantHash = new HashMap<String, Integer>();
		TextFile tf = new TextFile(listfile, TextFile.R);
		String ln = tf.readLine();
		ArrayList<String> variantnames = new ArrayList<String>();
		int ctr = 0;
		while (ln != null) {

			variantnames.add(ln);
			variantHash.put(ln, ctr);
			ctr++;

			ln = tf.readLine();
		}

		tf.close();

		double[] pvals = new double[variantHash.size()];
		AssociationFile a = new AssociationFile();
		ArrayList<AssociationResult> results = a.read(assocfile);

		for (AssociationResult result : results) {
			String variant = result.getSnp().getChromosome().getNumber() + "-" + result.getSnp().getStart() + "-" + result.getSnp().getName();
			Integer id = variantHash.get(variant);
			if (id != null) {
				pvals[id] = ZScores.pToZ(result.getPval());
				if (result.getBeta()[0] < 0) {
					pvals[id] = -pvals[id];
				}
			}
		}

		TextFile outf = new TextFile(out, TextFile.W);

		for (int i = 0; i < pvals.length; i++) {
			outf.writeln(variantnames.get(i) + "\t" + pvals[i]);
		}
		outf.close();

	}

	public void filterCorrelationMatrix(String caviarout, String correlationmatrix, String out) throws IOException {


		TextFile tf2 = new TextFile(caviarout, TextFile.R);
		ArrayList<Boolean> include = new ArrayList<>();
		String[] elems2 = tf2.readLineElems(TextFile.tab);
		TextFile out1 = new TextFile(out + "-zscores.txt", TextFile.W);
		TextFile out2 = new TextFile(out + "-variantindex.txt", TextFile.W);
		int ctr = 1;
		while (elems2 != null) {
			double d = Double.parseDouble(elems2[1]);
			if (d == 0) {
				include.add(false);
			} else {
				include.add(true);
				out1.writeln(Strings.concat(elems2, Strings.tab));
				out2.writeln(ctr + "\t" + elems2[0]);
			}
			ctr++;
			elems2 = tf2.readLineElems(TextFile.tab);
		}
		out1.close();
		out2.close();
		tf2.close();

		boolean[] includecol = new boolean[include.size()];
		for (int i = 0; i < include.size(); i++) {
			includecol[i] = include.get(i);
		}


		TextFile tf = new TextFile(correlationmatrix, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr2 = 0;
		TextFile out3 = new TextFile(out + "-matrix.txt", TextFile.W);
		while (elems != null) {
			boolean inc = include.get(ctr2);
			if (inc) {
				out3.writeln(Strings.concat(elems, includecol, Strings.tab));
			}
			elems = tf.readLineElems(TextFile.tab);
			ctr2++;
		}
		out3.close();
		tf.close();


	}

	public void reintroduce(String assocfile, String caviaroutput, String out) throws IOException {

	}


}
