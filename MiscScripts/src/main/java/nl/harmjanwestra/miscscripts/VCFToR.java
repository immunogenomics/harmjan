package nl.harmjanwestra.miscscripts;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.math.matrix2.DoubleMatrixDataset;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by hwestra on 3/23/17.
 */
public class VCFToR {

	public static void main(String[] args) {
		String vcf = "/Data/tmp/sh2b3fix/test.vcf";
		String covariates = "/Data/tmp/sh2b3fix/covarmerged.txtmergedCovariates.txt";
		String pheno = "/Data/tmp/sh2b3fix/covarmerged.txtmergeddisease.txt";
		String out = "/Data/tmp/sh2b3fix/testr";

		try {
			TextFile tf = new TextFile(vcf, TextFile.R);
			VCFGenotypeData d = new VCFGenotypeData(vcf);
			ArrayList<String> samples = d.getSamples();
			d.close();

			System.out.println(samples.size() + " samples");
			String ln = tf.readLine();
			ArrayList<VCFVariant> vars = new ArrayList<>();
			while (ln != null) {
				if (!ln.startsWith("#")) {
					VCFVariant v = new VCFVariant(ln);
					vars.add(v);
				}
				ln = tf.readLine();
			}
			tf.close();

			System.out.println(vars.size() + " variants");
			DoubleMatrixDataset ds = DoubleMatrixDataset.loadDoubleData(covariates);
			TextFile tf2 = new TextFile(pheno, TextFile.R);
			String[] elems = tf2.readLineElems(Strings.tab);
			HashMap<String, String> phenos = new HashMap<>();
			while (elems != null) {
				if (elems.length == 2) {
					phenos.put(elems[0], elems[1]);
				}
				elems = tf2.readLineElems(Strings.tab);
			}
			tf2.close();
			System.out.println(phenos.size() + " phenotypes");

//			TextFile tfout = new TextFile(out + "covarsonly.txt", TextFile.W);
//
//			tfout.close();


			String header = "sample\tpheno\tgt";
			ArrayList<String> covs = (ArrayList<String>) ds.getColObjects();
			for (int q = 0; q < ds.columns(); q++) {
				header += "\t" + covs.get(q);
			}

			for (VCFVariant v : vars) {
				TextFile outf = new TextFile(out + "-" + v.getId() + ".txt", TextFile.W);
				outf.writeln(header);

				DoubleMatrix2D dosages = v.getDosagesAsMatrix2D();

				for (int s = 0; s < samples.size(); s++) {
					double dosage = dosages.getQuick(s, 0);
					if (dosage != -1) {
						String sample = samples.get(s);
						Integer id = (Integer) ds.getHashRows().get(sample);
						String lnout = sample + "\t" + (Double.parseDouble(phenos.get(sample)) - 1) + "\t" + dosage;
						for (int q = 0; q < ds.columns(); q++) {
							lnout += "\t" + ds.getMatrix().getQuick(id, q);
						}
						outf.writeln(lnout);
					}
				}
				outf.close();
			}


		} catch (IOException e) {

		}
	}
}
