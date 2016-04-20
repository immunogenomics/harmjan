package nl.harmjanwestra.gwas;

import nl.harmjanwestra.gwas.CLI.CountVariantsOptions;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by Harm-Jan on 04/20/16.
 */
public class CountVariants {

	private final CountVariantsOptions options;
	public CountVariants(CountVariantsOptions options) throws IOException {
		this.options = options;
		count();
	}

	public void count() throws IOException {

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> bedRegions = reader.readAsList(options.bedfile);

		TextFile listIn = new TextFile(options.input, TextFile.R);
		String[] listArr = listIn.readAsArray();
		listIn.close();

		int nrVariantsTotal = 0;
		int nrVariantsInRegions = 0;
		int nrVariantsInRegionsWithProperImpQual = 0;
		int nrVariantsInRegionsWithProperMAF = 0;

		for (String file : listArr) {
			TextFile tf = new TextFile(file, TextFile.R);

			tf.readLine();
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				// SNP     Chr     Pos     ImputationQual  MAF     OverlapOK       MAFOk   ImpQualOK

				// rs77716349      1       197752167       0.0     null    true    null    false
				Boolean overlapok = false;
				Boolean mafOK = Boolean.parseBoolean(elems[6]);
				Boolean impQualOK = Boolean.parseBoolean(elems[7]);

				Feature f = new Feature();
				f.setChromosome(Chromosome.parseChr(elems[1]));
				Integer pos = Integer.parseInt(elems[2]);
				f.setStart(pos);
				f.setStop(pos);

				for (Feature r : bedRegions) {
					if (r.overlaps(f)) {
						overlapok = true;
						break;
					}
				}

				nrVariantsTotal++;
				if (overlapok) {
					nrVariantsInRegions++;
					if (impQualOK) {
						nrVariantsInRegionsWithProperImpQual++;
						if (mafOK) {
							nrVariantsInRegionsWithProperMAF++;
						}
					}
				}


				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
		}


		System.out.println(nrVariantsTotal + "\tTotal variants");
		System.out.println(nrVariantsInRegions + "\tVariants in regions");
		System.out.println(nrVariantsInRegionsWithProperImpQual + "\tVariants with proper impqual");
		System.out.println(nrVariantsInRegionsWithProperMAF + "\tVariant with proper impqual and maf");
	}

	

}
