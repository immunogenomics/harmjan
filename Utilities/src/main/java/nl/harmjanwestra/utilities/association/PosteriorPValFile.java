package nl.harmjanwestra.utilities.association;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 12/3/15.
 */
public class PosteriorPValFile {

	public ArrayList<AssociationResult> readVariantPValues(String pvaluefile, Feature region) throws IOException {

		// Chr	Pos	Id	CombinedId	Beta	Se	PVal	Posterior
		TextFile tf = new TextFile(pvaluefile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		ArrayList<AssociationResult> output = new ArrayList<AssociationResult>();
		while (elems != null) {

			if (elems.length > 7) {
				AssociationResult result = new AssociationResult();


				Chromosome chr = Chromosome.parseChr(elems[0]);
				Integer pos = Integer.parseInt(elems[1]);
				String id = elems[2];

				Double pval = Double.parseDouble(elems[6]);
				Double posterior = Double.parseDouble(elems[7]);

				SNPFeature snp = new SNPFeature(chr, pos, pos);
				snp.setName(id);
				result.setSnp(snp);
				result.setPosterior(posterior);
				result.setPval(pval);


				output.add(result);

			}


			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return output;
	}
}
