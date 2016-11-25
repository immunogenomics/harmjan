import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashSet;

/**
 * Created by hwestra on 11/25/16.
 */
public class Main {

	public static void main(String[] args) {

	}

	public void run(String synfile, String genelist) throws IOException {
		HashSet<String> query = null;

		TextFile tf = new TextFile(genelist, TextFile.R);
		String ln = tf.readLine();
		while (ln != null) {
			query.add(ln.trim());
			ln = tf.readLine();
		}
		tf.close();

		TextFile tf2 = new TextFile(synfile, TextFile.R);

		String ln2 = tf2.readLine();
		while (ln2 != null) {

			/*
			#tax_id GeneID
			Symbol  LocusTag        Synonyms        dbXrefs chromosome      map_location    description     type_of_gene    Symbol_from_nomenclature_authority      Full_name_from_nomenclature_authority   Nomenclature_status     Other_designations      Modification_date
			 */

			ln2 = tf2.readLine();
		}
		tf2.close();

	}
}
