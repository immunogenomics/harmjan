package nl.harmjanwestra.miscscripts;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.Set;

/**
 * Created by hwestra on 2/5/15.
 */
public class FilterSamples {
	public static void main(String[] args) {
		String f1 = "/Data/Projects/2014-FR-Reseq/outliers.txt";
		String f2 = "/Data/Projects/2014-FR-Reseq/SamplesToUse.list";
		String f3 = "/Data/Projects/2014-FR-Reseq/SamplesToUseOutliers.list";
		try {
			TextFile tf = new TextFile(f1, TextFile.R);
			Set<String> selected = tf.readAsSet(0, TextFile.tab);
			tf.close();

			TextFile tf2 = new TextFile(f2, TextFile.R);
			TextFile tf3 = new TextFile(f3, TextFile.W);

			String ln = tf2.readLine();
			while (ln != null) {

				String[] sampleNameElems = ln.split("/");
				String bamfile = sampleNameElems[sampleNameElems.length - 1];
				System.out.println(bamfile);
				String[] bamfileElems = bamfile.split("\\.");
				String sampleName = bamfileElems[0];
				if (selected.contains(sampleName)) {
					tf3.writeln(ln);
				}
				ln = tf2.readLine();
			}
			tf3.close();
			tf2.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
