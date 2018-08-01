package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by hwestra on 1/13/15.
 */
public class RewriteSampleNames {

	public static void main(String[] args) {
		try {
			System.out.println("Rewriting");
			TextFile tf = new TextFile("/Data/Projects/2014-FR-Reseq/RASequencingSamples.txt", TextFile.R);
			TextFile tfout = new TextFile("/Data/Projects/2014-FR-Reseq/RASequencingSamples-ReWrite.txt", TextFile.W);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {

				Integer s = Integer.parseInt(elems[1]);
				tfout.writeln(elems[0] + "\t" + s);

				elems = tf.readLineElems(TextFile.tab);
			}

			tf.close();
			tfout.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
