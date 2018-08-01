package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 6/30/16.
 */
public class CompareLists {

	public static void main(String[] args) {
		try {

			CompareLists c = new CompareLists();
			c.run("/Data/tmp/2016-06-29-quals/samples/RA-assoc0.3-COSMO-chr8-samplelist.txt",
					"/Data/tmp/2016-06-29-quals/samples/T1D-assoc0.3-COSMO-chr8-samplelist.txt");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String file1, String file2) throws IOException {
		TextFile tf = new TextFile(file1, TextFile.R);
		ArrayList<String> samples1 = tf.readAsArrayList();
		tf.close();

		TextFile tf2 = new TextFile(file2, TextFile.R);
		ArrayList<String> samples2 = tf2.readAsArrayList();
		tf.close();

		HashSet<String> samples1hash = new HashSet<String>();
		samples1hash.addAll(samples1);

		int shared = 0;
		for (String s : samples2) {
			if (samples1hash.contains(s)) {
				shared++;
			}
		}
		System.out.println(shared);
	}
}
