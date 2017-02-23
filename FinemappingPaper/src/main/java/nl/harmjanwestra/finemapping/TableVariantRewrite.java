package nl.harmjanwestra.finemapping;

import nl.harmjanwestra.utilities.enums.Chromosome;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 2/21/17.
 */
public class TableVariantRewrite {

	public static void main(String[] args) {

		String[] tables = new String[]{
				"",
				"",
				"",
		};

		String reference = "";
		int[][] columns = new int[][]{
				new int[]{},
				new int[]{},
				new int[]{},
		};

		TableVariantRewrite v = new TableVariantRewrite();
		for (int i = 0; i < tables.length; i++) {
			String table = tables[i];
			int[] column = columns[i];

			try {
				v.run(table, table + "-rewrite.txt", reference, column);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}


	}

	public static void run(String infile, String outfile, String reference, int[] columns) throws IOException {


		HashSet<String> variantPositions = new HashSet<String>();
		HashSet<Chromosome> chrs = new HashSet<>();

		TextFile tf = new TextFile(infile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			for (int i : columns) {
				String posStr = elems[i];
				String[] strElems = posStr.split("_");
				Chromosome chrobj = Chromosome.parseChr(strElems[0]);
				String finalPosStr = chrobj.toString() + "_" + strElems[1];
				variantPositions.add(finalPosStr);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		HashMap<String, String> strToStr = new HashMap<String, String>();
		for (Chromosome chr : chrs) {
			String refFile = reference.replace("CHR", "" + chr.getNumber());
			TextFile ref = new TextFile(refFile, TextFile.R);
			String ln = ref.readLine();
			while (ln != null) {
				String substr = ln.substring(0, 300);
				String[] substrElems = Strings.tab.split(substr);
				String finalPosStr = chr.toString() + "_" + substrElems[1];
				if (variantPositions.contains(finalPosStr)) {
					if (strToStr.containsKey(finalPosStr)) {
						System.out.println("Multiple variants at position:\t" + finalPosStr + "\tnew: " + substrElems[2] + "\told:" + strToStr.get(finalPosStr));
					} else {
						strToStr.put(finalPosStr, substrElems[2]);
					}
				}
			}
			ref.close();
		}

		TextFile tfout = new TextFile(outfile, TextFile.W);
		tf.open();
		tfout.writeln(tf.readLine());
		elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			for (int i : columns) {
				String posStr = elems[i];
				String[] strElems = posStr.split("_");
				Chromosome chrobj = Chromosome.parseChr(strElems[0]);
				String finalPosStr = chrobj.toString() + "_" + strElems[1];

				String newId = strToStr.get(finalPosStr);
				if (newId == null) {
					System.out.println("ID for " + finalPosStr + " could not be updated");
				} else {
					elems[i] = newId;
				}
			}
			tfout.writeln(Strings.concat(elems, TextFile.tab));
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		tfout.close();


	}

}
