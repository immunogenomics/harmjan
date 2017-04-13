package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

/**
 * Created by Harm-Jan on 02/29/16.
 */
public class EidRename {
	public static void main(String[] args) {

		EidRename r = new EidRename();
		r.chromhmmfilerename();


	}

	public void broshifterfilerename() {

		String key = "D:\\tmp\\2016-02-28\\bro\\ra\\EIDlegend.txt";
		String file = "D:\\tmp\\2016-02-28\\bro\\ra\\ChromHMM-summit150bp-Overall.txt";
		String out = "D:\\tmp\\2016-02-28\\bro\\ra\\ChromHMM-summit150bp-Overall-rewriteMDSMatrix.txt";

		int col = 15;

		try {
			TextFile in1 = new TextFile(key, TextFile.R);
			String[] elems = in1.readLineElems(TextFile.tab);
			HashMap<String, String> map = new HashMap<String, String>();
			while (elems != null) {
				map.put(elems[0], elems[1]);
				elems = in1.readLineElems(TextFile.tab);
			}

			in1.close();

			TextFile outf = new TextFile(out, TextFile.W);

			TextFile in2 = new TextFile(file, TextFile.R);
			outf.writeln(in2.readLine());
			elems = in2.readLineElems(TextFile.tab);

			while (elems != null) {
				String id = elems[col];

				String[] idelems = id.split("_");
				String[] idelems2 = id.split("-");

				String newid = map.get(idelems[0]) + "-" + idelems2[1];
				elems[col] = newid;
				outf.writeln(Strings.concat(elems, Strings.tab));


				elems = in2.readLineElems(TextFile.tab);
			}

			in1.close();

		} catch (IOException e) {


		}

	}

	public void chromhmmfilerename() {
		String list = "/Data/Epigenetics/chromhmm/allChromHMM.txt";
		String key = "/Data/Epigenetics/chromhmm/EIDlegend.txt";
		String out = "/Data/Epigenetics/chromhmm/allChromHMM-rewriteMDSMatrix.txt";
		String path = "/medpop/srlab/external-data/ENCODE/ChromHMM/imputed12marks/jointModel/final/split/";
		String removethis = "_25_imputed12marks_mnemonics-";

		try {
			TextFile in1 = new TextFile(key, TextFile.R);
			String[] elems = in1.readLineElems(TextFile.tab);
			HashMap<String, String> map = new HashMap<String, String>();
			while (elems != null) {
				map.put(elems[0], elems[1]);
				elems = in1.readLineElems(TextFile.tab);
			}

			TextFile outf = new TextFile(out, TextFile.W);

			TextFile in2 = new TextFile(list, TextFile.R);
			String ln = in2.readLine();
			elems = in2.readLineElems(TextFile.tab);

			while (ln != null) {
				String filename = ln.replaceAll(path, "");

				String[] felems = filename.split("_");
				String eid = felems[0];

				String name = map.get(eid);
				String mnemonic = felems[felems.length - 1];
				mnemonic = mnemonic.replaceAll("\\.bed\\.gz", "");
				mnemonic = mnemonic.replaceAll("mnemonics-", "");


				outf.writeln(name + "-" + mnemonic + "\t" + ln);


				ln = in2.readLine();
			}

			outf.close();
			in2.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
