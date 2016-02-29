package nl.harmjanwestra.miscscripts;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

/**
 * Created by Harm-Jan on 02/29/16.
 */
public class EidRename {
	public static void main(String[] args) {

		String key = "D:\\tmp\\2016-02-28\\bro\\ra\\EIDlegend.txt";
		String file = "D:\\tmp\\2016-02-28\\bro\\ra\\ChromHMM-summit150bp-Overall.txt";
		String out = "D:\\tmp\\2016-02-28\\bro\\ra\\ChromHMM-summit150bp-Overall-rewrite.txt";

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
}
