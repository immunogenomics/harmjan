package nl.harmjanwestra.txtr;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by Harm-Jan on 02/01/16.
 */
public class Main {

	public static void main(String[] args) {
		if(args.length<3){
			System.out.println("filein fileout nrlns");
		}
		try {
			Main main = new Main();
			main.split(args[0], args[1], Integer.parseInt(args[2]));
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void split(String file, String fileout, int lns) throws IOException {

		TextFile in = new TextFile(file, TextFile.R);
		int ctr = 0;
		int ctr2 = 1;
		String ln = in.readLine();
		TextFile out = new TextFile(fileout + "-" + ctr2, TextFile.W);
		while (ln != null) {
			ln = in.readLine();
			out.writeln(ln);
			ctr++;
			if (ctr % lns == 0) {
				out.close();
				out = new TextFile(fileout + "-" + ctr2, TextFile.W);
				ctr2++;
			}
		}
		out.close();

		in.close();
	}

}
