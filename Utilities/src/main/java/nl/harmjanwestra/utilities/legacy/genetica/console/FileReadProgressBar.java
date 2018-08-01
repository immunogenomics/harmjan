package nl.harmjanwestra.utilities.legacy.genetica.console;

/**
 * Created by hwestra on 12/15/16.
 */
public class FileReadProgressBar {

	boolean printctr2;
	int ctr = 0;
	int itr = 0;
	int printevery = 1000;
	String message = "lines read";
	int ctr2;

	public FileReadProgressBar(int printEvery, String message) {
		this.printevery = printEvery;
		this.message = message;
	}

	public FileReadProgressBar(int printEvery, boolean printctr2) {
		this.printevery = printEvery;
		this.printctr2 = printctr2;
	}

	public void iteratectr2() {
		ctr2++;
	}

	public void iterate() {
		ctr++;
		if (ctr % printevery == 0) {

			String out = "\r";
			switch (itr) {
				case 0:
					out += "|";
					itr++;
					break;
				case 1:
					out += "/";
					itr++;
					break;
				case 2:
					out += "-";
					itr++;
					break;
				case 3:
					out += "\\";
					itr = 0;
					break;
			}
			if (printctr2) {
				out += " - " + ctr2 + " / " + ctr + " " + message;
			} else {
				out += " - " + ctr + " " + message;
			}

			System.out.print(out);
		}
	}

	public void close() {
		System.out.print("\n");
	}


}
