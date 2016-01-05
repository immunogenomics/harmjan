package nl.harmjanwestra.broshifter;


import nl.harmjanwestra.assoc.CLI.LRTestOptions;
import nl.harmjanwestra.assoc.LRTest;
import nl.harmjanwestra.broshifter.CLI.BroShifterOptions;
import nl.harmjanwestra.broshifter.CLI.MainOptions;

import java.io.IOException;

/**
 * Created by hwestra on 11/23/15.
 */
public class Main {


	public static void main(String[] args) {

		try {
			MainOptions options = new MainOptions(args);
			if (options.mode.equals(MainOptions.MODE.NA)) {
				System.out.println("Please specify a mode");
			} else if (options.mode.equals(MainOptions.MODE.BROSHIFTER)) {
				new BroShifter(new BroShifterOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.POSTERIORPVAL)) {

			} else if (options.mode.equals(MainOptions.MODE.PLOT)) {

			} else if (options.mode.equals(MainOptions.MODE.MERGE)) {

			} else if (options.mode.equals(MainOptions.MODE.ASSOC)) {
				LRTest test = new LRTest(new LRTestOptions(args));
				test.test();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

	}


}
