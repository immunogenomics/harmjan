package nl.harmjanwestra.broshifter;


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
			System.out.println(options.mode);
			if (options.mode.equals(MainOptions.MODE.NA)) {
				options.printHelp();
			} else if (options.mode.equals(MainOptions.MODE.BROSHIFTER)) {
				BroShifter bs = new BroShifter(new BroShifterOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.POSTERIORPVAL)) {

			} else if (options.mode.equals(MainOptions.MODE.PLOT)) {

			} else if (options.mode.equals(MainOptions.MODE.MERGE)) {

			} else if (options.mode.equals(MainOptions.MODE.ASSOC)) {

			}
		} catch (IOException e){
			e.printStackTrace();
		}

	}


}
