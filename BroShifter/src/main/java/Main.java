import CLI.BroShifterOptions;
import CLI.MainOptions;

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
				System.out.println("meh");
				BroShifter bs = new BroShifter(new BroShifterOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.POSTERIORPVAL)) {

			} else if (options.mode.equals(MainOptions.MODE.PLOT)) {

			} else if (options.mode.equals(MainOptions.MODE.MERGE)) {

			}
		} catch (IOException e){
			e.printStackTrace();
		}

	}


}
