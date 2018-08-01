package nl.harmjanwestra.ngs;

import nl.harmjanwestra.utilities.plink.PedAndMapFunctions;

import java.io.IOException;

/**
 * Created by hwestra on 9/28/15.
 */
public class ReplacePEDHeader {

	public static void main(String[] args) {

		try {


			if (args.length < 2) {
				System.out.println("USAGE: ped fam");
			} else {

				String ped = args[0];
				String fam = args[1];

				PedAndMapFunctions p = new PedAndMapFunctions();
				p.updatePedWithFAM(ped, fam);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

	}


}
