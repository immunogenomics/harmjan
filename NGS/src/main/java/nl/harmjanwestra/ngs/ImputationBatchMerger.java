package nl.harmjanwestra.ngs;

import java.io.IOException;

/**
 * Created by hwestra on 10/13/15.
 */
public class ImputationBatchMerger {

	public static void main(String[] args) {
		PlayGround g = new PlayGround();
		if (args.length < 3) {
			System.out.println("Usage: prefix outfile nrbatches chr");
		} else {
			String prefix = args[0];
			String outfilename = args[1];
			Integer nrBatches = Integer.parseInt(args[2]);
			Integer chr = Integer.parseInt(args[3]);
			try {
				g.mergeImputedBatches(prefix, outfilename, nrBatches, chr);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

}
