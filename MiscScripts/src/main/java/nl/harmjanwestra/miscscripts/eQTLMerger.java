package nl.harmjanwestra.miscscripts;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by hwestra on 5/22/15.
 */
public class eQTLMerger {

	public static void main(String[] args) {
		eQTLMerger m = new eQTLMerger();

	}

	public eQTLMerger() {

		try {
			for (int perm = 0; perm < 11; perm++) {

				String outfilename = "/broad/hptmp/hwestra/eqtlmeta/output/merged/PermutedEQTLsPermutationRound" + perm + ".txt.gz";
				if (perm == 0) {
					outfilename = "/broad/hptmp/hwestra/eqtlmeta/output/merged/eQTLs.txt.gz";
				}
				TextFile mergedOut = new TextFile(outfilename, TextFile.W);

				System.out.println("out: " + outfilename);

				for (int run = 1; run < 11; run++) {



					String tfName = "/broad/hptmp/hwestra/eqtlmeta/output/snpset-x0" + run + "/PermutedEQTLsPermutationRound" + perm + ".txt.gz";
					if(perm == 0){
						tfName = "/broad/hptmp/hwestra/eqtlmeta/output/snpset-x0" + run + "/eQTLs.txt.gz";
					}

					TextFile tfin = new TextFile(tfName, TextFile.R);
					System.out.println("in: " + tfName);

					String header = tfin.readLine();

					if (run == 1) {
						mergedOut.writeln(header);
					}

					String ln = tfin.readLine();
					while (ln != null) {
						mergedOut.writeln(ln);
						ln = tfin.readLine();
					}
					tfin.close();

				}
				System.out.println();

				mergedOut.close();
			}

		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
