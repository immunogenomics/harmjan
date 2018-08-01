package nl.harmjanwestra.ngs;

import nl.harmjanwestra.utilities.genotypes.GenotypeTools;
import nl.harmjanwestra.utilities.plink.PedAndMapFunctions;
import nl.harmjanwestra.utilities.vcf.VCFFunctions;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;

import java.io.IOException;

/**
 * Created by hwestra on 10/13/15.
 */
public class HJGen {
	public static void main(String[] args) {

		PedAndMapFunctions pedAndMapFunctions = new PedAndMapFunctions();

		try {

			String plinkDataset = "/Sync/OneDrive/Genotypes/genotypes";
			GenotypeTools t = new GenotypeTools();
			String map = plinkDataset + ".map";
//			String refmap = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC-RsIDsFromT1DStudy.map";
//			String chrupdate = plinkDataset + "-raciChrNames.map";
//			pedAndMapFunctions.rewriteMapFileChromosomeNames(refmap, map, chrupdate);

//			String rsnameref = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC-RsIDsFromT1DStudy.map";
//			String rsupdate = plinkDataset + "-raciChrNames-updatedRS.map";
//			pedAndMapFunctions.updateMapFileRsIdsUsingMapFile(rsnameref, chrupdate, rsupdate);
//		Gpio.delete(chrupdate);

//			String rsupdatededup = plinkDataset + "-raciChrNames-updatedRS-dedup.map";
//			pedAndMapFunctions.deduplicateMAP(rsupdate, rsupdatededup);
//		Gpio.delete(rsupdate);
//
			String bedout = plinkDataset + "-raciChrNames-dedup.bed";
			pedAndMapFunctions.rewriteMapToBed(map, bedout);
//

			// liftover
			String lifted = plinkDataset + "-lifted.bed";
			String unlifted = plinkDataset + "-unlifted.bed";
//

			ProcessBuilder pb = new ProcessBuilder("/Data/Projects/2014-FR-Reseq/ImmunoChip/liftOver", bedout,
					"/Data/Projects/2014-FR-Reseq/ImmunoChip/hg18ToHg19.over.chain.gz", lifted, unlifted);
//			t.run(pb);
			System.out.println("Lifted over: " + lifted + " | " + unlifted);

			String hg19map = plinkDataset + "-hg19.map";
			pedAndMapFunctions.convertPostLiftOverMAP(map, hg19map, lifted);

			String hg19mapupd = plinkDataset + "-hg19-updRS.map";
			String dbsnpvcf = "/Data/Projects/2014-FR-Reseq/GATKResources/dbsnp_138.b37.vcf.gz";
			pedAndMapFunctions.updateRSNames(dbsnpvcf, hg19map, hg19mapupd);

			String hg19mapdedup = plinkDataset + "-hg19-updRS-dedup.map";
			pedAndMapFunctions.deduplicateMAP(hg19mapupd, hg19mapdedup);


			Gpio.copyFile(plinkDataset + ".map", plinkDataset + ".hg18map");
			Gpio.copyFile(hg19mapdedup, plinkDataset + ".map");

			VCFFunctions vcfFunctions = new VCFFunctions();

			String ped = plinkDataset;
			String vcf = plinkDataset + ".vcf";
			vcfFunctions.convertPEDToVCF(ped, vcf, false);

			Gpio.copyFile(plinkDataset + ".hg18map", plinkDataset + ".map");

			vcfFunctions.splitPerChromosome(vcf, plinkDataset);

		} catch (IOException e) {
			e.printStackTrace();

		}
	}
}
