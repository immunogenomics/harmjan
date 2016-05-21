package nl.harmjanwestra.vcfutils.plots;

/**
 * Created by hwestra on 5/20/16.
 */
public class PlotterImpQual {

	public static void main(String[] args) {


	}

	public void run() {

		String[] files = new String[]{
				"/Data/tmp/2016-05-20/T1D-EUR-stats.vcf.gz",
				"/Data/tmp/2016-05-20/T1D-COSMO-stats.vcf.gz",
				"/Data/tmp/2016-05-20/T1D-HRC-COSMO.vcf.gz"
		};
		String[] labels = new String[]{"EUR", "COSMO", "HRC-COSMO"};
		String variantsOnIC = "/Data/tmp/2016-05-20/T1D-recode-stats.vcf.gz";


		// plot 1: impqual scatterplot
		// plot 2: boxplot of maf vs impqual
		// plot 3: ??



	}
}
