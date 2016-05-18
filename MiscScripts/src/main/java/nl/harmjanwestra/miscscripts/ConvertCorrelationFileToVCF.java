package nl.harmjanwestra.miscscripts;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by hwestra on 5/4/16.
 */
public class ConvertCorrelationFileToVCF {

	public static void main(String[] args) {


		try {

			ConvertCorrelationFileToVCF c = new ConvertCorrelationFileToVCF();
//			String input = "/Data/tmp/2016-05-04/T1D-cosmo-merged.txt";
//			String output = "/Data/tmp/2016-05-04/T1D-cosmo-accurracy-stats.vcf.gz";
//			c.run(input, output);
//
//			input = "/Data/tmp/2016-05-04/T1D-HRC-merged.txt";
//			output = "/Data/tmp/2016-05-04/T1D-HRC-accurracy-stats.vcf.gz";
//			c.run(input, output);

			String input = "/Data/tmp/2016-05-18/T1D-COSMO-merged.txt";
			String output = "/Data/tmp/2016-05-18/T1D-COSMO-merged.vcf.gz";
			c.run(input, output);

			input = "/Data/tmp/2016-05-18/T1D-EUR-merged.txt";
			output = "/Data/tmp/2016-05-18/T1D-EUR-merged.vcf.gz";
			c.run(input, output);

			input = "/Data/tmp/2016-05-18/T1D-HRC-COSMO-merged.txt";
			output = "/Data/tmp/2016-05-18/T1D-HRC-COSMO-merged.vcf.gz";
			c.run(input, output);

		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void run(String in, String out) throws IOException {

		TextFile inf = new TextFile(in, TextFile.R);
		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");

		String[] elems = inf.readLineElems(TextFile.tab);
		while (elems != null) {
			// 1_114301335_rs61817589  T       C,T     0.2981157469717362      1.0     T       C,T     0.25443640329136513     1.0     1       439     1.0     1.0     null    1.0
			String var = elems[0];
			String alleles = elems[6];
			String info = elems[12];
			String[] varElems = var.split("_");
			String[] alleleElems = alleles.split(",");


			String minor = elems[5];
			String maf = elems[7];
			Double mafd = Double.parseDouble(maf);
			if (!minor.equals(alleleElems[0])) {
				mafd = 1 - mafd;
			}

			String infoStr = "INFO=" + info + ";AF=" + mafd;

			String output = varElems[0] + "\t" + varElems[1] + "\t" + varElems[2] + "\t" + alleleElems[0] + "\t" + alleleElems[1] + "\t.\t.\t" + infoStr;
			outf.writeln(output);

			elems = inf.readLineElems(TextFile.tab);
		}

		inf.close();
		outf.close();

	}

}
