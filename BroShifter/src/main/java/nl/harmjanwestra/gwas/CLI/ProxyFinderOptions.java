package nl.harmjanwestra.gwas.CLI;

import org.apache.commons.cli.*;

/**
 * Created by hwestra on 3/22/16.
 */
public class ProxyFinderOptions {
	public String tabixrefprefix;
	public int windowsize = 1000000;
	public double threshold = 0.8;
	public String snpfile;
	public String output;
	public int nrthreads = 1;

	private static final Options OPTIONS;

	static {
		OPTIONS = new Options();
		Option option = Option.builder().longOpt("proxy").build();
		OPTIONS.addOption(option);

		option = Option.builder("r")
				.longOpt("ref")
				.hasArg()
				.desc("Prefix for tabix file [format /path/to/chr$i.vcf.gz]")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("w")
				.longOpt("windowsize")
				.hasArg()
				.desc("Window size [default 1000000]")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("t")
				.longOpt("threshold")
				.hasArg()
				.desc("R-squared threshold [default: 0.8]")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("i")
				.longOpt("snps")
				.hasArg()
				.desc("SNP file")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("o")
				.longOpt("out")
				.hasArg()
				.desc("Output file")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("threads")
				.hasArg()
				.desc("Nr of threads")
				.build();
		OPTIONS.addOption(option);

	}

	public ProxyFinderOptions(String[] args) {

		boolean run = true;
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);


			if (cmd.hasOption("r")) {
				tabixrefprefix = cmd.getOptionValue("r");
			} else {
				System.out.println("Provide reference");
				run = false;
			}

			if (cmd.hasOption("w")) {
				windowsize = Integer.parseInt(cmd.getOptionValue("w"));
			}

			if (cmd.hasOption("t")) {
				threshold = Double.parseDouble(cmd.getOptionValue("t"));
			}

			if (cmd.hasOption("i")) {
				snpfile = cmd.getOptionValue("i");
			} else {
				System.out.println("Provide input");
				run = false;
			}

			if (cmd.hasOption("o")) {
				output = cmd.getOptionValue("o");
			} else {
				System.out.println("Provide output");
				run = false;
			}

			if (cmd.hasOption("threads")) {
				nrthreads = Integer.parseInt(cmd.getOptionValue("threads"));
			}


		} catch (ParseException e) {
			e.printStackTrace();


		}

		if (!run) {
			printHelp();
			System.exit(-1);
		}
	}

	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		System.exit(-1);
	}

}
