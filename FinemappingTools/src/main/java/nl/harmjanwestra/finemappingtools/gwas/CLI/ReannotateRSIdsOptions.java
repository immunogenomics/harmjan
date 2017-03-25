package nl.harmjanwestra.finemappingtools.gwas.CLI;

import org.apache.commons.cli.*;

/**
 * Created by hwestra on 8/10/16.
 */
public class ReannotateRSIdsOptions {


	private static final Options OPTIONS;

	static {
		OPTIONS = new Options();
		Option option = Option.builder().longOpt("updaters").build();
		OPTIONS.addOption(option);

		option = Option.builder("i")
				.hasArg()
				.desc("Input (can be comma separated)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("r")
				.hasArg()
				.desc("dbSNP VCF file")
				.build();
		OPTIONS.addOption(option);
	}

	public String assoc;
	public String vcf;

	public ReannotateRSIdsOptions(String[] args) {

		boolean run = true;
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);

			if (cmd.hasOption("i")) {
				assoc = cmd.getOptionValue("i");
			} else {
				run = false;
			}
			if (cmd.hasOption("r")) {
				vcf = cmd.getOptionValue("r");
			} else {
				run = false;
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
