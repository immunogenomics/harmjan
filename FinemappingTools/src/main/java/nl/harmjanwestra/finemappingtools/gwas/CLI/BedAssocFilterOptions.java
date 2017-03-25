package nl.harmjanwestra.finemappingtools.gwas.CLI;

import org.apache.commons.cli.*;

/**
 * Created by hwestra on 2/24/16.
 */
public class BedAssocFilterOptions {
	public String outfile;
	public String regionFile;
	public String assocFile;
	public String thresholdfile;
	public double threshold = 10E-6;

	private static final Options OPTIONS;

	static {
		OPTIONS = new Options();
		Option option = Option.builder().longOpt("bedfilter").build();
		OPTIONS.addOption(option);

		option = Option.builder("b")
				.hasArg()
				.desc("Bed path with regions")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("p")
				.desc("Print top associations per region")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("i")
				.hasArg()
				.desc("Association path")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("t")
				.hasArg()
				.desc("Threshold (default: 10-6)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("q")
				.hasArg()
				.desc("Thresholds per locus")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("o")
				.hasArg()
				.desc("Output")
				.build();
		OPTIONS.addOption(option);
	}

	public boolean printTopAssocPerRegion;

	public BedAssocFilterOptions(String[] args) {

		boolean run = true;
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);


			if (cmd.hasOption("b")) {
				this.regionFile = cmd.getOptionValue("b");
			} else {
				run = false;
			}

			if (cmd.hasOption("i")) {
				this.assocFile = cmd.getOptionValue("i");
			} else {
				run = false;
			}
			if (cmd.hasOption("p")) {
				printTopAssocPerRegion = true;
			}

			if (cmd.hasOption("q")) {
				thresholdfile = cmd.getOptionValue("q");
			}

			if (cmd.hasOption("t")) {
				try {
					this.threshold = Double.parseDouble(cmd.getOptionValue("t"));
				} catch (NumberFormatException e) {
					System.out.println(cmd.getOptionValue("t") + " is not a double for -t");
					run = false;
				}
			}

			if (cmd.hasOption("o")) {
				this.outfile = cmd.getOptionValue("o");
			} else {
				run = false;
			}

			if (!run) {
				printHelp();
				System.exit(-1);
			}

		} catch (ParseException e) {
			e.printStackTrace();
		}
	}

	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		System.exit(-1);
	}
}
