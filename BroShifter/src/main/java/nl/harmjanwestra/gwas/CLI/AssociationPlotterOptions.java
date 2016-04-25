package nl.harmjanwestra.gwas.CLI;

import org.apache.commons.cli.*;

/**
 * Created by hwestra on 11/23/15.
 */
public class AssociationPlotterOptions {
	private static final Options OPTIONS;

	static {
		OPTIONS = new Options();
		Option option = Option.builder().longOpt("plotposteriors").build();
		OPTIONS.addOption(option);

		option = Option.builder("i")
				.hasArg()
				.desc("Input files (or comma separated list)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("n")
				.hasArg()
				.desc("Input file names (comma separated list)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("a")
				.hasArg()
				.desc("GTF Annotation file")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("r")
				.hasArg()
				.desc("Region list")
				.build();
		OPTIONS.addOption(option);


		option = Option.builder("o")
				.hasArg()
				.desc("Output prefix")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("s")
				.hasArg()
				.desc("Sequenced regions file")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("m")
				.hasArg()
				.desc("Max p")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("p")
				.desc("Plot posteriors")
				.build();
		OPTIONS.addOption(option);
	}

	private Double maxp;

	public Double getMaxp() {
		return maxp;
	}

	String associationFiles;
	String associationFileNames;
	String annotationfile;
	String bedregionfile;
	String outputprefix;
	String sequencedRegionsFile;
	boolean plotPosterior;

	public String getAssociationFiles() {
		return associationFiles;
	}

	public String getAssociationFileNames() {
		return associationFileNames;
	}

	public String getAnnotationfile() {
		return annotationfile;
	}

	public String getBedregionfile() {
		return bedregionfile;
	}

	public String getOutputprefix() {
		return outputprefix;
	}

	public String getSequencedRegionsFile() {
		return sequencedRegionsFile;
	}

	public boolean isPlotPosterior() {
		return plotPosterior;
	}

	public AssociationPlotterOptions(String[] args) {

		boolean run = true;

		CommandLineParser parser = new DefaultParser();
		final CommandLine cmd;
		try {
			cmd = parser.parse(OPTIONS, args, false);
			if (cmd.hasOption("i")) {
				associationFiles = cmd.getOptionValue("i");
			} else {
				System.out.println("Please provide input: -i");
				run = false;
			}

			if (cmd.hasOption("n")) {
				associationFileNames = cmd.getOptionValue("n");
			} else {
				System.out.println("Please provide input: -n");
				run = false;
			}

			if (cmd.hasOption("a")) {
				annotationfile = cmd.getOptionValue("a");
			} else {
				System.out.println("Please provide input: -a");
				run = false;
			}

			if (cmd.hasOption("r")) {
				bedregionfile = cmd.getOptionValue("r");
			} else {
				System.out.println("Please provide input: -r");
				run = false;
			}

			if (cmd.hasOption("o")) {
				outputprefix = cmd.getOptionValue("o");
			} else {
				System.out.println("Please provide input: -o");
				run = false;
			}

			if (cmd.hasOption("s")) {
				sequencedRegionsFile = cmd.getOptionValue("s");
			} else {
				System.out.println("Please provide input: -s");
				run = false;
			}

			if (cmd.hasOption("m")) {
				maxp = Double.parseDouble(cmd.getOptionValue("m"));
			}

			if (cmd.hasOption("p")) {
				plotPosterior = true;
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
