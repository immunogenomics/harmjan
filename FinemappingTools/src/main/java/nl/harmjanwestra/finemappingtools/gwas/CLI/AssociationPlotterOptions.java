package nl.harmjanwestra.finemappingtools.gwas.CLI;

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
				.desc("Input path names (comma separated list)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("b")
				.hasArg()
				.desc("Threshold for credible sets (default: 0.9)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("a")
				.hasArg()
				.desc("GTF Annotation path")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("r")
				.hasArg()
				.desc("Region list")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("threshold")
				.hasArg()
				.desc("Single significance threshold for all regions")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("thresholds")
				.hasArg()
				.desc("File with bonferroni threshold per region")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("maxpvalues")
				.hasArg()
				.desc("File with max Pvalue per region")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("ldprefix")
				.hasArg()
				.desc("Tabix ld prefix")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("ldlimit")
				.hasArg()
				.desc("Tabix ld: limit to a list of samples")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("o")
				.hasArg()
				.desc("Output prefix")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("s")
				.hasArg()
				.desc("Sequenced regions path")
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

		option = Option.builder("i2")
				.hasArg()
				.desc("File template for conditional files. The text String ITER will be replaced by an iteration number.")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder().longOpt("thresholdconditional")
				.hasArg()
				.desc("Significance threshold for conditional results")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder().longOpt("thresholdsconditional")
				.hasArg()
				.desc("Significance thresholds for conditional results per region (supply file)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder().longOpt("nriters")
				.hasArg()
				.desc("Number of iterations to plot for conditional analysis")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder().longOpt("maxpvaluesconditional")
				.hasArg()
				.desc("Max p-values for each plot in the conditional analysis")
				.build();
		OPTIONS.addOption(option);


	}

	String associationFiles;
	String associationFileNames;
	String annotationfile;
	String bedregionfile;
	String outputprefix;
	String sequencedRegionsFile;
	boolean plotPosterior;
	private double maxpconditional;
	private Double maxp;
	private double credibleSetThreshold = 0.9;
	private String significanceThresholdFile;
	private String significanceConditionalThresholdFile;
	private String LDPrefix;
	private String LDLimit;
	private double defaultSignificance = 5E-8;
	private String maxPvalueFile;
	private String maxPvalueFileConditional;
	private double defaultSignificanceConditional = 5E-8;
	private int nrIters = 1;
	private String conditionalFiles;

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

			if (cmd.hasOption("i2")) {
				conditionalFiles = cmd.getOptionValue("i2");
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

			if (cmd.hasOption("ldprefix")) {
				LDPrefix = cmd.getOptionValue("ldprefix");
			}

			if (cmd.hasOption("ldlimit")) {
				LDLimit = cmd.getOptionValue("ldlimit");
			}

			if (cmd.hasOption("maxpvalues")) {
				maxPvalueFile = cmd.getOptionValue("maxpvalues");
			}

			if (cmd.hasOption("maxpvaluesconditional")) {
				maxPvalueFileConditional = cmd.getOptionValue("maxpvaluesconditional");
			}


			if (cmd.hasOption("r")) {
				bedregionfile = cmd.getOptionValue("r");
			} else {
				System.out.println("Please provide input: -r");
				run = false;
			}

			if (cmd.hasOption("b")) {
				credibleSetThreshold = Double.parseDouble(cmd.getOptionValue("b"));
			}

			if (cmd.hasOption("o")) {
				outputprefix = cmd.getOptionValue("o");
			} else {
				System.out.println("Please provide input: -o");
				run = false;
			}

			if (cmd.hasOption("s")) {
				sequencedRegionsFile = cmd.getOptionValue("s");
			}

			if (cmd.hasOption("threshold")) {
				defaultSignificance = Double.parseDouble(cmd.getOptionValue("threshold"));
			}

			if (cmd.hasOption("thresholdconditional")) {
				defaultSignificanceConditional = Double.parseDouble(cmd.getOptionValue("thresholdconditional"));
			}


			if (cmd.hasOption("thresholds")) {
				significanceThresholdFile = cmd.getOptionValue("thresholds");
			}

			if (cmd.hasOption("thresholdsconditional")) {
				significanceConditionalThresholdFile = cmd.getOptionValue("thresholdsconditional");
			}

			if (cmd.hasOption("m")) {
				maxp = Double.parseDouble(cmd.getOptionValue("m"));
			}

			if (cmd.hasOption("nriters")) {
				nrIters = Integer.parseInt(cmd.getOptionValue("nriters"));
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

	public Double getMaxp() {
		return maxp;
	}

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

	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		System.exit(-1);
	}

	public double getCredibleSetThreshold() {
		return credibleSetThreshold;
	}

	public String getSignificanceThresholdFile() {

		return significanceThresholdFile;
	}

	public String getLDPrefix() {
		return LDPrefix;
	}

	public String getLDLimit() {
		return LDLimit;
	}

	public double getDefaultSignificance() {
		return defaultSignificance;
	}

	public String getMaxPvalueFile() {
		return maxPvalueFile;
	}

	public String getSignificanceConditionalThresholdFile() {
		return significanceConditionalThresholdFile;
	}

	public String getMaxConditionalPvalueFile() {
		return maxPvalueFileConditional;
	}

	public double getDefaultSignificanceConditional() {
		return defaultSignificanceConditional;
	}

	public int getNrIters() {
		return nrIters;
	}

	public String getConditionalFiles() {
		return conditionalFiles;
	}
}
