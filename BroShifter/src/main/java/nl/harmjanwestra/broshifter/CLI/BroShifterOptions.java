package nl.harmjanwestra.broshifter.CLI;


import org.apache.commons.cli.*;
import umcg.genetica.text.Strings;


/**
 * Created by hwestra on 11/23/15.
 */
public class BroShifterOptions {

	private static final Options OPTIONS;

	static {
		OPTIONS = new Options();
		Option option = Option.builder().longOpt("broshifter").build();
		OPTIONS.addOption(option);
		option = Option.builder().longOpt("plotoverlap").build();
		OPTIONS.addOption(option);

		option = Option.builder("r")
				.hasArg()
				.desc("Input regions in bedfile format")
				.longOpt("regions")

				.build();
		OPTIONS.addOption(option);

		option = Option.builder("i")
				.hasArg()
				.desc("Input posterior p-value file (gwas-file format)")
				.longOpt("posteriors")

				.build();
		OPTIONS.addOption(option);

		option = Option.builder("a")
				.hasArg()
				.desc("List of annotations to test. One line per annotation file. Annotations can be in .xls or .bed file, and may be gzipped")
				.longOpt("annotations")

				.build();
		OPTIONS.addOption(option);

		option = Option.builder("p")
				.hasArg()
				.desc("Number of permutations to run [default = 10000]")
				.longOpt("iterations")
				.build();
		OPTIONS.addOption(option);


		option = Option.builder("o")
				.hasArg()
				.desc("Output file location")
				.longOpt("out")

				.build();
		OPTIONS.addOption(option);

		option = Option.builder("c")
				.desc("Run conditional analysis?")
				.longOpt("conditional")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("s")
				.desc("Use peak centers to test for enrichment")
				.longOpt("summit")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("e")
				.hasArg()
				.desc("Number of bases to extend around peak center [default = 150]. If using a weight, this acts as the max distance from peak center.")
				.longOpt("peakextend")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Trim regions to fit around variants. [default = false]")
				.longOpt("trim")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("x")
				.hasArg()
				.desc("After trimming, extend the region with a fixed number of basepairs [default = 100]")
				.longOpt("regionextend")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("t")
				.hasArg()
				.desc("Number of threads to use [default = 1]")
				.longOpt("threads")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Force the use of the set number of threads")
				.longOpt("force-threads")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Gene annotation GTF")
				.longOpt("gtf")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Annotation to plot is a continuous trait (bedgraph).")
				.longOpt("continuous")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Make overlap matrix")
				.longOpt("matrix")
				.build();
		OPTIONS.addOption(option);


		option = Option.builder("w")
				.desc("Use a weight for distance. [none|linear|sqrt|inverse|exp|hoverd]. Default: none")
				.longOpt("weight")
				.hasArg()
				.build();
		OPTIONS.addOption(option);
	}

	public String regionFile;
	public String posteriorFile;
	public String listOfAnnotations;
	public int nrIterations = 10000;
	public String outfile;
	public boolean conditional = false;
	public boolean usePeakCenter = false;
	public int bpToExtendAnnotation = 150;
	public boolean trimRegions = false;
	public int defaultRegionExtend = 100;
	public int nrThreads = 1;
	public double credibleSetThreshold = 0.95;
	public DISTANCEWEIGHT distanceweight = DISTANCEWEIGHT.NONE;
	public String geneAnnotationFile;
	// TODO: allow code to determine mean + sd annotation size
	// set some parameters using those stats
	public int maxAllowedDistance = 150;
	public boolean overlapmatrix;
	public boolean continuous = false;

	public BroShifterOptions(String[] args) {
		try {

			System.out.println(Strings.concat(args, Strings.comma));

			System.out.println("Parsing options");
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);

			boolean run = true;

			if (cmd.hasOption("regions")) {
				regionFile = cmd.getOptionValue("regions");
			} else {
				System.out.println("Path to regions file not provided");
				run = false;
			}

			if (cmd.hasOption("posteriors")) {
				posteriorFile = cmd.getOptionValue("posteriors");
			} else {
				System.out.println("Path to posteriors not provided");
				run = false;
			}

			if (cmd.hasOption("gtf")) {
				geneAnnotationFile = cmd.getOptionValue("gtf");
			}

			if (cmd.hasOption("continuous")) {
				continuous = true;
			}

			if (cmd.hasOption("annotations")) {
				listOfAnnotations = cmd.getOptionValue("annotations");
			} else {
				System.out.println("Path to annotations file not provided");
				run = false;
			}

			if (cmd.hasOption("out")) {
				outfile = cmd.getOptionValue("out");
			} else {
				System.out.println("Path to out file not provided");
				run = false;
			}

			if (cmd.hasOption("iterations")) {
				try {
					nrIterations = Integer.parseInt(cmd.getOptionValue("iterations"));
				} catch (NumberFormatException e) {
					System.out.println(cmd.getOptionValue("iterations") + " provided for number of iterations is not an integer");
					run = false;
				}
			}

			if (cmd.hasOption("conditional")) {
				conditional = true;
			}

			if (cmd.hasOption("matrix")) {
				overlapmatrix = true;
			}

			if (cmd.hasOption("summit")) {
				usePeakCenter = true;
			}

			if (cmd.hasOption("peakextend")) {
				try {
					bpToExtendAnnotation = Integer.parseInt(cmd.getOptionValue("peakextend"));

				} catch (NumberFormatException e) {
					System.out.println(cmd.getOptionValue("peakextend") + " provided for peak extension is not an integer");
					run = false;
				}
			}
			this.maxAllowedDistance = bpToExtendAnnotation;

			if (cmd.hasOption("trim")) {
				trimRegions = true;
			}

			if (cmd.hasOption("regionextend")) {
				try {
					defaultRegionExtend = Integer.parseInt(cmd.getOptionValue("regionextend"));
				} catch (NumberFormatException e) {
					System.out.println(cmd.getOptionValue("regionextend") + " provided for region extension is not an integer");
					run = false;
				}
			}


			if (cmd.hasOption("threads")) {
				try {
					nrThreads = Integer.parseInt(cmd.getOptionValue("threads"));
				} catch (NumberFormatException e) {
					System.out.println(cmd.getOptionValue("threads") + " provided for number of threads is not an integer");
					run = false;
				}
				int cores = Runtime.getRuntime().availableProcessors();

				boolean forcethreads = false;
				if (cmd.hasOption("force-threads")) {
					forcethreads = true;
				}

				if (nrThreads > cores && !forcethreads) {
					System.out.println("Warning: " + nrThreads + " is a bigger number of threads than there are cores in this machine. Are you sure? Use --force-threads to force");
					run = false;
				}
			}

			if (cmd.hasOption("weight")) {
				String val = cmd.getOptionValue("weight");

				// none|linear|sqrt|inverse|exp
				if (val.equals("none")) {
					distanceweight = DISTANCEWEIGHT.NONE;
				} else if (val.equals("linear")) {
					distanceweight = DISTANCEWEIGHT.LINEAR;
				} else if (val.equals("sqrt")) {
					distanceweight = DISTANCEWEIGHT.SQUAREROOT;
				} else if (val.equals("inverse")) {
					distanceweight = DISTANCEWEIGHT.INVERSE;
				} else if (val.equals("exp")) {
					distanceweight = DISTANCEWEIGHT.EXPONENT;
				} else if (val.equals("hoverd")) {
					distanceweight = DISTANCEWEIGHT.HEIGHTOVERDISTANCE;
				}
			}

			if (!run) {
				printHelp();
			}
		} catch (ParseException ex) {
			System.out.println(ex.getMessage());
			;
			printHelp();
		}
	}

	public BroShifterOptions() {

	}


	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		System.exit(-1);
	}


	public enum HEIGHTWEIGHT {
		NONE,
		LINEAR,
		SQUAREROOT
	}

	public enum DISTANCEWEIGHT {
		NONE,
		LINEAR,
		INVERSE,
		SQUAREROOT,
		EXPONENT,
		HEIGHTOVERDISTANCE
	}

}
