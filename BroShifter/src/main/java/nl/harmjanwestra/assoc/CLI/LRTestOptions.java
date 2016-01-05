package nl.harmjanwestra.assoc.CLI;

import nl.harmjanwestra.broshifter.CLI.MainOptions;
import org.apache.commons.cli.*;
import umcg.genetica.text.Strings;

/**
 * Created by hwestra on 1/4/16.
 */
public class LRTestOptions {

	private static final Options OPTIONS;

	static {
		OPTIONS = MainOptions.OPTIONS;

		Option option = Option.builder("i")
				.hasArg()
				.desc("Input VCF")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("o")
				.hasArg()
				.desc("Output prefix")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("d")
				.hasArg()
				.desc("Disease status file")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("c")
				.hasArg()
				.desc("Covariate file")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("e")
				.hasArg()
				.desc("List of samples to exclude")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("q")
				.hasArg()
				.desc("Imputation quality threshold (default = 0.8)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("n")
				.hasArg()
				.desc("Observe allele at least n times (default == 5)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("t")
				.hasArg()
				.desc("Number of threads")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("f")
				.hasArg()
				.desc("FAM file")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("s")
				.hasArg()
				.desc("Limit to list of SNPs")
				.build();
		OPTIONS.addOption(option);


		option = Option.builder()
				.longOpt("maxiter")
				.hasArg()
				.desc("Run iterative analysis until this number of iters")
				.build();
		OPTIONS.addOption(option);


		option = Option.builder()
				.longOpt("includecov")
				.hasArg()
				.desc("Use only these covariates when testing (point to file with list)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("onlyimputed")
				.hasArg()
				.desc("Test only variants that were imputed")
				.build();
		OPTIONS.addOption(option);
	}

	private String vcf;
	private String outputdir;
	private String diseaseStatusFile;
	private String covariateFile;
	private String samplesToExclude;
	private boolean testVariantsWithoutImputationQuality;
	private double imputationqualitythreshold = 0.8;
	private int minNObservedAllele;
	private String famfile;
	private String snpLimitFile;
	private int maxIter = 1;
	private double mafthresholdD = 0.005;
	private String covariatesToInclude;

	public LRTestOptions(String[] args) {

		boolean run = true;
		try {

			System.out.println(Strings.concat(args, Strings.comma));

			System.out.println("Parsing options");
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);

			if (cmd.hasOption("i")) {
				vcf = cmd.getOptionValue("i");
			} else {
				System.out.println("Please provide input: -i");
				run = false;
			}

			if (cmd.hasOption("o")) {
				outputdir = cmd.getOptionValue("o");
			} else {
				System.out.println("Please provide output: -o");
				run = false;
			}

			if (cmd.hasOption("d")) {
				diseaseStatusFile = cmd.getOptionValue("d");
			} else {
				System.out.println("Please provide input: -d");
				run = false;
			}

			if (cmd.hasOption("c")) {
				covariateFile = cmd.getOptionValue("c");
			} else {
				System.out.println("Please provide input: -c");
				run = false;
			}

			if (cmd.hasOption("e")) {
				samplesToExclude = cmd.getOptionValue("e");
			}

			if (cmd.hasOption("q")) {
				String impq = cmd.getOptionValue("q");
				try {
					imputationqualitythreshold = Double.parseDouble(impq);
					if (imputationqualitythreshold > 1 || imputationqualitythreshold < 0) {
						System.out.println("Please value between 0 and 1 for -q");
						run = false;
					}
				} catch (NumberFormatException e) {
					System.out.println(impq + " is not a decimal/double for option -q");
					run = false;
				}
			}

			if (cmd.hasOption("n")) {
				String nStr = cmd.getOptionValue("n");
				try {
					minNObservedAllele = Integer.parseInt(nStr);
				} catch (NumberFormatException e) {
					System.out.println(nStr + " is not a integer for option -n");
					run = false;
				}
			}

			if (cmd.hasOption("f")) {
				famfile = cmd.getOptionValue("f");
			} else {
				System.out.println("Please provide input: -f");
				run = false;
			}

			if (cmd.hasOption("s")) {
				snpLimitFile = cmd.getOptionValue("s");
			}
			// testVariantsWithoutImputationQuality -- onlyimputed
			// maxiter = maxiter
			// covariatesToInclude = includecov

			if (cmd.hasOption("onlyimputed")) {
				testVariantsWithoutImputationQuality = false;
			} else {
				testVariantsWithoutImputationQuality = true;
			}

			if (cmd.hasOption("maxiter")) {
				String nStr = cmd.getOptionValue("maxiter");
				try {
					maxIter = Integer.parseInt(nStr);
				} catch (NumberFormatException e) {
					System.out.println(nStr + " is not an integer for option --maxiter");
				}
			}

			if (cmd.hasOption("includecov")) {
				covariatesToInclude = cmd.getOptionValue("includecov");
			}


		} catch (ParseException e) {
			e.printStackTrace();

		}

		if (!run) {
			printHelp();
			System.exit(-1);
		}

	}

	public String getVcf() {
		return vcf;
	}

	public String getOutputdir() {
		return outputdir;
	}

	public String getDiseaseStatusFile() {
		return diseaseStatusFile;
	}

	public String getCovariateFile() {
		return covariateFile;
	}

	public String getSamplesToExclude() {
		return samplesToExclude;
	}

	public boolean isTestVariantsWithoutImputationQuality() {
		return testVariantsWithoutImputationQuality;
	}

	public double getImputationqualitythreshold() {
		return imputationqualitythreshold;
	}

	public int getMinNObservedAllele() {
		return minNObservedAllele;
	}

	public String getFamfile() {
		return famfile;
	}

	public String getSnpLimitFile() {
		return snpLimitFile;
	}

	public String getCovariatesToInclude() {
		return covariatesToInclude;
	}

	public int getMaxIter() {
		return maxIter;
	}

	public double getMafthresholdD() {
		return mafthresholdD;
	}

	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		System.exit(-1);
	}
}
