package nl.harmjanwestra.assoc.CLI;

import nl.harmjanwestra.broshifter.CLI.MainOptions;
import org.apache.commons.cli.*;
import umcg.genetica.text.Strings;

/**
 * Created by hwestra on 1/4/16.
 */
public class LRTestOptions {

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

	public int getThreadnum() {
		return threadnum;
	}

	public String getFamfile() {
		return famfile;
	}

	public String getSnpLimitFile() {
		return snpLimitFile;
	}

	public static Options getOPTIONS() {
		return OPTIONS;
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

	private String vcf;
	private String outputdir;
	private String diseaseStatusFile;
	private String covariateFile;
	private String samplesToExclude;
	private boolean testVariantsWithoutImputationQuality;
	private double imputationqualitythreshold;
	private int minNObservedAllele;
	private int threadnum;
	private String famfile;
	private String snpLimitFile;
	private int maxIter = 5;
	private double mafthresholdD = 0.005;
	private String covariatesToInclude;

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
				.desc("Output dir")
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
				.desc("Imputation quality threshold")
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
	}

	public LRTestOptions(String[] args) {

		try {

			System.out.println(Strings.concat(args, Strings.comma));

			System.out.println("Parsing options");
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);

			if (cmd.hasOption("i")) {
				vcf = cmd.getOptionValue("i");
			} else {
				vcf = null;
			}

			if (cmd.hasOption("o")) {
				outputdir = cmd.getOptionValue("o");
			}

			if (cmd.hasOption("d")) {
				diseaseStatusFile = cmd.getOptionValue("d");
			}

			if (cmd.hasOption("c")) {
				covariateFile = cmd.getOptionValue("c");
			}

			if (cmd.hasOption("e")) {
				samplesToExclude = cmd.getOptionValue("e");
			}

			if (cmd.hasOption("q")) {
				String impq = cmd.getOptionValue("q");
				try {
					imputationqualitythreshold = Double.parseDouble(impq);
				} catch (NumberFormatException e) {
					System.out.println(impq + " is not a number/double for option -q");
				}
			}

			if (cmd.hasOption("n")) {
				String nStr = cmd.getOptionValue("n");
				try {
					minNObservedAllele = Integer.parseInt(nStr);
				} catch (NumberFormatException e) {
					System.out.println(nStr + " is not a number/double for option -n");
				}
			}

			if (cmd.hasOption("t")) {
				String nStr = cmd.getOptionValue("t");
				try {
					threadnum = Integer.parseInt(nStr);
				} catch (NumberFormatException e) {
					System.out.println(nStr + " is not a number/double for option -t");
				}
			}

			if (cmd.hasOption("f")) {
				famfile = cmd.getOptionValue("f");
			}

			if (cmd.hasOption("s")) {
				snpLimitFile = cmd.getOptionValue("s");
			}
		} catch (ParseException e) {
			e.printStackTrace();

		}

	}
}
