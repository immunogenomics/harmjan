package nl.harmjanwestra.gwas.CLI;

import org.apache.commons.cli.*;

/**
 * Created by hwestra on 1/4/16.
 */
public class LRTestOptions {

	private static final Options OPTIONS;

	static {
		OPTIONS = new Options();
		Option option = Option.builder().longOpt("gwas").build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("conditions")
				.hasArg()
				.desc("Run second iteration conditional on these variants, in specified order (semicolon separated). Format: chrx-start-rsname;chrx2-start2-rsname2")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("exhaustive")
				.desc("Run pairwise conditional analysisType (on all pairs in regions); pair vs null model that includes only covariates")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("conditional")
				.desc("Run single variant conditional analysisType (on all pairs in regions); null model includes conditional variant")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("haplotype")
				.desc("Determine haplotypes from variants and perform haplotype based test (warning: can only handle a limited number of haplotypes!)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("i")
				.hasArg()
				.desc("Input VCF")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("maf")
				.hasArg()
				.desc("MAF threshold (default: 0.005)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("hapfreq")
				.hasArg()
				.desc("Haplotype frequency threshold (default: 0.01)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("hwep")
				.hasArg()
				.desc("HWEP threshold (default: 10-5)")
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

		option = Option.builder("r")
				.hasArg()
				.desc("Test within these regions")
				.argName(".bed file")
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
				.longOpt("genotypeprobthreshold")
				.desc("Genotype probability threshold for haplotype test")
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
				.longOpt("snplimit")
				.desc("Limit to list of SNPs")
				.build();
		OPTIONS.addOption(option);


		option = Option.builder()
				.longOpt("maxiter")
				.hasArg()
				.desc("Run iterative analysisType until this number of iters")
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
				.desc("Test only variants that were imputed")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("nomissing")
				.desc("Assume there is no missing data (aka calculate the null model only once).")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("testmultiallelicindepently")
				.desc("Test each allele in a multiple allelic variant independently. Does not split variant.")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("splitmultiallelic")
				.desc("Split multi-allelic variants before testing.")
				.build();
		OPTIONS.addOption(option);
	}

	public boolean assumeNoMissingData;
	public boolean testMultiAllelicVariantsIndependently;
	public boolean splitMultiAllelicIntoMultipleVariants;
	private ANALYSIS analysisType = ANALYSIS.NORMAL;
	private String conditional;
	private double HWEPThreshold = 1E-5;
	private boolean haplotypeAnalysis;
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
	private String bedfile;
	private int nrThreads = 1;
	private double haplotypeFrequencyThreshold = 0.01;
	private double genotypeProbThreshold = 0;

	public LRTestOptions(String[] args) {

		boolean run = true;
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);

			if (cmd.hasOption("i")) {
				vcf = cmd.getOptionValue("i");
			} else {
				System.out.println("Please provide input: -i");
				run = false;
			}

			if (cmd.hasOption("r")) {
				bedfile = cmd.getOptionValue("r");
			}

			if (cmd.hasOption("nomissing")) {
				assumeNoMissingData = true;
			}
			if (cmd.hasOption("testmultiallelicindepently")) {
				testMultiAllelicVariantsIndependently = true;
			}

			if (cmd.hasOption("splitmultiallelic")) {
				splitMultiAllelicIntoMultipleVariants = true;
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

			if (cmd.hasOption("maf")) {
				mafthresholdD = Double.parseDouble(cmd.getOptionValue("maf"));
			}

			if (cmd.hasOption("hwep")) {
				HWEPThreshold = Double.parseDouble(cmd.getOptionValue("hwep"));
			}

			if (cmd.hasOption("conditions")) {
				this.conditional = cmd.getOptionValue("conditional");
				this.maxIter = 1;
			}

			if (cmd.hasOption("e")) {
				samplesToExclude = cmd.getOptionValue("e");
			}

			if (cmd.hasOption("genotypeprobthreshold")) {
				genotypeProbThreshold = Double.parseDouble(cmd.getOptionValue("genotypeprobthreshold"));
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

			if (cmd.hasOption("t")) {
				String nStr = cmd.getOptionValue("t");
				try {
					nrThreads = Integer.parseInt(nStr);
				} catch (NumberFormatException e) {
					System.out.println(nStr + " is not a integer for option -t");
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

			if (cmd.hasOption("exhaustive")) {
				analysisType = ANALYSIS.EXHAUSTIVE;
			} else if (cmd.hasOption("conditional")) {
				analysisType = ANALYSIS.CONDITIONAL;
			} else if (cmd.hasOption("haplotype")) {
				analysisType = ANALYSIS.HAPLOTYPE;
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

			if (cmd.hasOption("hapfreq")) {
				haplotypeFrequencyThreshold = Double.parseDouble(cmd.getOptionValue("hapfreq"));
			}


		} catch (ParseException e) {
			e.printStackTrace();

		}

		if (!run) {
			printHelp();
			System.exit(-1);
		}

	}

	public double getHaplotypeFrequencyThreshold() {
		return haplotypeFrequencyThreshold;
	}

	public String getConditional() {
		return conditional;
	}

	public int getNrThreads() {
		return nrThreads;
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

	public void setMaxIter(int maxIter) {
		this.maxIter = maxIter;
	}

	public double getMafthresholdD() {
		return mafthresholdD;
	}

	public String getBedfile() {
		return bedfile;
	}

	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		System.exit(-1);
	}

	public double getHWEPThreshold() {
		return HWEPThreshold;
	}

	public ANALYSIS getAnalysisType() {
		return analysisType;
	}

	public double getGenotypeProbThreshold() {
		return genotypeProbThreshold;
	}

	public void setSplitMultiAllelic(boolean splitMultiAllelic) {
		this.splitMultiAllelicIntoMultipleVariants = splitMultiAllelic;
	}

	public enum ANALYSIS {
		CONDITIONAL,
		EXHAUSTIVE,
		HAPLOTYPE,
		NORMAL
	}


}
