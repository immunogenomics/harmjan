package nl.harmjanwestra.vcfutils;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.oxford.HapSample;
import nl.harmjanwestra.utilities.vcf.VCFFunctions;
import org.apache.commons.cli.*;

import java.io.IOException;
import java.util.Random;

/**
 * Created by hwestra on 11/24/15.
 */
public class Main {

	private static Options OPTIONS;

	static {
		OPTIONS = new Options();

		Option option;

		option = Option.builder()
				.desc("Batch splitter")
				.longOpt("splitbatch")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Filter for other alleles than ones containing A,C,T, or G")
				.longOpt("filteralleles")
				.build();
		OPTIONS.addOption(option);


		option = Option.builder()
				.desc("Replace sample names")
				.longOpt("samplereplace")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Strip the info field")
				.longOpt("stripinfo")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Replace sample names from plink VCF")
				.longOpt("samplereplaceplink")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Proxy finder")
				.longOpt("proxy")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Filter overlapping variants")
				.longOpt("filtervariantoverlap")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Concatenate variants for shared samples")
				.longOpt("concat")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Concatenate correlation files")
				.longOpt("concatcorrelation")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Concatenate variants for shared samples for a list of files")
				.longOpt("concatlist")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Sample filter (using a list of samples)")
				.longOpt("filtersample")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Sample filter (overlap between two VCF files)")
				.longOpt("filtersampleoverlap")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Detect genotyping sample mix-ups")
				.longOpt("mixups")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Plot r-squared")
				.longOpt("plotrsq")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Plot correlations")
				.longOpt("plotcor")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Update rs names in a VCF")
				.longOpt("updaters")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("DB SNP VCF file")
				.longOpt("dbsnp")
				.hasArg()
				.build();
		OPTIONS.addOption(option);


		option = Option.builder()
				.desc("Filter variants on maf,callrate,genotypingquals and allelic depth")
				.longOpt("filtergenotype")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Filter for bed regions")
				.longOpt("filterregions")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Reintroduce unimputed variants (can also be used to mergecheese two VCF files)")
				.longOpt("reintroduce")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Merge two VCF files (with non-overlapping samples) as a big chunk of swiss cheese; with holes that is.")
				.longOpt("mergecheese")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Merge overlapping variants only. Do strand checks etc.")
				.longOpt("merge")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Match variants to reference VCF")
				.longOpt("match")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Match variants to reference VCF")
				.longOpt("matchusingsummarystats")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Correlate imputation output (between two files) ")
				.longOpt("correlatevcf")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Generate Random seed")
				.longOpt("rng")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Dedup (very naive)")
				.longOpt("dedup")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Replace missing codes")
				.longOpt("rmc")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Convert Hap/sample")
				.longOpt("hapsample")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Variant sampler")
				.longOpt("sample")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Variant correlation matrix (within a region)")
				.longOpt("ldmatrix")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Use dprime (in stead of rsquared) for --ldmatrix")
				.longOpt("dprime")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Use correlation (in stead of dprime or rsquared) for --ldmatrix")
				.longOpt("usecorrelation")
				.build();
		OPTIONS.addOption(option);


		option = Option.builder()
				.desc("Chr splitter")
				.longOpt("splitchr")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Merge imputation batches")
				.longOpt("mergeimputation")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Merge imputation batches")
				.longOpt("checkimputation")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Make pseudo controls")
				.longOpt("pseudo")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Combine correlation results")
				.longOpt("combinecorrelation")
				.build();
		OPTIONS.addOption(option);


		option = Option.builder()
				.desc("RSquared threshold")
				.argName("double")
				.longOpt("rsquared")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Chromosome")
				.argName("int")
				.longOpt("chr")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Callrate threshold")
				.argName("double")
				.longOpt("callrate")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Genotype qual threshold")
				.argName("double")
				.longOpt("gqual")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Read depth threshold")
				.argName("double")
				.longOpt("readdepth")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Allelic balance")
				.argName("double")
				.longOpt("allelicbalance")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Maf threshold")
				.argName("double")
				.longOpt("maf")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("i")
				.desc("Input VCF")
				.argName("file")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Reference name")
				.argName("string")
				.longOpt("refname")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Number of iterations")
				.argName("int")
				.longOpt("nriter")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Number of batches")
				.argName("int")
				.longOpt("nrbatches")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("i2")
				.desc("Input VCF #2")
				.argName("file")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("f")
				.argName("file")
				.hasArg()
				.desc("Input family information (PED or FAM file)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("b")
				.hasArg()
				.desc("Region BED file")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("l")
				.hasArg()
				.desc("List file")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("r")
				.hasArg()
				.desc("Reference VCF")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("o")
				.hasArg()
				.desc("Output prefix")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("s")
				.hasArg()
				.desc("Random number generator seed")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("countoverlap")
				.desc("Count overlapping variants between two vcfs")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.hasArg()
				.longOpt("rsquared")
				.desc("RSquared threshold")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("listsamples")
				.desc("List samples in VCF")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("list")
				.desc("Specify that input is a list.")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("exclude")
				.desc("Switch mode of filter")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.hasArg()
				.longOpt("separator")
				.desc("VCF Genotype allele separator")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("skipheader")
				.desc("Skip printing the header (for --ldmatrix)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.hasArg()
				.longOpt("window")
				.desc("Window size")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.hasArg()
				.longOpt("threads")
				.desc("Number of threads to use (for --proxy)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.hasArg()
				.longOpt("tabix")
				.desc("Tabix file name prefix")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("p")
				.argName("double")
				.hasArg()
				.desc("Sampling percentage for sampling")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("n")
				.hasArg()
				.argName("int")
				.desc("Batch size for splitting")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("vcfsort")
				.hasArg()
				.argName("path")
				.desc("Location of vcfsort")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("mac")
				.desc("Is machine a mac machine (default: no)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("keepoverlapping")
				.desc("Keep overlapping variants or samples (for --merge or --filtersample; default: no)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("updateimputationquals")
				.desc("Correlate imputation probs with best guesses as a measure for imputation quality.")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("summarize")
				.desc("Calculate summary statistics. ")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("summarize2")
				.desc("Calculate summary statistics (2). ")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("summarize3")
				.desc("Calculate AF, AN and AC for each variant in input. Output as VCF, without sample data.")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("summarize3compare")
				.desc("Compare summarized VCF files.")
				.build();
		OPTIONS.addOption(option);

	}

	public static void main(String[] args) {

		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, true);

			boolean run = true;
			String ref = null;
			String input = null;
			String out = null;
			double perc = 1;


			if (cmd.hasOption("rng")) {
				Random random = new Random();

				System.out.println("Here's your random number :)");
				long seed = random.nextLong();
				System.out.println(seed);
				System.exit(-1);
			}


			if (!cmd.hasOption("i")) {
				run = false;
				System.out.println("Please provide input vcf");
			} else {
				input = cmd.getOptionValue("i");
			}

			if (cmd.hasOption("listsamples") && run) {
				VCFListSamples s = new VCFListSamples();
				s.run(input);
				System.exit(-1);
			}

			if (!cmd.hasOption("o")) {
				System.out.println("Please provide output prefix");
				run = false;
			} else {
				out = cmd.getOptionValue("o");
			}

//			if (cmd.hasOption("proxy")) {
//
//				ProxyFinder finder = new ProxyFinder();
//
//				String tabix = null;
//				int windowsize = 1000000;
//				double threshold = 0.8;
//				int nrthreads = 1;
//
//
//				if (cmd.hasOption("rsquared")) {
//					threshold = Double.parseDouble(cmd.getOptionValue("rsquared"));
//				}
//
//				if (cmd.hasOption("window")) {
//					windowsize = Integer.parseInt(cmd.getOptionValue("window"));
//				}
//
//				if (cmd.hasOption("threads")) {
//					nrthreads = Integer.parseInt(cmd.getOptionValue("threads"));
//				}
//
//				if (cmd.hasOption("tabix")) {
//					tabix = cmd.getOptionValue("tabix");
//				} else {
//					run = false;
//				}
//
//				if (run) {
//					finder.find(tabix, windowsize, threshold, input, out, nrthreads);
//				}
//
//			} else
			if (cmd.hasOption("stripinfo")) {
				if (cmd.hasOption("i")) {
					StripInfo s = new StripInfo();
					s.strip(cmd.getOptionValue("i"));
				}
			} else if (cmd.hasOption("countoverlap")) {

				CompareOverlappingVariants q = new CompareOverlappingVariants();
				if (cmd.hasOption("i2")) {
					q.run(input, cmd.getOptionValue("i2"));
				} else {
					System.out.println("use -i2 with --countoverlap");
				}

			} else if (cmd.hasOption("updateimputationquals")) {

				VCFCorrelator correlator = new VCFCorrelator();
				boolean infoscore = true;
				correlator.updateVCFInfoScore(input, out, infoscore);

			} else if (cmd.hasOption("filtersampleoverlap")) {

				VCFSampleFilter filter = new VCFSampleFilter();
				if (cmd.hasOption("i2")) {


					boolean keep = false;
					if (cmd.hasOption("keepoverlapping")) {
						keep = true;
					}

					filter.filteroverlapping(input, cmd.getOptionValue("i2"), out, keep);

				} else {
					System.out.println("use -i2 with --filtersampleoverlap");
				}

			} else if (cmd.hasOption("concatcorrelation")) {

				CorrelationResultCombiner combiner = new CorrelationResultCombiner();
				String[] files = input.split(",");
				combiner.concat(files, out);

			} else if (cmd.hasOption("filteralleles")) {

				VCFFilter filter = new VCFFilter();

				filter.filterNonACTGVariants(input, out);


			} else if (cmd.hasOption("concat")) {
				VCFMerger merger = new VCFMerger();
				if (cmd.hasOption("i2")) {

					merger.concatenate(input, cmd.getOptionValue("i2"), out);

				} else {
					System.out.println("Use -i2 with --concat");
				}
			} else if (cmd.hasOption("concatlist")) {
				VCFMerger merger = new VCFMerger();

				merger.concatenate(input, out);

			} else if (cmd.hasOption("filtervariantoverlap")) {
				VCFRemoveOverlappingVariants t = new VCFRemoveOverlappingVariants();
				if (cmd.hasOption("i2")) {

					t.remove(input, cmd.getOptionValue("i2"), out);

				} else {
					System.out.println("Use -i2 for --filtervariantoverlap");
				}
			} else if (cmd.hasOption("mixups")) {

				MixupTest test = new MixupTest();

				if (cmd.hasOption("i2")) {
					test.test(input, cmd.getOptionValue("i2"), out);
				} else {
					System.out.println("Use -i2 for --mixups");
				}


			} else if (cmd.hasOption("samplereplace")) {

				if (cmd.hasOption("l")) {
					VCFSampleNameReplace vpl = new VCFSampleNameReplace();
					vpl.replace(cmd.getOptionValue("l"), input, out);
				} else {
					System.out.println("Use -l with --samplereplace");
				}

			} else if (cmd.hasOption("samplereplaceplink")) {
				if (cmd.hasOption("f")) {
					VCFSampleNameReplace vpl = new VCFSampleNameReplace();
					vpl.replacePlinkDups(cmd.getOptionValue("f"), input, out);
				} else {
					System.out.println("Use -f with --samplereplaceplink");
				}

			} else if (cmd.hasOption("filtersample")) {

				VCFSampleFilter filter = new VCFSampleFilter();
				String samplefile = null;
				if (cmd.hasOption("l")) {
					samplefile = cmd.getOptionValue("l");
					boolean keep = false;
					if (cmd.hasOption("keepoverlapping")) {
						keep = true;
					}
					filter.filter(input, out, samplefile, keep);
				} else {
					System.out.println("Use params -i -o and -l");
				}


			} else if (cmd.hasOption("updaters")) {

				VCFVariantRSNameUpdater updatr = new VCFVariantRSNameUpdater();
				if (cmd.hasOption("dbsnp")) {
					updatr.updateRSNames(cmd.getOptionValue("dbsnp"), input, out);
				} else {
					System.out.println("Use --dbsnp, -i and -o with --updaters");
				}


			} else if (cmd.hasOption("filtergenotype")) {


				VCFFilter filter = new VCFFilter();
				Integer gqual = null;
				Integer readdepth = null;
				Double allelicBalance = null;
				double maf = 0;
				double callrate = 0;
				boolean onlyautosomes = false;
				if (cmd.hasOption("readdepth")) {
					readdepth = Integer.parseInt(cmd.getOptionValue("readdepth"));
				}
				if (cmd.hasOption("gqual")) {
					gqual = Integer.parseInt(cmd.getOptionValue("gqual"));
				}

				if (cmd.hasOption("allelicbalance")) {
					allelicBalance = Double.parseDouble("allelicbalance");
				}

				if (cmd.hasOption("maf")) {
					maf = Double.parseDouble(cmd.getOptionValue("maf"));
				}

				if (cmd.hasOption("callrate")) {
					callrate = Double.parseDouble(cmd.getOptionValue("callrate"));
				}

				if (cmd.hasOption("autosomes")) {
					onlyautosomes = true;
				}


				filter.filter(input, out, maf, callrate, readdepth, gqual, allelicBalance, onlyautosomes);

			} else if (cmd.hasOption("plotrsq")) {

				RSquaredPlot plot = new RSquaredPlot();


				double threshold = 0.8;
				if (cmd.hasOption("rsquared")) {
					threshold = Double.parseDouble(cmd.getOptionValue("rsquared"));
				}

				double mafthreshold = 0.0;
				if (cmd.hasOption("maf")) {
					mafthreshold = Double.parseDouble(cmd.getOptionValue("maf"));
				}

				String bedfile = null;
				if (cmd.hasOption("b")) {
					bedfile = cmd.getOptionValue("b");
				}
				plot.plot(input, out, threshold, mafthreshold, bedfile);


			} else if (cmd.hasOption("plotcor")) {
				CorrelationResultPlotter plotter = new CorrelationResultPlotter();


				if (cmd.hasOption("refname")) {
					plotter.run(input, cmd.getOptionValue("refname"), out);
				} else {
					plotter.run(input, out);
				}
			} else if (cmd.hasOption("ldmatrix")) {
				VariantLDMatrix mat = new VariantLDMatrix();
				if (cmd.hasOption("b")) {

					boolean dprime = false;
					boolean usecorrelation = false;
					boolean printheader = true;
					if (cmd.hasOption("usecorrelation")) {
						usecorrelation = true;
					} else if (cmd.hasOption("dprime")) {
						dprime = true;
					}
					if (cmd.hasOption("skipheader")) {
						printheader = false;
					}

					mat.correlate(input, cmd.getOptionValue("b"), out, printheader, dprime, usecorrelation);

				} else {
					System.out.println("Please use -b with --ldmatrix");
				}
			} else if (cmd.hasOption("correlatevcf")) {

				String v2 = null;
				String li = null;

				if (cmd.hasOption("i2")) {
					v2 = cmd.getOptionValue("i2");
				} else {
					System.err.println("Provide -i2 with --correlatevcf");
					run = false;
				}

				if (cmd.hasOption("l")) {
					li = cmd.getOptionValue("l");
				}
//				else {
//					System.err.println("Provide -l with --correlatevcf");
//					run = false;
//				}
				if (run) {
					VCFCorrelator c = new VCFCorrelator();

					c.run(input, v2, li, out);

				} else {
					printHelp();
				}
			} else if (cmd.hasOption("summarize")) {
				VCFFunctions f = new VCFFunctions();

				f.summarizeVCF(input, out);

			} else if (cmd.hasOption("summarize2")) {
				VCFFunctions f = new VCFFunctions();

				f.determineVCFSummaryStatistics(input, out);

			} else if (cmd.hasOption("summarize3")) {
				VCFVariantStats stats = new VCFVariantStats();
				stats.run(input, out);

			} else if (cmd.hasOption("summarize3compare")) {
				VCFVariantStats stats = new VCFVariantStats();
				if (cmd.hasOption("b") && cmd.hasOption("i2")) {
					boolean filterExcludes = false;
					if (cmd.hasOption("exclude")) {
						filterExcludes = true;
					}
					String filter = null;
					if (cmd.hasOption("l")) {
						filter = cmd.getOptionValue("l");
					}
					stats.compare(input, cmd.getOptionValue("i2"), out, cmd.getOptionValue("b"), filter, filterExcludes);
				} else {
					System.out.println("Provide -b and -i2 additional to -i and -o. [and optionally -l for a list of variants]");
				}
			} else if (cmd.hasOption("rmc")) {
				VCFFunctions f = new VCFFunctions();

				f.replaceMissingCodes(input, out);

			} else if (cmd.hasOption("reintroduce")) {
				VCFMerger merger = new VCFMerger();
				boolean linux = true;
				String v2 = "";
				String vcfsort = "";

				if (cmd.hasOption("mac")) {
					linux = false;
				}
				if (cmd.hasOption("i2")) {
					v2 = cmd.getOptionValue("i2");
				} else {
					System.err.println("Provide -i2 with --reintroduce");
					run = false;
				}
				if (cmd.hasOption("vcfsort")) {
					vcfsort = cmd.getOptionValue("vcfsort");
				} else {
					System.out.println("Please supply --vcfsort");
				}
				if (run) {

					merger.reintroducteNonImputedVariants(input, v2, out, linux, vcfsort);

				} else {
					printHelp();
				}
			} else if (cmd.hasOption("dedup")) {
				VCFFunctions f = new VCFFunctions();

				f.removedups(input, out);

			} else if (cmd.hasOption("combinecorrelation")) {
				CorrelationResultCombiner combine = new CorrelationResultCombiner();

				String refname = null;
				int nriter = 0;
				int nrbatches = 0;
				run = true;
				if (cmd.hasOption("refname")) {
					refname = cmd.getOptionValue("refname");
				} else {
					run = false;
					System.out.println("Provide --refname");
				}

				if (cmd.hasOption("nriter")) {
					nriter = Integer.parseInt(cmd.getOptionValue("nriter"));
				} else {
					run = false;
					System.out.println("Provide --nriter");
				}

				if (cmd.hasOption("nrbatches")) {
					nrbatches = Integer.parseInt(cmd.getOptionValue("nrbatches"));
				} else {
					run = false;
					System.out.println("Provide --nrbatches");
				}

				if (run) {

					combine.run(input, out, refname, nriter, nrbatches);

				}
			} else if (cmd.hasOption("merge")) {
				VCFMerger merger = new VCFMerger();
				boolean linux = true;
				if (cmd.hasOption("mac")) {
					linux = false;
				}

				int chrint = 1;
				if (cmd.hasOption("chr")) {
					chrint = Integer.parseInt(cmd.getOptionValue("chr"));
				} else {
					System.out.println("Please use --chr with --merge");
					run = false;
				}

				String vcfsort = null;
				if (cmd.hasOption("vcfsort")) {
					vcfsort = cmd.getOptionValue("vcfsort");
				} else {
					System.out.println("Please use --vcfsort with --merge");
					run = false;
				}

				String vcf2 = null;
				if (cmd.hasOption("i2")) {
					vcf2 = cmd.getOptionValue("i2");
				} else {
					System.out.println("Please use -i2 with --merge");
					run = false;
				}

				String sep = "/";
				if (cmd.hasOption("separator")) {
					sep = cmd.getOptionValue("separator");
				}

				boolean keepoverlapping = false;
				if (cmd.hasOption("keepoverlapping")) {
					keepoverlapping = true;
				}

				if (run) {
					merger.mergeAndIntersect(linux, chrint, vcfsort, vcf2, input, out, keepoverlapping, sep);
				}
			} else if (cmd.hasOption("match")) {
				boolean linux = true;
				if (cmd.hasOption("mac")) {
					linux = false;
				}

//				int chrint = 1;
//				if (cmd.hasOption("chr")) {
//					chrint = Integer.parseInt(cmd.getOptionValue("chr"));
//				} else {
//					System.out.println("Please use --chr with --match");
//					run = false;
//				}

				String vcfsort = null;
				if (cmd.hasOption("vcfsort")) {
					vcfsort = cmd.getOptionValue("vcfsort");
				} else {
					System.out.println("Please use --vcfsort with --match");
					run = false;
				}

				String vcf2 = null;
				if (cmd.hasOption("i2")) {
					vcf2 = cmd.getOptionValue("i2");
				} else {
					System.out.println("Please use -i2 with --march");
					run = false;
				}

				String sep = "/";
				if (cmd.hasOption("separator")) {
					sep = cmd.getOptionValue("separator");
				}

				if (run) {
					VCFMatcher matcher = new VCFMatcher();
					matcher.run(vcf2, input, out, linux, vcfsort);
				}
			} else if (cmd.hasOption("matchusingsummarystats")) {
				boolean linux = true;
				if (cmd.hasOption("mac")) {
					linux = false;
				}

//				int chrint = 1;
//				if (cmd.hasOption("chr")) {
//					chrint = Integer.parseInt(cmd.getOptionValue("chr"));
//				} else {
//					System.out.println("Please use --chr with --match");
//					run = false;
//				}

				String vcfsort = null;
				if (cmd.hasOption("vcfsort")) {
					vcfsort = cmd.getOptionValue("vcfsort");
				} else {
					System.out.println("Please use --vcfsort with --match");
					run = false;
				}

				String vcf2 = null;
				if (cmd.hasOption("i2")) {
					vcf2 = cmd.getOptionValue("i2");
				} else {
					System.out.println("Please use -i2 with --march");
					run = false;
				}

				String sep = "/";
				if (cmd.hasOption("separator")) {
					sep = cmd.getOptionValue("separator");
				}

				if (run) {
					VCFMatcher matcher = new VCFMatcher();
					matcher.runReferenceOnlyHasVariantPositions(vcf2, input, out, linux, vcfsort);
				}
			} else if (cmd.hasOption("mergeimputation")) {
				VCFMerger merger = new VCFMerger();
				int nrbatches = 0;
				if (cmd.hasOption("nrbatches")) {
					nrbatches = Integer.parseInt(cmd.getOptionValue("nrbatches"));
				} else {
					run = false;
					System.out.println("Provide --nrbatches");
				}

				String variantList = null;
				if (cmd.hasOption("l")) {
					variantList = cmd.getOptionValue("l");
				}

				if (run) {
					merger.mergeImputationBatches(input, out, variantList, nrbatches);
				}
			} else if (cmd.hasOption("checkimputation")) {
				VCFMerger merger = new VCFMerger();
				int nrbatches = 0;
				if (cmd.hasOption("nrbatches")) {
					nrbatches = Integer.parseInt(cmd.getOptionValue("nrbatches"));
				} else {
					run = false;
					System.out.println("Provide --nrbatches");
				}

				String variantList = null;
				if (cmd.hasOption("l")) {
					variantList = cmd.getOptionValue("l");
				}

				if (run) {
					merger.checkImputationBatches(input, out, variantList, nrbatches);
				}
			} else if (cmd.hasOption("mergecheese")) {
				VCFMerger merger = new VCFMerger();
				String vcf1 = input;
				System.out.println("Merging cheese :)");
				if (cmd.hasOption("i2")) {
					String vcf2 = cmd.getOptionValue("i2");
					try {
						merger.merge(vcf1, vcf2, out);
					} catch (IOException e) {
						e.printStackTrace();
					}
				} else {
					System.out.println("Please provide -i2 for --mergecheese");
					printHelp();
				}
			} else if (cmd.hasOption("hapsample")) {
				HapSample s = new HapSample();
				boolean linux = true;

				if (cmd.hasOption("mac")) {
					linux = false;
				}
				if (cmd.hasOption("vcfsort")) {
					try {
						s.covertHapsSampleToVCF(input, out, linux, cmd.getOptionValue("vcfsort"));
					} catch (IOException e) {
						e.printStackTrace();
					}
				} else {
					System.out.println("Please supply --vcfsort");
					printHelp();
				}
			} else if (cmd.hasOption("splitbatch")) {

				int n = 0;
				String ped = null;
				if (cmd.hasOption("f")) {
					ped = cmd.getOptionValue("f");
				} else {
					System.out.println("Please provide pedigree data (ped or fam file)");
					run = false;
				}

				if (cmd.hasOption("n")) {
					String dstr = cmd.getOptionValue("n");
					try {
						n = Integer.parseInt(dstr);
						if (n <= 0) {
							System.out.println("Value " + dstr + " provided should be larger than 0");
							run = false;
						}
					} catch (NumberFormatException e) {
						System.out.println("Value " + dstr + " provided for -n is not an int");
						run = false;
					}
				} else {
					run = false;
					System.out.println("Please provide batch size");
				}
				Random random = new Random();

				long seed = random.nextLong();
				if (cmd.hasOption("s")) {
					seed = Long.parseLong(cmd.getOptionValue("s"));
					System.out.println("Using preset seed: " + seed);
				}

				System.out.println(run);
				if (run) {

					BatchSplitter f = new BatchSplitter();
					System.out.println(input);
					System.out.println(ped);
					System.out.println(out);
					System.out.println(seed);
					f.splitVCFOverRandomBatches(input, ped, out, n, seed);

				} else {
					printHelp();
				}


			} else if (cmd.hasOption("sample")) {
				if (cmd.hasOption("p")) {
					String dstr = cmd.getOptionValue("p");
					try {
						perc = Double.parseDouble(dstr);
						if (perc >= 1 || perc <= 0) {
							System.out.println("Value " + dstr + " provided should be between 0 and 1");
							run = false;
						}
					} catch (NumberFormatException e) {
						System.out.println("Value " + dstr + " provided for -p is not a double.");
						run = false;
					}
				} else {
					run = false;
					System.out.println("Please provide sampling percentage");
				}

				if (cmd.hasOption("r")) {
					ref = cmd.getOptionValue("r");
				}


				if (run) {

					VariantSampler sampler = new VariantSampler();
					sampler.sample(input, ref, perc, out);

				} else {
					printHelp();
				}
			} else if (cmd.hasOption("splitchr")) {
				VCFFunctions f = new VCFFunctions();

				f.splitPerChromosome(input, out);

			} else if (cmd.hasOption("filterregions")) {

				if (cmd.hasOption("b")) {

					String chr = null;
					if (cmd.hasOption("chr")) {
						chr = cmd.getOptionValue("chr");
					}

					boolean list = false;
					if (cmd.hasOption("list")) {
						list = true;
					}

					VCFBedFilter f = new VCFBedFilter();
					f.filter(input, out, cmd.getOptionValue("b"), chr, list);

				} else {
					System.out.println("Please provide -b for --filter");
					printHelp();
				}

			} else if (cmd.hasOption("pseudo")) {

				String fam = null;

				if (cmd.hasOption("f")) {
					fam = cmd.getOptionValue("f");
				} else {
					System.out.println("Please provide -f for --pseudo");
					printHelp();
					System.exit(0);
				}


				PseudoControls c = new PseudoControls();
				c.make(input, out, fam, out + ".fam");


			} else {
				printHelp();
			}


		} catch (ParseException ex) {
			printHelp();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
	}


	public static void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		System.exit(-1);
	}

}
