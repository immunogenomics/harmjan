package nl.harmjanwestra.miscscripts;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 1/22/15.
 */
public class GATKBashScriptWriter {

	public static void main(String[] args) {
//		String sampleFile = "/Data/Projects/2014-FR-Reseq/sampleList.txt";
//		String genotypeJobsOutput = "/Data/Projects/2014-FR-Reseq/SequencingBashFiles/";
//		String gvcfdir = "/broad/hptmp/hwestra/variantcalling/variantcalls/gvcf/";
//		String combineGVCFJobsOut = "/Data/Projects/2014-FR-Reseq/SequencingBashFilesGVCFCombine/";
//		String tmpDir = "/dev/shm/";
//
//		createGenotypeCallingBashFiles(sampleFile, genotypeJobsOutput);
//
//
//		createCombineGVCFBashFiles(sampleFile, gvcfdir, tmpDir, combineGVCFJobsOut);
//
//
//		String regenotypeJobOut = "/Data/Projects/2014-FR-Reseq/SequencingBashFilesGVCFGenotype/";
//		makeGVCFGenotypeCommand(gvcfdir, regenotypeJobOut);

		String list = "/Data/Projects/2014-FR-Reseq/2015-finalRun/BQSR/allBamFiles.list";
		String bqsrtabledir = "/medpop/srlab/hwestra/fr-reseq/3BQSR/";
		String outputdir = "/broad/hptmp/hwestra/variantcalling/2015-03-30/BQSR/";
		String jobdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/BQSR/jobsDetermineBQSR/";

		try {
//			bqsrApply(list, dir, output);
//			Gpio.createDir(output);
//			 reGenotypeJobs(list, output, dir);
//			reGenotypeGVCFJobs(list, dir, output);

			list = "/Data/Projects/2014-FR-Reseq/2015-finalRun/Jobs/IndelRealign/2015-03-30-AllSelectedBAMFiles.txt";
			outputdir = "/broad/hptmp/hwestra/variantcalling/2015-03-30/IndelRealigner/";
			jobdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/Jobs/IndelRealign/jobs/";
//			indelrealign(list, jobdir, outputdir);

			jobdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/Jobs/BQSR/jobsDetermineBQSR/";
			bqsr(list, bqsrtabledir, jobdir);

			jobdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/Jobs/BQSR/jobsDetermineBQSRRecal/";
			Gpio.createDir(jobdir);
			bqsrRecal(list, bqsrtabledir, jobdir);
//			Gpio.createDir(jobdir);
//			bqsrApply(list, bqsrtabledir, outputdir, jobdir);


//			list = "/Data/Projects/2014-FR-Reseq/2015-finalRun/Regenotyping/BQSRBamFiles.txt";
//			outputdir = "/broad/hptmp/hwestra/variantcalling/2015-03-30/Regenotyping/";
//			jobdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/Regenotyping/jobs/";
//			Gpio.createDir(jobdir);
//			reGenotypeJobs(list, jobdir, outputdir);
//
//			list = "/Data/Projects/2014-FR-Reseq/2015-finalRun/Regenotyping/gvcfFiles.txt";
//			jobdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/Regenotyping/jobsgvcf/";
//			outputdir = "/broad/hptmp/hwestra/variantcalling/2015-03-30/MergedGenotypes/";
//			reGenotypeGVCFJobs(list, outputdir, jobdir);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void indelrealign(String fileList, String jobdir, String output) throws IOException {

		TextFile tf = new TextFile(fileList, TextFile.R);
		String[] files = tf.readAsArray();
		tf.close();


		for (String file : files) {
			String[] fileelems = file.split("/");
			String[] sampleNameElems = fileelems[fileelems.length - 1].split("\\.");
			String sampleName = sampleNameElems[0];

			TextFile tf2 = new TextFile(jobdir + sampleName + ".sh", TextFile.W);
			String outStr = "nice -n 20 java -Xmx4g -jar ~/Tools/GATK3.3/GenomeAnalysisTK.jar \\\n" +
					"       -T IndelRealigner \\\n" +
					"       -R /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta \\\n" +
					"       -I " + file + " \\\n" +
					"       -o " + output + sampleName + ".bam \\\n" +
					"       -targetIntervals /medpop/srlab/hwestra/fr-reseq/2015-03-30-FinalRun/2015-03-30-intervalListWithKnownIndels.intervals \\\n" +
					"       -known /humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_phase1.indels.b37.vcf \\\n" +
					"       -known /humgen/gsa-hpprojects/GATK/bundle/current/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \\\n" +
					"       ";
			tf2.writeln(outStr);
			tf2.close();
		}
	}

	public static void bqsr(String fileList, String outputdir, String jobdir) throws IOException {
		TextFile tf = new TextFile(fileList, TextFile.R);
		String[] files = tf.readAsArray();
		tf.close();


		for (String file : files) {

			String[] fileelems = file.split("/");
			String[] sampleNameElems = fileelems[fileelems.length - 1].split("\\.");
			String sampleName = sampleNameElems[0];

			TextFile out = new TextFile(jobdir + sampleName + ".sh", TextFile.W);
			String outStr = "nice -n 20 java -Xmx4g -jar /broad/hptmp/hwestra/variantcalling/GATK3.3/GenomeAnalysisTK.jar \\\n" +
					"   -T BaseRecalibrator \\\n" +
					"   -R /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta \\\n" +
					"   -I " + file + "\\\n" +
					"   -knownSites /humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_138.b37.vcf \\\n" +
					"   -knownSites /humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_phase1.indels.b37.vcf \\\n" +
					"   -knownSites /humgen/gsa-hpprojects/GATK/bundle/current/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \\\n" +
					"   -L /medpop/srlab/hwestra/fr-reseq/2015-03-30-FinalRun/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed \\\n" +
					"   -o " + outputdir + sampleName + ".bqsrTable \\\n" +
					"   \n" +
					"\n";
			out.writeln(outStr);
			out.close();
		}


	}


	public static void bqsrRecal(String fileList, String outputdir, String jobdir) throws IOException {
		TextFile tf = new TextFile(fileList, TextFile.R);
		String[] files = tf.readAsArray();
		tf.close();


		for (String file : files) {

			String[] fileelems = file.split("/");
			String[] sampleNameElems = fileelems[fileelems.length - 1].split("\\.");
			String sampleName = sampleNameElems[0];

			TextFile out = new TextFile(jobdir + sampleName + ".sh", TextFile.W);
			String outStr = "nice -n 20 java -Xmx4g -jar /broad/hptmp/hwestra/variantcalling/GATK3.3/GenomeAnalysisTK.jar \\\n" +
					"   -T BaseRecalibrator \\\n" +
					"   -R /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta \\\n" +
					"   -I " + file + "\\\n" +
					"   -knownSites /humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_138.b37.vcf \\\n" +
					"   -knownSites /humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_phase1.indels.b37.vcf \\\n" +
					"   -knownSites /humgen/gsa-hpprojects/GATK/bundle/current/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \\\n" +
					"   -L /medpop/srlab/hwestra/fr-reseq/2015-03-30-FinalRun/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed \\\n" +
					"   -BQSR " + outputdir + sampleName + ".bqsrTable \\\n" +
					"   -o " + outputdir + sampleName + "-after.bqsrTable \\\n" +
					"   \n" +
					"\n";
			outStr += "nice -n 20 java -Xmx4g -jar /broad/hptmp/hwestra/variantcalling/GATK3.3/GenomeAnalysisTK.jar \\\n" +
					"   -T AnalyzeCovariates \\\n" +
					"   -R /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta \\\n" +
					"   -before " + outputdir + sampleName + ".bqsrTable \\\n" +
					"   -after " + outputdir + sampleName + "-after.bqsrTable \\\n" +
					"   -plots " + outputdir + sampleName + "-cal.pdf \\\n" +
					"   \n" +
					"\n";
			out.writeln(outStr);
			out.close();


			/*
			java	–jar	GenomeAnalysisTK.jar	–T	AnalyzeCovariates	\

–R	human.fasta	\

–before	recal.table	\

–aSer	aSer_recal.table	\

–plots	recal_plots.pdf
			 */

		}


	}

	public static void bqsrApply(String fileList, String bqsrtabledir, String outputdir, String joboutputdir) throws IOException {
		TextFile tf = new TextFile(fileList, TextFile.R);
		String[] files = tf.readAsArray();
		tf.close();


		for (String file : files) {

			String[] fileelems = file.split("/");
			String[] sampleNameElems = fileelems[fileelems.length - 1].split("\\.");
			String sampleName = sampleNameElems[0];

			TextFile out = new TextFile(joboutputdir + sampleName + ".sh", TextFile.W);
			String outStr = "nice -n 20 java -Xmx2g -jar /broad/hptmp/hwestra/variantcalling/GATK3.3/GenomeAnalysisTK.jar \\\n" +
					"   -T PrintReads \\\n" +
					"   -R /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta \\\n" +
					"   -I " + file + "\\\n" +
					"   -BQSR " + bqsrtabledir + sampleName + ".bqsrTable \\\n" +
					"   -o " + outputdir + sampleName + "-bqsr.bam \\\n" +
					"   1> " + outputdir + sampleName + "-bqsr.log \\\n" +
					"   2> " + outputdir + sampleName + "-bqsr.err \\\n" +
					"\n";
			out.writeln(outStr);
			out.close();
		}


	}

	public static void reGenotypeJobs(String fileList, String jobOutput, String outdir) throws IOException {
		TextFile tf = new TextFile(fileList, TextFile.R);
		String[] files = tf.readAsArray();
		tf.close();

		for (String file : files) {

			String[] fileelems = file.split("/");
			String[] sampleNameElems = fileelems[fileelems.length - 1].split("\\.");
			String sampleName = sampleNameElems[0];

			String out = "nice -n 15 java -Xmx5g -jar /broad/hptmp/hwestra/variantcalling/GATK3.3/GenomeAnalysisTK.jar -T HaplotypeCaller \\\n" +
					"                -R /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta \\\n" +
					"                -I " + file + " \\\n" +
					"                -o " + outdir + sampleName + ".vcf \\\n" +
					"                --dbsnp /humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_138.b37.vcf \\\n" +
					"                -minPruning 3 \\\n" +
					"                --emitRefConfidence GVCF \\\n" +
					"                --variant_index_type LINEAR \\\n" +
					"                --variant_index_parameter 128000 \\\n" +
					"                -L /medpop/srlab/hwestra/fr-reseq/2015-03-30-FinalRun/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed \\\n" +
					"                &> " + outdir + sampleName + ".log \\\n" +
					"";
			TextFile jobOut = new TextFile(jobOutput + sampleName + ".sh", TextFile.W);
			jobOut.writeln(out);
			jobOut.close();
		}


	}

	public static void reGenotypeGVCFJobs(String filelist, String outdir, String joboutput) throws IOException {
		TextFile tf = new TextFile(filelist, TextFile.R);
		String[] files = tf.readAsArray();
		tf.close();

		String out = "nice -n 15 java -Xmx5g -jar /broad/hptmp/hwestra/variantcalling/GATK3.3/GenomeAnalysisTK.jar \\\n" +
				"   -T GenotypeGVCFs \\\n" +
				"   -R /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta \\\n" +
				"   --dbsnp /humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_138.b37.vcf \\\n" +
				"   -L /broad/hptmp/hwestra/variantcalling/2015-01-07-targetRegions.bed \\\n" +
				"   -o " + outdir + "/merged.vcf \\\n";
		for (String file : files) {
			out += "   --variant " + file + " \\\n";
		}


		TextFile outfi = new TextFile(joboutput + "genotypegvcf.sh", TextFile.W);
		outfi.writeln(out);
		outfi.close();


	}


	// regenotype gvcf files
	public static void makeGVCFGenotypeCommand(String gvcfdir, String regenotypeJobOut) {
		String out = "nice -n 15 java -Xmx5g -jar /broad/hptmp/hwestra/variantcalling/GATK3.3/GenomeAnalysisTK.jar \\\n" +
				"   -T GenotypeGVCFs \\\n" +
				"   -R /dev/shm/human_g1k_v37.fasta \\\n" +
				"   --dbsnp /humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_138.b37.vcf \\\n" +
				"   -L /broad/hptmp/hwestra/variantcalling/2015-01-07-targetRegions.bed \\\n" +
				"   -o /broad/hptmp/hwestra/variantcalling/variantcalls/merged.vcf \\\n";

		for (int i = 1; i < 24; i++) {
			String chr = "" + i;
			if (i == 23) {
				chr = "X";
			}

			out += "   --variant " + gvcfdir + "merged-" + chr + ".gvcf \\\n";
		}
		try {

			TextFile outfi = new TextFile(regenotypeJobOut + "genotype.sh", TextFile.W);
			outfi.writeln(out);
			outfi.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// combine GVCF files per chromosome
	public static void createCombineGVCFBashFiles(String sampleFile, String gvcfdir, String tmpDir, String bashoutdir) {
		try {
			TextFile tfIn = new TextFile(sampleFile, TextFile.R);

			ArrayList<String> samples = new ArrayList<String>();
			String ln = tfIn.readLine();
			while (ln != null) {
				if (ln.trim().length() > 0) {
					samples.add(ln);
				}
				ln = tfIn.readLine();
			}
			tfIn.close();
			System.out.println(samples.size() + " samples loaded.");

			for (int i = 1; i < 24; i++) {
				String chr = "" + i;
				if (i == 23) {
					chr = "X";
				}

				TextFile tf = new TextFile(bashoutdir + chr + ".sh", TextFile.W);
				String gatkCommand = "nice -n 15 java -Xmx5g -jar /broad/hptmp/hwestra/variantcalling/GATK3.3/GenomeAnalysisTK.jar \\\n" +
						"   -R /dev/shm/human_g1k_v37.fasta \\\n" +
						"   -T CombineGVCFs \\\n" +
						"   -L /broad/hptmp/hwestra/variantcalling/2015-01-07-targetRegions.bed \\\n" +
						"   -o " + gvcfdir + "merged-" + chr + ".gvcf \\\n";


				String gvcfStr = "";
				String copyCommand = "";
				String delCommand = "";
				for (String sample : samples) {
					String origSample = sample;
					String parsedSample = sample.replaceAll("/", "-");
					String sampleGVCFIn = gvcfdir + parsedSample + "/" + chr + ".vcf";
					String sampleGVCFOut = tmpDir + parsedSample + "-" + chr + ".vcf";

					copyCommand += "cp " + sampleGVCFIn + " " + sampleGVCFOut + "\n";
					delCommand += "rm " + sampleGVCFOut + "\n";
					gvcfStr += "   --variant " + sampleGVCFOut + " \\\n";
				}

				gatkCommand += gvcfStr;
				tf.writeln(copyCommand);
				tf.writeln(gatkCommand);
				tf.writeln(delCommand);
				tf.close();
			}

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// split chromosome per sample, and generate GVCF files
	public static void createGenotypeCallingBashFiles(String sampleFile, String outdir) {

		try {
			TextFile tfIn = new TextFile(sampleFile, TextFile.R);

			ArrayList<String> samples = new ArrayList<String>();
			String ln = tfIn.readLine();
			while (ln != null) {
				if (ln.trim().length() > 0) {
					samples.add(ln);
				}
				ln = tfIn.readLine();
			}
			tfIn.close();
			System.out.println(samples.size() + " samples loaded.");

			for (String sample : samples) {


				String origSample = sample;
				String parsedSample = sample.replaceAll("/", "-");
				TextFile bashOut = new TextFile(outdir + parsedSample + ".sh", TextFile.W);
				String script = "#!/bin/bash\n" +
						"\n" +
						"mkdir -p /broad/hptmp/hwestra/variantcalling/variantcalls/gvcf/" + parsedSample + "/ \n" +
						"mkdir -p /dev/shm/hwestra/ \n" +
						"\n" +
						"for i in {1..22}\n" +
						"do\n" +
						"samtools view -b -h -r " + origSample + " /broad/hptmp/hwestra/variantcalling/perchr/$i.bam -o /dev/shm/hwestra/" + parsedSample + "-$i.bam\n" +
						"samtools index /dev/shm/hwestra/" + parsedSample + "-$i.bam /dev/shm/hwestra/" + parsedSample + "-$i.bam.bai \n" +
						"\n" +
						"nice -n 15 java -Xmx5g -jar /broad/hptmp/hwestra/variantcalling/GATK3.3/GenomeAnalysisTK.jar -T HaplotypeCaller \\\n" +
						"                -R /dev/shm/human_g1k_v37.fasta \\\n" +
						"                -I /dev/shm/hwestra/" + parsedSample + "-$i.bam \\\n" +
						"                -o /broad/hptmp/hwestra/variantcalling/variantcalls/gvcf/" + parsedSample + "/$i.vcf \\\n" +
						"                --dbsnp /humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_138.b37.vcf \\\n" +
						"                -stand_call_conf 30 \\\n" +
						"                -stand_emit_conf 10 \\\n" +
						"                -minPruning 3 \\\n" +
						"                --sample_name " + origSample + " \\\n" +
						"                --emitRefConfidence GVCF \\\n" +
						"                --variant_index_type LINEAR \\\n" +
						"                --variant_index_parameter 128000 \\\n" +
						"                -L /broad/hptmp/hwestra/variantcalling/perchr/$i.bed \\\n" +
						"                &> /broad/hptmp/hwestra/variantcalling/variantcalls/gvcf/" + parsedSample + "/$i-log.txt\n" +
						"rm /dev/shm/hwestra/" + parsedSample + "-$i.* \n" +
						"done\n" +
						"\n" +
						"samtools view -b -h -r " + origSample + " /broad/hptmp/hwestra/variantcalling/perchr/X.bam -o /dev/shm/hwestra/" + parsedSample + "-X.bam\n" +
						"samtools index /dev/shm/hwestra/" + parsedSample + "-X.bam /dev/shm/hwestra/" + parsedSample + "-X.bam.bai \n" +
						"nice -n 15 java -Xmx5g -jar /broad/hptmp/hwestra/variantcalling/GATK3.3/GenomeAnalysisTK.jar -T HaplotypeCaller \\\n" +
						"                -R /dev/shm/human_g1k_v37.fasta \\\n" +
						"                -I /dev/shm/hwestra/" + parsedSample + "-X.bam \\\n" +
						"                -o /broad/hptmp/hwestra/variantcalling/variantcalls/gvcf/" + parsedSample + "/X.vcf \\\n" +
						"                --dbsnp /humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_138.b37.vcf \\\n" +
						"                -stand_call_conf 30 \\\n" +
						"                -stand_emit_conf 10 \\\n" +
						"                -minPruning 3 \\\n" +
						"                --sample_name " + origSample + " \\\n" +
						"                --emitRefConfidence GVCF \\\n" +
						"                --variant_index_type LINEAR \\\n" +
						"                --variant_index_parameter 128000 \\\n" +
						"                -L /broad/hptmp/hwestra/variantcalling/perchr/X.bed \\\n" +
						"                &> /broad/hptmp/hwestra/variantcalling/variantcalls/gvcf/" + parsedSample + "/X-log.txt\n" +
						"";
				bashOut.writeln(script);
				bashOut.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
