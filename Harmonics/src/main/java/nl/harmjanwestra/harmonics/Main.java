/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.harmonics;

import nl.harmjanwestra.harmonics.posthoc.PeakMerge;

/**
 * @author hwestra
 */
public class Main {
	public static void main(String[] args) {
		try {


			Harmonics e = new Harmonics();
//
//
//			String bedfile = "/Data/tmp/atac/greenLeafRepresentativeRegion.bed";
//			String[] bamfilelocs = new String[]{
//					"/Data/tmp/atac/bwa/20141202-0-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/bwa/20141202-1-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/bwa/20141202-2-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/bwa/20141202-4-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/bwa/20141202-8-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/bwa/20141202-12-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/bwa/20141202-24-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/bwa/20141202-48-hrs-sorted-dedup.bam"
//			};
//			String[] sampleNames = new String[]{
//					"0hrs",
//					"1hrs",
//					"2hrs",
//					"4hrs",
//					"8hrs",
//					"12hrs",
//					"24hrs",
//					"48hrs",
//			};
//
//			String outputdir = "/Data/tmp/atac/output/";
//
//			int posshift = 4;
//			int negshift = -5;
//			int distsize = 50;
//			int smoothingwindow = 5;
//
//
//			e.run(bedfile, bamfilelocs, sampleNames, posshift, negshift, distsize, smoothingwindow, outputdir);
//
//			bamfilelocs = new String[]{
//					"/Data/tmp/atac/greenleaf/bwa/SRR891268-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/greenleaf/bwa/SRR891269-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/greenleaf/bwa/SRR891270-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/greenleaf/bwa/SRR891271-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/greenleaf/bwa/SRR891272-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/greenleaf/bwa/SRR891273-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/greenleaf/bwa/SRR891274-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/greenleaf/bwa/SRR891275-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/greenleaf/bwa/SRR891276-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/greenleaf/bwa/SRR891277-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/greenleaf/bwa/SRR891278-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/greenleaf/bwa/SRR891279-hrs-sorted-dedup.bam",
//					"/Data/tmp/atac/greenleaf/bwa/SRR891280-hrs-sorted-dedup.bam",
//			};
//			sampleNames = new String[]{
//					"SRR891268",
//					"SRR891269",
//					"SRR891270",
//					"SRR891271",
//					"SRR891272",
//					"SRR891273",
//					"SRR891274",
//					"SRR891275",
//					"SRR891276",
//					"SRR891277",
//					"SRR891278",
//					"SRR891279",
//					"SRR891280",
//			};
//
//			outputdir = "/Data/tmp/atac/output/greenleaf";
//
//			e.run(bedfile, bamfilelocs, sampleNames, posshift, negshift, distsize, smoothingwindow, outputdir);



			PeakMerge merger = new PeakMerge();



		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
