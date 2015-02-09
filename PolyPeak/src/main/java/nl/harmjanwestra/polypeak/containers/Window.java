/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.polypeak.containers;

import java.util.ArrayList;

/**
 * @author hwestra
 */
public class Window {

	private int[][] coverage;
	private final int start;
	private final int stop;
	private final Sequence chr;
	private int maxCoverage = 0;
	private int minCoverage = 0;
	private boolean stranded;
	private double posMean = 0;
	private double negMean = 0;
	private double posVariance = 0;
	private double negVariance = 0;
	private int nrReadsNeg = 0;
	private int nrReadsPos = 0;
	private double insertSizeVariance = 0;

	private double insertSizeMean = 0;
	private int nrRecordsPassingFilter;
	private int nrRecordsTotal;
	private ArrayList<Integer> insertSizes;
	private ArrayList<Integer> readLengths;
	private double readLengthMean;
	private double readLengthVariance;

	public Window(int start, int stop, Sequence chr) {
		this.start = start;
		this.stop = stop;
		this.chr = chr;
	}

	public int[][] getCoverage() {
		return coverage;
	}

	public int getStart() {
		return start;
	}

	public int getStop() {
		return stop;
	}

	public Sequence getChr() {
		return chr;
	}

	public int getMaxCoverage() {
		return maxCoverage;
	}

	public int getMinCoverage() {
		return minCoverage;
	}

	public void setCoverage(int[][] coverage) {
		this.coverage = coverage;
	}

	public void setStranded(boolean stranded) {
		this.stranded = stranded;
	}

	public boolean isStranded() {
		return stranded;
	}

	public int getPosNrReads() {
		return nrReadsPos;
	}

	public int getNegNrReads() {
		return nrReadsNeg;
	}

	public double getInsertSizeMean() {
		return insertSizeMean;
	}

	public double getInsertSizeVariance() {
		return insertSizeVariance;
	}

	public void setMaxCoverage(int maxCoverage) {
		this.maxCoverage = maxCoverage;
	}

	public void setMinCoverage(int minCoverage) {
		this.minCoverage = minCoverage;
	}

	public void setPosMean(double posMean) {
		this.posMean = posMean;
	}

	public double getPosMean() {
		return posMean;
	}

	public void setNegMean(double negMean) {
		this.negMean = negMean;
	}

	public double getNegMean() {
		return negMean;
	}

	public void setPosVariance(double posVariance) {
		this.posVariance = posVariance;
	}

	public double getPosVariance() {
		return posVariance;
	}


	public void setNegVariance(double negVariance) {
		this.negVariance = negVariance;
	}

	public double getNegVariance() {
		return negVariance;
	}

	public void setNrReadsNeg(int nrReadsNeg) {
		this.nrReadsNeg = nrReadsNeg;
	}

	public void setNrReadsPos(int nrReadsPos) {
		this.nrReadsPos = nrReadsPos;
	}

	@Override
	public String toString() {
		return chr.getName() + ":" + start + "-" + stop;
	}

	public void setNrRecordsPassingFilter(int nrRecordsPassingFilter) {
		this.nrRecordsPassingFilter = nrRecordsPassingFilter;
	}

	public int getNrRecordsPassingFilter() {
		return nrRecordsPassingFilter;
	}

	public void setNrRecordsTotal(int nrRecordsTotal) {
		this.nrRecordsTotal = nrRecordsTotal;
	}

	public int getNrRecordsTotal() {
		return nrRecordsTotal;
	}

	public void setInsertSizes(ArrayList<Integer> insertSizes) {
		this.insertSizes = insertSizes;
	}

	public void setReadLengths(ArrayList<Integer> readLengths) {
		this.readLengths = readLengths;
	}

	public ArrayList<Integer> getInsertSizes() {
		return insertSizes;
	}

	public ArrayList<Integer> getReadLengths() {
		return readLengths;
	}

	public void setInsertSizeMean(double insertSizeMean) {
		this.insertSizeMean = insertSizeMean;
	}

	public void setInsertSizeVariance(double insertSizeVariance) {
		this.insertSizeVariance = insertSizeVariance;
	}

	public void setReadLengthMean(double readLengthMean) {
		this.readLengthMean = readLengthMean;
	}

	public void setReadLengthVariance(double readLengthVariance) {
		this.readLengthVariance = readLengthVariance;
	}

	public double getReadLengthMean() {
		return readLengthMean;
	}

	public double getReadLengthVariance() {
		return readLengthVariance;
	}
}
