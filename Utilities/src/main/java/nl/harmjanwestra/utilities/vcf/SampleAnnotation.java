package nl.harmjanwestra.utilities.vcf;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.enums.Gender;

import java.util.ArrayList;

/**
 * Created by hwestra on 9/7/16.
 */
public class SampleAnnotation {
	private DiseaseStatus[][] sampleDiseaseStatus;
	private Gender[] individualGender;
	private String[] sampleName;
	private ArrayList<String> covariateNames;

	public DoubleMatrix2D getCovariates() {
		return covariates;
	}

	public void setCovariates(DoubleMatrix2D covariates) {
		this.covariates = covariates;
	}

	private DoubleMatrix2D covariates;

	public DiseaseStatus[][] getSampleDiseaseStatus() {
		return sampleDiseaseStatus;
	}

	public void setSampleDiseaseStatus(DiseaseStatus[][] sampleDiseaseStatus) {
		this.sampleDiseaseStatus = sampleDiseaseStatus;
	}

	public void setSampleDiseaseStatus(DiseaseStatus[] sampleDiseaseStatus) {
		this.sampleDiseaseStatus = new DiseaseStatus[sampleDiseaseStatus.length][1];
		for (int q = 0; q < sampleDiseaseStatus.length; q++) {
			this.sampleDiseaseStatus[q][0] = sampleDiseaseStatus[q];
		}
	}

	public Gender[] getIndividualGender() {
		return individualGender;
	}

	public void setIndividualGender(Gender[] individualGender) {
		this.individualGender = individualGender;
	}

	public String[] getSampleName() {
		return sampleName;
	}

	public void setSampleName(String[] sampleName) {
		this.sampleName = sampleName;
	}

	public void setCovariateNames(ArrayList<String> covariateNames) {
		this.covariateNames = covariateNames;
	}
}
