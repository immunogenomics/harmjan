package nl.harmjanwestra.utilities.vcf;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import nl.harmjanwestra.utilities.legacy.genetica.io.BinaryFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 5/18/16.
 */
public class VCFGenotypeDataBinary extends BinaryFile {

	private String[] sampleNames;
	private boolean isPhased;

	public VCFGenotypeDataBinary(String loc, boolean mode) throws IOException {
		super(loc, mode);
		if (mode == BinaryFile.R) {
			readHeader();
		}
	}

	public VCFGenotypeDataBinary(String loc, boolean mode, int buffersize) throws IOException {
		super(loc, mode, buffersize);
		if (mode == BinaryFile.R) {
			readHeader();
		}
	}

	public VCFVariant readVariant() throws IOException {

		String chr = this.is.readUTF();
		Integer pos = this.is.readInt();
		String id = this.is.readUTF();
		String alleleStr = this.is.readUTF();
		String info = this.is.readUTF();
		byte nrAlleles = this.is.readByte();
		boolean hasProbs = this.is.readBoolean();

		DoubleMatrix2D alleles = new DenseDoubleMatrix2D(sampleNames.length, nrAlleles);
//		// read all
//		for (int row = 0; row < alleles.rows(); row++) {
//			byte[] rowB = new byte[sampleNames.length];
//			this.is.read(rowB);
//			for (int col = 0; col < rowB.length; col++) {
//				alleles.setQuick(row, col, rowB[col]);
//			}
//		}
		DoubleMatrix2D dosages = null;
//		if (hasProbs) {
//			dosages = new DenseDoubleMatrix2D(sampleNames.length, nrAlleles);
//			for (int col = 0; col < dosages.columns(); col++) {
//				byte[] b = new byte[sampleNames.length];
//				this.is.read(b);
//				for (int row = 0; row < dosages.rows(); row++) {
//					double d = (127d / 2) / b[row];
//					dosages.set(row, col, d);
//				}
//			}
//		}

		return new VCFVariant(chr, pos, id, alleleStr, info, alleles, dosages, null);

	}

	public String[] getSampleNames() {
		return sampleNames;
	}

	public void writeVariant(VCFVariant variant) throws IOException {

		this.os.writeUTF(variant.getChr().toString());
		this.os.writeInt(variant.getPos());
		this.os.writeUTF(variant.getId());
		this.os.writeUTF(Strings.concat(variant.getAlleles(), Strings.comma));
		this.os.writeUTF(variant.getInfoString());
		this.os.writeByte(variant.getAlleles().length);
		this.os.writeBoolean(variant.hasImputationDosages());

		DoubleMatrix2D alleles = variant.getGenotypeAllelesAsMatrix2D(); // format: [nrSamples][alleles]
		for (int row = 0; row < alleles.rows(); row++) {
			byte[] balleles = new byte[alleles.columns()];
			for (int col = 0; col < alleles.columns(); col++) {
				balleles[col] = (byte) alleles.getQuick(row, col);
			}
			this.os.write(balleles);
		}

		if (variant.hasImputationDosages()) {
			DoubleMatrix2D dosages = variant.getGenotypeDosagesAsMatrix2D(); // format: [nrSamples][alleles];

			// flip for more efficient storage.
			for (int col = 0; col < dosages.columns(); col++) {
				byte[] outb = new byte[dosages.rows()];
				for (int row = 0; row < dosages.rows(); row++) {
					double dosage = dosages.getQuick(row, col);
					// convert to byte..
					// dosage varies between 0 and 2
					// byte varies between -128 (reserved for -1) and 127. lets use -1 as missing and 0-127 for the range
					byte bd = -1;
					if (dosage != -1) {
						bd = (byte) (Math.ceil((double) 127 / 2 * dosage));
					}
					outb[row] = bd;
				}
				this.os.write(outb);
			}
		}
	}

	public void setIsPhased(boolean b) {
		isPhased = b;
	}

	public void writeHeader(ArrayList<String> samples) throws IOException {
		writeHeader(samples.toArray(new String[0]));
	}

	public void writeHeader(String[] samples) throws IOException {
		this.os.writeByte(1);
		this.os.writeBoolean(isPhased);
		this.os.writeInt(samples.length);
		for (int s = 0; s < samples.length; s++) {
			this.os.writeChars(samples[s]);
		}
	}

	public void readHeader() throws IOException {
		byte magic = this.is.readByte();
		isPhased = this.is.readBoolean();
		Integer nrSamples = this.is.readInt();
		this.sampleNames = new String[nrSamples];
		for (int s = 0; s < nrSamples; s++) {
			String sampleName = this.is.readUTF();
			sampleNames[s] = sampleName;
		}
	}

	public void write(String ln) throws IOException {
		if (!ln.startsWith("#")) {
			VCFVariant variant = new VCFVariant(ln);
			writeVariant(variant);
		}


	}
}
