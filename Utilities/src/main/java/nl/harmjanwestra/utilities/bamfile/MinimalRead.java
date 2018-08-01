package nl.harmjanwestra.utilities.bamfile;

import htsjdk.samtools.Cigar;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.Strand;

import java.util.Objects;

/**
 * Created by hwestra on 7/2/15.
 */
public class MinimalRead {

	private Chromosome chr;
	private String name;
	private Cigar cigar;
	private int end;
	private int start;
	private MinimalRead mate;
	private Strand strand;

	public MinimalRead(String name) {
		this.name = name;
	}

	public MinimalRead(String name, Chromosome chr, int start, int end, Cigar cigar, Strand strand) {
		this.name = name;
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.cigar = cigar;
		this.strand = strand;
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (!(o instanceof MinimalRead)) return false;
		MinimalRead that = (MinimalRead) o;
		return Objects.equals(name, that.name);
	}

	@Override
	public int hashCode() {
		return Objects.hash(name);
	}

	public void setChr(Chromosome chr) {

		this.chr = chr;
	}

	public void setName(String name) {
		this.name = name;
	}

	public void setCigar(Cigar cigar) {
		this.cigar = cigar;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public void setMate(MinimalRead mate) {
		this.mate = mate;
	}

	public Chromosome getChr() {
		return chr;
	}

	public String getName() {
		return name;
	}

	public Cigar getCigar() {
		return cigar;
	}

	public int getEnd() {
		return end;
	}

	public int getStart() {
		return start;
	}

	public MinimalRead getMate() {
		return mate;
	}

	public MinimalRead copy() {
		MinimalRead read = new MinimalRead(name, chr, start, end, cigar, strand);
		return read;
	}

	public void setStrand(Strand strand) {
		this.strand = strand;
	}

	public Strand getStrand() {

		return strand;
	}
}
