import java.util.Date;

/**
 * Created by hwestra on 9/19/15.
 */
public class Relationship {
	Date start;
	Date end;
	Individual individual1;
	Individual individual2;
	KIND kind;

	public Relationship(Individual ind1, Individual ind2, Date start, Date end, KIND kind) {
		this.individual1 = ind1;
		this.individual2 = ind2;
		this.start = start;
		this.end = end;
		this.kind = kind;
	}

	public enum KIND {
		MARRIED,
		NOTMARRIED,
		ADOPTION
	}

	public Date getStart() {
		return start;
	}

	public Date getEnd() {
		return end;
	}

	public Individual getIndividual1() {
		return individual1;
	}

	public Individual getIndividual2() {
		return individual2;
	}

	public KIND getKind() {
		return kind;
	}
}
