import java.util.Comparator;

/**
 * Created by hwestra on 9/20/15.
 */
public class IndividualComparator implements Comparator<Individual> {
	public int compare(Individual o1, Individual o2) {
		if(o1.getDateofbirth() == null){
			return 1;
		}
		if(o2.getDateofbirth() == null){
			return -1;
		}
		if(o1.getDateofbirth().before(o2.getDateofbirth())){
			return -1;
		} else if(o2.getDateofbirth().after(o2.getDateofbirth())){
			return 1;
		}
		return 0;
	}
}
