import java.util.ArrayList;
import java.util.Objects;

/**
 * Created by hwestra on 9/19/15.
 */
public class Location {

	String locationName;
	String countryName;

	ArrayList<Individual> individualsAtLocation;

	public Location(String location, String country) {
		this.locationName = location;
		this.countryName = country;
		this.individualsAtLocation = new ArrayList<Individual>();
	}

	public String getLocationName() {
		return locationName;
	}

	public String getCountryName() {
		return countryName;
	}

	public ArrayList<Individual> getIndividualsAtLocation() {
		return individualsAtLocation;
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (!(o instanceof Location)) return false;
		Location location = (Location) o;
		return Objects.equals(locationName, location.locationName) &&
				Objects.equals(countryName, location.countryName);
	}

	@Override
	public int hashCode() {
		return Objects.hash(locationName, countryName);
	}
}
