import java.util.ArrayList;
import java.util.Date;

/**
 * Created by hwestra on 9/19/15.
 */
public class Individual {


	private int XPixel;

	public void setXPixel(int XPixel) {
		this.XPixel = XPixel;
	}
	public int getXPixel(){
		return XPixel;
	}

	public enum Gender {
		MALE,
		FEMALE,
		UNKNOWN;

		public static Gender parseGender(String gstr) {
			if (gstr.toLowerCase().equals("m")) {
				return MALE;
			} else if (gstr.toLowerCase().equals("v")) {
				return FEMALE;
			}
			return UNKNOWN;
		}
	}

	ArrayList<Individual> children;
	ArrayList<Relationship> relationships;

	int id;
	Gender gender;
	String firstName;
	String lastName;
	Date dateofbirth;
	Date dateofdeath;
	Location placeOfBirth;
	Location placeOfDeath;
	Location placeOfLiving;
	Individual father;
	Individual mother;


	int x = -1;
	int y = -1;
	int familyId = -1;


	public int getFamilyId() {
		return familyId;
	}

	public void setFamilyId(int familyId) {
		this.familyId = familyId;
	}


	public int getX() {
		return x;
	}

	public void setX(int x) {
		this.x = x;
	}

	public int getY() {
		return y;
	}

	public void setY(int y) {
		this.y = y;
	}

	public Individual getFather() {
		return father;
	}

	public void setFather(Individual father) {
		this.father = father;
	}

	public Individual getMother() {
		return mother;
	}

	public void setMother(Individual mother) {
		this.mother = mother;
	}

	public Individual(int id,
					  String firstname,
					  String lastname,
					  Gender gender,
					  Date dayofbirth,
					  Date dayofdeath,
					  Location placeOfBirth,
					  Location placeOfDeath,
					  Location placeOfLiving
	) {
		this.children = new ArrayList<Individual>();
		this.relationships = new ArrayList<Relationship>();
		this.id = id;
		this.firstName = firstname;
		this.lastName = lastname;
		this.gender = gender;
		this.dateofbirth = dayofbirth;
		this.dateofdeath = dayofdeath;
		this.placeOfBirth = placeOfBirth;
		this.placeOfDeath = placeOfDeath;
		this.placeOfLiving = placeOfLiving;
	}

	public void addChild(Individual child) {
		children.add(child);
	}

	public void addRelationship(Relationship relationship) {
		relationships.add(relationship);
	}

	public ArrayList<Individual> getChildren() {
		return children;
	}

	public void setChildren(ArrayList<Individual> children) {
		this.children = children;
	}

	public ArrayList<Relationship> getRelationships() {
		return relationships;
	}

	public void setRelationships(ArrayList<Relationship> relationships) {
		this.relationships = relationships;
	}

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	public Gender getGender() {
		return gender;
	}

	public void setGender(Gender gender) {
		this.gender = gender;
	}

	public String getFirstName() {
		return firstName;
	}

	public void setFirstName(String firstName) {
		this.firstName = firstName;
	}

	public String getLastName() {
		return lastName;
	}

	public void setLastName(String lastName) {
		this.lastName = lastName;
	}

	public Date getDateofbirth() {
		return dateofbirth;
	}

	public void setDateofbirth(Date dateofbirth) {
		this.dateofbirth = dateofbirth;
	}

	public Date getDateofdeath() {
		return dateofdeath;
	}

	public void setDateofdeath(Date dateofdeath) {
		this.dateofdeath = dateofdeath;
	}

	public Location getPlaceOfBirth() {
		return placeOfBirth;
	}

	public void setPlaceOfBirth(Location placeOfBirth) {
		this.placeOfBirth = placeOfBirth;
	}

	public Location getPlaceOfDeath() {
		return placeOfDeath;
	}

	public void setPlaceOfDeath(Location placeOfDeath) {
		this.placeOfDeath = placeOfDeath;
	}

	public Location getPlaceOfLiving() {
		return placeOfLiving;
	}

	public void setPlaceOfLiving(Location placeOfLiving) {
		this.placeOfLiving = placeOfLiving;
	}

	@Override
	public String toString() {
		return
				"id=" + id +
						", firstName='" + firstName + '\'' +
						", lastName='" + lastName + '\'' +
						", children='" + children.size();
	}
}
