import com.itextpdf.text.DocumentException;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * Created by hwestra on 9/19/15.
 */
public class FamilyTree {

	public static void main(String[] args) {
		FamilyTree f = new FamilyTree();
		try {
			f.run("/Sync/OneDrive/genealogie.txt");
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void run(String treeInfo) throws IOException {

		ArrayList<Individual> individuals = parseTreeInfo(treeInfo);
		System.out.println(individuals.size() + " individuals loaded");


		// iterate tree, find persons without parents


		ArrayList<Date> allDates = new ArrayList<Date>();
		HashMap<Date, ArrayList<Individual>> individualsPerDate = new HashMap<Date, ArrayList<Individual>>();
		for (int i = 0; i < individuals.size(); i++) {
			Individual ind = individuals.get(i);
			Date dob = ind.getDateofbirth();
			if (dob != null) {
				if (!individualsPerDate.containsKey(dob)) {
					allDates.add(dob);
				}

				ArrayList<Individual> individualsForDate = individualsPerDate.get(dob);
				if (individualsForDate == null) {
					individualsForDate = new ArrayList<Individual>();
				}
				individualsForDate.add(ind);
				individualsPerDate.put(dob, individualsForDate);
			}
		}

		Collections.sort(allDates);

		Date start = allDates.get(0);
		ArrayList<Individual> rootInds = individualsPerDate.get(start);
		System.out.println("Root of the tree: " + allDates.get(0));
		for (int i = 0; i < rootInds.size(); i++) {
			System.out.println(rootInds.get(i).getFirstName() + "\t" + rootInds.get(i).getLastName() + "\t" + rootInds.get(i).getChildren().size() + " children");

		}

		Individual root = rootInds.get(0);

		int y = 0;
		int maxWidth = 0;
		boolean run = true;
		Individual current = root;
		ArrayList<Individual> individualsAtYLevel = rootInds;
		HashMap<Integer, HashSet<Individual>> individualsPerYLevel = new HashMap<Integer, HashSet<Individual>>();

		while (run) {
			// iterate all individuals at this level
			HashSet<Individual> remainingParentsAtThisLevel = new HashSet<Individual>();
			HashSet<Individual> childrenSet = new HashSet<Individual>();
			HashSet<Individual> knownParents = new HashSet<Individual>();
			knownParents.addAll(individualsAtYLevel);

			HashSet<Individual> visitedParents = new HashSet<Individual>();

			for (int i = 0; i < individualsAtYLevel.size(); i++) {
				// count the children
				Individual parent = individualsAtYLevel.get(i);
				visitedParents.add(parent);

				// System.out.println(parent.getId() + "\t" + parent.getFirstName() + "\t" + parent.getLastName() + "\t" + parent.getChildren().size() + " children");

				ArrayList<Individual> children = parent.getChildren();
				for (int j = 0; j < children.size(); j++) {
					Individual child = children.get(j);
					childrenSet.add(child);
					Individual mother = child.getMother();
					Individual father = child.getFather();
					visitedParents.add(parent);
					if (mother != null) {
						if (!mother.equals(parent)) {
							// we're not descending from the mother
							remainingParentsAtThisLevel.add(mother);
						}
					}

					if (father != null) {
						if (!father.equals(parent)) {
							// we're not descending from the father
							remainingParentsAtThisLevel.add(father);
						}
					}
				}
			}

			// System.out.println(remainingParentsAtThisLevel.size());

			while (!remainingParentsAtThisLevel.isEmpty()) {
				// we have some stray parents at this level. Check whether they have children
				HashSet<Individual> tmpRemainingParentsAtThisLevel = new HashSet<Individual>();
				for (Individual parent : remainingParentsAtThisLevel) {
					if (!visitedParents.contains(parent)) {
						visitedParents.add(parent);
						// System.out.println(parent.getId() + "\t" + parent.getFirstName() + "\t" + parent.getLastName() + "\t" + parent.getChildren().size() + " children");
						ArrayList<Individual> children = parent.getChildren();
						for (int j = 0; j < children.size(); j++) {
							Individual child = children.get(j);
							childrenSet.add(child);
							Individual mother = child.getMother();
							Individual father = child.getFather();
							if (mother != null) {
								if (!mother.equals(parent) && !visitedParents.contains(parent)) {
									// we're not descending from the mother
									tmpRemainingParentsAtThisLevel.add(mother);
								}
							}

							if (father != null) {
								if (!father.equals(parent) && !visitedParents.contains(parent)) {
									// we're not descending from the father
									tmpRemainingParentsAtThisLevel.add(father);
								}
							}
						}
					}
				}
				remainingParentsAtThisLevel = tmpRemainingParentsAtThisLevel;

			}

//			System.out.println(visitedParents.size() + " parents at level " + y);
//			System.out.println(childrenSet.size() + " children at level " + (y + 1));


			HashSet<Individual> inds = individualsPerYLevel.get(y);
			if (inds == null) {
				inds = new HashSet<Individual>();
			}
			inds.addAll(visitedParents);
			individualsPerYLevel.put(y, inds);
			inds = individualsPerYLevel.get(y + 1);
			if (inds == null) {
				inds = new HashSet<Individual>();
			}
			inds.addAll(childrenSet);
			individualsPerYLevel.put(y + 1, inds);

			individualsAtYLevel = new ArrayList<Individual>();
			individualsAtYLevel.addAll(childrenSet);

			if (childrenSet.isEmpty()) {
				run = false;
			}

			y++;
		}


		HashSet<Individual> allplacedIndividuals = new HashSet<Individual>();
		HashMap<Individual, Integer> individualToYLevel = new HashMap<Individual, Integer>();
		for (int i = 0; i < y; i++) {

			HashSet<Individual> inds = individualsPerYLevel.get(i);
			if (inds == null) {
				inds = new HashSet<Individual>();
			}

			for (Individual ind : inds) {
				individualToYLevel.put(ind, i);
			}
			allplacedIndividuals.addAll(inds);
			System.out.println(inds.size() + " individuals at level: " + y);
		}

		ArrayList<Individual> individualsNotPlaced = new ArrayList<Individual>();
		for (Individual ind : individuals) {
			if (!allplacedIndividuals.contains(ind)) {
				individualsNotPlaced.add(ind);
			}
		}


		System.out.println(individualsNotPlaced.size() + " not located properly: ");
		System.out.println("linking these people via relationships");

		for (Individual ind : individualsNotPlaced) {
			Individual father = ind.getFather();
			String fatherStr = "unknown";
			Individual mother = ind.getMother();
			String motherStr = "unknown";

			if (father != null) {
				fatherStr = father.toString();
			}
			if (mother != null) {
				motherStr = mother.toString();
			}


			ArrayList<Relationship> relationships = ind.getRelationships();

			boolean recovered = false;
			for (Relationship r : relationships) {
				Individual ind1 = r.individual1;
				Individual ind2 = r.individual2;
				if (ind.equals(ind1)) {
					Integer ylevel = individualToYLevel.get(ind2);
					if (ylevel != null) {
						individualsPerYLevel.get(ylevel).add(ind);
						individualToYLevel.put(ind, ylevel);
						recovered = true;
					}

				}
				if (ind.equals(ind2)) {
					Integer ylevel = individualToYLevel.get(ind1);
					if (ylevel != null) {
						individualsPerYLevel.get(ylevel).add(ind);
						individualToYLevel.put(ind, ylevel);
						recovered = true;
					}
				}
			}

			if (!recovered) {
				ArrayList<Individual> children = ind.getChildren();
				for (Individual child : children) {
					Integer ylevel = individualToYLevel.get(child);
					if (ylevel != null) {
						individualsPerYLevel.get(ylevel - 1).add(ind);
						individualToYLevel.put(ind, ylevel - 1);
						recovered = true;
					}
				}
				if (!recovered) {
					System.out.println("individual could not be placed: " + ind.toString() + "\t" + fatherStr + "\t" + motherStr);
				}
			}
		}


		int maxNrInds = 0;
		int total = 0;
		for (int i = 0; i < y; i++) {

			HashSet<Individual> inds = individualsPerYLevel.get(i);
			if (inds == null) {
				inds = new HashSet<Individual>();
			}
			total += inds.size();
			if (inds.size() > maxNrInds) {
				maxNrInds = inds.size();
			}
			System.out.println(inds.size() + " individuals at level: " + i);
		}
		System.out.println(total + " individuals");
		System.out.println(maxNrInds + " max size of generation");

		// give each individual an x-position


		// iterate the individuals
		// reorder them so that they are underneath their parents

		ArrayList<ArrayList<Individual>> finalOrder = new ArrayList<ArrayList<Individual>>();
		int maxNrFamilies = 0;
		for (int ylevel = 0; ylevel < y; ylevel++) {
			HashSet<Individual> individualsOnYLevel = individualsPerYLevel.get(ylevel);

			// order by birthdate
			ArrayList<Individual> orderedIndividuals = orderIndividuals(individualsOnYLevel);

			HashSet<Individual> visitedIndividuals = new HashSet<Individual>();
			if (ylevel == 0) {
				ArrayList<Individual> finalOrderAtLevel = new ArrayList<Individual>();
				int x = 0;
				int familyId = 0;
				for (int i = 0; i < orderedIndividuals.size(); i++) {
					Individual ind = orderedIndividuals.get(i);

					if (!visitedIndividuals.contains(ind)) {
						ind.setY(ylevel);
						ind.setX(x);
						ind.setFamilyId(familyId);
						x++;
						visitedIndividuals.add(ind);
						finalOrderAtLevel.add(ind);

						// find partner(s)
						ArrayList<Relationship> relationships = ind.getRelationships();
						for (int j = 0; j < relationships.size(); j++) {
							Relationship r = relationships.get(j);

							Individual ind1 = r.getIndividual1();
							Individual ind2 = r.getIndividual2();
							Individual partner = ind1;
							if (ind1.equals(ind)) {
								partner = ind2;
							}
							if (!visitedIndividuals.contains(partner)) {
								partner.setY(ylevel);
								partner.setX(x);
								partner.setFamilyId(familyId);
								x++;
								visitedIndividuals.add(partner);
								finalOrderAtLevel.add(partner);
							}
						}
						familyId++;

					}
				}
				finalOrder.add(finalOrderAtLevel);
			} else {
				// get nrOfIndividualsAt X


				// find the individuals that can be linked to previous generations
				HashMap<Integer, ArrayList<Individual>> families = new HashMap<Integer, ArrayList<Individual>>();
				for (int i = 0; i < orderedIndividuals.size(); i++) {
					Individual ind = orderedIndividuals.get(i);
					if (!visitedIndividuals.contains(ind)) {
						visitedIndividuals.add(ind);
						Individual father = ind.getFather();
						Individual mother = ind.getMother();
						Integer previousGenFamId = null;
						if (father != null) {
							int famId = father.getFamilyId();
							if (previousGenFamId == null || famId > previousGenFamId) {
								previousGenFamId = famId;
							}
						}

						if (mother != null) {
							int famId = mother.getFamilyId();
							if (previousGenFamId == null || famId > previousGenFamId) {
								previousGenFamId = famId;
							}
						}

						if (previousGenFamId != null) {
							ArrayList<Individual> members = families.get(previousGenFamId);
							if (members == null) {
								members = new ArrayList<Individual>();
							}
							members.add(ind);
							// put the partners next to the individual

							ArrayList<Relationship> relationships = ind.getRelationships();
							for (Relationship r : relationships) {
								Individual ind1 = r.getIndividual1();
								Individual ind2 = r.getIndividual2();
								Individual partner = ind1;
								if (ind1.equals(ind)) {
									partner = ind2;
								}
								if (!visitedIndividuals.contains(partner)) {
									// ind2 is the partner
									members.add(partner);
									visitedIndividuals.add(partner);
								}
							}
							families.put(previousGenFamId, members);
						}
					}
				}

				System.out.println();
				System.out.println(families.size() + " families at level: " + (ylevel - 1));
				if (families.size() > maxNrFamilies) {

				}

				Set<Integer> keys = families.keySet();
				ArrayList<Integer> keyList = new ArrayList<Integer>();
				for (Integer i : keys) {
					keyList.add(i);
				}
				Collections.sort(keyList);
				HashSet<Individual> individualsInFamilies = new HashSet<Individual>();
				int familyId = 0;
				int x = 0;
				ArrayList<Individual> finalOrderAtLevel = new ArrayList<Individual>();
				for (int q = 0; q < keyList.size(); q++) {
					// get the individuals per family, split them up in new families
					ArrayList<Individual> familyInds = families.get(keyList.get(q));
					System.out.println(familyInds.size() + " individuals in family " + keyList.get(q) + " at level: " + ylevel);
					for (int z = 0; z < familyInds.size(); z++) {
						Individual ind = familyInds.get(z);
						if (!individualsInFamilies.contains(ind)) {

							ind.setY(ylevel);
							ind.setX(x);
							ind.setFamilyId(familyId);
							x++;
							individualsInFamilies.add(ind);
							finalOrderAtLevel.add(ind);
							// find partner(s)
							ArrayList<Relationship> relationships = ind.getRelationships();
							for (int j = 0; j < relationships.size(); j++) {
								Relationship r = relationships.get(j);

								Individual ind1 = r.getIndividual1();
								Individual ind2 = r.getIndividual2();
								Individual partner = ind1;
								if (ind1.equals(ind)) {
									partner = ind2;
								}
								if (!individualsInFamilies.contains(partner)) {
									partner.setY(ylevel);
									partner.setX(x);
									partner.setFamilyId(familyId);
									x++;
									individualsInFamilies.add(partner);
									finalOrderAtLevel.add(partner);
								}
							}
							familyId++;

						}
					}
				}

				System.out.println(individualsInFamilies.size() + " individuals placed out of " + orderedIndividuals.size());
				// place the unplaced
				int unplaced = 0;
				for (int i = 0; i < orderedIndividuals.size(); i++) {
					Individual ind = orderedIndividuals.get(i);
					if (!individualsInFamilies.contains(ind)) {
						ind.setY(ylevel);
						ind.setX(x);
						ind.setFamilyId(familyId);
						x++;
						visitedIndividuals.add(ind);
						finalOrderAtLevel.add(ind);
						// find partner(s)
						ArrayList<Relationship> relationships = ind.getRelationships();
						for (int j = 0; j < relationships.size(); j++) {
							Relationship r = relationships.get(j);

							Individual ind1 = r.getIndividual1();
							Individual ind2 = r.getIndividual2();
							Individual partner = ind1;
							if (ind1.equals(ind)) {
								partner = ind2;
							}
							if (!individualsInFamilies.contains(partner)) {
								partner.setY(ylevel);
								partner.setX(x);
								partner.setFamilyId(familyId);
								x++;
								individualsInFamilies.add(partner);
								finalOrderAtLevel.add(partner);
							}
						}
						familyId++;
						unplaced++;

					}
				}
				System.out.println(unplaced + " unplaced individuals at level " + ylevel);
				finalOrder.add(finalOrderAtLevel);
			}
		}

		int boxwidth = 50;
		int boxmargin = 20;
		int width = (maxNrInds * boxwidth) + ((maxNrInds - 1) * boxmargin);

		int halfWidth = width / 2;
		System.out.println("max width: " + width);

		System.out.println();
		// now we can finally draw some stuff..
		int ylevelWithMaxNrInds = 0;
		maxNrInds = 0;
		for (int ylevel = 0; ylevel < finalOrder.size(); ylevel++) {
			ArrayList<Individual> inds = finalOrder.get(ylevel);

			if (inds.size() > maxNrInds) {
				maxNrInds = inds.size();
				ylevelWithMaxNrInds = ylevel;
			}
		}

		// determine position for each individual at max ylevel
		ArrayList<Individual> indsAtMaxYlevel = finalOrder.get(ylevelWithMaxNrInds);
		for (int i = 0; i < indsAtMaxYlevel.size(); i++) {
			Individual ind = indsAtMaxYlevel.get(i);
			int position = (i * boxwidth) + (i * boxmargin);
			ind.setXPixel(position);
		}

		// go up
		// place parent relative to children
		for (int ylevel = ylevelWithMaxNrInds - 1; ylevel > -1; ylevel--) {

			ArrayList<Individual> inds = finalOrder.get(ylevel);

			int ylevelxwidth = ((inds.size() - 1) * boxmargin) + (inds.size() * boxwidth);
			int ylevelmidpoint = halfWidth - (ylevelxwidth / 2);

			HashSet<Individual> individualsVisited = new HashSet<Individual>();
			for (int i = 0; i < inds.size(); i++) {
				Individual ind = inds.get(i);

				if (!individualsVisited.contains(ind)) {
					individualsVisited.add(ind);
					// get partners
					ArrayList<Individual> partnerlist = getPartners(ind);
					partnerlist = orderIndividuals(partnerlist);
					// now we know the width of this family


					ArrayList<Individual> children = ind.getChildren();
					if (children == null || children.isEmpty()) {
						// person has no children
						// check whether any partner has children
						int minX = Integer.MAX_VALUE;
						int maxX = -Integer.MAX_VALUE;
						HashSet<Individual> visitedChildren = new HashSet<Individual>();
						for (int p = 0; p < partnerlist.size(); p++) {
							Individual partner = partnerlist.get(p);
							ArrayList<Individual> partnerchildren = partner.getChildren();
							for (int c = 0; c < partnerchildren.size(); c++) {
								Individual child = partnerchildren.get(c);
								if (!visitedChildren.contains(child)) {
									int pixel = child.getXPixel();

									if (pixel > maxX) {
										maxX = pixel;
									}
									if (pixel < minX) {
										minX = pixel;
									}
									visitedChildren.add(child);
								}
							}
						}

						if (visitedChildren.isEmpty()) {
							// no children for individual
							// draw person right next to the last known x position

							int indposition = ylevelmidpoint - (ylevelxwidth / 2) + (i * boxwidth) + (i * boxmargin);
							ind.setXPixel(indposition);
							for (int p = 0; p < partnerlist.size(); p++) {
								Individual partner = partnerlist.get(p);
								int position = indposition + (((p + 1) * boxwidth) + ((p + 1) * boxmargin));
								partner.setXPixel(position);
								individualsVisited.add(partner);
							}
						} else {
							// draw individual and partners relative to midpoint of children
							// determine midpoint
							int midpoint = (maxX - minX) / 2;
							// determine width of family
							int nrIndsInFam = partnerlist.size() + 1;
							int widthOfFam = (nrIndsInFam * boxwidth) + ((nrIndsInFam - 1) * boxmargin);
							int startX = midpoint - (widthOfFam / 2);

							ind.setXPixel(startX);
							for (int p = 0; p < partnerlist.size(); p++) {
								Individual partner = partnerlist.get(p);
								int position = startX + (((p + 1) * boxwidth) + ((p + 1) * boxmargin));
								partner.setXPixel(position);
								individualsVisited.add(partner);
							}
						}


					} else {
						// get x position of children
						int minX = Integer.MAX_VALUE;
						int maxX = -Integer.MAX_VALUE;


						HashSet<Individual> visitedChildren = new HashSet<Individual>();
						for (int c = 0; c < children.size(); c++) {
							Individual child = children.get(c);
							int pixel = child.getXPixel();

							if (pixel > maxX) {
								maxX = pixel;
							}
							if (pixel < minX) {
								minX = pixel;
							}
							visitedChildren.add(child);
						}

						// check whether any of the partners had any unaccounted for children
						for (int p = 0; p < partnerlist.size(); p++) {
							Individual partner = partnerlist.get(p);
							ArrayList<Individual> partnerchildren = partner.getChildren();
							for (int c = 0; c < partnerchildren.size(); c++) {
								Individual child = partnerchildren.get(c);
								if (!visitedChildren.contains(child)) {
									int pixel = child.getXPixel();

									if (pixel > maxX) {
										maxX = pixel;
									}
									if (pixel < minX) {
										minX = pixel;
									}
									visitedChildren.add(child);
								}
							}
						}

						// determine midpoint
						int midpoint = (maxX - minX) / 2;
						// determine width of family
						int nrIndsInFam = partnerlist.size() + 1;
						int widthOfFam = (nrIndsInFam * boxwidth) + ((nrIndsInFam - 1) * boxmargin);
						int startX = midpoint - (widthOfFam / 2);

						ind.setXPixel(startX);
						for (int p = 0; p < partnerlist.size(); p++) {
							Individual partner = partnerlist.get(p);
							int position = startX + (((p + 1) * boxwidth) + ((p + 1) * boxmargin));
							partner.setXPixel(position);
							individualsVisited.add(partner);
						}
					}
				}


//				int position = (i * boxwidth) + (i * margin);
//				ind.setXPixel(position);
			}
		}

		// go down
		// place children relative to parents
		for (int ylevel = ylevelWithMaxNrInds + 1; ylevel < finalOrder.size(); ylevel++) {
			ArrayList<Individual> inds = finalOrder.get(ylevel);
			int ylevelxwidth = ((inds.size() - 1) * boxmargin) + (inds.size() * boxwidth);
			int ylevelmidpoint = halfWidth - (ylevelxwidth / 2);
			HashSet<Individual> individualsVisited = new HashSet<Individual>();

			// group individuals and their partners based on parents
			for (int i = 0; i < inds.size(); i++) {
				Individual ind = inds.get(i);
				if (!individualsVisited.contains(ind)) {
// get brothers and sisters
					// including half brothers and sisters

					HashSet<Individual> siblings = new HashSet<Individual>();
					Individual father = ind.getFather();
					Individual mother = ind.getMother();
					if (father != null) {
						siblings.addAll(getChildrenFromParent(father));
					}
					if (mother != null) {
						siblings.addAll(getChildrenFromParent(mother));
					}

					if (!siblings.contains(ind)) {
						siblings.add(ind);
					}

					ArrayList<Individual> allChildren = new ArrayList<Individual>();
					// add all the partners as well
					for (Individual sibling : siblings) {
						ArrayList<Individual> partners = getPartners(sibling);
						allChildren.addAll(partners);
						allChildren.add(sibling);
					}

					// order by birthdate
					allChildren = orderIndividuals(allChildren);

					// get lowest and highest X for parent
					int minX = Integer.MAX_VALUE;
					int maxX = -Integer.MAX_VALUE;
					for (int j = 0; j < allChildren.size(); j++) {
						Individual child = allChildren.get(j);
						Individual childFather = child.getFather();
						Individual childMother = child.getMother();
						if (childFather != null) {
							if (childFather.getXPixel() > maxX) {
								maxX = childFather.getXPixel();
							} else if (childFather.getXPixel() < minX) {
								minX = childFather.getXPixel();
							}
						}
						if (childMother != null) {

							if (childMother.getXPixel() > maxX) {
								maxX = childMother.getXPixel();
							} else if (childMother.getXPixel() < minX) {
								minX = childMother.getXPixel();
							}
						}
					}

					// determine midpoint
					int midpoint = (maxX - minX) / 2;
					// determine width of family
					int nrIndsInFam = allChildren.size();
					int widthOfFam = (nrIndsInFam * boxwidth) + ((nrIndsInFam - 1) * boxmargin);
					int startX = midpoint - (widthOfFam / 2);

					// position the children
					for (int p = 0; p < allChildren.size(); p++) {
						Individual child = allChildren.get(p);
						int position = startX + ((p * boxwidth) + (p * boxmargin));
						child.setXPixel(position);
						individualsVisited.add(child);
					}
				}


			}

		}

		// all individuals should now have an X and an Y coordinate..

		int yspacing = 100;
		int height = y * yspacing;
		height += yspacing;
		try {
			Plot plot = new Plot(treeInfo + ".pdf", width, height);

			for (int ylevel = 0; ylevel < y; ylevel++) {
				int pixelY = ylevel * yspacing;
				ArrayList<Individual> inds = finalOrder.get(ylevel);
				for (int i = 0; i < inds.size(); i++) {
					Individual ind = inds.get(i);
					Individual father = ind.getFather();
					Individual mother = ind.getMother();
					if (father != null) {
						// draw line to dad
					}
					if (mother != null) {
						// draw line to mom
					}

					plot.drawInd(ind, pixelY, boxwidth);


				}
			}

			plot.close();
		} catch (DocumentException e) {
			e.printStackTrace();
		}

	}

	private ArrayList<Individual> getPartners(Individual individual) {
		ArrayList<Relationship> relationships = individual.getRelationships();
		HashSet<Individual> partnerset = new HashSet<Individual>();

		for (Relationship r : relationships) {
			Individual ind1 = r.getIndividual1();
			if (ind1.equals(individual)) {
				partnerset.add(r.getIndividual2());
			} else {
				partnerset.add(r.getIndividual1());
			}
		}

		ArrayList<Individual> partners = new ArrayList<Individual>();
		for (Individual i : partnerset) {
			partners.add(i);
		}
		return partners;
	}

	private HashSet<Individual> getChildrenFromParent(Individual parent) {

		HashSet<Individual> childrenOut = new HashSet<Individual>();
		ArrayList<Individual> children = parent.getChildren();
		childrenOut.addAll(children);
		ArrayList<Relationship> relationships = parent.getRelationships();

		for (int r = 0; r < relationships.size(); r++) {
			Individual partner = relationships.get(r).getIndividual1();
			if (partner.equals(parent)) {
				partner = relationships.get(r).getIndividual2();
			}

			ArrayList<Individual> childrenOfPartner = partner.getChildren();
			childrenOut.addAll(childrenOfPartner);

		}

		return childrenOut;
	}

	private ArrayList<Individual> orderIndividuals(HashSet<Individual> inds) {
		ArrayList<Individual> individuals = new ArrayList<Individual>();
		individuals.addAll(inds);
		return orderIndividuals(individuals);
	}

	private ArrayList<Individual> orderIndividuals(ArrayList<Individual> individuals) {
		ArrayList<Date> allDates = new ArrayList<Date>();

		HashMap<Date, ArrayList<Individual>> individualsPerDate = new HashMap<Date, ArrayList<Individual>>();
		for (int i = 0; i < individuals.size(); i++) {
			Individual ind = individuals.get(i);
			Date dob = ind.getDateofbirth();

			if (dob == null) {
				dob = new Date(System.currentTimeMillis());
			}
			if (!individualsPerDate.containsKey(dob)) {
				allDates.add(dob);
			}

			ArrayList<Individual> individualsForDate = individualsPerDate.get(dob);
			if (individualsForDate == null) {
				individualsForDate = new ArrayList<Individual>();
			}
			individualsForDate.add(ind);
			individualsPerDate.put(dob, individualsForDate);


		}

		Collections.sort(allDates);

		ArrayList<Individual> ordered = new ArrayList<Individual>();
		for (int d = 0; d < allDates.size(); d++) {
			ArrayList<Individual> indsOnDate = individualsPerDate.get(allDates.get(d));

			ordered.addAll(indsOnDate);

		}

		return ordered;
	}


	private ArrayList<Individual> parseTreeInfo(String file) throws IOException {

	/*
	Nummer
	Naam
	Voornaam
	Geslacht
	Geboren
	Geboorteplaats
	Land
	Gestorven
	Sterfplaats
	Vader
	Moeder
	Gehuwd met
	Huwelijksdatum
	Woonplaats
	Land
	Gescheiden
	Generatie
	 */

		ArrayList<Individual> individuals = new ArrayList<Individual>();

		TextFile tf = new TextFile(file, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(Strings.tab);

		while (elems != null) {
			Integer id = Integer.parseInt(elems[0]);


			String name = elems[1];
			String firstname = elems[2];
			String gender = elems[3];

			Date dodStr = null;
			Date dobStr = null;
			Location placeOfBirth = null;
			Location placeOfDeath = null;
			Location placeOfLiving = null;
			if (elems.length > 4) {
				dobStr = formatDate(elems[4]);
			}
			if (elems.length > 8) {
				String placeOfBirthStr = elems[5];
				String countryOfBirth = elems[6];


				placeOfBirth = new Location(placeOfBirthStr, countryOfBirth);

				dodStr = formatDate(elems[7]);
				String placeOfDeathStr = elems[8];
				placeOfDeath = new Location(placeOfDeathStr, null);


				if (elems.length > 14) {
					String placeOfLivingStr = elems[13];
					String countryOfLivingStr = elems[14];
					placeOfLiving = new Location(placeOfLivingStr, countryOfLivingStr);
				}


			} else {
				System.err.println("No relationship info: " + Strings.concat(elems, Strings.tab));
			}

			System.out.println("Adding: " + firstname + "\t" + name + "\t" + dobStr);


			Individual i = new Individual(
					id,
					firstname,
					name,
					Individual.Gender.parseGender(gender),
					dobStr,
					dodStr,
					placeOfBirth,
					placeOfDeath,
					placeOfLiving

			);
			individuals.add(i);
			elems = tf.readLineElems(Strings.tab);
		}

		tf.close();


		// index
		HashMap<Integer, Individual> intToInd = new HashMap<Integer, Individual>();
		for (Individual i : individuals) {
			intToInd.put(i.getId(), i);
		}

		// add relationships and children
		tf.open();
		tf.readLine();
		elems = tf.readLineElems(Strings.tab);

		while (elems != null) {
			Integer id = Integer.parseInt(elems[0]);

			Individual ind = intToInd.get(id);

			Individual father = null;
			Individual mother = null;

			if (elems.length > 9) {
				try {
					Integer fid = Integer.parseInt(elems[9]);
					Individual fathInd = intToInd.get(fid);
					father = fathInd;
					ind.setFather(fathInd);
					fathInd.addChild(ind);
				} catch (NumberFormatException e) {

				}
			}
			if (elems.length > 10) {
				try {
					Integer mid = Integer.parseInt(elems[10]);
					Individual mothInd = intToInd.get(mid);
					mother = mothInd;
					ind.setMother(mothInd);
					mothInd.addChild(ind);
				} catch (NumberFormatException e) {

				}
			}

			if (father != null && mother != null) {
				// add relationship based on child
				Relationship r = new Relationship(father, mother, null, null, Relationship.KIND.NOTMARRIED);
				father.addRelationship(r);
				mother.addRelationship(r);
			}


			if (elems.length > 12) {
				String spouses = elems[11];
				if (spouses.trim().length() > 0) {
					String relationshipdates = elems[12];
					String[] spouseelems = spouses.split(",");
					String[] relationShipDateElems = relationshipdates.split(",");
					int q = 0;
					for (String spouse : spouseelems) {
						Integer iid = Integer.parseInt(spouse);

						Individual iiiiid = intToInd.get(iid);
						Date iidd = null;
						if (relationShipDateElems.length == spouseelems.length) {
							iidd = formatDate(relationShipDateElems[q]);
						}
						Relationship s = new Relationship(ind, iiiiid, iidd, null, Relationship.KIND.MARRIED);
						ind.addRelationship(s);
						q++;
					}
				}
			}


			elems = tf.readLineElems(Strings.tab);
		}

		tf.close();


		return individuals;
	}

	private Date formatDate(String s) {

		String[] elems = s.split("-");
		if (elems.length > 1) {
			if (elems[elems.length - 1].length() == 2) {
				elems[elems.length - 1] = "19" + elems[elems.length - 1];
			}
			SimpleDateFormat formatter = new SimpleDateFormat("dd-MM-yyyy");
			try {
				return formatter.parse(s);
			} catch (ParseException e) {
				e.printStackTrace();
			}

		}
		return null;


	}

}
