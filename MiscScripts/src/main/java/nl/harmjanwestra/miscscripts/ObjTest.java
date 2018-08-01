package nl.harmjanwestra.miscscripts;

/**
 * Created by hwestra on 6/13/16.
 */
public class ObjTest {

	public static void main(String[] args) {
		ObjTest t = new ObjTest();
		t.run();
	}

	public void run() {

		int[] array = new int[5];
		Obj obj = new Obj(array);
		System.out.println(obj.getSize());
		array = new int[10];
		System.out.println(obj.getSize());

	}


	public class Obj {
		int[] array;

		public Obj(int[] arr) {
			this.array = arr;
		}

		public int getSize() {
			return array.length;
		}
	}
}
