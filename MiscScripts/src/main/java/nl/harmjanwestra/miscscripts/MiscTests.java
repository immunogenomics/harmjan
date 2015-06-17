/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package nl.harmjanwestra.miscscripts;

import java.util.regex.Pattern;

/**
 * @author hwestra
 */
public class MiscTests {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
//        int lastI = 0;
////        for(int i=0;i<123456789; i+=100000){
////            System.out.println(i);
////            lastI = i;
////        }
////        System.out.println(lastI);
////        int a = 0x1;
////        for(int i=0;i<10; i++){
////
//////            System.out.println((a << 8)+"\t"+Integer.toBinaryString(a << 8)+"\t"+Integer.bitCount(a<<8)+"\t"+Integer.);
////        }
////        int v = 1<<5;
//        DecimalFormat f = new DecimalFormat("#.##");
//
//        double a = 0.999;
//        double b = 0.899;
//        double c = 99.999;
//        double d = 999.99;
//        System.out.println(f.format(a));
//        System.out.println(f.format(b));
//        System.out.println(f.format(c));

		String test = "split|this|string";
		Pattern p = Pattern.compile("\\|");
		String[] elems = p.split(test);
		System.out.println(elems.length);
		p = Pattern.compile("|");
		elems = p.split(test);
		System.out.println(elems.length);

	}


}
