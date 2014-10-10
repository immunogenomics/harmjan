/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.ngs.containers;

/**
 *
 * @author hwestra
 */
public class ReadGroup {
    public String name;
    public int[] readsPerChr = new int[26];
    public int[] readsPairedPerChr = new int[26];
    
    public int[][] mapqPerChr = new int[26][100];
    public int[] dupsPerChr = new int[26];
    public int[][] baseQualPerPos = new int[100][60]; 
    int nrReads = 0;

    
}
