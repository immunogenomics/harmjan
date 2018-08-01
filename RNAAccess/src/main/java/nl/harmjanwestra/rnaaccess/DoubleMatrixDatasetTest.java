/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package nl.harmjanwestra.rnaaccess;

import java.io.IOException;
import java.util.Arrays;
import nl.harmjanwestra.utilities.legacy.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author Harm-Jan
 */
public class DoubleMatrixDatasetTest {
    public static void main(String[] args) throws Exception{
        try{
            DoubleMatrixDataset ds = new DoubleMatrixDataset();
            double[][] matrix=  new double[10][50];
            String[] rowNames = new String[10];
            String[] colNames = new String[50];
            for(int t=0;t<10;t++){
                rowNames[t] = "row"+t;
                for(int q=0;q<50;q++){
                    colNames[q] = "col"+q;
                    matrix[t][q] = Math.random();
                }
            }
            ds.setMatrix(matrix);
            ds.setRowObjects(Arrays.asList(rowNames));
            ds.setColObjects(Arrays.asList(colNames));
            ds.save("c:\\work\\test.txt");
            
        } catch (IOException e){
            
        }
    }
}
