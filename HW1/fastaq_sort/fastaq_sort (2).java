/*Title: BINF 8500 HW1
 * Purpose: Implement the quick sort algorithm to sort entries in a fastq file alphabetically by the nucleotide sequence of each read.
 * Author: Jiani Chen
 * To run the code, fastq file directory need to type in the main method.
 */


import java.io.*;
import java.util.Scanner;
public class fastaq_sort {
	// Quick Sort method
	
	public static void quickSort(String[][]list) {
		quickSort(list,0,list.length - 1);
	}
	
	//quick sort helper
	public static void quickSort(String[][] list, int first, int last) {
		if (last>first) {
			int pivotIndex=partition(list,first,last);
			quickSort(list,first,pivotIndex-1);
			quickSort(list,pivotIndex+1,last);
						
		}
	}
	
	public static int partition(String[][] list, int first, int last) {
		String pivot =list[first][1];//choose the first element as pivot
		String pivot0 =list[first][0];
		String pivot3 =list[first][3];
		int low =first +1;
		int high =last;
		while(high > low) {
			while (low<=high && list[low][1].compareTo(pivot)<=0) {
				low++;
				
			}
			while (low<=high && list[high][1].compareTo(pivot)>0) {
				high--;
			}
			if (high>low) {
				String temp =list [high][1];
				list[high][1]=list[low][1];
				list[low][1]=temp;
				
				//also swap the first line and the forth line
				String temp1 =list [high][0];
				list[high][0]=list[low][0];
				list[low][0]=temp1;
				
				String temp4 =list [high][3];
				list[high][3]=list[low][3];
				list[low][3]=temp4;
			}
		}
		while (high>first && list[high][1].compareTo(pivot)>0) {
			high--;
		}
		if (list[high][1].compareTo(pivot)<0) {
			list[first][1]=list[high][1];
			list[high][1]=pivot;
			
			list[first][0]=list[high][0];
			list[high][0]=pivot0;
			
			list[first][2]=list[high][2];
			list[high][3]=pivot3;
						
			return high;
		
		}
		else {
			return first;
		}
	}
	
	
	
	

	public static void main(String[] args) throws FileNotFoundException{
		//creat input and output file
		File text = new File("sample1M.fastq");
		Scanner scnr = new Scanner(text);
		PrintWriter outFile = new PrintWriter("sort.txt");
		
		// read the fastq file into 2D array`
		String[][]fasq= new String[1000000][4];// declare the length of array according to the read
		int r=0;// start from the first read
		while(scnr.hasNextLine()) {
			fasq[r][0]=scnr.nextLine();
			fasq[r][1]=scnr.nextLine();
			fasq[r][2]=scnr.nextLine();
			fasq[r][3]=scnr.nextLine();
			
			r++;
		}
		// sort the file using quicksort method and write into the output file

		quickSort(fasq);
		
		for (int read =0; read<fasq.length;read++) {
			for(int line=0;line<4;line++) {
				outFile.println(fasq[read][line]);
				
			}
		}
		
			outFile.print("\n");
			
		
	
	outFile.close();
	scnr.close();

}

	}


