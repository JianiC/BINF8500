import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;

public class cluster {
	public static void main(String[] args) throws FileNotFoundException{	
		File text = new File("Bacteria+Archaea.txt");
		Scanner scnr = new Scanner(text);
		PrintWriter outFile = new PrintWriter("k-means.txt");
		outFile.println("k\tWCSS\tAIC\tBIC");
		// read the fastq file into 2D array`
				int r = 0;
				String[]temp1=new String[131];
				while(scnr.hasNextLine()) {
					temp1[r]=scnr.nextLine();
					r++;
				}
		//discard first row
				String[]temp1_2=new String[130];
				for(int n=0;n<130;n++) {
					temp1_2[n]=temp1[n+1];
				}
								
		// covert1d array into 2d:spilt a list of string into array
				String[][]temp2=new String[temp1_2.length] [21];
				
				for(int i=0; i<temp1_2.length;i++) {
					for(int j=0;j<temp1_2.length;j++) {
						temp2[i]=temp1_2[i].split("\t",21);						
					}					
				}
			
		//discrad the first column
				String[][] temp2_2=new String[temp2.length][temp2[0].length-1];
				
				for(int i=0; i<temp2_2.length;i++) {
					for(int j=0;j<temp2_2[0].length;j++) {
						temp2_2[i][j]=temp2[i][j+1];				
					}					
				}
				
	
		//convert string to double
				double content[][]=new double[temp2_2.length][temp2_2[0].length];
				for(int i=0; i<temp2_2.length;i++) {
					for(int j=0;j<temp2_2[0].length;j++) {
						content[i][j]=Double.parseDouble(temp2_2[i][j]);						
					}
				}
								
		//content nomalize
				double newcontent[][]=new double [content.length][content[0].length];
				newcontent=normalize(content);
				
				
				for(int k=1;k<15;k++) {
					
				double WCSS[] =new double[300];
				//get min with 300 times
				for(int t=0;t<300;t++) {
				
				double c[][]=new double[k][content[0].length];
				c=random(k,newcontent);
				
				// get min_distance
				double tempmin;
				
				double[] min_d= new double [newcontent.length];
				int []cluster=new int [newcontent.length];
				
				// initialize min_distance as distance to c[0]
				for(int p_id=0;p_id<newcontent.length;p_id++) {					
					min_d[p_id]=distance(p_id,0,newcontent,c);
					}
							
				for(int p_id=0;p_id<newcontent.length;p_id++) {
					for(int c_id=0;c_id<c.length;c_id++) {
						tempmin=distance(p_id,c_id,newcontent,c);
						if(tempmin<=min_d[p_id]) {
							min_d[p_id]=tempmin;
							cluster[p_id]=c_id;
					}					
				}
			}				
		
				double [][]c_sum=new double[k][20];
				for(int i=0;i<c_sum.length;i++) {
					for(int j=0;j<20;j++) {
						c_sum[i][j]=0;
					}					
				}
				
				
				double [][]c_mean=new double[k][30];				
				double[] newmin_d= new double [newcontent.length];			
				int []newcluster=new int[newcontent.length];
				double new_tempmin;				
				Boolean stopchange= false;
				
				
				//get minimum WCSS with 300 times
		
				
				while(stopchange==false){ 
					
					for(int i=0;i<k;i++) {
						for (int j=0; j<newcontent[0].length;j++) {
							double sum=0;
							int count=0;
							for(int x=0; x<newcontent.length;x++) {
								if(cluster[x]==i) {
									sum+=newcontent[x][j];
									count++;
								}
							}
							c_mean[i][j]=sum/count;
						}
					}
					for(int p_id=0;p_id<newcontent.length;p_id++) {
						newmin_d[p_id]=distance(p_id,0,newcontent,c_mean);
						}
						
					for(int p_id=0;p_id<newcontent.length;p_id++) {
						for(int c_id=0;c_id<c.length;c_id++) {
							new_tempmin=distance(p_id,c_id,newcontent,c_mean);
							if(new_tempmin<=newmin_d[p_id]) {
								newmin_d[p_id]=new_tempmin;
								newcluster[p_id]=c_id;
							}												
						}	
					
					}		
					
				/*	
					System.out.println("\nnewcluster assignment");
					for(int i=0;i<cluster.length;i++) {
						System.out.print(newcluster[i]);
					}
					System.out.println();
					
					*/
					stopchange=Arrays.equals(newcluster,cluster);																
					cluster=newcluster.clone();
	
				}
				
					
				for(int i=0;i<newcontent.length;i++) {
					WCSS[t]+=newmin_d[i];
				   }
				//System.out.println(WCSS[t]);				
				}
				// get minmum WCSS[t]
				
				double min_WCSS=WCSS[0];
				for (int i=0; i<300;i++) {
					if (WCSS [i] < min_WCSS) {
						min_WCSS = WCSS[i];
					}
				}
				
				
				//calculate  AIC/BIC				
				double AIC=2*k*content[0].length+min_WCSS;
				double BIC=Math.log(content.length)*k*content[0].length+min_WCSS;
				
				System.out.println(k+"\t"+min_WCSS+"\t"+AIC+"\t"+BIC);				
				outFile.println(k+"\t"+min_WCSS+"\t"+AIC+"\t"+BIC);	
			
				}
				scnr.close();
				outFile.close();		
			
	}
	
	
	
	
	
	
	
	
	

	public static double[][] normalize(double[][]array) {
		double var[]=new double[array.length];
		double std[]=new double[array.length];
	
// initialize sum=0 for each dimension
		double column_mean[]= new double[array[0].length];
		double psum[]=new double[array[0].length];
		for(int i=0;i<array[0].length;i++) {
			psum[i]=0;
			var[i]=0;
		}
		

		for(int i=0;i<array[0].length;i++) {			
			for(int j=0;j<array.length;j++) {
			psum[i]=psum[i]+array[j][i];
		}
		column_mean[i]=psum[i]/array.length;			
	}
	
		for(int i=0;i<array[0].length;i++) {
			for(int j=0;j<array.length;j++) {
			var[i]= var[i]+(array[j][i]-column_mean[i])*(array[j][i]-column_mean[i]);			
		}
		std[i]=Math.sqrt(var[i]/array.length);
	}
		double normalize[][]=new double[array.length][array[0].length];
		for(int i=0;i<array.length;i++) {
			for(int j=0;j<array[0].length;j++) {
				normalize[i][j]=(array[i][j]-column_mean[j])/std[j];			
		}			
	}
	return normalize;
}


//k means clustering helper
	
//random to create a list filled with point index
	public static  double[][]random ( int k,double[][]array) {
		int id[]=new int [k];
		Random random=new Random();
		for(int i=0;i<k;i++) {
			int randomInt=random.nextInt(array.length);
			id[i]=randomInt;
	}
	
//create initial centroid list using the index list
	double point[][]=new double[k][array[0].length];
	for(int i=0;i<k;i++) {
		for(int j=0;j<array[0].length;j++) {
			point[i][j]=array[id[i]][j];				
		}			
	}
	return point;
}
	
//assignment helper to get distance 
	public static double distance(int p_id, int c_id,double[][]p,double [][]c) {
		double distance =0;
		for(int j=0;j<p[0].length;j++) {
			distance=distance+ (p[p_id][j]-c[c_id][j])*(p[p_id][j]-c[c_id][j]);	
	}
	return distance;
}

}
