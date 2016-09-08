/* Author: Riku Katainen @ University of Helsinki
 * 
 * Tumor Genomics Group (http://research.med.helsinki.fi/gsb/aaltonen/) 	
 * Contact: riku.katainen@helsinki.fi / help@baseplayer.fi
 * 
 * LICENSE: 
 * 
 * GNU AFFERO GENERAL PUBLIC LICENSE
 * Version 3, 19 November 2007
 *  
 */
package base.BasePlayer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;



public class Launcher  {	
	static String maindir;
	static boolean fromMain = false, firstStart = false;
//	private static long timer;	
	static String memlimit = "4G", defaultDir = "demo", line, ctrldir = "demo", defaultAnnotation ="", defaultGenome = "";
	static String gerpfile;	
	static ArrayList<String> config = new ArrayList<String>();
	public static String trackDir = "";
	public static String projectDir= "", downloadDir = "";
	
	public static void main(String[] args) {
		try {
		
	//	Logo.main(argsit);
	//	 Logo.frame.dispatchEvent(new WindowEvent(Logo.frame, WindowEvent.WINDOW_CLOSING));
	//	timer = System.currentTimeMillis();
		maindir = new File(Launcher.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParent().replace("%20", " ");	
		
		  if(new File(maindir +"/config.txt").exists()) {
			
			  BufferedReader fileReader = new BufferedReader(new FileReader(maindir +"/config.txt"));				 
			 
			  
			  while((line = fileReader.readLine()) != null) {
				  config.add(line);
				  if(line.startsWith("Memory")) {
					//  if(!line.contains(memlimit)) {
						  memlimit = line.substring(line.indexOf("=")+1).replace(" ", "");
					 // }
				  }
				  else if(line.startsWith("DefaultDir")) {
				//	  if(!line.contains(defaultDir)) {
						  defaultDir = line.substring(line.indexOf("=")+1).replace(" ", "");
				//	  }				  
				  }
				  else if(line.startsWith("DefaultControlDir")) {
					  ctrldir = line.substring(line.indexOf("=")+1).replace(" ", "");
				  }
				  else if(line.startsWith("GERP")) {
					  gerpfile = line.substring(line.indexOf("=")+1).replace(" ", "");
				  }
				  else if(line.startsWith("DefaultGenome")) {
					  defaultGenome = line.substring(line.indexOf("=")+1).replace(" ", "");
				  }
				  else if(line.startsWith("DefaultGenes")) {
					  defaultAnnotation = line.substring(line.indexOf("=")+1).replace(" ", "");
				  }
				  else if(line.startsWith("DefaultProject")) {
					  projectDir = line.substring(line.indexOf("=")+1).replace(" ", "");
				  }
				  else if(line.startsWith("DefaultTrack")) {
					  trackDir = line.substring(line.indexOf("=")+1).replace(" ", "");
				  }
				  else if(line.startsWith("DownloadDir")) {
					  downloadDir = line.substring(line.indexOf("=")+1).replace(" ", "");
				  }
				  else if(line.startsWith("FirstStart")) {
					  if(line.substring(line.indexOf("=")+1).replace(" ", "").contains("rue")) {
						  firstStart = true;
					  }					 
				  }
			  }
			  if(!new File(defaultDir).exists()) {
				  defaultDir = "";
			  }
			  fileReader.close();	
			  if(!fromMain) {
				  
				  ProcessBuilder pb;
				 try {
					
					  pb = new ProcessBuilder("java","-Xmx" +memlimit, "-Dprism.verbose=true", "-Dsun.java2d.d3d=false","-jar", maindir +"/BasePlayer.jar");				
				     pb.start();
			//      Logo.frame.dispatchEvent(new WindowEvent(Logo.frame, WindowEvent.WINDOW_CLOSING));
				 /*    
				     
					 BufferedReader read = new BufferedReader(new InputStreamReader(process.getErrorStream()));
					
					 if(read.readLine() !=null) {
					
						 process.destroy();
						
						 read.close();
						 pb = new ProcessBuilder("java","-Xmx1200m","-Dprism.verbose=true",  "-Dsun.java2d.d3d=false","-jar", maindir +"/BasePlayer.jar");				
						 pb.start();
						 
					      
					 }
					 else {						
						 read.close();					      
					 }
					*/
				 }
				 catch(Exception e) {
					 e.printStackTrace();
					
				 }
				
		    	/* while(System.currentTimeMillis() - timer < 3000) {			    		 
		    	 }			    	
		    	Logo.frame.dispatchEvent(new WindowEvent(Logo.frame, WindowEvent.WINDOW_CLOSING));
			      */
			      
			     
			  }
			  else {
			//	  Main.path = defaultDir;
			//	  Control.path = ctrldir;
			  }
		  }
		  else {
			  if(!fromMain) {
				  ProcessBuilder pb;
					 try {
						  pb = new ProcessBuilder("java","-Xmx4G", "-Dprism.verbose=true", "-Dsun.java2d.d3d=false","-jar", maindir +"/BasePlayer.jar");				
					     pb.start();
						
						
					 }
					 catch(Exception e) {
						 e.printStackTrace();
						
					 }			      
			     
			  }
		  }		  
			
		}
		catch(Exception e) {
		//	JOptionPane.showMessageDialog(Main.chromDraw, e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
			//e.printStackTrace();
			
			if(!fromMain) {
				try {
			//		ProcessBuilder pb;					
			//		pb = new ProcessBuilder("java","-Xmx1200m","-Dprism.verbose=true",  "-Dsun.java2d.d3d=false", "-jar", maindir +"/BasePlayer.jar");					
			//		pb.start();					
				}
				catch(Exception ex) {
					ex.printStackTrace();
				}
			 }
			// System.exit(0);
		}
		
	
	}	
	
	/*public static class ProcMon implements Runnable {

		  private final Process _proc;
		  private volatile boolean _complete;

		  public boolean isComplete() { return _complete; }

		  public void run() {
			  try {
				 
				  _proc.waitFor();
				 
				  _complete = true;
			  }
			  catch(Exception e) {
				  e.printStackTrace();
			  }
		  }
		 
		  public ProcMon (Process proc) {
			 
		    this._proc = proc;
		  
		    Thread t = new Thread();
		    t.start();
		   
		  }
		}*/

	
}
