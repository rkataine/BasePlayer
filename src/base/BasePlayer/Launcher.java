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
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;

public class Launcher  {	
	static String maindir;
	static boolean fromMain = false, firstStart = false;
	static File configfile;
//	private static long timer;	
	static String memlimit = "1200m", defaultDir = "demo", line, ctrldir = "demo", defaultAnnotation ="", defaultGenome = "", defaultSaveDir = "", genomeDir = "", fontSize = "", backColor ="";
	static int wallpaperIndex = 0, alphaValue = 0;
	static String proxyHost ="", proxyPort ="", proxyType = "";
	static boolean isProxy = false;
	static String gerpfile;	
	static ArrayList<String> config = new ArrayList<String>();
	public static String trackDir = "";
	public static String projectDir= "", downloadDir = "";
	static String lineSeparator;
	static Boolean exampledata = false;
	public static void main(String[] args) {
		try {
				lineSeparator =  System.getProperty("line.separator");
		//	Logo.main(argsit);
		//	 Logo.frame.dispatchEvent(new WindowEvent(Logo.frame, WindowEvent.WINDOW_CLOSING));
		//	timer = System.currentTimeMillis();
		
				maindir = new File(Launcher.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParent().replace("%20", " ");	
				if(new File(maindir +"/demo").exists()) {
					exampledata = true;
				}
				
				BufferedReader fileReader = null;
				//  if(new File(maindir +"/config.txt").exists()) {
				String userHome = System.getProperty("user.home");
				try {
					if(!new File(userHome +"/config_baseplayer.txt").exists()) {
						File conf = new File(userHome +"/config_baseplayer.txt");
						createConf(conf);
						configfile = conf;
					}
					else {
						configfile = new File(userHome +"/config_baseplayer.txt");
					}
				}
				catch(Exception e) {
					try {
						if(!new File(maindir +"/config_baseplayer.txt").exists()) {
							File conf = new File(maindir +"/config_baseplayer.txt");
							createConf(conf);
							configfile = conf;
						}
						else {
							configfile = new File(maindir +"/config_baseplayer.txt");
						}
					}
					catch(Exception ex) {
						if(!fromMain) {
							  ProcessBuilder pb;
								 try {
									 pb = new ProcessBuilder("java","-Xmx1200m", "-Dprism.verbose=true", "-Dsun.java2d.d3d=false","-jar", maindir +"/BasePlayer.jar");				
								     pb.start();						
								 }
								 catch(Exception exc) {
									 e.printStackTrace();									
								 }						     
						  }
						else {
							Main.showError("Could not read/create config-file.\nPlease, add writing permissions to your BasePlayer folder.", "Error");
							return;
						}
					}
				}
			
			  fileReader = new BufferedReader(new FileReader(configfile));			 
			  
			  while((line = fileReader.readLine()) != null) {
				  config.add(line);
				  if(line.startsWith("Memory")) {
					//  if(!line.contains(memlimit)) {
						  memlimit = line.substring(line.indexOf("=")+1).replace(" ", "");
					 // }
				  }
				  else if(line.startsWith("DefaultDir")) {
				//	  if(!line.contains(defaultDir)) {
						  defaultDir = line.substring(line.indexOf("=")+1).trim();
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
				  else if(line.startsWith("backColor")) {
					  backColor = line.substring(line.indexOf("=")+1).replace(" ", "");
				  }
				  else if(line.startsWith("wallpaper")) {
					  try {
						
						  wallpaperIndex = Integer.parseInt(line.substring(line.indexOf("=")+1).replace(" ", ""));
						  
					  }
					  catch(Exception e) {
						  
					  }
				  }
				  else if(line.startsWith("alphaValue")) {
					  try {
						  alphaValue = Integer.parseInt(line.substring(line.indexOf("=")+1).replace(" ", ""));
					  }
					  catch(Exception e) {
						  
					  }
				  }
				  else if(line.startsWith("DefaultSave")) {
					  defaultSaveDir = line.substring(line.indexOf("=")+1).replace(" ", "");
				  }
				  else if(line.startsWith("genomeDir")) {
					  genomeDir = line.substring(line.indexOf("=")+1).replace(" ", "");
				  }
				  else if(line.startsWith("fontSize")) {					  
					  fontSize = line.substring(line.indexOf("=")+1).replace(" ", "");
				  }
				  else if(line.startsWith("proxyHost")) {					  
					  proxyHost = line.substring(line.indexOf("=")+1).replace(" ", "");
				  }
				  else if(line.startsWith("proxyPort")) {					  
					  proxyPort = line.substring(line.indexOf("=")+1).replace(" ", "");
				  }
				  else if(line.startsWith("proxyType")) {					  
					  proxyType = line.substring(line.indexOf("=")+1).replace(" ", "");
				  }
				  else if(line.startsWith("isProxy")) {		
					  if(line.contains("rue")) {
						  isProxy = true;
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
		//  }
		/*  else {
			  if(!fromMain) {
				  ProcessBuilder pb;
					 try {
						 pb = new ProcessBuilder("java","-Xmx1200m", "-Dprism.verbose=true", "-Dsun.java2d.d3d=false","-jar", maindir +"/BasePlayer.jar");				
					     pb.start();						
					 }
					 catch(Exception e) {
						 e.printStackTrace();
						
					 }			      
			     
			  }
		  }		*/  
			
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
	public static void createConf(File conf) {
		try {		
			String[] args = {}; 
			firstStart = true;
			BufferedWriter writer = new BufferedWriter(new FileWriter(conf));
			writer.write("#Config file for BasePlayer" +lineSeparator);
			writer.write("MemoryLimit=1200m"+lineSeparator);
			
			if(exampledata) {
				
				writer.write("DefaultDir=demo/samples"+lineSeparator);
				writer.write("DefaultControlDir=demo/controls"+lineSeparator);
				writer.write("DefaultTrackDir=demo/tracks"+lineSeparator);
			}
			else {
				writer.write("DefaultDir="+lineSeparator);
				writer.write("DefaultControlDir="+lineSeparator);
				writer.write("DefaultTrackDir="+lineSeparator);
			}			
			
			writer.write("DefaultGenome="+lineSeparator);
			writer.write("DefaultGenes="+lineSeparator);			
			writer.write("DefaultProjectDir="+lineSeparator);
			writer.write("DownloadDir="+lineSeparator);
			writer.write("FirstStart=true"+lineSeparator);
			writer.close();
		}
		catch(Exception e) {
			e.printStackTrace();
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
