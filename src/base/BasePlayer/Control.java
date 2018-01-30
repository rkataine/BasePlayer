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
import htsjdk.tribble.readers.TabixReader;

import javax.swing.*;

import java.io.*;
import java.util.ArrayList;
import java.util.Map.Entry;
import java.util.zip.GZIPInputStream;
public class Control {
	
    static String path = Launcher.ctrldir;
	static String offsetString ="";
	static boolean hold = false;
	static int row = 1;		
	static int runIndex = -1;
	static boolean isCancel;	
	private static TabixReader tabixreader;	
	//private 
	static ControlData controlData = new ControlData();
	private static int index;
	private static int refindex;
	private static int allelenumber;
	private static String[] infosplit;
	private static String[] coverages;	
	
public Control() {	
	    
}

static void applyControl() {
	try {
	//	String[] controlvars = null;
		String teststring = null;
		
		for(int c = 0; c<controlData.fileArray.size(); c++) {	
			if(!Main.drawCanvas.loading) {
	    		 break;
	    	 }
			if(!controlData.fileArray.get(c).controlOn || controlData.fileArray.get(c).controlled) {
				continue;
			}
		
			 Main.drawCanvas.loadbarAll = (int)((c/(double)(controlData.fileArray.size()))*100);
			 try {
				 tabixreader = new TabixReader(controlData.fileArray.get(c).getTabixFile());
				 
			 }
			 catch(Exception tab) {
				 Main.showError("Could not find " +controlData.fileArray.get(c).getTabixFile(), "Note");
				 //JOptionPane.showMessageDialog(Main.chromDraw, "Could not find " +controlData.fileArray.get(c).getTabixFile(), "Note", JOptionPane.INFORMATION_MESSAGE);
			 }
			 Main.drawCanvas.loadingtext = "Applying control " +controlData.fileArray.get(c).getName();				
			
			 TabixReader.Iterator iterator=null; 		
		    
			 VarNode current = FileRead.head.getNext();	     
			if (current == null) {
				return;
			}
		     try {
		    	
		    	 
		    	iterator = tabixreader.query(Main.chromosomeDropdown.getSelectedItem().toString() +":" +(current.getPosition()-1));
		    	teststring = iterator.next();
		    	
		    	if(teststring == null) {	
		    		
	    			iterator = tabixreader.query("chr" +Main.chromosomeDropdown.getSelectedItem().toString() +":" +(current.getPosition()-1));
	    			teststring = iterator.next();
	    			 
		    		if(teststring == null) {
		    			 try {
					    	 if(Main.chromosomeDropdown.getSelectedItem().toString().equals("X")) {
					    		
					    		 iterator = tabixreader.query("23:" +(current.getPosition()-1));
					    		 
					    	 }
					    	 else if(Main.chromosomeDropdown.getSelectedItem().toString().equals("Y")) {
					    		 iterator = tabixreader.query("24:" +(current.getPosition()-1));
					    	 }
					    	 else if(Main.chromosomeDropdown.getSelectedItem().toString().contains("M")) {
					    		 iterator = tabixreader.query("25:" +(current.getPosition()-1));
					    	 }
					    	 teststring = iterator.next();
					    		if(teststring == null) {
					    			ErrorLog.addError("Chromosome " +Main.chromosomeDropdown.getSelectedItem().toString() +" not found in " +controlData.fileArray.get(c).getName());
				    		    	 
					    		}
		    			 }
		    			 catch(Exception exc) {
		    				 ErrorLog.addError("Chromosome " +Main.chromosomeDropdown.getSelectedItem().toString() +" not found in " +controlData.fileArray.get(c).getName());
		    		    	 
		    			 }
		    		}
		    	}
		     }
		     catch(Exception e) {	    
		    	 
		    	 ErrorLog.addError(e.getStackTrace());		    		
		     }
		     
		     if (iterator == null || teststring == null) {
		    	continue; 
		     }
		    /* if (1==1) {
		    	 String line = "";
		    	 String[] split;
		    	 int position;
		    	 while(line !=null) {	 
		    		 line = iterator.next();
			    	 
			    	 if(line == null) {
			    		 break;
			    	 }
			    	 split = line.split("\t");	    	
			    	 position = Integer.parseInt(split[1])-1;
			    	 System.out.println(position);		    		
		    	 }
		    	 return;
		     }
		     */
		     if(controlData.fileArray.get(c).getTabixFile().endsWith(".vcf.gz")) {	
		    	 if(!controlData.fileArray.get(c).remOverlaps.isSelected()) {
		    		 useVCFstrict(iterator, current, c, teststring);				  
		    	 }
		    	 else {
		    		 useVCFoverlap(iterator, current, c, teststring);
		    	 }
			 }
		     else {
		    	 useCTRL(iterator, current, c);		    	 
		     }
		     current = null;
		     controlData.fileArray.get(c).controlled = true;
			 
		}
		
		/*if(VariantHandler.commonSlider.getValue() > 1 && VariantHandler.clusterSize > 0) {
	    	Main.drawCanvas.calcClusters(FileRead.head);
		}*/
		Draw.calculateVars = true;
		Draw.updatevars = true;
		Main.drawCanvas.repaint();
		VariantHandler.table.repaint();
	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
	}
	 
	
}

static void useCTRL(TabixReader.Iterator iterator, VarNode current, int c) {
	try {
	String line, controlvar;
	String[] split;
	int position, count = 0;
	ArrayList<SampleNode> samplelist = null;
	
	   while((line = iterator.next()) !=null) {	 
	    	
	    	 split = line.split("\t");
	    	
	    	 
	    	 if(split[2].length() > 1 && !split[2].contains(",")) {
	    		 
	    		 position = Integer.parseInt(split[1])-1;
	    		 
	    		 if(split[2].startsWith("-")) {
	    			 if(split[2].substring(1).length() == 2) {
	    				 controlvar = "del" +split[2].substring(2);
	    			 }
	    			 else {
	    				 controlvar = "del" +(split[2].length()-2);
	    			 }
	    		 }
	    		 else {
	    			
	    			 if(split[2].substring(1).length() == 1) {
	    				 controlvar = "ins" +split[2].substring(1);
	    			 }
	    			 else {
	    				 controlvar = "ins" +(split[2].length()-1);
	    			 }		    			 
	    		 }		    			
	    	 }
	    	 else {
	    		 position = (Integer.parseInt(split[1])-1);
	    		 controlvar = split[2];
	    	 }		    	 
	    	  
	    	 while(current != null && current.getPosition() < position) {
	    		
	    		 current = current.getNext();		    		 
	    	 }
	    	 if(current == null) {
					break;
			 }
	    	 if(position < current.getPosition()) {
	    		 continue;
	    	 }
	    	 if(current.getPosition() == position) { 
	    		 current.controlled = true;
	    		
	    		 for(int i = 0; i<current.vars.size(); i++) {
	    			 if(split[2].contains(",")) {					    				
	    				samplelist = current.vars.get(i).getValue();		    				
	    				 samplelist.add(new SampleNode(Integer.parseInt(split[3])*(Integer.parseInt(split[4])+1),controlData.fileArray.get(c).varcount,controlData.fileArray.get(c)));		    				
	    				 continue;
	    			 }
	    			 else if(current.vars.get(i).getKey().equals(controlvar)) {		    			
	    				 samplelist = current.vars.get(i).getValue();		
	    				 if(samplelist.get(samplelist.size()-1).getControlSample() != null && samplelist.get(samplelist.size()-1).getControlSample().equals(controlData.fileArray.get(c)) ) {
	    					 samplelist.get(samplelist.size()-1).alleles += Integer.parseInt(split[3])*(Integer.parseInt(split[4])+1);
	    				 }
	    				 else {
	    					 samplelist.add(new SampleNode(Integer.parseInt(split[3])*(Integer.parseInt(split[4])+1),controlData.fileArray.get(c).varcount,controlData.fileArray.get(c)));		    				
	    		    			
	    				 }
	    					 continue;
	    			 }
	    		 }	    		 
	    	//	 current = current.getNext();
	    		 
	    	 }
	    	 
			if(count % 1000 == 0) {
				 Main.drawCanvas.loadBarSample = (int)((current.getPosition()/(double)(Main.drawCanvas.splits.get(0).chromEnd))*100);     
				 Draw.updatevars = true;
				 Main.drawCanvas.repaint();
			}
			/*if(current == null) {
				break;
			}*/
					count++;
		}
	   if(samplelist != null) {
			 samplelist = null;
		}
	   current = null;
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}

static void useVCFoverlap(TabixReader.Iterator iterator, VarNode current, int c, String firstline) {
	try {
		String line = firstline, controlvar = "";
		String[] split, templist = null;
		int position, count = 0, ref= 0, alt = 0, endindex=-1, acendindex =-1;
		ArrayList<SampleNode> samplelist = null;
		int alleles = 0, baselength = 0;		 
		
		   while(line !=null) {	    	 
		    	 if(!Main.drawCanvas.loading) {
		    		 break;
		    	 }
		    	 split = line.split("\t");	    	
		    	 baselength = MethodLibrary.getControlBaseLength(split[3], split[4], 0);
		    	
		    	 position = Integer.parseInt(split[1])-1;
		    	 while(current != null && current.getPosition() < position) {	    		
		    		 current = current.getNext();		    		 
		    	 }
		    	 if(current == null) {
					break;
				 }
		    	 if(position+baselength < current.getPosition()) {
		    		 line = iterator.next();
		    		 continue;
		    	 }	    	
		    	
		    	 if(current.getPosition() >= position && current.getPosition() <= position+baselength) { 
		    		 if(!split[4].contains(",")) {
			    		 controlvar = FileRead.getVariant(split[3], split[4]);		    	
			    	 }
			    	 else {
			    		 templist = split[4].split(",");
			    		
			    		 for(int i = 0; i<templist.length; i++) {
			    			 templist[i] = FileRead.getVariant(split[3], templist[i]);	
			    			
			    		 }	    		
			    	 }
		    		 current.controlled = true;
		    		 index = split[7].indexOf(";AC=")+4;
		    		 if(index < 0) {
		    			 index = split[7].indexOf("AC=")+3;
		    		 }
		    		 acendindex = split[7].indexOf(";", index);
		    		 if(acendindex == -1) {
		    			 acendindex = split[7].length();
		    		 }
		    		 refindex = split[7].indexOf(";AN=")+4;
		    		
		    		 endindex = split[7].indexOf(";", refindex);
		    		 if(endindex == -1) {
		    			 endindex = split[7].length();
		    		 }
		    		 
		    		 for(int i = 0; i<current.vars.size(); i++) {
		    			 count++;
		    			 if(split[4].contains(",")) {
		    				 
		    				 for(int t = 0; t<templist.length; t++) {
		    					 if(current.vars.get(i).getKey().equals(templist[t]) || baselength > 1) {
		    						
		    						 if(controlData.fileArray.get(c).varcount > 2) {
		    							 alleles = Integer.parseInt(split[7].substring(index, acendindex).split(",")[t]);	    							 					    		
		 					    		 allelenumber = Integer.parseInt(split[7].substring(refindex, endindex));
		 					    	
		    						 }
		    						 else {
		    							 infosplit = split[split.length-1].split(":");		    							 
		    							 coverages = infosplit[split[8].indexOf("AD")/3].split(",");
		    							 ref = Integer.parseInt(coverages[0]);
		    							 alt = Integer.parseInt(coverages[1]);
		    							 
		    							 alleles = alt;				 	
		    							 allelenumber = ref+alt;
		    						 }
		    						 if(alleles > 0) {
		    							 
			    	    				 samplelist = current.vars.get(i).getValue();		    				
			    	    				 samplelist.add(new SampleNode(alleles,allelenumber,controlData.fileArray.get(c)));	
		    						 }
		    	    				
		    					 }	    					
		    	    		 }
		    				 
		    				 continue;
		    			 }
		    			 else if(current.vars.get(i).getKey().equals(controlvar) || baselength > 1) {
		    				 
	    					 if(controlData.fileArray.get(c).varcount > 2) {
	    						 allelenumber = Integer.parseInt(split[7].substring(refindex, endindex));
	    						 alleles = Integer.parseInt(split[7].substring(index, acendindex));  					 
	    					 
	    					 }
	    					 else {
	    						 infosplit = split[split.length-1].split(":");
								 if(infosplit[0].length() > 2 ) { 
	    							 if(infosplit[0].charAt(0) != infosplit[0].charAt(2)) {
	    								 alleles = 1;
	    							 }
	    							 else {
	    								 alleles = 2;
	    							 }
							 	 }
								 else {
									 alleles = 0;
								 }
								 allelenumber = 2;
								/* coverages = infosplit[split[8].indexOf("AD")/3].split(",");
								 ref = Integer.parseInt(coverages[0]);
								 alt = Integer.parseInt(coverages[1]);
								 
								 alleles = alt;				 	
								 allelenumber = ref+alt;
	    					 */
	    					 }
		    				
		    				 if(alleles > 0) {
		    				    samplelist = current.vars.get(i).getValue();		   
		    				 	samplelist.add(new SampleNode(alleles,allelenumber,controlData.fileArray.get(c)));	    				
		    				 }
		    				 continue;
		    			 }
		    		 }	    		 
		    	 }
		    	 
				if(count % 10 == 0) {
					 Main.drawCanvas.loadBarSample = (int)((current.getPosition()/(double)(Main.drawCanvas.splits.get(0).chromEnd))*100);     
					 Draw.updatevars = true;
					 Main.drawCanvas.repaint();
				}
					try {
						line = iterator.next();
					}
					catch(htsjdk.samtools.FileTruncatedException e) {
						
					}
			}
		   if(samplelist != null) {
				 samplelist = null;
			}
		   current = null;
		}
		catch(Exception e) {
			
			e.printStackTrace();
		}
}

static void useVCFstrict(TabixReader.Iterator iterator, VarNode current, int c, String firstline) {
	try {
	String line = firstline, controlvar = "";
	String[] split, templist = null;
	int position, count = 0, endindex = -1, acendindex = -1;
	ArrayList<SampleNode> samplelist = null;
	int alleles = 0;
	 
	
	   while(line !=null) {	    	 
	    	 if(!Main.drawCanvas.loading) {
	    		 break;
	    	 }
	    	 split = line.split("\t");	    	
	    	 position = Integer.parseInt(split[1])-1;
	    	 while(current != null && current.getPosition() < position) {	    		
	    		 current = current.getNext();		    		 
	    	 }
	    	 if(current == null) {
				break;
			 }
	    	 if(position < current.getPosition()) {
	    		 line = iterator.next();
	    		 continue;
	    	 }	    	
	    	
	    	 if(current.getPosition() == position) { 
	    		 if(!split[4].contains(",")) {
		    		 controlvar = FileRead.getVariant(split[3], split[4]);		    	
		    	 }
		    	 else {
		    		 templist = split[4].split(",");
		    		
		    		 for(int i = 0; i<templist.length; i++) {
		    			 templist[i] = FileRead.getVariant(split[3], templist[i]);		 
		    			 
		    		 }	    		
		    	 }
	    		 current.controlled = true;
	    		
	    		 index = split[7].indexOf(";AC=")+4;
	    		 if(index < 0) {
	    			 index = split[7].indexOf("AC=")+3;
	    		 }
	    		 acendindex = split[7].indexOf(";", index);
	    		
	    		 if(acendindex == -1) {
	    			 acendindex = split[7].length();
	    		 }
	    		 refindex = split[7].indexOf(";AN=")+4;
	    		 
	    		/* if(refindex == 3) {
	    			 refindex = split[7].indexOf("AN=")+3;
	    			 
	    		 }*/
	    		 endindex = split[7].indexOf(";", refindex);
	    		 if(endindex == -1) {
	    			 endindex = split[7].length();
	    		 }
	    		 for(int i = 0; i<current.vars.size(); i++) {
	    			 count++;
	    			 if(split[4].contains(",")) {
	    				 
	    				 for(int t = 0; t<templist.length; t++) {
	    					 
	    					 if(current.vars.get(i).getKey().equals(templist[t])) {
	    						 if(controlData.fileArray.get(c).varcount > 2) {
	    							 alleles = Integer.parseInt(split[7].substring(index, acendindex).split(",")[t]);
	    							
	    						//	 System.out.println(position +" " +alleles); 
	    							try {
	 					    		 allelenumber = Integer.parseInt(split[7].substring(refindex, endindex));
	    							}
	    							catch(Exception e) {
	    								System.out.println(line);
	    							}
	    					//		 allelenumber = Integer.parseInt(split[7].substring(refindex, split[7].indexOf(";", refindex)));
	    							
	 					    		
	    						 }
	    						 else {
	    						//	 infosplit = split[split.length-1].split(":");
	    							 if(infosplit[0].length() > 2 ) { 
		    							 if(infosplit[0].charAt(0) != infosplit[0].charAt(2)) {
		    								 alleles = 1;
		    							 }
		    							 else {
		    								 alleles = 2;
		    							 }
	    						 	 }
	    							 else {
	    								 System.out.println(line);
	    								 alleles = 0;
	    							 }
	    						/*	 coverages = infosplit[split[8].indexOf("AD")/3].split(",");
	    							
	    							 ref = Integer.parseInt(coverages[0]);
	    							
	    							 alt = Integer.parseInt(coverages[1]);
	    							 */
	    						//	 alleles = alt;				 	
	    						//	 allelenumber = ref+alt;
	    							 allelenumber = 2;
	    						 }
	    						
	    						 if(alleles > 0) {
	    							samplelist = current.vars.get(i).getValue();		    				
	    	    				 	samplelist.add(new SampleNode(alleles,allelenumber,controlData.fileArray.get(c)));	
	    						 }
	    	    				// break;
	    					 }	    					
	    	    		 }
	    				 
	    				 continue;
	    			 }
	    			 else if(current.vars.get(i).getKey().equals(controlvar)) {
	    				 
    					 if(controlData.fileArray.get(c).varcount > 2) {
    						 allelenumber = Integer.parseInt(split[7].substring(refindex, endindex));
    						 alleles = Integer.parseInt(split[7].substring(index, acendindex));    					 
    					 
    					 }
    					 else {
    						 infosplit = split[split.length-1].split(":");
							 if(infosplit[0].length() > 2 ) { 
    							 if(infosplit[0].charAt(0) != infosplit[0].charAt(2)) {
    								 alleles = 1;
    							 }
    							 else {
    								 alleles = 2;
    							 }
						 	 }
							 else {
								 alleles = 0;
							 }
							 allelenumber = 2;
							/* coverages = infosplit[split[8].indexOf("AD")/3].split(",");
							 ref = Integer.parseInt(coverages[0]);
							 alt = Integer.parseInt(coverages[1]);
							 
							 alleles = alt;				 	
							 allelenumber = ref+alt;
    					 */
    					 }
	    				
	    				 if(alleles > 0) {
	    					 samplelist = current.vars.get(i).getValue();	    				
	    				 	 samplelist.add(new SampleNode(alleles,allelenumber,controlData.fileArray.get(c)));	    
	    				 }
	    				 continue;
	    			 }
	    		 }	    		 
	    	 }
	    	 
			if(count % 10 == 0) {
				 Main.drawCanvas.loadBarSample = (int)((current.getPosition()/(double)(Main.drawCanvas.splits.get(0).chromEnd))*100);     
				 Draw.updatevars = true;
				 Main.drawCanvas.repaint();
			}
				try {
				line = iterator.next();
				}
				catch(htsjdk.samtools.FileTruncatedException e) {
					
				}
		}
	   if(samplelist != null) {
			 samplelist = null;
		}
	   current = null;
	}
	catch(Exception e) {
		
		e.printStackTrace();
	}
}
public static class runner extends SwingWorker<String, Object> {	
		
	protected String doInBackground() {
		
		Main.drawCanvas.loading("Controlling");
		applyControl();
		Main.drawCanvas.ready("Controlling");
		return "";
	}
	
}

static void addFiles(File[] filestemp) {
//  	int calc = 0;
  	try {
  	
  	int samplecount = 0;
  	if(Main.trackPane.isVisible() && Control.controlData.fileArray.size() == 0) {
  		Main.trackPane.setDividerLocation(Main.varpane.getDividerLocation());
  	}
	
	for(int i=0;i<filestemp.length; i++) {
		if(filestemp[i].isDirectory()) {
    		for(int j=0;j< filestemp[i].listFiles().length; j++) {
    			if(filestemp[i].listFiles()[j].length() == 0) {
    				continue;
    			}
    			if(filestemp[i].listFiles()[j].getAbsolutePath().matches(".*.vcf") || filestemp[i].listFiles()[j].getAbsolutePath().matches(".*.ctrl")) {
    				ControlFile addSample = new ControlFile(filestemp[i].listFiles()[j].getName(), (short)j, filestemp[i].listFiles()[j].getCanonicalPath());
    				controlData.fileArray.add(addSample);
    			
    				
    				if(filestemp[i].listFiles()[j].getName().contains(".ctrl")) {
    					Object[] addi = {row, filestemp[i].listFiles()[j].getName(), new Boolean(true)};
    					 Main.controlDraw.trackDivider.add(0.0);
    					
    					try {
    						samplecount = 0;
    						BufferedReader reader = new BufferedReader(new FileReader(filestemp[i].listFiles()[j]));
    						String line = reader.readLine();
    						while(line.startsWith("#")) {
    							if(line.startsWith("##Chr")) {
    								break;
    							}
    							if(line.contains("#SampleCount")) {
    								samplecount = Integer.parseInt(line.split("=")[1]);
    								controlData.total+= samplecount;
    							}
    							else {
    								samplecount++;
    								controlData.total++;
    							}
 //   							controlData.sampleArray.add(line.substring(1));
    				//			Object[] add = {row, line.substring(1), ""};
    							
    						
    							row++;
    							line = reader.readLine();
    						}
    						reader.close();					
    				//		controlsamples.put(addi[1].toString(), samplecount);
    						addi[1] = ""+samplecount +"_"+addi[1];
    	    				
    					}
    					catch(Exception e){
    						ErrorLog.addError(e.getStackTrace());
    						e.printStackTrace();
    					}
    				}
    				else {   	
    					row++; 
    					samplecount++;
    					controlData.total++;    				
    				}    				
    			}    			
    		}   
		}
		else {	
			ControlFile addSample = new ControlFile(filestemp[i].getName(), (short)(controlData.fileArray.size()), filestemp[i].getCanonicalPath());			
			controlData.fileArray.add(addSample);		
			MethodLibrary.addHeaderColumns(addSample);
			Main.controlDraw.trackDivider.add(0.0);
		/*	VariantHandler.table.addRowGeneheader(addSample);	
			VariantHandler.table.addRowGeneheader("OR");
			 if(VariantHandler.tables.size() > 0) {
				  for(int t = 0; t < VariantHandler.tables.size(); t++) {
					  VariantHandler.tables.get(t).addRowGeneheader(Control.controlData.fileArray.get(i));
					  VariantHandler.tables.get(t).addRowGeneheader("OR");
					  VariantHandler.tables.get(t).repaint();
				  }
			  }
			  VariantHandler.clusterTable.addColumnGeneheader(addSample);
			  VariantHandler.clusterTable.addColumnGeneheader("OR");
			  VariantHandler.clusterTable.repaint();
			  VariantHandler.table.repaint();
			*/  
			if(filestemp[i].getName().endsWith(".ctrl")) {
				
				Object[] addi = { "", filestemp[i].getName(), new Boolean(true) };
				try {
					samplecount = 0;
					BufferedReader reader = new BufferedReader(new FileReader(filestemp[i]));
					String line = reader.readLine();
					while(line.startsWith("#")) {
						if(line.startsWith("##Chr")) {
							break;
						}
						if(line.contains("#SampleCount")) {
							samplecount = Integer.parseInt(line.split("=")[1]);
							controlData.total += samplecount*2;
						}
						else {
							samplecount++;
							controlData.total+=2;
						}
//						controlData.sampleArray.add(line);
			//			Object[] add = {row, line, new Boolean(true)};
						
						row++;
						line = reader.readLine();
					}
					reader.close();					
		//			controlsamples.put(addi[1].toString(), samplecount);
					addi[1] = ""+samplecount +"_"+addi[1];
    				
				}
				catch(Exception e){
					ErrorLog.addError(e.getStackTrace());
					e.printStackTrace();
				}
			}
			else if(filestemp[i].getName().endsWith(".ctrl.gz")) {
		//		Object[] addi = { "", filestemp[i].getName(), new Boolean(true) };
		//		model.addRow(addi);
			
				try {
					samplecount = 0;
					GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(filestemp[i]));
					BufferedReader reader = new BufferedReader(new InputStreamReader(gzip));
				
					String line = reader.readLine();
					while(line.startsWith("#")) {
						if(line.startsWith("##")) {
							break;
						}
						if(line.startsWith("#SampleCount")) {
							samplecount = Integer.parseInt(line.split("=")[1]);
							controlData.total += samplecount*2;
						}
						else {
							samplecount++;
							controlData.total+=2;
						}						
						
						line = reader.readLine();
					}
					addSample.varcount = samplecount*2;
					reader.close();					
		//			controlsamples.put(addi[1].toString(), samplecount);
					controlData.sampleCount +=samplecount;
				//	addi[1] = ""+samplecount +"_"+addi[1];
    		//		model.fireTableDataChanged();
    				
    				
				}
				catch(Exception e){
					ErrorLog.addError(e.getStackTrace());
					e.printStackTrace();
				}
				
				
			}
			else if(filestemp[i].getName().endsWith(".vcf.gz")) {
			
				try {
					samplecount = 0;
					GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(filestemp[i]));
					BufferedReader reader = new BufferedReader(new InputStreamReader(gzip));				
					String line = reader.readLine(), population = "";
					int idindex = 0, dotindex = 0;
					while(line.startsWith("#")) {
						if(line.contains("##INFO")) {
							if(line.contains("ID=")) {
								idindex = line.indexOf("ID=");
								dotindex = line.indexOf(",");
								if(idindex > 0 && dotindex > 0) {
									if(line.substring(idindex, dotindex).contains("AC")) {
										population = line.substring(idindex+3, dotindex).replace("AC", "").replace("_", "");
										if(population.length() == 0) {
											population = "ALL";
										}
										
										
									}
								}								
							}
						}
						else if(line.startsWith("#CHROM")) {
							if(line.split("\t").length -9 < 0) {
								samplecount = 0;
							}
							else {
								samplecount = line.split("\t").length -9;
							}
							controlData.total += samplecount*2;
						}					
						
						line = reader.readLine();
					}
					addSample.varcount = samplecount*2;
					if(samplecount == 0) {
						int count =0, maxnumber = 0, endindex = -1;
						String[] split = line.split("\t");
						if(line.contains("AN=")) {
							while(count < 1000) {
								
								refindex = split[7].indexOf(";AN=")+4;
					    		if(refindex == -1) {
					    			refindex = split[7].indexOf("AN=")+3;
					    		}
					    		endindex = split[7].indexOf(";", refindex);
					    		if(endindex == -1) {
					    			endindex = split[7].length();					    			
					    		}				
					    		
					    		allelenumber = Integer.parseInt(split[7].substring(refindex, endindex));
					    		
					    		if(allelenumber > maxnumber) {
					    			maxnumber = allelenumber;
					    		}
								count++;
								line = reader.readLine();
								if(line != null) {
									split = line.split("\t");
								}
								else {
									break;
								}
							}
							
						}
						addSample.varcount = maxnumber;
						samplecount = maxnumber/2;
					}
				
					else if(samplecount == 1) {
						addSample.varcount = 2;
						samplecount = 1;
					}
				
					reader.close();					
					gzip.close();
					controlData.sampleCount +=samplecount;
					
				}
				catch(Exception e){
					ErrorLog.addError(e.getStackTrace());
					e.printStackTrace();
				}				
			}	
			
		}
		 
	}
  
	
 if(controlData.fileArray.size() > 0) {	
	 for(int i = 0 ; i<Main.controlDraw.trackDivider.size(); i++) {
		 Main.controlDraw.trackDivider.set(i, ((i+1)*(Main.varpane.getDividerLocation()/(double)Main.controlDraw.trackDivider.size())));
	 }	
	 if(!Main.trackPane.isVisible()) {
		  Main.trackPane.setVisible(true);	
		  Main.varpane.setDividerSize(2);		  
		  Main.varpane.setDividerLocation(controlData.fileArray.size()*80);		  
	 }
	 else {
		  
		  Main.varpane.setDividerLocation(Main.varpane.getDividerLocation()+80);
		  
		  if(Main.bedScroll.isVisible()) {
			  Main.trackPane.setDividerSize(2);
		  }
	 }
	 
	 Main.controlScroll.setVisible(true);
	 Main.controlDraw.setVisible(true);
	 hold = false;
 }
 
 }  
  	catch(Exception e) {
  		ErrorLog.addError(e.getStackTrace());
  		 e.printStackTrace();
  	 }
}
static void clearControls() {
	Entry<String, ArrayList<SampleNode>> entry;	
	 VarNode current = FileRead.head.getNext();
	 while(current != null) {
		 if(current.controlled) {
			 for(int v = 0; v<current.vars.size(); v++) {
					entry = current.vars.get(v);
					
					for(int i = entry.getValue().size()-1; i> -1;i-- ) {
						if(entry.getValue().get(i).alleles == null) {
							
							break;
						}
						
						entry.getValue().remove(i);
					
					}							
			 }
			 current.controlled = false;
		 }
		
		 current = current.getNext();
	 }
	 controlData.fileArray.clear();
//	 controlData.sampleArray.clear();	
	 controlData.total = 0;
	 controlData.sampleCount = 0;
	 entry = null;	
	 current = null;
}

static void dismissControls(VarNode head) {
	VarNode current = head;
	Entry<String, ArrayList<SampleNode>> entry;
	while(current != null) {
		 if(current.controlled) {
			 for(int v = 0; v<current.vars.size(); v++) {
					entry = current.vars.get(v);
					
					for(int i = entry.getValue().size()-1; i> -1;i-- ) {
						if(entry.getValue().get(i).alleles == null) {
							break;
						}
						entry.getValue().remove(i);
						
					}							
			 }
			 current.controlled = false;
		 }
		
		 current = current.getNext();
	 }
	for(int i = 0 ; i<controlData.fileArray.size(); i++) {
		controlData.fileArray.get(i).controlled = false;
	}
	entry = null;
	current = null;
	head = null;
}
static void dismissControl(VarNode head, ControlFile sample) {
	VarNode current = head;
	Entry<String, ArrayList<SampleNode>> entry;
	sample.controlOn = false;
	Main.controlDraw.repaint();
	sample.controlled = false;
	while(current != null) {
		 if(current.controlled) {
			 for(int v = 0; v<current.vars.size(); v++) {
					entry = current.vars.get(v);
					
					for(int i = entry.getValue().size()-1; i> -1;i-- ) {
						if(entry.getValue().get(i).alleles == null) {
							break;
						}
						if(entry.getValue().get(i).getControlSample().equals(sample)) {
							entry.getValue().remove(i);
						}
					}							
			 }
			
		 }
		
		 current = current.getNext();
	 }
	for(int i = 0 ; i<controlData.fileArray.size(); i++) {
		controlData.fileArray.get(i).controlled = false;
	}
	
	Control.applyControl();
	
	entry = null;
	current = null;
	head = null;
}
static void control(VarNode head) {
	try {
	//	Control.controlData.alleleThreshold = Double.parseDouble(allelebox.getText());
	}
	catch(Exception ec) {
		ec.printStackTrace();
		ErrorLog.addError(ec.getStackTrace());
	}
	try {
		
		runner runner = new runner();
	    runner.execute();
	     
	}
	catch(Exception ex) {
		ex.printStackTrace();
		ErrorLog.addError(ex.getStackTrace());
	}
}


static class MyFilter extends javax.swing.filechooser.FileFilter {
	public boolean accept(File file) { 
			if (file.isDirectory()) {
				return true;
		    }			
			
			else if(file.getName().matches(".*.ctrl.gz")) {
				return true;
			}
		/*	else if(file.getName().matches(".*.vcf")) {
				return true;
			}*/
			else if(file.getName().matches(".*.vcf.gz")) {
				return true;
			}
			else {
				return false;
			}		
	}
	
	public String getDescription() { return "*.ctrl.gz, *.vcf.gz"; }
	}


	

}


