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
	private static final long serialVersionUID = 1L;
	
   
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
				 JOptionPane.showMessageDialog(Main.chromDraw, "Could not find " +controlData.fileArray.get(c).getTabixFile(), "Note", JOptionPane.INFORMATION_MESSAGE);
			 }
			 Main.drawCanvas.loadingtext = "Applying control " +controlData.fileArray.get(c).getName();				
			
			 TabixReader.Iterator iterator=null; 		
		    
			 VarNode current = FileRead.head.getNext();	     
		     
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
		     if(controlData.fileArray.get(c).getTabixFile().endsWith(".vcf.gz")) {
		    	
		    	useVCF(iterator, current, c, teststring);		    	 
				  
			 }
		     else {
		    	 useCTRL(iterator, current, c);
		    	 
		    	 
		     }
		     current = null;
		     controlData.fileArray.get(c).controlled = true;
			 
		}
		if(VariantHandler.commonSlider.getValue() > 1 && !VariantHandler.clusterBox.getText().equals("0") && !VariantHandler.clusterBox.getText().equals("")) {
	    	Main.drawCanvas.calcClusters(FileRead.head);
		}
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
static void useVCF(TabixReader.Iterator iterator, VarNode current, int c, String firstline) {
	try {
	String line = firstline, controlvar = "";
	String[] split, templist = null;
	int position, count = 0;
	ArrayList<SampleNode> samplelist = null;
	int alleles = 0;
	
	   while(line !=null) {	 
	    	
	    	 split = line.split("\t");	    	
	    	 position = Integer.parseInt(split[1])-1;
	    	 if(!Main.drawCanvas.loading) {
	    		 break;
	    	 }
	    	 
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
	    			 index = split[7].indexOf("\tAC=")+4;
	    		 }
	    		 refindex = split[7].indexOf(";AN=")+4;
	    		 
	    		 for(int i = 0; i<current.vars.size(); i++) {
	    			 count++;
	    			 if(split[4].contains(",")) {
	    				 
	    				 for(int t = 0; t<templist.length; t++) {
	    					 if(current.vars.get(i).getKey().equals(templist[t])) {
	    						 if(controlData.fileArray.get(c).varcount > 2) {
	    						 alleles = Integer.parseInt(split[7].substring(index, split[7].indexOf(";", index)).split(",")[t]);
	    						 allelenumber = Integer.parseInt(split[7].substring(refindex, split[7].indexOf(";", refindex)));
	    						 }
	    						 else {
	    							 alleles = 2;
	    							 allelenumber = 2;
	    						 }
	    	    				 samplelist = current.vars.get(i).getValue();		    				
	    	    				 samplelist.add(new SampleNode(alleles,allelenumber,controlData.fileArray.get(c)));	
	    	    				 break;
	    					 }	    					
	    	    		 }
	    				 
	    				 continue;
	    			 }
	    			 else if(current.vars.get(i).getKey().equals(controlvar)) {
	    				 allelenumber = Integer.parseInt(split[7].substring(refindex, split[7].indexOf(";", refindex)));
	    				 if(controlData.fileArray.get(c).varcount > 2) {
	    					 alleles = Integer.parseInt(split[7].substring(index, split[7].indexOf(";", index)));
	    				 }
	    				 else {
	    					 alleles = 2;
	    				 }
	    				 
	    				 samplelist = current.vars.get(i).getValue();		   
	    				
	    				 samplelist.add(new SampleNode(alleles,allelenumber,controlData.fileArray.get(c)));	    				
	    				 continue;
	    			 }
	    		 }	    		 
	    	 }
	    	 
			if(count % 10 == 0) {
				 Main.drawCanvas.loadBarSample = (int)((current.getPosition()/(double)(Main.drawCanvas.splits.get(0).chromEnd))*100);     
				 Draw.updatevars = true;
				 Main.drawCanvas.repaint();
			}
				
				line = iterator.next();
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
			VariantHandler.table.addRowGeneheader("AF: " +filestemp[i].getName());	
			VariantHandler.table.addRowGeneheader("OR");
			 if(VariantHandler.tables.size() > 0) {
				  for(int t = 0; t < VariantHandler.tables.size(); t++) {
					  VariantHandler.tables.get(t).addRowGeneheader("AF: " +Control.controlData.fileArray.get(i).getName());
					  VariantHandler.tables.get(t).addRowGeneheader("OR");
					  VariantHandler.tables.get(t).repaint();
				  }
			  }
			  VariantHandler.clusterTable.addRowGeneheader("AF: " +Control.controlData.fileArray.get(i).getName());
			  VariantHandler.clusterTable.addRowGeneheader("OR");
			  VariantHandler.clusterTable.repaint();
			  VariantHandler.table.repaint();
			  
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
					String line = reader.readLine();
					
					while(line.startsWith("#")) {
					
						if(line.startsWith("#CHROM")) {
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
						int count =0, maxnumber = 0;
						String[] split = line.split("\t");
						if(line.contains("AN=")) {
							while(count < 1000) {
								
								refindex = split[7].indexOf(";AN=")+4;
					    		 
					    		
					    		allelenumber = Integer.parseInt(split[7].substring(refindex, split[7].indexOf(";", refindex)));
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
	
	 if(!Main.trackPane.isVisible()) {
		  Main.trackPane.setVisible(true);
	//	  Main.trackPane.setDividerSize(3);
		  Main.varpane.setDividerSize(4);	  
		//  Main.varpane.setResizeWeight(0.1);
		  
		  Main.varpane.setDividerLocation(controlData.fileArray.size()*64);
		  
	 }
	 else {
		  Main.trackPane.setDividerLocation(0.5);		 
		  Main.varpane.setDividerLocation(Main.varpane.getDividerLocation()+100);
		  if(Main.bedScroll.isVisible()) {
			  Main.trackPane.setDividerSize(4);
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
void clearControls() {
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
	entry = null;
	current = null;
	head = null;
}
static void control(VarNode head) {
	try {
	//	Control.controlData.alleleThreshold = Double.parseDouble(allelebox.getText());
	}
	catch(Exception ec) {
	//	allelebox.setText("0.0");
	//	Control.controlData.alleleThreshold = 0.0;
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


