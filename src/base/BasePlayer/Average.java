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

//import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
//import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.tribble.readers.TabixReader;

import java.awt.Component;
import java.awt.FileDialog;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.*;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;






import javax.swing.table.TableCellRenderer;

import java.awt.Color;

import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.SwingWorker;
import javax.swing.table.DefaultTableModel;


@SuppressWarnings("deprecation")
public class Average extends JPanel implements TableCellRenderer, MouseListener  {
	private static final long serialVersionUID = 1L;
	
    static JFrame frame;    	
	//UI	
	//Labels	
    static java.util.List<int[]> mergeVector = null;
	//Buttons etc.
    static JButton openbutton = new JButton("Open BED-file");
	static JButton output = new JButton("Write output");
	static int startpos=0,endpos=0,startchrom=0, endchrom=0; 
	static JButton execute = new JButton("Execute");
	JButton cancel = new JButton("Cancel");
	static JTextArea info = new JTextArea();
	static java.util.List<java.util.List<String>> outVector =Collections.synchronizedList(new ArrayList<java.util.List<String>>());
	
	JLabel fileLabel = new JLabel("Open BED-file to calculate average coverage");
	static JPanel panel = new JPanel(new GridBagLayout());
	static Object[] headers = {"Sample", "Average coverage", "Average mapping quality", "Soft clip rate", "Zero quality rate", "Covered (%) (Coverage : Percentage)", "Status"}; 
    static Object[][] data = {}; 
    static DefaultTableModel model = new DefaultTableModel(data, headers); 
    static JCheckBox smallregion = new JCheckBox("Small region?");
	static JTable table = new JTable(model);
	JScrollPane infoScroll = new JScrollPane(table,  
    JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,  
    JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
	
	public static void main(String[] args) {				 
		 javax.swing.SwingUtilities.invokeLater(new Runnable() {
	            public void run() {	            
	            		createAndShowGUI();	            	
	            }
	        });			
	}
	public static void setSamples() {
		outVector.clear();
		model.setRowCount(0);
		if(outVector.size() == 0) {
			for(int i = 0; i<Main.samples; i++) {
				if(Main.drawCanvas.sampleList.get(i).samFile != null) {
					
					java.util.List<String> addVector = Collections.synchronizedList(new ArrayList<String>());	
					outVector.add(addVector);
					Object[] row =  {Main.drawCanvas.sampleList.get(i).getName()};
					model.addRow(row);
					
				}
			}
		}
		else {
			outVector.clear();
			model.getDataVector().clear();
			for(int i = 0; i<Main.samples; i++) {
				if(Main.drawCanvas.sampleList.get(i).samFile != null) {					
					java.util.List<String> addVector = Collections.synchronizedList(new ArrayList<String>());	
					outVector.add(addVector);
					Object[] row =  {Main.drawCanvas.sampleList.get(i).getName()};
					model.addRow(row);					
				}
			}
		}		
	}
	
public Average() {
	super(new GridBagLayout());		
	
	panel.setBackground(Draw.sidecolor);
	
	output.setEnabled(false);
	execute.setEnabled(false);
	cancel.setEnabled(false);
	
	GridBagConstraints c = new GridBagConstraints();
	c.anchor = GridBagConstraints.WEST;
	c.insets = new Insets(5, 5, 2, 5);
	c.gridx = 0;
	c.gridy = 0;
	c.gridwidth = 1;  
	panel.add(openbutton, c);
 //   c.gridy = 1;
   
    c.gridx = 1;
    panel.add(execute, c);
    c.gridx = 2;
    panel.add(cancel,c);
    c.gridx = 3;
    panel.add(output,c);   
    c.gridx = 4;
    panel.add(smallregion,c);
 //   c.gridx = 0;
    
   // c.gridy = 1;
  //  c.gridx = 0;
  //  panel.add(fileLabel, c);
    c.gridy = 1;
    c.gridx = 0;
    c.gridwidth =5;
    panel.add(fileLabel, c);
   
    
    c = new GridBagConstraints(
    		0,2, // position
    	    5,1, // size
    	    1.0,1.0, // fill ratio
    	    GridBagConstraints.CENTER, GridBagConstraints.BOTH, // position inside the cell
    	    new Insets(2,2,2,2),0,0);    
  //  c.gridwidth =4;
    table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
    infoScroll.getViewport().setBackground(Color.white);
    infoScroll.setPreferredSize(infoScroll.getPreferredSize());
//    infoScroll.setMinimumSize(infoScroll.getPreferredSize());
    table.setAutoCreateRowSorter(true);
          
    

    panel.add(infoScroll, c);
      
    add(panel, c);
    openbutton.addMouseListener(this);
    output.addMouseListener(this);
    execute.addMouseListener(this);
    cancel.addMouseListener(this);
    setFonts(Main.menuFont);
    
}

static void setFonts(Font menuFont) {
	for(int i = 0 ; i<Average.panel.getComponentCount(); i++) {
		Average.panel.getComponent(i).setFont(menuFont);
	}
	table.getTableHeader().setFont(menuFont);
	table.setFont(menuFont);
	table.setRowHeight(menuFont.getSize()+4);
}

public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
	//throw new UnsupportedOperationException("Not supported yet.");
	Component comp = getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
	
	if (table.getValueAt(row, column).toString().contains("Running")) {
		comp.setForeground(Color.red);
	} else {
		comp.setForeground(Color.YELLOW);
	}
	return (comp);
}
	
static class MyFilter extends javax.swing.filechooser.FileFilter {
	public boolean accept(File file) { 
			if (file.isDirectory()) {
				return true;
		    }	
			if (file.getName().contains(".bed") && !file.getName().contains(".tbi")) {
				return true;
			}					
			
			else {
				return false;
			}	
	} 
	
	public String getDescription() { return "*.bed"; }
} 



public static class runner extends SwingWorker<String, Object> {
	
	
	protected String doInBackground() {		
		if(smallregion.isSelected()) {
			calcAverageRegions();
		}
		else {
			calcAverage();
		}
		Average.execute.setEnabled(true);
		Average.output.setEnabled(true);
		Average.openbutton.setEnabled(true);
		return "";
	}
	
}
static void calcAverageRegions() {
	
	try {		
		 java.util.List<int[]> mergeVector = Average.mergeVector;
	 
	  int zeros = 0;
	 
	  int start = 0, end = 0;
	 		
	 
	  int stat = 0;
	  
	  String status = "";
	  String percent = "";
	  int[] coverageWidth = new int[10];
	  FileRead reader = new FileRead();
	  int[] coverageArray = {};
	  Reads readClass;
	  
	  Iterator<SAMRecord> iterator;
	 // SamReader inputSam;
	  for(int s=0;s<Main.drawCanvas.sampleList.size();s++) {
		 
		  if (Main.cancel) {
				break;
			}
		  if(Main.drawCanvas.sampleList.get(s).samFile == null) {
			 continue; 
		  }
		  //SamReader inputSam = SamReaderFactory.make().open(Main.drawCanvas.sampleList.get(s).samFile);
		  long result = 0, qualities = 0, cliplengths = 0, readlengths = 0, reads =0;
			
		  int currentLength = 0;
		  readClass = Main.drawCanvas.sampleList.get(s).getreadHash().get(Main.drawCanvas.splits.get(0));
		  coverageWidth = new int[10];
		  for(int i = 0; i<10; i++) {
			  coverageWidth[i] = 0;
		  }
		  boolean first = true;
		
		  for(int c = 0; c<mergeVector.size(); c++) {
			  			  
			  //inputSam = new SAMFileReader(Main.drawCanvas.sampleList.get(s).samFile);
			  //Iterator<SAMRecord> iterator = inputSam.queryOverlapping(Main.chromosomeDropdown.getItemAt(mergeVector.get(c)[0]), mergeVector.get(c)[1], mergeVector.get(c)[2]);
			  iterator = reader.getBamIterator(readClass, Main.chromosomeDropdown.getItemAt(mergeVector.get(c)[0]), mergeVector.get(c)[1], mergeVector.get(c)[2]);
			  SAMRecord samRecord = null;
			  
			  coverageArray = new int[mergeVector.get(c)[2]-mergeVector.get(c)[1]];
				for(int i = 0 ; i<coverageArray.length; i++) {
					coverageArray[i] = 0;
				}
				currentLength +=mergeVector.get(c)[2]-mergeVector.get(c)[1];
				if (Main.cancel) {
					break;
				}
				
				
			  while(iterator.hasNext()) {
				  
				  try {
					  samRecord = iterator.next();
					 
				  }
				  catch(Exception ex) {
					  continue;
				  }
				  if (Main.cancel) {
						break;
					}
				  	if(samRecord == null) {
				  		continue;
				  	}
				  	 
				  	end = samRecord.getAlignmentEnd();
				  	
					if(end == 0) {
						continue;
					}	
					
					start = samRecord.getAlignmentStart();
				//	if(start > mergeVector.get(c)[2]) {
						 
				//		if(!first) {
							
							
				//		}
						
						
				//	}
					
						
						
					
					reads++;
					qualities += samRecord.getMappingQuality();
					cliplengths += (samRecord.getUnclippedEnd()-samRecord.getUnclippedStart())-(end-start);
					readlengths += (samRecord.getUnclippedEnd()-samRecord.getUnclippedStart());				
					if(samRecord.getMappingQuality() == 0) {
						zeros++;  
						continue;
					}
					try {
					for(int i = start; i<end; i++ ) {
						if(i<mergeVector.get(c)[1]) {
							continue;
						}
						if(i > mergeVector.get(c)[2]-1) {
							
							break;
						}
						
						coverageArray[i-mergeVector.get(c)[1]]++;
					}
					}
					catch(Exception e) {
						ErrorLog.addError(e.getStackTrace());
						e.printStackTrace();
						continue;
					}
						
				   }
		  		
		  	//	inputSam.close();			
		  		for(int i = 0; i<coverageArray.length; i++) {
					result+=coverageArray[i];
					if(coverageArray[i] > 0) {
						if(coverageArray[i] >= 10) {
							for(int j = 0 ; j<10; j++) {
								coverageWidth[j]++;
							}									
						}
						else {
							for(int j = coverageArray[i]-1; j>=0; j--) {
								coverageWidth[j]++;
							}
						}
					}
				}
			
				if(stat != (int)(c/(double)mergeVector.size()*100)) {		
					if(currentLength != 0) {
			  			Average.model.setValueAt(""+MethodLibrary.round((result/(double)currentLength),2), s, 1);	
			  		}				  		
			  		percent = "1 : "+(int)((coverageWidth[0]/(double)currentLength)*100)
			  				+" | 2 : "+(int)((coverageWidth[1]/(double)currentLength)*100)
			  				+" | 3 : "+(int)((coverageWidth[2]/(double)currentLength)*100)
			  				+" | 4 : "+(int)((coverageWidth[3]/(double)currentLength)*100)
			  				+" | 5 : "+(int)((coverageWidth[4]/(double)currentLength)*100)
			  				+" | 6 : "+(int)((coverageWidth[5]/(double)currentLength)*100)
			  				+" | 7 : "+(int)((coverageWidth[6]/(double)currentLength)*100)
			  				+" | 8 : "+(int)((coverageWidth[7]/(double)currentLength)*100)
			  				+" | 9 : "+(int)((coverageWidth[8]/(double)currentLength)*100)
			  				+" | 10+ : "+(int)((coverageWidth[9]/(double)currentLength)*100);
			  		
			  		Average.model.setValueAt(MethodLibrary.round(qualities/(double)reads, 2), s, 2);
			  		Average.model.setValueAt(MethodLibrary.round(cliplengths/(double)readlengths,4), s, 3);
			  		Average.model.setValueAt(MethodLibrary.round(zeros/(double)reads,4), s, 4);
			  		
			  		Average.model.setValueAt(percent +"%", s, 5);
			  		stat = (int)(c/(double)mergeVector.size()*100);
			  		status = ""+c/(double)mergeVector.size()*100;						  	
			  		Average.model.setValueAt("Running " +status.substring(0,status.indexOf("."))+"%", s, 6);
				}
		  }
		  Average.model.setValueAt("100%", s, 6);		
		  
	  }
		}
		catch(Exception ex) {
			ErrorLog.addError(ex.getStackTrace());			
			ex.printStackTrace();
		}

}
static void calcAverage() {
	
	try {		
	  java.util.List<int[]> mergeVector = Average.mergeVector;
	  ArrayList<int[]> coverageMatrix = new ArrayList<int[]>();
	  for(int i = 0 ; i< mergeVector.size(); i++) {
		  int[] coverageArray = new int[mergeVector.get(i)[2]-mergeVector.get(i)[1]];
		 
		  coverageMatrix.add(coverageArray);
	  }
	 
	  int zeros = 0;
	  long result = 0, qualities = 0, cliplengths = 0, readlengths = 0, reads =0;
	
	  int currentLength = 0;
	  int start = 0, end = 0;
	  Reads readClass;
	  FileRead reader = new FileRead();
	  Iterator<SAMRecord> iterator;
	 
	  int pointer = 0, stat = 0, extrapointer = 0;
	  
	  String status = "";
	  String percent = "";
	  int[] coverageWidth = new int[10];
	  boolean first = true;
	  SamReader inputSam = null;
	  //int[] coverageArray = {};
	  for(int s=0;s<Main.drawCanvas.sampleList.size();s++) {
		 
		  if (Main.cancel) {
				break;
			}
		  if(Main.drawCanvas.sampleList.get(s).samFile == null) {
			 continue; 
		  }
		  //SamReader inputSam = SamReaderFactory.make().open(Main.drawCanvas.sampleList.get(s).samFile);
		   //SAMFileReader inputSam = new SAMFileReader(Main.drawCanvas.sampleList.get(s).samFile);
		  for(int i = 0 ; i< coverageMatrix.size(); i++) {			  
			  for(int j = 0 ; j<coverageMatrix.get(i).length; j++) {
				  coverageMatrix.get(i)[j] = 0;
			  }			 
		  }
		  pointer = 0;
		  extrapointer = 0;
		  qualities = 0;
		  readClass = Main.drawCanvas.sampleList.get(s).getreadHash().get(Main.drawCanvas.splits.get(0));
		  currentLength = 0;
		  result = 0;
		  zeros = 0;
		  reads = 0;
		  cliplengths = 0;
		  readlengths = 0;
		  coverageWidth = new int[10];
		  HashMap<String, Integer> chrmap = null;
		  //Iterator<SAMRecord> iterator = null;
		  //System.out.println(Main.chromosomeDropdown.getItemAt(startchrom) +" " +Average.startpos +" " +Average.endpos);
		  if(startchrom == endchrom) {
			  //iterator = inputSam.queryOverlapping(Main.chromosomeDropdown.getItemAt(startchrom), Average.startpos,  Average.endpos);
			  iterator = reader.getBamIterator(readClass, Main.chromosomeDropdown.getItemAt(startchrom),  Average.startpos, Average.endpos);
				 
		  }
		  else {
			  inputSam = SamReaderFactory.make().open(Main.drawCanvas.sampleList.get(s).samFile);
			  //inputSam = new SAMFileReader(Main.drawCanvas.sampleList.get(s).samFile);
			  chrmap = MethodLibrary.mapChrnameToIndex(inputSam.getFileHeader());
			  iterator = inputSam.iterator();
		  }
		 
		  SAMRecord samRecord = null;
		  for(int i = 0; i<10; i++) {
			  coverageWidth[i] = 0;
		  }
		  while(iterator.hasNext()) {
		 
			  try {
				  samRecord = iterator.next();
			  }
			  catch(Exception ex) {
				  ex.printStackTrace();
				  continue;
			  }
			  if (Main.cancel) {
					break;
				}
			  	if(samRecord == null) {
			  		continue;
			  	}
			 
			  	end = samRecord.getAlignmentEnd();
			  	
				if(end == 0) {
					continue;
				}	
				if(pointer > mergeVector.size()-1) {
					break;
				}
				
				if(chrmap.get(samRecord.getReferenceName()) == null) {
					continue;
				}
				 
				if(chrmap.get(samRecord.getReferenceName()) != mergeVector.get(pointer)[0]) {
					pointer = 0;
					int chrtest = chrmap.get(samRecord.getReferenceName());
					try {
						while(chrtest != mergeVector.get(pointer)[0] ) {
							pointer++;
							
						}
					}
					catch(Exception e) {
						break;
					}
					first = true;
				}
				
				/*if(samRecord.getReferenceIndex() < mergeVector.get(pointer)[0]) {
					
					continue;
				}
				if(mergeVector.get(pointer)[0] < samRecord.getReferenceIndex()) {
					
					pointer++;
					if(pointer > mergeVector.size()-1) {
						break;
					}	
					first = true;
				}*/
				if(end < mergeVector.get(pointer)[1]) {
					
					continue;
				}
				
				start = samRecord.getAlignmentStart();
				if(start > mergeVector.get(pointer)[2]) {
					
					if(!first) {
						
						for(int i = 0; i<coverageMatrix.get(pointer).length; i++) {
							result+=coverageMatrix.get(pointer)[i];
							if(coverageMatrix.get(pointer)[i] > 0) {
								if(coverageMatrix.get(pointer)[i] >= 10) {
									for(int j = 0 ; j<10; j++) {
										coverageWidth[j]++;
									}									
								}
								else {
									for(int j = coverageMatrix.get(pointer)[i]-1; j>=0; j--) {
										coverageWidth[j]++;
									}
								}
							}
						}
					
						if(stat != (int)(pointer/(double)mergeVector.size()*100)) {		
							if(currentLength != 0) {
					  			Average.model.setValueAt(""+MethodLibrary.round((result/(double)currentLength),2), s, 1);	
					  		}				  		
					  		percent = "1 : "+(int)((coverageWidth[0]/(double)currentLength)*100)
					  				+" | 2 : "+(int)((coverageWidth[1]/(double)currentLength)*100)
					  				+" | 3 : "+(int)((coverageWidth[2]/(double)currentLength)*100)
					  				+" | 4 : "+(int)((coverageWidth[3]/(double)currentLength)*100)
					  				+" | 5 : "+(int)((coverageWidth[4]/(double)currentLength)*100)
					  				+" | 6 : "+(int)((coverageWidth[5]/(double)currentLength)*100)
					  				+" | 7 : "+(int)((coverageWidth[6]/(double)currentLength)*100)
					  				+" | 8 : "+(int)((coverageWidth[7]/(double)currentLength)*100)
					  				+" | 9 : "+(int)((coverageWidth[8]/(double)currentLength)*100)
					  				+" | 10+ : "+(int)((coverageWidth[9]/(double)currentLength)*100);
					  		
					  		Average.model.setValueAt(MethodLibrary.round(qualities/(double)reads, 2), s, 2);
					  		Average.model.setValueAt(MethodLibrary.round(cliplengths/(double)readlengths,4), s, 3);
					  		Average.model.setValueAt(MethodLibrary.round(zeros/(double)reads,4), s, 4);
					  		
					  		Average.model.setValueAt(percent +"%", s, 5);
					  		stat = (int)(pointer/(double)mergeVector.size()*100);
					  		status = ""+pointer/(double)mergeVector.size()*100;						  	
					  		Average.model.setValueAt("Running " +status.substring(0,status.indexOf("."))+"%", s, 6);
						}
					}
					pointer++;
					if(pointer > mergeVector.size()-1) {
						break;
					}	
					first = true;
				}
				if(first) {
					
					
					currentLength +=mergeVector.get(pointer)[2]-mergeVector.get(pointer)[1];
					first = false;
				}
				reads++;
				qualities += samRecord.getMappingQuality();
				cliplengths += (samRecord.getUnclippedEnd()-samRecord.getUnclippedStart())-(end-start);
				readlengths += (samRecord.getUnclippedEnd()-samRecord.getUnclippedStart());				
				if(samRecord.getMappingQuality() == 0) {
					zeros++;  
					continue;
				}
				try {
					extrapointer = pointer;
					while(extrapointer < mergeVector.size() && end > mergeVector.get(extrapointer)[1]) {
						for(int i = start; i<end; i++ ) {
							if(i<mergeVector.get(extrapointer)[1]) {
								continue;
							}
							if(i > mergeVector.get(extrapointer)[2]-1) {								
								break;
							}
							
							coverageMatrix.get(extrapointer)[i-mergeVector.get(extrapointer)[1]]++;
						}
						extrapointer++;
					}
				}
				catch(Exception e) {
					ErrorLog.addError(e.getStackTrace());
					e.printStackTrace();
					continue;
				}
					
			    }
		  		Average.model.setValueAt("100%", s, 6);		
		  		inputSam.close();			
							  
		  }
		  
		}
		catch(Exception ex) {
			ErrorLog.addError(ex.getStackTrace());			
			ex.printStackTrace();
		}

}
static void calcAverageOld() {
	
	try {		
		 java.util.List<int[]> mergeVector = Average.mergeVector;
	 
	  int zeros = 0;
	  long result = 0, qualities = 0, cliplengths = 0, readlengths = 0, reads =0;
	
	  int currentLength = 0;
	  int start = 0, end = 0;
	  Reads readClass;
	  FileRead reader = new FileRead();
	  Iterator<SAMRecord> iterator;
	 
	  int pointer = 0, stat = 0;
	  
	  String status = "";
	  String percent = "";
	  int[] coverageWidth = new int[10];
	  boolean first = true;
	  SamReader inputSam = null;
	  int[] coverageArray = {};
	  for(int s=0;s<Main.drawCanvas.sampleList.size();s++) {
		 
		  if (Main.cancel) {
				break;
			}
		  if(Main.drawCanvas.sampleList.get(s).samFile == null) {
			 continue; 
		  }
		  //SamReader inputSam = SamReaderFactory.make().open(Main.drawCanvas.sampleList.get(s).samFile);
		   //SAMFileReader inputSam = new SAMFileReader(Main.drawCanvas.sampleList.get(s).samFile);
		  pointer = 0;
		  qualities = 0;
		  readClass = Main.drawCanvas.sampleList.get(s).getreadHash().get(Main.drawCanvas.splits.get(0));
		  currentLength = 0;
		  result = 0;
		  zeros = 0;
		  reads = 0;
		  cliplengths = 0;
		  readlengths = 0;
		  coverageWidth = new int[10];
		  HashMap<String, Integer> chrmap = null;
		  //Iterator<SAMRecord> iterator = null;
		  //System.out.println(Main.chromosomeDropdown.getItemAt(startchrom) +" " +Average.startpos +" " +Average.endpos);
		  if(startchrom == endchrom) {
			  //iterator = inputSam.queryOverlapping(Main.chromosomeDropdown.getItemAt(startchrom), Average.startpos,  Average.endpos);
			  iterator = reader.getBamIterator(readClass, Main.chromosomeDropdown.getItemAt(startchrom),  Average.startpos, Average.endpos);
				 
		  }
		  else {
			  inputSam = SamReaderFactory.make().open(Main.drawCanvas.sampleList.get(s).samFile);
			  //inputSam = new SAMFileReader(Main.drawCanvas.sampleList.get(s).samFile);
			  chrmap = MethodLibrary.mapChrnameToIndex(inputSam.getFileHeader());
			  iterator = inputSam.iterator();
		  }
		 
		  SAMRecord samRecord = null;
		  for(int i = 0; i<10; i++) {
			  coverageWidth[i] = 0;
		  }
		  while(iterator.hasNext()) {
		 
			  try {
				  samRecord = iterator.next();
			  }
			  catch(Exception ex) {
				  ex.printStackTrace();
				  continue;
			  }
			  if (Main.cancel) {
					break;
				}
			  	if(samRecord == null) {
			  		continue;
			  	}
			 
			  	end = samRecord.getAlignmentEnd();
			  	
				if(end == 0) {
					continue;
				}	
				if(pointer > mergeVector.size()-1) {
					break;
				}
				
				if(chrmap.get(samRecord.getReferenceName()) == null) {
					continue;
				}
				 
				if(chrmap.get(samRecord.getReferenceName()) != mergeVector.get(pointer)[0]) {
					pointer = 0;
					int chrtest = chrmap.get(samRecord.getReferenceName());
					try {
						while(chrtest != mergeVector.get(pointer)[0] ) {
							pointer++;
							
						}
					}
					catch(Exception e) {
						break;
					}
					first = true;
				}
				
				/*if(samRecord.getReferenceIndex() < mergeVector.get(pointer)[0]) {
					
					continue;
				}
				if(mergeVector.get(pointer)[0] < samRecord.getReferenceIndex()) {
					
					pointer++;
					if(pointer > mergeVector.size()-1) {
						break;
					}	
					first = true;
				}*/
				if(end < mergeVector.get(pointer)[1]) {
					
					continue;
				}
				
				start = samRecord.getAlignmentStart();
				if(start > mergeVector.get(pointer)[2]) {
					
					if(!first) {
						
						for(int i = 0; i<coverageArray.length; i++) {
							result+=coverageArray[i];
							if(coverageArray[i] > 0) {
								if(coverageArray[i] >= 10) {
									for(int j = 0 ; j<10; j++) {
										coverageWidth[j]++;
									}									
								}
								else {
									for(int j = coverageArray[i]-1; j>=0; j--) {
										coverageWidth[j]++;
									}
								}
							}
						}
					
						if(stat != (int)(pointer/(double)mergeVector.size()*100)) {		
							if(currentLength != 0) {
					  			Average.model.setValueAt(""+MethodLibrary.round((result/(double)currentLength),2), s, 1);	
					  		}				  		
					  		percent = "1 : "+(int)((coverageWidth[0]/(double)currentLength)*100)
					  				+" | 2 : "+(int)((coverageWidth[1]/(double)currentLength)*100)
					  				+" | 3 : "+(int)((coverageWidth[2]/(double)currentLength)*100)
					  				+" | 4 : "+(int)((coverageWidth[3]/(double)currentLength)*100)
					  				+" | 5 : "+(int)((coverageWidth[4]/(double)currentLength)*100)
					  				+" | 6 : "+(int)((coverageWidth[5]/(double)currentLength)*100)
					  				+" | 7 : "+(int)((coverageWidth[6]/(double)currentLength)*100)
					  				+" | 8 : "+(int)((coverageWidth[7]/(double)currentLength)*100)
					  				+" | 9 : "+(int)((coverageWidth[8]/(double)currentLength)*100)
					  				+" | 10+ : "+(int)((coverageWidth[9]/(double)currentLength)*100);
					  		
					  		Average.model.setValueAt(MethodLibrary.round(qualities/(double)reads, 2), s, 2);
					  		Average.model.setValueAt(MethodLibrary.round(cliplengths/(double)readlengths,4), s, 3);
					  		Average.model.setValueAt(MethodLibrary.round(zeros/(double)reads,4), s, 4);
					  		
					  		Average.model.setValueAt(percent +"%", s, 5);
					  		stat = (int)(pointer/(double)mergeVector.size()*100);
					  		status = ""+pointer/(double)mergeVector.size()*100;						  	
					  		Average.model.setValueAt("Running " +status.substring(0,status.indexOf("."))+"%", s, 6);
						}
					}
					pointer++;
					if(pointer > mergeVector.size()-1) {
						break;
					}	
					first = true;
				}
				if(first) {
					
					coverageArray = new int[mergeVector.get(pointer)[2]-mergeVector.get(pointer)[1]];
					for(int i = 0 ; i<coverageArray.length; i++) {
						coverageArray[i] = 0;
					}
					currentLength +=mergeVector.get(pointer)[2]-mergeVector.get(pointer)[1];
					first = false;
				}
				reads++;
				qualities += samRecord.getMappingQuality();
				cliplengths += (samRecord.getUnclippedEnd()-samRecord.getUnclippedStart())-(end-start);
				readlengths += (samRecord.getUnclippedEnd()-samRecord.getUnclippedStart());				
				if(samRecord.getMappingQuality() == 0) {
					zeros++;  
					continue;
				}
				try {
				for(int i = start; i<end; i++ ) {
					if(i<mergeVector.get(pointer)[1]) {
						continue;
					}
					if(i > mergeVector.get(pointer)[2]-1) {
						
						break;
					}
					
					coverageArray[i-mergeVector.get(pointer)[1]]++;
				}
				}
				catch(Exception e) {
					ErrorLog.addError(e.getStackTrace());
					e.printStackTrace();
					continue;
				}
					
			    }
		  		Average.model.setValueAt("100%", s, 6);		
		  		inputSam.close();			
							  
		  }
		  
		}
		catch(Exception ex) {
			ErrorLog.addError(ex.getStackTrace());			
			ex.printStackTrace();
		}

}
static void createAndShowGUI() {
	
	frame = new JFrame("Coverage calculator");    	
	JFrame.setDefaultLookAndFeelDecorated(true);
    frame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
  
    frame.setResizable(true);
    
    JComponent newContentPane = new Average();
    newContentPane.setOpaque(true); 
    frame.setContentPane(newContentPane);

    frame.pack();
 
}

@Override
public void mouseClicked(MouseEvent e) {
	if (e.getSource() == openbutton) {
		 try {
			
			 output.setEnabled(false);
			 execute.setEnabled(false);
	         cancel.setEnabled(false);
	    	 //String path = new java.io.File(".").getCanonicalPath();
	    	
	    	 //path = Main.trackDir;	
	    	 
	    	 TabixReader tabixReader = null;
	    	 BufferedReader bedfilereader = null;
	    	 /* JFileChooser chooser = new JFileChooser(path);
	    	  MyFilter bedFilter = new MyFilter();
	    	  chooser.setMultiSelectionEnabled(true);
	    	  chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
	    	  chooser.setAcceptAllFileFilterUsed(false);
	    	  chooser.addChoosableFileFilter(bedFilter); 	
	    	  
	    	  chooser.setDialogTitle("Add region file");
	          int returnVal = chooser.showOpenDialog((Component)this.getParent());
	          */
	          FileDialog fc = new FileDialog(frame, "Choose region file", FileDialog.LOAD);
	    	  fc.setDirectory(Main.trackDir);
	    	  fc.setFile("*.bed");
	    	  fc.setFilenameFilter(new FilenameFilter() {
		  			public boolean accept(File dir, String name) {
				        return name.endsWith(".bed") || name.endsWith(".bed.gz");
				     }
				 });
	    	  fc.setMultipleMode(false);
	    	  fc.setVisible(true);
	    	  String openfile = fc.getFile();			
	      
	       if (openfile != null) {
	    	   File bedfile = new File(fc.getDirectory() +"/" +openfile);
	    	   	 
	           Boolean tabix = false;
	          
	         // if (returnVal == JFileChooser.APPROVE_OPTION) {
	        	 
	        	//  File bedfile = chooser.getSelectedFile(); 
	        	 // Main.trackDir = chooser.getSelectedFile().getParent();	        	 
	        	  String line;
	        	  String[] split;	        	  	        	  
	        	  if(bedfile.getAbsolutePath().endsWith(".tbi")) {
	        		  Main.showError("Oops, you opened the index file.", "Note");
	        		  return;
	        	  }
	              if(bedfile.getAbsolutePath().endsWith(".gz")) {
	            	  tabix = true;
	            	  tabixReader = new TabixReader(bedfile.getCanonicalPath());	            	 
	              }            	 
	              
	              else {	            	  
	            	bedfilereader = new BufferedReader(new FileReader(bedfile));   			
	      			
	             }
	              
	              java.util.List<String[]> tempVector = Collections.synchronizedList(new ArrayList<String[]>());
	              
	      		while(true) {	      				
	      			
	      			if(tabix) {
	      				line = tabixReader.readLine();
	      			}
	      			else {
	      				line = bedfilereader.readLine();
	      			}
	      			
	      			if(line == null) {
	      				break;
	      			}
	      				      			
	      			if(line.startsWith("#")) {
	    				
	    				continue;
	    			}
	      			split = line.split("\\t+");	      			
	    			
	    			if(split.length < 3) {
	    				continue;
	    			}
	    			if(split[0].equals("chr")) {				 
	    				continue;	    				 
	    			}
	    		
	    			if(!split[0].matches("(chr)?\\d+")) {
	    				if(split[0].contains("X")) {
	    					split[0] = ""+Main.chromModel.getIndexOf("X");
	    				}
	    				else if(split[0].contains("Y")) {
	    					split[0] =""+Main.chromModel.getIndexOf("Y");
	    			
	    				}
	    				else if(split[0].matches("(chr)?MT?")) {
	    					if(Main.chromModel.getIndexOf("MT") > -1) {
	    						split[0] = ""+Main.chromModel.getIndexOf("MT");
	    					}
	    					else if(Main.chromModel.getIndexOf("M") > -1){
	    						split[0] = ""+Main.chromModel.getIndexOf("M");
	    					}
	    					
	    				}
	    				else {
	    					continue;
	    				}
	    			}
	    			else {
	    				if(split[0].contains("hr")) {
	    					split[0] = split[0].substring(3); 
	    				}
	    			}		
	      			split[0] = ""+(Integer.parseInt(split[0])-1);
	    			String[] add = { split[0], split[1], split[2]};		
	    			
					tempVector.add(add);
				 		
	      		
	      		}
	      		if(bedfilereader != null) { 
	      			bedfilereader.close();
	      		}
	      		
	      		Collections.sort(tempVector, new MethodLibrary.BEDSorter());
	      		Average.mergeVector = Collections.synchronizedList(new ArrayList<int[]>());
	      		int[] adderi = { 0, 0, 0 };   
	      		Average.mergeVector.add(adderi);		
	      		int index,start=0, end=0, chrom=0, size = 0; 
	      		boolean first  = true;
	      		
	      		for(int i = 0; i<tempVector.size(); i++) {
	      			
	      			chrom = Integer.parseInt(tempVector.get(i)[0]);
	      			
	      			start = Integer.parseInt(tempVector.get(i)[1]);
	      			if(first) {
	      				startchrom = chrom;
	      				startpos = start;
	      				first = false;
	      			}
	     			end = Integer.parseInt(tempVector.get(i)[2]);
	     			index = Average.mergeVector.size()-1;         		
	     			
	     				if(start > Average.mergeVector.get(index)[2]) {  
	     					int[] add = {chrom, start, end};
	     					Average.mergeVector.add(add);
	     					size += (end-start);
	     					
	     				}
	     				else {
	     					
	     					while(start <= Average.mergeVector.get(index)[1] && chrom == Average.mergeVector.get(index)[0]) {         						
	     						
	     						if(start <= Average.mergeVector.get(index)[1] && end >= Average.mergeVector.get(index)[1]) {
	     							if(end > Average.mergeVector.get(index)[2]) {         								
	     								Average.mergeVector.remove(index);
	     								
	     							}
	     							else {
	     								end = Average.mergeVector.get(index)[2];
	     								Average.mergeVector.remove(index);         								
	     							}
	     						}         						
	     						
	     						index--;
	     					}
	     					
	     					if(start < Average.mergeVector.get(index)[2] ) {         						
	     						if(end > Average.mergeVector.get(index)[2]) {
	     							Average.mergeVector.get(index)[2] = end;
	     						}        				         						
	     					}
	     					else {         						
	         					int[] add = {chrom, start, end};
	         					
	         					Average.mergeVector.add(index+1, add);
	         					size += (end-start);
	     					}
	     					
	     				}
	     				if(chrom != Average.mergeVector.get(index)[0]) {
	     					int[] add = {chrom, start, end};
	     					size += (end-start);
        					Average.mergeVector.add(index+1, add);
	     				}
	      			
	      		}	      	
  				endchrom = chrom;
  				endpos = end;      			
	      		Average.mergeVector.remove(0);
	      		output.setEnabled(true);
	        	execute.setEnabled(true);
	            cancel.setEnabled(true);	     
	            fileLabel.setText(""+bedfile.getName() +" | region size: " +size);	           
		      }	 
		 }
		 catch(Exception ex) {
			 ex.printStackTrace();
		 }
		
	}
	else if (e.getSource() == execute) {
		  execute.setEnabled(false);
		  openbutton.setEnabled(false);
		  output.setEnabled(false);
		  Main.cancel = false;		 
	   	  runner run = new runner();
	   	  model.setValueAt("0%", 0, 3);
	   	
	   	  run.execute();		
	}
	
	else if (e.getSource() == cancel) {
		
			Main.cancel = true;
		
		
	}
	
	else if(e.getSource() == output) {
		  try {
			  String outfilename = "";
		  	/*	
		  		String path ="";
		 
	    		path = new java.io.File(".").getCanonicalPath();
	    		    
	    	  JFileChooser chooser = new JFileChooser(path);
	    	  chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
	          int returnVal = chooser.showSaveDialog((Component)this.getParent());        
	          
	    	  
	         if (returnVal == JFileChooser.APPROVE_OPTION)
	      	{
	        	 outfilename = chooser.getSelectedFile().getAbsolutePath();
	      	}
	         */
	       BufferedWriter output;
	       FileDialog fs = new FileDialog(frame, "Save output file to...", FileDialog.SAVE);
 		    fs.setDirectory(Main.savedir);
 		    fs.setFile("*.tsv");
 		    fs.setVisible(true);			        	 
 		    
 		    String filename = fs.getFile();
	        
		    if(filename != null) {   		    	  
	    	
		    	outfilename = fs.getDirectory() +"/" +filename;   	
   		 		
				Main.savedir = fs.getDirectory();
   	 			Main.writeToConfig("DefaultSaveDir=" +Main.savedir);    
	    	  
    	  		File outfile = new File(outfilename);
    	    	output = new BufferedWriter(new FileWriter(outfile));
    	    	output.write("#" +fileLabel.getText() +"\n");
    	    	System.out.println(fileLabel.getText());
    	    	output.write("#Sample\tAverage coverage\tAverage mapping quality\tSoft clip rate\tZero quality rate\tCovered (%) (Coverage : Percentage)\n");
    	  		for(int j = 0; j< model.getRowCount(); j++) {
    	  			if(model.getValueAt(j, 1) != null) {
    	  				output.write(model.getValueAt(j, 0).toString() +"\t" +model.getValueAt(j, 1).toString() +"\t" +model.getValueAt(j, 2).toString() +"\t" +model.getValueAt(j, 3).toString() +"\t" +model.getValueAt(j, 4).toString() +"\t" +model.getValueAt(j, 5).toString()  +"\n");
    	  			}
    	  		}
    	  		output.close();
    	  	}	    	  
	    	  	
	      
		  }
		  catch(Exception ex) {
			
			  ex.printStackTrace();
		  }		  
	}		
	
}

@Override
public void mousePressed(MouseEvent e) {

	
}

@Override
public void mouseReleased(MouseEvent e) {
	
	
}

@Override
public void mouseEntered(MouseEvent e) {

	
}

@Override
public void mouseExited(MouseEvent e) {

	
}

}
