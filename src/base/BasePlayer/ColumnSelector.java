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
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.io.BufferedReader;

import java.io.FileInputStream;

import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;

import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTable;
import javax.swing.table.TableColumn;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;

public class ColumnSelector extends JPanel implements ActionListener, ComponentListener {
	private static final long serialVersionUID = 1L;
	JFrame frame = new JFrame("Column selector");    
	JPanel panel = new JPanel(new GridBagLayout());
	JTable table;	
	JLabel info = new JLabel("Select correct columns for values (if TSV-file, select at least chromosome and start values)");
	BedTrack track = null;
	JButton save = new JButton("OK");
	boolean found = false, tableSet = false;
	
	void setWindow() {
		
		try {
			JFrame.setDefaultLookAndFeelDecorated(false);
			 if(Main.ref == null) {
				 frame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE); 
				 frame.setVisible(true);
			 }
			 else {
				 frame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE); 
				 frame.setVisible(false);
			 }
			  
			 JComponent newContentPane = this;
			 
			 frame.setContentPane(newContentPane);
			 frame.setResizable(true);    
			
			 frame.setMinimumSize(new Dimension(600,300));	
			 frame.addComponentListener(this);
			 if(frame.isVisible()) {
				 createTable();
				
			 }
		
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
	    
	}
	
	public ColumnSelector() {
		super(new GridBagLayout());		
		setWindow();
	}
	public ColumnSelector(BedTrack track) {
		super(new GridBagLayout());		
		this.track = track;
		
		setWindow();
	}
	void createTable() {
		try {
		tableSet = true;
		GZIPInputStream gzip = null;
		String line;
		String[] split;
		BufferedReader reader;
		reader = null;
		String[][] data = null;
		if(this.track.file == null) {
			try {
			SeekableStream stream = SeekableStreamFactory.getInstance().getStreamFor(this.track.url);
			TabixReaderMod tabixReader = new TabixReaderMod(track.url.toString(), track.index.toString(), stream);
			int count = 0;
			String[] columns = null;
			
		  while((line = tabixReader.readLine()) != null) {	    	
		    			    	
		    	if(line.startsWith("#")) {		    		
		    		continue;
		    	}
		    	split = line.split("\\s+");
		    	if(count == 0) {
		    		data = new String[11][split.length];
		    		columns = new String[split.length];
		    		
		    		for(int i = 0; i<split.length; i++) {
			    		data[0][i] = "Select";
			    		if(track.chromcolumn != null && track.chromcolumn == i) {
			    			data[0][i] = "Chromosome";				    			
			    		}
			    		else if(track.startcolumn != null && track.startcolumn == i) {
			    			data[0][i] = "Start";				    			
			    		}
			    		else if(track.endcolumn != null && track.endcolumn == i) {
			    			data[0][i] = "End";				    			
			    		}
			    		else if(track.namecolumn != null && track.namecolumn == i) {
			    			data[0][i] = "Name";				    			
			    		}
			    		else if(track.valuecolumn != null && track.valuecolumn == i) {
			    			data[0][i] = "Value";				    			
			    		}
			    		else if(track.strandcolumn != null && track.strandcolumn == i) {
			    			data[0][i] = "Strand";				    			
			    		}
			    		else if(track.basecolumn != null && track.basecolumn == i) {
			    			data[0][i] = "Base";				    			
			    		}
			    	}
		    	}
		    	
		    	for(int i = 0; i<split.length; i++) {
		    		data[count+1][i] = split[i];
		    	}		
		    	
		    	count++;
		    	if(count > 9) {
		    		break;
		    	}		    
		  }
		  tabixReader.close();
		  stream.close();
		  for(int i = 0 ; i<columns.length; i++) {
			  columns[i] = "";
		  }
		
		panel.setBackground(Draw.sidecolor);
		table = new JTable(data, columns);
		for(int i = 0 ; i<columns.length; i++) {
			TableColumn column = table.getColumnModel().getColumn(i);			
			JComboBox<String> comboBox = new JComboBox<String>();
			comboBox.addItem("None");
			comboBox.addItem("Chromosome");
			comboBox.addItem("Start");
			comboBox.addItem("End");
			comboBox.addItem("Name");
			comboBox.addItem("Value");	
			comboBox.addItem("Strand");
			comboBox.addItem("Base");		
			column.setCellEditor(new DefaultCellEditor(comboBox));
			
			if(track.chromcolumn != null && track.chromcolumn == i) {
				comboBox.setSelectedItem("Chromosome");					
			}
		}
							
			
								
				
			}
			catch(Exception e) {
				e.printStackTrace();
				
				track.getHead().putNext(null);						
			}
		}
		else {
			gzip = new GZIPInputStream(new FileInputStream(this.track.file));			
			reader = new BufferedReader(new InputStreamReader(gzip));	
		
		
		
		int count = 0;
		String[] columns = null;
		
	  while((line = reader.readLine()) != null) {	    	
	    			    	
	    	if(line.startsWith("#")) {		    		
	    		continue;
	    	}
	    	split = line.split("\\s+");
	    	if(count == 0) {
	    		data = new String[11][split.length];
	    		columns = new String[split.length];
	    		
	    		for(int i = 0; i<split.length; i++) {
		    		data[0][i] = "Select";
		    		if(track.chromcolumn != null && track.chromcolumn == i) {
		    			data[0][i] = "Chromosome";				    			
		    		}
		    		else if(track.startcolumn != null && track.startcolumn == i) {
		    			data[0][i] = "Start";				    			
		    		}
		    		else if(track.endcolumn != null && track.endcolumn == i) {
		    			data[0][i] = "End";				    			
		    		}
		    		else if(track.namecolumn != null && track.namecolumn == i) {
		    			data[0][i] = "Name";				    			
		    		}
		    		else if(track.valuecolumn != null && track.valuecolumn == i) {
		    			data[0][i] = "Value";				    			
		    		}
		    		else if(track.strandcolumn != null && track.strandcolumn == i) {
		    			data[0][i] = "Strand";				    			
		    		}
		    		else if(track.basecolumn != null && track.basecolumn == i) {
		    			data[0][i] = "Base";				    			
		    		}
		    	}
	    	}
	    	
	    	for(int i = 0; i<split.length; i++) {
	    		data[count+1][i] = split[i];
	    	}		
	    	
	    	count++;
	    	if(count > 9) {
	    		break;
	    	}		    
	  }
	  
	  for(int i = 0 ; i<columns.length; i++) {
		  columns[i] = "";
	  }
	  reader.close(); 	
	panel.setBackground(Draw.sidecolor);
	table = new JTable(data, columns);
	for(int i = 0 ; i<columns.length; i++) {
		TableColumn column = table.getColumnModel().getColumn(i);			
		JComboBox<String> comboBox = new JComboBox<String>();
		comboBox.addItem("None");
		comboBox.addItem("Chromosome");
		comboBox.addItem("Start");
		comboBox.addItem("End");
		comboBox.addItem("Name");
		comboBox.addItem("Value");	
		comboBox.addItem("Strand");
		comboBox.addItem("Base");		
		column.setCellEditor(new DefaultCellEditor(comboBox));
		
		if(track.chromcolumn != null && track.chromcolumn == i) {
			comboBox.setSelectedItem("Chromosome");					
		}
	}
		}
	
	GridBagConstraints c = new GridBagConstraints();	
	c.gridx=0;
	c.gridy=0;
	table.setRowHeight(0, 30);		
	panel.add(save,c);
	c.gridx = 1;
	panel.add(info);
	//info.setForeground(Color.white);
	c.gridx=0;
	save.addActionListener(this);
	save.setPreferredSize(Main.buttonDimension);
	save.setMinimumSize(Main.buttonDimension);
	c = new GridBagConstraints(
	0,2, // position
    4,1, // size
    1.0,1.0, // fill ratio
    GridBagConstraints.CENTER, GridBagConstraints.BOTH, // position inside the cell
    new Insets(2,2,2,2),0,0);	
	table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
	table.setVisible(true);
	panel.add(table, c);
	add(panel,c);
	frame.addWindowListener(new java.awt.event.WindowAdapter() {
	    @Override
	    public void windowClosing(java.awt.event.WindowEvent windowEvent) {
	    	if(track.file.getName().endsWith(".tsv.gz")) {
				if(track.chromcolumn == null || track.startcolumn == null) {
					info.setForeground(Color.red);
					info.revalidate();
				}
				else {
					info.setForeground(Color.black);
					frame.dispose();
				}
			}
	    	else {
				frame.dispose();
			}	    	
	    }
	});
	 frame.pack();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	void setFonts(Font menuFont) {
		for(int i = 0 ; i<this.panel.getComponentCount(); i++) {
			this.panel.getComponent(i).setFont(menuFont);
		}
		if(table != null) {
			table.getTableHeader().setFont(menuFont);
			table.setFont(menuFont);
			table.setRowHeight(menuFont.getSize()+4);
		}
		frame.pack();
		
	}
	public static void main(String[] args) {
	//	BedTrack track = new BedTrack(new File("X:/cg8/BED/CADD/whole_genome_SNVs.tsv.gz"),0);
	//	 ColumnSelector select = new ColumnSelector(track);
		 /*
		 javax.swing.SwingUtilities.invokeLater(new Runnable() {
	           public void run() {
	           
	           		createAndShowGUI();
	           	
	           }
	       });
		*/
	}
	
	@Override
	public void actionPerformed(ActionEvent event ) {
		if(event.getSource() == save) {
			if(track != null) {
				for(int i = 0; i<table.getColumnCount(); i++) {
					String column = table.getValueAt(0, i).toString();
					
					if(column.contains("Chr")) {
						track.chromcolumn = i;
					}
					else if(column.contains("Sta")) {
						track.startcolumn = i;
					}
					else if(column.contains("End")) {
						track.endcolumn = i;
					}
					else if(column.contains("Name")) {
						track.namecolumn = i;
					}
					else if(column.contains("Value")) {
						track.valuecolumn = i;
					}
					else if(column.contains("Str")) {
						track.strandcolumn = i;
					}
					else if(column.contains("Base")) {
						track.basecolumn = i;
					}
					
				}	
				if(track.chromcolumn != null) {
					if(table.getValueAt(0, track.chromcolumn).toString().contains("Non")) {
						track.chromcolumn = null;
					}
				}
				if(track.startcolumn != null) {
					if(table.getValueAt(0, track.startcolumn).toString().contains("Non")) {
						track.startcolumn = null;
					}
				}
				if(track.endcolumn != null) {
					if(table.getValueAt(0, track.endcolumn).toString().contains("Non")) {
						track.endcolumn = null;
					}
				}
				if(track.namecolumn != null) {
					if(table.getValueAt(0, track.namecolumn).toString().contains("Non")) {
						track.namecolumn = null;
					}
				}
				if(track.valuecolumn != null) {
					if(table.getValueAt(0, track.valuecolumn).toString().contains("Non")) {
						track.valuecolumn = null;
					}
				}
				if(track.strandcolumn != null) {
					if(table.getValueAt(0, track.strandcolumn).toString().contains("Non")) {
						track.strandcolumn = null;
					}
				}
				if(track.basecolumn != null) {
					if(table.getValueAt(0, track.basecolumn).toString().contains("Non")) {
						track.basecolumn = null;
					}
				}
				
			}
			if((track.file != null && track.file.getName().endsWith(".tsv.gz")) || (track.url != null && track.url.toString().toLowerCase().endsWith(".tsv.gz"))) {
				if(track.chromcolumn == null || track.startcolumn == null) {
					info.setForeground(Color.red);
					info.revalidate();
				}
				else {
					info.setForeground(Color.black);
					FileRead.setBedTrack(track);
					frame.dispose();
				}
			}
			else {
				FileRead.setBedTrack(track);
				frame.dispose();
			}
		}		
	}
	@Override
	public void componentResized(ComponentEvent e) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void componentMoved(ComponentEvent e) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void componentShown(ComponentEvent e) {
		//frame.setLocation(Main.frame.getLocation().x, Main.frame.getLocation().y);
		frame.setLocation(Main.frame.getLocationOnScreen().x+Main.frame.getWidth()/2 - this.frame.getWidth()/2, Main.frame.getLocationOnScreen().y+Main.frame.getHeight()/6);
		if(!tableSet) {
			createTable();
		}
		
		// TODO Auto-generated method stub
		
	}
	@Override
	public void componentHidden(ComponentEvent e) {
		// TODO Auto-generated method stub
		
	}
}
