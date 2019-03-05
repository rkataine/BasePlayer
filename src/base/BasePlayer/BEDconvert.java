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
import java.awt.Component;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import base.BasePlayer.ExternalSort;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import java.util.HashMap;
import java.util.zip.GZIPInputStream;
import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JComboBox;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.SwingWorker;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;
import java.awt.FileDialog;
import javax.swing.table.TableColumnModel;


public class BEDconvert extends JPanel implements ActionListener, MouseListener {
	private static final long serialVersionUID = 1L;
	
	JFrame frame = new JFrame("BED converter");    
	JPanel panel = new JPanel(new GridBagLayout());
	
	static JTable table;	
	JLabel info = new JLabel("Open tab separated files for BED conversion.");
	Object[][] data = {};
	static HashMap<Integer, Boolean> editables = new HashMap<Integer, Boolean>();
	static HashMap<String, Integer> columnValues = new HashMap<String, Integer>();
	static JScrollPane tableScroll = new JScrollPane();
	static TableColumnModel columnModel;
	static File infile = null;
	 //String[] columns = null, selectcols = null;
	static DefaultTableModel tablemodel;
    String[][] selectdata = {};
    //static ArrayList<String[]> datarows;
    static String[] items = { "Select","Chromosome","Start","End","Name","Score","Strand" };
   // String[] header;
    Object[] headers = new Object[]{};     
    //String[] selheaders = new String[]{"select", "select", "select", "select"};    
    ArrayList<String> header = new ArrayList<String>();
    boolean changing = false;
    ArrayList<JComboBox<String>> comboboxes = new ArrayList<JComboBox<String>>();
    JTableHeader tableHeader;
   	//DefaultTableModel selectmodel = new DefaultTableModel(data, selheaders); 
   	JButton open = new JButton("Open");
   	//JButton add = new JButton("Add column");
	JButton write = new JButton("Convert");
   	boolean found = false, tableSet = false;
   // DefaultTableModel model = new DefaultTableModel(data, headers); 
  
    ActionListener comboActionListener = new ActionListener() {		    	

		public void actionPerformed(ActionEvent actionEvent) {
    		
    	  try {
    	
    	  if (actionEvent.getActionCommand() == "comboBoxChanged") {    	
    		 
    		 if(columnValues.size() == 0) {
    			 return;
    		 }
    		
    		 if(changing) {
    			 return;
    		 }
    		
    		  JComboBox<String> comboBox =  (JComboBox<String>)actionEvent.getSource();
    		 /* if(columnValues.get(comboBox.getSelectedItem()) != -1) {
    			  System.out.println(columnValues.get(comboBox.getSelectedItem()));
    			  JComboBox box = (JComboBox)selectortable.getCellEditor(0, columnValues.get(comboBox.getSelectedItem())).getTableCellEditorComponent(selectortable, 0, false, 0, columnValues.get(comboBox.getSelectedItem()));
    			  columnValues.put(comboBox.getSelectedItem().toString(), -1);
    			  box.setSelectedIndex(0);
    			  box.revalidate();
    		  }*/
    		
    		  checkHeadSelect(comboBox);
    		  checkHeaders();
    		 /* for(int i = 0; i<comboboxes.size(); i++) {    	
    			 
    			  columnValues.put(comboboxes.get(i).getSelectedItem().toString(), i);
    			  if(comboboxes.get(i).getSelectedItem().toString().equals("Editable")) {
    				  editables.put(i, true);
    			  }
    			  
    		  }		
    		  for(int i = 0; i < table.getColumnModel().getColumnCount(); i++) {  		    	    		    	
    		    	
    		    	EditableHeaderTableColumn col =(EditableHeaderTableColumn) table.getColumnModel().getColumn(i);
    		 	    if(col.getHeaderRenderer() == null) {
    		 	    	editables.put(i, true);
    		 	    }    		    	
    		  }
    		  */
    		
    	  }
    	  }
    	  catch(Exception e) {
    		  e.printStackTrace();
    	  }
		}
    };
    
	static int addedItems = 0;
    void checkHeadSelect(JComboBox<String> comboBox) {
    	 for(int i = 0; i<comboboxes.size(); i++) {
			  
			  
			if(comboboxes.get(i).equals(comboBox)) {
				  continue;
			  }
			  if(comboboxes.get(i).getSelectedItem().toString().equals("Editable")) {
				  continue;
			  }
			  if(comboboxes.get(i).getSelectedItem().toString().equals(comboBox.getSelectedItem().toString())) {    				  
				  changing = true;    				 
				 
				  comboboxes.get(i).setSelectedIndex(0);
				  EditableHeaderTableColumn col =(EditableHeaderTableColumn) table.getColumnModel().getColumn(i);
				  col.setHeaderValue("Select");
				 // comboboxes.get(i).setSelectedIndex(0);
				  //selectortable.setValueAt("select", 0, i);
				  comboboxes.get(i).revalidate();
				  comboboxes.get(i).repaint();
				  table.revalidate();
				  table.repaint();
				  changing = false;
				  break;    				  
			  }
		  }
    }
	void checkHeaders() {
    	columnValues.put("Select", -1);	
    	columnValues.put("Chromosome", -1);
    	columnValues.put("Start", -1);
    	columnValues.put("End", -1);
    	columnValues.put("Name", -1);
    	columnValues.put("Score", -1);
    	columnValues.put("Strand", -1);
		// editables.clear();
		  
    	 for(int i = 0; i < table.getColumnModel().getColumnCount(); i++) {  		    	    		    	
		    	
	    	EditableHeaderTableColumn col =(EditableHeaderTableColumn) table.getColumnModel().getColumn(i);	 	    	
 	    	columnValues.put(col.getHeaderValue().toString(), i);		 	    	
		 	
		  }    
    	 table.revalidate();
    	 tableScroll.revalidate();
    	 table.repaint();
    }
   /* DefaultTableModel model = new DefaultTableModel(data, headers) {

     
		private static final long serialVersionUID = 1L;

		@Override
        public boolean isCellEditable(int row, int column) {
           //all cells false
           return false;
        }
    };*/
   
	
	void setWindow() {
		try {
			JFrame.setDefaultLookAndFeelDecorated(false);
			//frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
			 if(VariantHandler.aminoCount == null) {				
				 frame.setVisible(true);
				 frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
			 }
			 else {				
				 frame.setVisible(false);
			 }
			 
			// JComponent newContentPane = this;			 
			 frame.setContentPane(this);
			 frame.setResizable(true);			
			 frame.setMinimumSize(new Dimension(600,300));			
			 createTable();
		
		}
		catch(Exception e) {
			e.printStackTrace();
		}	    
	}
	
	public BEDconvert() {
		super(new GridBagLayout());		
		
		setWindow();
	}
	void createTable(File file) {
		try {
			
			addedItems = 0;
		    int headercount = 4, rows=0;		   
		    BufferedReader reader = null;
		    GZIPInputStream gzip = null;
			  String line;
			  tablemodel = new DefaultTableModel(data,new Object[]{"","","",""});
			  table = new JTable(tablemodel) {
			    	
					private static final long serialVersionUID = 1L;
					

					@Override
					public boolean isCellEditable(int rowIndex, int columnIndex) {				
						
						return false;
					}
			    };
			  
			    
			    table.addMouseListener(this);
			   // columns = new String[tablelength];	
			    columnModel = table.getColumnModel();
			    
			    table.setTableHeader(new EditableHeader(columnModel));
			    
			  if(file.getName().endsWith(".gz")) {
				  try {
					  gzip = new GZIPInputStream(new FileInputStream(file));
					  reader = new BufferedReader(new InputStreamReader(gzip));		
				  }
				  catch(Exception e) {
					  e.printStackTrace();
				  }
			  }
			  else {
				  reader = new BufferedReader(new FileReader(file));	
			  }
		   
		    comboboxes.clear();
		   
		    tablemodel.setRowCount(0);
		    int rowlimit = 20;
		    boolean first = true;
		    header.clear();
			while((line = reader.readLine()) != null) {
				if(line.startsWith("!")) {
					continue;
				}
				if(line.startsWith("#") || line.startsWith("track")) {
					if(!line.startsWith("#")) {
						line = "#" +line;
					}
					header.add(line);
					continue;
				}
				if(first) {
					if(header.size() > 0) {
						String[] row = header.get(header.size()-1).split("\t");
						if(row.length > headercount) {
							for(int i=0;i<row.length-headercount; i++) {
								tablemodel.addColumn("new");
							}
							headercount = row.length;					
						}		    
						tablemodel.addRow(row);
						
					}
					first = false;
				}
				if(rows > rowlimit) {
					//Main.showError("Showing " +rowlimit +" rows.", "Note", BEDconvert.tableScroll);
					break;
				}
				String[] row = line.split("\t");
				rows++;
				
				//datarows.add(row);
				if(row.length > headercount) {
					
					for(int i=0;i<row.length-headercount; i++) {
						tablemodel.addColumn("new");
					}
					headercount = row.length;					
				}		    
				
				tablemodel.addRow(row);
			}
			reader.close();
			
			if(rows == 0) {
				Main.showError("Could not find any rows in the file.", "Error");
				return;
			}
			/*table = new JTable(rows,headercount) {
		    	
				private static final long serialVersionUID = 1L;
				

				@Override
				public boolean isCellEditable(int rowIndex, int columnIndex) {
					
					if(editables.containsKey(columnIndex)) {
						return true;
					}
					return false;
				}
		    };
		   */
		    	//table.addMouseListener(tableMouseListener);
			   // columns = new String[tablelength];	
			    columnModel = table.getColumnModel();
			    table.setTableHeader(new EditableHeader(columnModel));
			    //table.addMouseListener(this);
			    //frame.addMouseListener(this);
			    String[] items = { "Select","Chromosome","Start","End","Name","Score", "Strand" };					   
			    ComboRenderer renderer = new ComboRenderer(items);
			    ColorColumnRenderer columnRenderer = new ColorColumnRenderer(new Color(150,255,150), table);
			    for(int i = 0; i < columnModel.getColumnCount(); i++) {
			    	
			    	JComboBox<String> combo = new JComboBox<String>();
			    	for (int c = 0; c < items.length; c++) {
			   	      combo.addItem(items[c]);
			   	    }
			    	
			    	combo.addActionListener(comboActionListener);
			    	
			    	EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(i);
			    	col.setHeaderValue("Select");
			    	col.setCellRenderer(columnRenderer);
			 	    col.setHeaderRenderer(renderer);
			 	    col.setHeaderEditor(new DefaultCellEditor(combo));
			 	   
			 	    comboboxes.add(combo);
			 	   
			    }
			    
			   checkTable(table);
			   /* for(int i = 0; i<rows;i++) {
			    	for(int j=0;j<datarows.get(i).length;j++) {
			    		table.setValueAt(datarows.get(i)[j], i, j);
			    	}			    	
			    }*/
			    tableScroll.getViewport().removeAll();
				tableScroll.getViewport().add(table);
			/*String[] row = new String[headercount];
			for(int i = 0 ; i<headercount; i++ ) {
				row[i] = "select";
			}
			
			selectmodel.setRowCount(0);
			selectmodel.addRow(row);
			//selectortable = new JTable(selectdata,columns);
			//selectcols = new String[headercount];		
			/*
			for(int i = 0 ; i<selectcols.length; i++) {
				TableColumn column = selectortable.getColumnModel().getColumn(i);		
				JComboBox<String> comboBox = new JComboBox<String>();
				comboBox.addItem("select");
				comboBox.addItem("Position");
				comboBox.addItem("Chromosome");
				comboBox.addItem("Start");
				comboBox.addItem("End");				
				comboBox.addItem("Sample");			
				column.setCellEditor(new DefaultCellEditor(comboBox));		
			}
			 */
			 //table = new JTable(model);
			 table.revalidate();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	void checkTable(JTable table) {
		TableColumnModel columnModel = table.getColumnModel();
		if(table.getRowCount() < 2) {
			return;
		}
		int chrindex = -1, startindex = -1, endindex = -1, nameindex =-1, scoreindex = -1, strandindex = -1; 
		//SEARCH chrom
		for(int i = 0; i < columnModel.getColumnCount(); i++) {
			try {
	    	if(table.getValueAt(0, i).toString().contains("chr") || table.getValueAt(1, i).toString().contains("chr") || table.getValueAt(0, i).toString().contains("genoName")) {
	    		
	    		chrindex = i;
	    		EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(i);
	    		comboboxes.get(i).setSelectedItem("Chromosome");
	    		col.setHeaderValue(comboboxes.get(i).getSelectedItem().toString());
	    		checkHeadSelect(comboboxes.get(i));	    		
	    		break;
	    	}
			}
			catch(Exception e) {
				System.out.println(table.getValueAt(0, 0).toString());
				e.printStackTrace();
			}
	    }
		//SEARCH start
		
		
			for(int i = 0; i < columnModel.getColumnCount()-1; i++) {	 
				if(table.getValueAt(0, i+1) == null ||table.getValueAt(0, i) == null) {
					continue;
				}
		    	if(table.getValueAt(0, i).toString().toLowerCase().contains("start") && table.getValueAt(0, i+1).toString().toLowerCase().contains("end")) {
		    		try {
		    			Integer.parseInt(table.getValueAt(1, i).toString());
		    			Integer.parseInt(table.getValueAt(1, i+1).toString());
		    		}
		    		catch(Exception e) {
		    			continue;
		    		}
		    		startindex = i;
		    		endindex = i+1;
		    		EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(i);
		    		comboboxes.get(i).setSelectedItem("Start");
		    		comboboxes.get(i+1).setSelectedItem("End");
		    		col.setHeaderValue(comboboxes.get(i).getSelectedItem().toString());
		    		col =(EditableHeaderTableColumn) columnModel.getColumn(i+1);
		    		col.setHeaderValue(comboboxes.get(i+1).getSelectedItem().toString());		
		    		break;
		    	}
		    }
			if(startindex < 0) {
				if(columnModel.getColumnCount() > chrindex +2) {
					int start=-1, end =0;
					try {
		    			start = Integer.parseInt(table.getValueAt(1, chrindex+1).toString());
		    			end = Integer.parseInt(table.getValueAt(1, chrindex+2).toString());
		    		}
		    		catch(Exception e) {
		    			
		    		}
					if(start > -1) {
						if(end > start) {
							startindex = chrindex+1;
							endindex = chrindex+2;
							
				    		EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(startindex);
				    		comboboxes.get(startindex).setSelectedItem("Start");		    		
				    		col.setHeaderValue(comboboxes.get(startindex).getSelectedItem().toString());
				    		col =(EditableHeaderTableColumn) columnModel.getColumn(endindex);
				    		comboboxes.get(endindex).setSelectedItem("End");		    		
				    		col.setHeaderValue(comboboxes.get(endindex).getSelectedItem().toString());
						}
					}
				}
			}
			//SEARCH name
			for(int i = 0; i < columnModel.getColumnCount(); i++) {	 
				if(table.getValueAt(0, i) == null) {
					continue;
				}
		    	if(table.getValueAt(0, i).toString().toLowerCase().contains("name") && i != chrindex) {
		    		nameindex = i;
		    		EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(i);
		    		comboboxes.get(i).setSelectedItem("Name");		    		
		    		col.setHeaderValue(comboboxes.get(i).getSelectedItem().toString());
		    		break;
		    	}
			}
			if(nameindex < 0) {
				if(endindex == 2 && table.getColumnCount() > 3 && table.getValueAt(1,3) !=null) {
					
		    		nameindex = 3;
		    		EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(3);
		    		if(col.getHeaderValue().toString().equals("Select")) {
		    			comboboxes.get(3).setSelectedItem("Name");		    		
			    		col.setHeaderValue(comboboxes.get(3).getSelectedItem().toString());			    		
		    		} 		    	
				}				
			}
			for(int i = 0; i < columnModel.getColumnCount(); i++) {	    
				if(table.getValueAt(0, i) == null) {
					continue;
				}
		    	if(table.getValueAt(0, i).toString().toLowerCase().contains("score")) {
		    		try {
		    			Double.parseDouble(table.getValueAt(1, i).toString());
		    			
		    		}
		    		catch(Exception e) {
		    			continue;
		    		}
		    		
		    		scoreindex = i;
		    		EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(i);
		    		comboboxes.get(i).setSelectedItem("Score");		    		
		    		col.setHeaderValue(comboboxes.get(i).getSelectedItem().toString());
		    		break;
		    	}
			}
			
			for(int i = 0; i < columnModel.getColumnCount(); i++) {	    	
				if(table.getValueAt(0, i) == null) {
					continue;
				}
		    	if(table.getValueAt(0, i).toString().toLowerCase().contains("strand")) {		    		
		    		if(table.getValueAt(1, i).toString().equals("+") || table.getValueAt(1, i).toString().equals("-")) {
		    			strandindex = i;
			    		EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(i);
			    		comboboxes.get(i).setSelectedItem("Strand");		    		
			    		col.setHeaderValue(comboboxes.get(i).getSelectedItem().toString());
			    		break;
		    		}		    		
		    	}
			}
			if(strandindex < 0) {
				for(int i = 0; i < columnModel.getColumnCount(); i++) {	    	
					if(table.getValueAt(1, i) == null) {
						continue;
					}
		    		if(table.getValueAt(1, i).toString().equals("+") || table.getValueAt(1, i).toString().equals("-")) {
		    			strandindex = i;
			    		EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(i);
			    		comboboxes.get(i).setSelectedItem("Strand");		    		
			    		col.setHeaderValue(comboboxes.get(i).getSelectedItem().toString());
			    		break;
		    		}		    	
				}
			}
			if(scoreindex < 0 && chrindex > -1) {
				for(int i = 0; i < columnModel.getColumnCount(); i++) {	    	
					if(table.getValueAt(1, i) == null) {
						continue;
					}
					
					try {
		    			Double.parseDouble(table.getValueAt(1, i).toString());
		    			
		    		}
		    		catch(Exception e) {
		    			continue;
		    		}		    		
		    			
		    		EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(i);
		    		if(col.getHeaderValue().equals("Select")) {
		    			comboboxes.get(i).setSelectedItem("Score");		    		
			    		col.setHeaderValue(comboboxes.get(i).getSelectedItem().toString());
			    		break;
		    		}
			    	
		    				    	
				}
			}
			checkHeaders();
	}
	void createTable() {
		try {
		
	    tablemodel = new DefaultTableModel(data,new Object[]{"","","",""});
	    table = new JTable(tablemodel) {
	    	
			private static final long serialVersionUID = 1L;
			

			@Override
			public boolean isCellEditable(int rowIndex, int columnIndex) {				
				
				return false;
			}
	    };
	    for(int i = 0; i<40; i++) {
	    	tablemodel.addRow(new String[] {"","","",""});		   
	    }
	    
	    table.addMouseListener(this);
	   // columns = new String[tablelength];	
	    columnModel = table.getColumnModel();
	    
	    table.setTableHeader(new EditableHeader(columnModel));
	    //table.addMouseListener(this);
	    frame.addMouseListener(this);
	    
	    ComboRenderer renderer = new ComboRenderer(items);
	    ColorColumnRenderer columnRenderer = new ColorColumnRenderer(new Color(150,255,150), table);
	    for(int i = 0; i < columnModel.getColumnCount(); i++) {
	    	
	    	JComboBox<String> combo = new JComboBox<String>();
	    	for (int c = 0; c < items.length; c++) {
	   	      combo.addItem(items[c]);
	   	    }
	    	//combo.setName(""+i);
	    	combo.addActionListener(comboActionListener);
	    	
	    	EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(i);
	    	col.setCellRenderer(columnRenderer);
	 	    col.setHeaderValue(combo.getItemAt(0));
	 	    col.setHeaderRenderer(renderer);
	 	    col.setHeaderEditor(new DefaultCellEditor(combo));	 	   
	 	    comboboxes.add(combo);	 	   
	    }
	  
		panel.setBackground(Draw.sidecolor);		
		
		
	columnValues.put("Select", -1);	
	columnValues.put("Chromosome", -1);
	columnValues.put("Start", -1);
	columnValues.put("End", -1);
	columnValues.put("Name", -1);
	columnValues.put("Score", -1);
	columnValues.put("Strand", -1);
	GridBagConstraints c = new GridBagConstraints();	
	c.gridx=0;
	c.gridy=0;
	c.anchor = GridBagConstraints.NORTHWEST;
	//table.setRowHeight(0, 30);	
	//add.setToolTipText("Add column");
	//add.addActionListener(this);
	
	panel.add(open,c);
	c.gridx++;
	panel.add(write,c);
	c.gridx++;
	panel.add(info,c);
	c.gridx++;;
	//c.anchor = GridBagConstraints.EAST;
	//panel.add(add,c);
	
	c.gridx=0;
	open.addActionListener(this);
	write.addActionListener(this);
	write.setPreferredSize(Main.buttonDimension);
	write.setMinimumSize(Main.buttonDimension);
	open.setPreferredSize(Main.buttonDimension);
	open.setMinimumSize(Main.buttonDimension);
	
	c.gridy++;
	c.gridwidth = 4;
	c.weightx = 1.0;
	
	c.fill = GridBagConstraints.HORIZONTAL;
	table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
	table.setVisible(true);
	//selectortable.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
	//selectortable.setVisible(true);
	//panel.add(selectortable, c);
	//c.gridy++;
	c.weighty = 1.0;
	c.fill = GridBagConstraints.BOTH;
	tableScroll.getViewport().add(table);
	panel.add(tableScroll, c);
	add(panel,c);
	
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
		//	selectortable.getTableHeader().setFont(menuFont);
		//	selectortable.setFont(menuFont);
		}
		frame.pack();
		
	}
	class OutputRunner extends SwingWorker<String, Object> {
		File outfile;
		public OutputRunner(File outfile) {
			this.outfile = outfile;
		}
		protected String doInBackground() {
			if(Main.drawCanvas != null) {
				Main.drawCanvas.loading("Converting " +outfile.getName());
			}
			frame.setState(Frame.ICONIFIED);
			frame.setCursor( Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
			ExternalSort.processFile(infile, outfile);
			frame.setCursor( Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
			if(Main.drawCanvas != null) {
				Main.drawCanvas.ready("all");
			}
			frame.setState(Frame.NORMAL);
			 if (JOptionPane.showConfirmDialog(BEDconvert.tableScroll, "Open BED file?", "Open?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE) == JOptionPane.YES_OPTION){
		    	   File[] bedfile = new File[1];
		    	   bedfile[0] = outfile;
		    	   FileRead filereader = new FileRead(bedfile);		        	
		          filereader.readBED(bedfile);
		        }
			return "";
		}
	}
	public static void main(String[] args) {
		BEDconvert conv = new BEDconvert();
		
		conv.frame.setVisible(true);
		conv.frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

	}
	
	@Override
	public void actionPerformed(ActionEvent event ) {
		if(event.getSource() == write) {
			try {	  	    
				File savefile = null;
				FileDialog fs = new FileDialog(frame, "Save BED file to...", FileDialog.SAVE);
		  		  fs.setDirectory(Main.trackDir);
		  		  fs.setFile("*.bed.gz");
		  		  fs.setVisible(true);	         
		        	 
		        	 
		        	 while(true) {
		        		String filename = fs.getFile();
			        
					    if(filename != null) {   		    	  
				    	   savefile = new File(fs.getDirectory() +"/" +filename);
					       if(!Files.isWritable(Paths.get(savefile.getParent()))) {
					    	  Main.showError("No permission to write.", "Error", tableScroll);
					    	  continue;
					       }
					       //outfile = chooser.getSelectedFile();
					       if(!savefile.getName().endsWith(".bed.gz")) {
					    	   if(savefile.getName().endsWith(".bed")) {
					    		   savefile = new File(savefile.getCanonicalPath() +".gz");
					    	   }
					    	   else {
					    		   savefile = new File(savefile.getCanonicalPath() +".bed.gz");
					    	   }				    	   
					       }
						       break;
					      }
				    	 else {					    	 
				    		 break;					    	 
				    	 }		        	 
		         }
		         
		    	/* JFileChooser chooser = new JFileChooser(Main.trackDir);
		    
		    	 chooser.setDialogTitle("Save BED file in...");
		    	  File outfile = null;
		    	  while(true) {
		    		  int returnVal = chooser.showSaveDialog((Component)this.getParent());	           	  
		          
		        
				    	 if(returnVal == JFileChooser.APPROVE_OPTION) {   		    	  
					       
					       if(!Files.isWritable(Paths.get(chooser.getSelectedFile().getParent()))) {
					    	  Main.showError("No permission to write.", "Error", tableScroll);
					    	  continue;
					       }
					       outfile = chooser.getSelectedFile();
					       if(!outfile.getName().endsWith(".bed.gz")) {
					    	   if(outfile.getName().endsWith(".bed")) {
					    		   outfile = new File(outfile.getCanonicalPath() +".gz");
					    	   }
					    	   else {
					    		   outfile = new File(outfile.getCanonicalPath() +".bed.gz");
					    	   }				    	   
					       }
					       break;
				      }
				    	 if(returnVal == JFileChooser.CANCEL_OPTION) {   	
					    	  outfile = null;
					    	  
					    	  break;
					      }
		    	  }*/
		    	  if(savefile != null) {
				       if(columnValues.get("Chromosome") < 0) {
				    	 Main.showError("Select chromosome column.", "Note", BEDconvert.tableScroll);
					     return;
				       }
				       if(columnValues.get("Start") < 0) {
				    	   Main.showError("Select start position column.", "Note", BEDconvert.tableScroll);
						     return;
				       }
				       /*if(columnValues.get("End") < 0) {
				    	   Main.showError("Select end position column.", "Note", BEDconvert.tableScroll);
						     return;
				       }*/
				       ExternalSort.columns.clear();
				       ExternalSort.columns.add(-1);
				       ExternalSort.columns.add(-1);
				       ExternalSort.columns.add(-1);
				       ExternalSort.columns.add(-1);
				       ExternalSort.columns.add(-1);
				       ExternalSort.columns.add(-1);
				       
				       if(columnValues.get("Strand") > -1) {
				    	   ExternalSort.columns.set(5,columnValues.get("Strand"));
				       }
				       else {
				    	   ExternalSort.columns.remove(5);
				       }
				       if(columnValues.get("Score") > -1) {				    	   
				    	   ExternalSort.columns.set(4,columnValues.get("Score"));
				       }
				       else {
				    	   if(ExternalSort.columns.size() == 5) {
				    		   ExternalSort.columns.remove(4);
				    	   }
				       }
				       if(columnValues.get("Name") > -1) {				    	   
				    	   ExternalSort.columns.set(3,columnValues.get("Name"));
				       }
				       else {
				    	   if(ExternalSort.columns.size() == 4) {
				    		   ExternalSort.columns.remove(3);
				    	   }
				       }
				       ExternalSort.columns.set(0, columnValues.get("Chromosome"));
				       ExternalSort.columns.set(1, columnValues.get("Start"));
				       if(columnValues.get("End") >= 0) {
				    	   ExternalSort.columns.set(2, columnValues.get("End"));
				       }
				       else {
				    	   ExternalSort.columns.set(2, columnValues.get("Start"));
				       }
				       ExternalSort.chrindex = columnValues.get("Chromosome");
				       
				       ExternalSort.startindex = columnValues.get("Start");
				       Main.trackDir = savefile.getParent();
				       Main.writeToConfig("DefaultTrackDir=" +Main.trackDir);
				       OutputRunner runner = new OutputRunner(savefile);
				       runner.execute();
				      // MethodLibrary.blockCompressAndIndex(infile, outfile, true); 
				      
		    	  }
				}
				catch(Exception ex) {
					ex.printStackTrace();
				}
		}
		else if(event.getSource() == open) {
			FileDialog chooser = new FileDialog(frame, "Choose a tab separated file", FileDialog.LOAD);
	  		  chooser.setDirectory(Main.trackDir);
	  		  //chooser.setFile("*.tsv");
	  		  chooser.setVisible(true);
	  		  String filename = chooser.getFile();			
	        
	         if (filename != null) {
	        	 File addfile = new File(chooser.getDirectory() +"/" +filename);
	        	
	        	 if(addfile.exists()) {
	        		 infile = addfile;	 	        	
		        	 createTable(infile);
		        	 info.setText("Select correct header values for BED columns. File: " +infile.getName());
	        	 }
	        	 else {
	        		 Main.showError("File does not exists.", "Error", frame);
	        	 }
	        	 
	         }
	     //  JFileChooser chooser = new JFileChooser(Main.trackDir);	 
			  
		    	//  chooser.setMultiSelectionEnabled(false);
		    	 // chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
		    	//  chooser.setAcceptAllFileFilterUsed(false);	    	  
		    	  
		    	 
		    	 // chooser.setDialogTitle("Open file");
		    	 // int returnVal = chooser.showOpenDialog((Component)this.getParent());	         
		}
		/*else if(event.getSource() == add) {
			
			tablemodel.addColumn("NewColumn");			
			
			columnModel = table.getColumnModel();
		    table.setTableHeader(new EditableHeader(columnModel));
		    ComboRenderer renderer = new ComboRenderer(items);		  
		  
		    for(int i = 0; i<comboboxes.size(); i++) {
		    	EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(i);
		 	    col.setHeaderValue(comboboxes.get(i).getSelectedItem());
		 	    col.setHeaderRenderer(renderer);
		 	    col.setHeaderEditor(new DefaultCellEditor(comboboxes.get(i)));
		 	 
		    }
		    
		    checkHeaders();
			table.revalidate();
			
		}*/
	}	
	

		
		@Override
		public void mouseClicked(MouseEvent arg0) {
			// TODO Auto-generated method stub
			
		}

		@Override
		public void mouseEntered(MouseEvent arg0) {
			
			EditableHeader header = (EditableHeader)table.getTableHeader();
			header.removeEditor();
			// TODO Auto-generated method stub
			
		}

		@Override
		public void mouseExited(MouseEvent arg0) {
			
			// TODO Auto-generated method stub
			EditableHeader header = (EditableHeader)table.getTableHeader();
			header.removeEditor();
		}

		@Override
		public void mousePressed(MouseEvent arg0) {
			
			
		}

		@Override
		public void mouseReleased(MouseEvent arg0) {
			// TODO Auto-generated method stub
			
		}
	
}

