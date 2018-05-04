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
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.nio.file.Files;
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
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;

public class TableBrowser extends JPanel implements ActionListener, MouseListener {
	private static final long serialVersionUID = 1L;
	
	JFrame frame = new JFrame("Table Browser");    
	JPanel panel = new JPanel(new GridBagLayout());
	static JTable table;	
	JLabel info = new JLabel("Open TSV files.");
	Object[][] data = {};
	static HashMap<Integer, Boolean> editables = new HashMap<Integer, Boolean>();
	static HashMap<String, Integer> columnValues = new HashMap<String, Integer>();
	static JScrollPane tableScroll = new JScrollPane();
	static TableColumnModel columnModel;
	 //String[] columns = null, selectcols = null;
	static DefaultTableModel tablemodel;
    String[][] selectdata = {};
    //static ArrayList<String[]> datarows;
    static String[] items = { "Select","Position","Chromosome","Start","End","Sample","Gene","Editable" };
   // String[] header;
    Object[] headers = new Object[]{};     
    //String[] selheaders = new String[]{"select", "select", "select", "select"};    
    ArrayList<String> header = new ArrayList<String>();
    boolean changing = false;
    ArrayList<JComboBox> comboboxes = new ArrayList<JComboBox>();
    JTableHeader tableHeader;
   	//DefaultTableModel selectmodel = new DefaultTableModel(data, selheaders); 
   	JButton open = new JButton("Open");
   	//JButton add = new JButton("Add column");
	JButton write = new JButton("Write");
   	boolean found = false, tableSet = false;
   // DefaultTableModel model = new DefaultTableModel(data, headers); 
    MouseAdapter tableMouseListener = new MouseAdapter() {
    	
			@Override
			 public void mouseEntered(java.awt.event.MouseEvent evt) {
				EditableHeader header = (EditableHeader)table.getTableHeader();
				header.removeEditor();
			 }
			@Override
			
			 public void mouseClicked(java.awt.event.MouseEvent evt) {
				
				if(evt.getClickCount() == 2) {
					
					int row = table.rowAtPoint(evt.getPoint());
				    int col = table.columnAtPoint(evt.getPoint());
				   
				    EditableHeaderTableColumn column =(EditableHeaderTableColumn) table.getColumnModel().getColumn(col);
				    
				    if(column.getHeaderValue() != null && !column.getHeaderValue().toString().equals("Select") && !column.getHeaderValue().toString().equals("Editable")) {
				    	
				    	 
				    	 if (row >= 0 && col >= 0) {

					    	
				    		 if(column.getHeaderValue().toString().equals("Position")) {
				    			 String pos = (String)table.getValueAt(row, col);
				    			 if(!pos.contains(":") && columnValues.get("Chromosome") < 0) {
				    				 Main.showError("Please, select chromosome column.", "Note", TableBrowser.tableScroll);
				    				 return;
				    			 }
			    			 	 Main.searchField.setText(pos);
				    			 Main.search(pos);
				    		 }	
				    		 else if(column.getHeaderValue().toString().equals("Chromosome")) {
				    			 String pos = (String)table.getValueAt(row, col);
				    			 Main.chromosomeDropdown.setSelectedItem(pos.replace("chr", ""));
				    		 }
				    		 else if(column.getHeaderValue().toString().equals("Start")) {
				    			 
				    			 if(columnValues.get("Chromosome") < 0) {
				    				 Main.showError("Please, select chromosome column.", "Note", TableBrowser.tableScroll);
				    				 return;
				    			 }
				    			 if(columnValues.get("End") < 0) {
				    				 int chrom = columnValues.get("Chromosome");
				    				 String start = (String)table.getValueAt(row, col);
				    				 String search = (String)table.getValueAt(row, chrom) +":" +start;
				    				
				    				 Main.searchField.setText(search);
					    			 Main.search(search);
				    			 }
				    			 else if(columnValues.get("End") > -1) {
				    				 int chrom = columnValues.get("Chromosome");
				    				 String start = (String)table.getValueAt(row, col);
				    				 int end = columnValues.get("End");
				    				 
				    				 String search = (String)table.getValueAt(row, chrom) +":" +start +"-" +(String)table.getValueAt(row, end);
				    				
				    				 Main.searchField.setText(search);
					    			 Main.search(search);
				    			 }
				    		 }
				    		 else if(column.getHeaderValue().toString().equals("Gene")) {
				    			 String gene = (String)table.getValueAt(row, col);
				    			 Main.searchField.setText(gene);
				    			 Main.search(gene);
				    		 }
				    		 if(columnValues.get("Sample") > -1) {
						    		
					    		 int sample = columnValues.get("Sample");
					    		
					    		 Main.search("s " +(String)table.getValueAt(row, sample));
					    	 }
					    	 
						 }
					}
				    else {
				    	if(column.getHeaderValue() != null && column.getHeaderValue().toString().equals("Select")) {
				    		Main.showError("Please, select column type.", "Note", TableBrowser.tableScroll);
				    		return;
				    	}
				    }
				   
				}			    
			 }
			
    };
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
    void checkHeaders() {
    	 columnValues.put("Select", -1);
		 columnValues.put("Position", -1);
		 columnValues.put("Chromosome", -1);
		 columnValues.put("Start", -1);
		 columnValues.put("End", -1);
		 columnValues.put("Sample", -1);
		 columnValues.put("Gene", -1);
		 columnValues.put("Editable", -1);
		 editables.clear();
		
    	 for(int i = 0; i < table.getColumnModel().getColumnCount(); i++) {  		    	    		    	
    		 
		    	EditableHeaderTableColumn col =(EditableHeaderTableColumn) table.getColumnModel().getColumn(i);
		 	    if(col.getHeaderRenderer() == null) {
		 	    	
		 	    	editables.put(i, true);
		 	    }
		 	    else {
		 	    	
		 	    	if(col.getHeaderValue().toString().equals("Editable")) {
		 	    		
		 	    		editables.put(i, true);
		 	    	}
		 	    	else {
		 	    		columnValues.put(col.getHeaderValue().toString(), i);
		 	    	}
		 	    }
		  }
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
	
	public TableBrowser() {
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
			
			data = new Object[][]{};
			tablemodel = new DefaultTableModel(data,new Object[]{"","","",""});
			table = new JTable(tablemodel) {			    	
				private static final long serialVersionUID = 1L;
				@Override
				
				public boolean isCellEditable(int rowIndex, int columnIndex) {			
					if(editables.containsKey(columnIndex)) {
						return true;
					}
					return false;
				}
		    };			  
			    
		    table.addMouseListener(tableMouseListener);
		    columnModel = table.getColumnModel();			
		    for(int i = columnModel.getColumnCount()-1 ;i>=0; i--) {
		    	columnModel.removeColumn(columnModel.getColumn(i));
			}
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
		    header.clear();
		    int rowlimit = 100000;
		    boolean first = true;
			while((line = reader.readLine()) != null) {
				
				if(line.startsWith("#")) {
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
					Main.showError("Showing " +rowlimit +" rows.", "Note", TableBrowser.tableScroll);
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
			editables.clear();
			reader.close();
			tablemodel.addColumn("Editable");
			headercount++;
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
			    String[] items = { "Select","Position","Chromosome","Start","End","Sample","Gene", "Editable" };					   
			    ComboRenderer renderer = new ComboRenderer(items);
			    BrowserColorColumnRenderer columnRenderer = new BrowserColorColumnRenderer(new Color(150,255,150), Color.black);
			    
			    for(int i = 0; i < columnModel.getColumnCount()-1; i++) {
			    	
			    	JComboBox<String> combo = new JComboBox<String>();
			    	for (int c = 0; c < items.length; c++) {
			   	      combo.addItem(items[c]);
			   	    }
			    	
			    	combo.addActionListener(comboActionListener);
			    	
			    	EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(i);
			    	col.setCellRenderer(columnRenderer);
			 	    col.setHeaderValue(combo.getItemAt(0));
			 	    col.setHeaderRenderer(renderer);
			 	    col.setHeaderEditor(new DefaultCellEditor(combo));
			 	   
			 	    comboboxes.add(combo);
			 	   
			    }
			    int columnIndex = columnModel.getColumn(columnModel.getColumnCount()-1).getModelIndex();
			    EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(columnIndex);
			    col.setHeaderValue("Editable");
			   
		    	col.setCellRenderer(columnRenderer);
			    checkTable(table);
			    checkHeaders();
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
	
		changing = true;
		if(table.getValueAt(0, 0).toString().equals("#Sample")) {			
			
			for(int i = 0; i < columnModel.getColumnCount()-1; i++) {	    
				EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(i);
		    	if(table.getValueAt(0, i).toString().contains("#Sample")) {		    		
		    		comboboxes.get(i).setSelectedItem("Sample");
		    		col.setHeaderValue(comboboxes.get(i).getSelectedItem().toString());
		    	}
		    	else if(table.getValueAt(0, i).toString().contains("Gene")) {		    		
		    		comboboxes.get(i).setSelectedItem("Gene");
		    		col.setHeaderValue(comboboxes.get(i).getSelectedItem().toString());	    		
		    	}
		    	else if(table.getValueAt(0, i).toString().contains("Position")) {		    		
		    		comboboxes.get(i).setSelectedItem("Position");
		    		col.setHeaderValue(comboboxes.get(i).getSelectedItem().toString());	    		
		    	}		    	
		    	
		    }
		}
		changing = false;
		//checkHeaders();
	}
	void createTable() {
		try {
		
	    tablemodel = new DefaultTableModel(data,new Object[]{"","","","","Editable"});
	    table = new JTable(tablemodel) {
	    	
			private static final long serialVersionUID = 1L;
			

			@Override
			public boolean isCellEditable(int rowIndex, int columnIndex) {
				
				if(editables.containsKey(columnIndex)) {
					return true;
				}
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
	  
	    frame.addMouseListener(this);
	    BrowserColorColumnRenderer columnRenderer = new BrowserColorColumnRenderer(new Color(150,255,150), Color.black);
	    ComboRenderer renderer = new ComboRenderer(items);
	    for(int i = 0; i < columnModel.getColumnCount()-1; i++) {
	    	
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
	    checkHeaders();
		panel.setBackground(Draw.sidecolor);		
		table.addMouseListener(tableMouseListener); 
		
	columnValues.put("Select", -1);
	columnValues.put("Position", -1);
	columnValues.put("Chromosome", -1);
	columnValues.put("Start", -1);
	columnValues.put("End", -1);
	columnValues.put("Sample", -1);
	columnValues.put("Gene", -1);
	columnValues.put("Editable", -1);
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
	//c.gridx++;;
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
	public static void main(String[] args) {
		TableBrowser bro = new TableBrowser();
		
		bro.frame.setVisible(true);
		bro.frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

	}
	
	@Override
	public void actionPerformed(ActionEvent event ) {
		if(event.getSource() == write) {
			try {	  	    
				File savefile = null;
				FileDialog fs = new FileDialog(frame, "Save table as...", FileDialog.SAVE);
		  		 
		  		  fs.setFile("*.tsv");
		  		  fs.setVisible(true);	         
		        	 
		        	 
		        	 while(true) {
		        		String filename = fs.getFile();
			        
					    if(filename != null) {   		    	  
				    	   savefile = new File(fs.getDirectory() +"/" +filename);
				    	 
				    	 
					       if(!Files.isWritable(Paths.get(savefile.getParent()))) {
					    	  Main.showError("No permission to write.", "Error");
					    	  continue;
					       }
					  
					       if(!savefile.getName().contains(".")) {
					    	   savefile = new File(savefile.getAbsolutePath() +".tsv");
					       }				     
					       Main.savedir = fs.getDirectory();
				           Main.writeToConfig("DefaultSaveDir=" +Main.savedir);
					       BufferedWriter writer = new BufferedWriter(new FileWriter(savefile));
					       
					       for(int i = 0; i<header.size()-1; i++) {
					    	   writer.write(header +"\n");
					       }
					       int columns = table.getColumnCount();
					       
					       for(int i=0;i<table.getRowCount(); i++) {
					    	   for(int j = 0; j<columns; j++) {
					    		   if(j > 0) {
					    			   writer.write("\t");
					    		   }
					    		  
					    		   if(table.getValueAt(i, j) == null) {
					    			   writer.write("\t");
					    		   }
					    		   else {
					    			   writer.write(""+table.getValueAt(i, j));
					    		   }				    		  
					    	   }
					    	   writer.write("\n");
					       }
					       
					       writer.close();
						   break;
					      }
				    	 else {					    	 
				    		 break;					    	 
				    	 }		        	 
		         }
		         if(1==1) {
		        	 return;
		         }
		    	  JFileChooser chooser = new JFileChooser();
		    	  chooser.setDialogTitle("Save table as...");
		          int returnVal = chooser.showSaveDialog((Component)this.getParent());	           	  
			       
			      if(returnVal == JFileChooser.APPROVE_OPTION) {   		    	  
				       File outfile = chooser.getSelectedFile();	
				       if(!outfile.getName().contains(".")) {				    	   
				    	   outfile = new File(outfile.getCanonicalPath() +".tsv");				    	   
				       }
				       BufferedWriter writer = new BufferedWriter(new FileWriter(outfile));
				       
				       for(int i = 0; i<header.size()-1; i++) {
				    	   writer.write(header +"\n");
				       }
				       int columns = table.getColumnCount();
				       
				       for(int i=0;i<table.getRowCount(); i++) {
				    	   for(int j = 0; j<columns; j++) {
				    		   if(j > 0) {
				    			   writer.write("\t");
				    		   }
				    		  
				    		   if(table.getValueAt(i, j) == null) {
				    			   writer.write("\t");
				    		   }
				    		   else {
				    			   writer.write(""+table.getValueAt(i, j));
				    		   }				    		  
				    	   }
				    	   writer.write("\n");
				       }
				       
				       writer.close();
			      }
				}
				catch(Exception ex) {
					ex.printStackTrace();
				}
		}
		else if(event.getSource() == open) {
			FileDialog fs = new FileDialog(frame, "Choose a tab separated file", FileDialog.LOAD);
	  		  fs.setDirectory(Main.savedir);	  		
	  		  fs.setVisible(true);
	  		  String filename = fs.getFile();			
	       
	         if (filename != null) {
	        	 File addfile = new File(fs.getDirectory() +"/" +filename);
	        	 Main.savedir = fs.getDirectory();
	        	 Main.writeToConfig("DefaultSaveDir=" +Main.savedir);
	        	 if(addfile.exists()) {
	        		     	
		        	 createTable(addfile);
		        	 info.setText("Select correct header values for BED columns. File: " +addfile.getName());
	        	 }
	        	 else {
	        		 Main.showError("File does not exists.", "Error", frame);
	        	 }
	        	 
	         }
	         if(1==1) {
	        	 return;
	         }
			  JFileChooser chooser = new JFileChooser("");	 
			  
	    	  chooser.setMultiSelectionEnabled(false);
	    	  //chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
	    	  chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
	    	  chooser.setAcceptAllFileFilterUsed(false);	    	  
	    	  
	    	 
	    	  chooser.setDialogTitle("Open file");
	    	  //chooser.setPreferredSize(new Dimension((int)Main.screenSize.getWidth()/3, (int)Main.screenSize.getHeight()/3));
	          int returnVal = chooser.showOpenDialog((Component)this.getParent());	         
	        
	         if (returnVal == JFileChooser.APPROVE_OPTION) {	    
	        	  File file = chooser.getSelectedFile();
	        	
	        	  createTable(file);
	        	  info.setText("Select correct header values for columns.");
	         }
		
		}
		//else if(event.getSource() == add) {
			//TableColumnModel columnModel = table.getColumnModel();
			//EditableHeaderTableColumn col = new EditableHeaderTableColumn(true);
		/*	
			tablemodel.addColumn("NewColumn");			
			
			columnModel = table.getColumnModel();
		    table.setTableHeader(new EditableHeader(columnModel));
		    ComboRenderer renderer = new ComboRenderer(items);		  
		    BrowserColorColumnRenderer columnRenderer = new BrowserColorColumnRenderer(new Color(150,255,150), Color.black);
		    for(int i = 0; i<comboboxes.size(); i++) {
		    	EditableHeaderTableColumn col =(EditableHeaderTableColumn) columnModel.getColumn(i);
		 	    col.setHeaderValue(comboboxes.get(i).getSelectedItem());
		 	    col.setHeaderRenderer(renderer);
		 	    col.setHeaderEditor(new DefaultCellEditor(comboboxes.get(i)));
		 	    col.setCellRenderer(columnRenderer);
		    }
		    
		    checkHeaders();
			table.revalidate();
			*/
	//	}
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
class BrowserColorColumnRenderer extends DefaultTableCellRenderer
{
   /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
Color bkgndColor, fgndColor;
     
   public BrowserColorColumnRenderer(Color bkgnd, Color foregnd) {
      super();
      bkgndColor = bkgnd;
      fgndColor = foregnd;
   }
     
   public Component getTableCellRendererComponent
        (JTable table, Object value, boolean isSelected,boolean hasFocus, int row, int column)
   {
      Component cell = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
      EditableHeaderTableColumn col =(EditableHeaderTableColumn) table.getColumnModel().getColumn(column);	 	 
      if(!isSelected) {
	      if(!col.getHeaderValue().equals("Select")) {
	    	  cell.setBackground( bkgndColor );
	    	  
	      }
	      else {
	    	  cell.setBackground( Color.white );
	    	  
	      }
      }
      return cell;
   }
}
class ComboRenderer extends JComboBox implements TableCellRenderer {
	
    ComboRenderer(String[] items) {
      for (int i = 0; i < items.length; i++) {
        addItem(items[i]);
      }
    }

    public Component getTableCellRendererComponent(JTable table,
        Object value, boolean isSelected, boolean hasFocus, int row,
        int column) {
      setSelectedItem(value);
      return this;
    }
  }