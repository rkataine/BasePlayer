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
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;

import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JTree;
import javax.swing.SwingWorker;
import javax.swing.UIManager;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;



import javax.swing.table.DefaultTableModel;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;



public class AddGenome  extends JPanel implements ActionListener, MouseListener {

	private static final long serialVersionUID = 1L;
	static JFrame frame = null;
	static boolean annotation = false;
	static Color green = new Color(150,255,150);
	JLabel genomeFileText;
	JLabel annotationFileText;
	static JPopupMenu menu = new JPopupMenu();
	static JTextArea area = new JTextArea();
	static JScrollPane menuscroll = new JScrollPane();
	static JTextField genomeName;	
	static JLabel sizeError = new JLabel("Not enough space on storage.");
	JButton openRef, openAnno, add, checkUpdates;
	static JButton download, remove, getLinks, checkEnsembl;
	static int longestName = 0;
	static String userDir = new File(Main.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParent().replace("%20", " ");
	File genomeFile, annotationFile;
	static boolean ensemblfetch = false;
	static boolean downloading = false;
	static HashMap<String, URL[]> genomeHash = new HashMap<String, URL[]>();
	static ArrayList<String> removables = new ArrayList<String>();
	static Object[] headers = {"Ensembl genomes", "Size (MB)"}, remheaders = {"Installed genomes"}; 
    static String[][] data = {}, remdata = {}; 
    static HashMap<String, Integer[]> sizeHash = new HashMap<String, Integer[]>();
    static JPanel panel = new JPanel(new GridBagLayout());
    static DefaultMutableTreeNode root = new DefaultMutableTreeNode("Genomes");
 //   static JLabel welcomeLabel = new JLabel("Welcome to BasePlayer! \nPlease use the dropdown below to download reference genome and gene annotation from Ensembl.");
    static JScrollPane treescroll;
    static JTextField genomedirectory = new JTextField();
    static JTree tree;
    
    static DefaultMutableTreeNode selectedNode = null;
    static DefaultTreeModel treemodel;
    static ArrayList<String> organisms = new ArrayList<String>();
    static JScrollPane remscroll;
    static DefaultTableModel model = new DefaultTableModel(data, headers) {
	
		private static final long serialVersionUID = 1L;
		@Override		
	    public boolean isCellEditable(int row, int column) {	      
	       return false;
	    }
	};
	static DefaultTableModel remmodel = new DefaultTableModel(remdata, remheaders) {
		 
		private static final long serialVersionUID = 1L;
		@Override		
	    public boolean isCellEditable(int row, int column) {	      
	       return false;
	    }
	};
	static JTable genometable = new JTable(model);
	static JTable remtable = new JTable(remmodel);
	static ArrayList<String[]> news = new ArrayList<String[]>();
	
	static void checkGenomes() {
		
		//DefaultMutableTreeNode
		File checkdir = Main.genomeDir, checkAnnodir;
		File[] addDir, annodir;
		root.removeAllChildren();
		int counter = 0;
		for(int i = model.getRowCount() - 1; i >=0; i--) {			
		   model.removeRow(i); 
		}
		for(int i = remmodel.getRowCount() - 1; i >=0; i--) {
			remmodel.removeRow(i); 
		}	
		
		removables.clear();
		
		int currentlen = 0, length = 0;
		if(checkdir == null) {
			try {
				checkdir = new File(userDir +"/genomes");
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
  		addDir = checkdir.listFiles();
  		for(int f = 0; f<addDir.length; f++) { 
  		
  			try {
  				if(!addDir[f].isDirectory()) {
  					continue;
  				}
	  			currentlen = genometable.getFontMetrics(genometable.getFont()).stringWidth(addDir[f].getName());
				if(currentlen > length) {
					length = currentlen;				
				}	
				
				AddGenome.removables.add(addDir[f].getName());
	  			checkAnnodir = new File(checkdir +"/" +addDir[f].getName() +"/annotation/");
	  			annodir = checkAnnodir.listFiles();	  			
	  			DefaultMutableTreeNode genome = new DefaultMutableTreeNode(addDir[f].getName());
	  			root.add(genome);
	  			if(annodir == null) {
	  				counter++;
	  				genome.add(new DefaultMutableTreeNode("Add new annotation..."));
	  			}
	  			else {
	  				counter+=annodir.length+3;
		  			for(int a = 0; a<annodir.length; a++) {
		  				currentlen = genometable.getFontMetrics(genometable.getFont()).stringWidth(annodir[a].getName());
						if(currentlen > length) {
							length = currentlen;				
						}	
						genome.add(new DefaultMutableTreeNode(annodir[a].getName()));
		  			}
		  			genome.add(new DefaultMutableTreeNode("Add new annotation..."));
	  			}	  			
  			}
  			catch(Exception e) {
  				e.printStackTrace();
  			}
  		}	  		
  		if(counter == 0) {
  			counter = 1;
  			length = genometable.getFontMetrics(genometable.getFont()).stringWidth("Add new annotation...");
  		}
  		counter++;
  		root.add(new DefaultMutableTreeNode("Add new reference..."));  	
  		genometable.setAutoResizeMode(JTable.AUTO_RESIZE_LAST_COLUMN);
		for(int i = 0 ; i< organisms.size(); i++) {
			if(!AddGenome.removables.contains(organisms.get(i))) {
				Object[] row = {organisms.get(i), ""+sizeHash.get(organisms.get(i))[0]/1048576};			
				model.addRow(row);			
				currentlen = genometable.getFontMetrics(genometable.getFont()).stringWidth(organisms.get(i));
				if(currentlen > length) {
					length = currentlen;				
				}	
			}
		}
		
		
		AddGenome.longestName = length;  	
		
		if(genometable.getRowCount() > 15) {
			genometable.setPreferredScrollableViewportSize((new Dimension(AddGenome.longestName+20+Main.defaultFontSize*8,genometable.getRowHeight()*15)));
			genometable.setMinimumSize(new Dimension(AddGenome.longestName+20,genometable.getRowHeight()*15));
		}
		else {			
			genometable.setPreferredScrollableViewportSize((new Dimension(AddGenome.longestName+20+Main.defaultFontSize*8,genometable.getRowHeight()*(genometable.getRowCount()+1))));
			genometable.setMinimumSize(new Dimension(AddGenome.longestName+20,genometable.getRowHeight()*(genometable.getRowCount()+1)));
		}
		
		if(remtable.getRowCount() > 15) {
			remtable.setPreferredScrollableViewportSize((new Dimension(AddGenome.longestName+20+Main.defaultFontSize*8,remtable.getRowHeight()*15)));
		}
		else {			
			remtable.setPreferredScrollableViewportSize((new Dimension(AddGenome.longestName+20+Main.defaultFontSize*8,remtable.getRowHeight()*(remtable.getRowCount()+1))));
		}
		//Main.defaultFontSize = 12;
		genometable.getColumnModel().getColumn(0).setPreferredWidth(AddGenome.longestName+10);
		genometable.getColumnModel().getColumn(0).setMinWidth(AddGenome.longestName+10);
	//	genometable.getColumnModel().getColumn(1).setPreferredWidth(Main.defaultFontSize*8);
		DefaultTreeModel model = (DefaultTreeModel) tree.getModel();
		model.reload();
		int rowheight = tree.getRowHeight();
		if(rowheight < 1) {
			rowheight=Main.defaultFontSize+4;
		}
	
		treescroll.setPreferredSize(new Dimension(AddGenome.longestName+20+Main.defaultFontSize*8, counter*rowheight));
  		treescroll.setMinimumSize(new Dimension(AddGenome.longestName+20+Main.defaultFontSize*8, counter*rowheight));
  		remscroll.setPreferredSize(new Dimension(AddGenome.longestName+20+Main.defaultFontSize*8, counter*rowheight));
  		remscroll.setMinimumSize(new Dimension(AddGenome.longestName+20+Main.defaultFontSize*8, counter*rowheight));
		tree.expandPath(new TreePath(root)); 
		
		setFonts(Main.menuFont);
	}
	
	public static void setFonts(Font menuFont) {
		if(menuFont == null) {
			menuFont = new Font("SansSerif", Font.BOLD, Main.defaultFontSize);
		}
		for(int i = 0 ; i<panel.getComponentCount(); i++) {
			panel.getComponent(i).setFont(menuFont);			
		}
		AddGenome.genometable.getTableHeader().setFont(menuFont);
		AddGenome.genometable.setFont(menuFont);
		AddGenome.genometable.setRowHeight(menuFont.getSize()+2);		
		tree.setFont(menuFont);
		tree.setRowHeight(menuFont.getSize()+2);
		frame.pack();
		
	}
	
	public AddGenome() {
		super(new BorderLayout());	
		
		makeGenomes();		
		tree = new JTree(root);
		tree.getSelectionModel().setSelectionMode(TreeSelectionModel.DISCONTIGUOUS_TREE_SELECTION);
		sizeError.setForeground(Draw.redColor);
		sizeError.setVisible(true);
		treemodel = (DefaultTreeModel) tree.getModel();
		remscroll = new JScrollPane(remtable);
		tree.setCellRenderer(new DefaultTreeCellRenderer() {
			private static final long serialVersionUID = 1L;
			private Icon collapsedIcon = UIManager.getIcon("Tree.collapsedIcon");
			private Icon expandedIcon = UIManager.getIcon("Tree.expandedIcon");
		//	private Icon leafIcon = UIManager.getIcon("Tree.leafIcon");
			private Icon addIcon = UIManager.getIcon("Tree.closedIcon");
			
    //        private Icon saveIcon = UIManager.getIcon("OptionPane.informationIcon");
            @Override
            public Component getTreeCellRendererComponent(JTree tree, Object value, boolean selected, boolean expanded,boolean isLeaf, int row, boolean focused) {
                Component c = super.getTreeCellRendererComponent(tree, value,selected, expanded, isLeaf, row, focused);
              
                if(!isLeaf) {
                	
                	
	                	//setFont(getFont().deriveFont(Font.PLAIN));
		                if(expanded) {
		                	setIcon(expandedIcon);
		                }
		                else {                	
		                	setIcon(collapsedIcon);
		                }
                	
	             /*   if(((DefaultMutableTreeNode) value).getUserObject().toString().equals("Annotations")) {
                		this.setFocusable(false);
                		setFont(getFont().deriveFont(Font.BOLD));
                		setIcon(null);
                	}
                */  
                }
                else {
                	if(((DefaultMutableTreeNode) value).getUserObject().toString().equals("Annotations")) {
                		
                //		setFont(getFont().deriveFont(Font.PLAIN));
                		setIcon(null);
                	}
                	else if( ((DefaultMutableTreeNode) value).getUserObject().toString().startsWith("Add new")) {
                     	
                 //       setFont(getFont().deriveFont(Font.PLAIN));
                		
                		setIcon(addIcon);
                	}
                	else {
                //		setFont(getFont().deriveFont(Font.ITALIC));
                		setIcon(null);
                	//	setIcon(leafIcon);
                	}
                	
                }
                
                return c;
            }
        });
		tree.addMouseListener(this);
		tree.addTreeSelectionListener(new TreeSelectionListener() {	

		public void valueChanged(TreeSelectionEvent e) {
			try {
			    DefaultMutableTreeNode node = (DefaultMutableTreeNode) tree.getLastSelectedPathComponent();
			    
			    if (node == null)		   
			    	return;
			  
			    selectedNode = node;
			    
			    if(node.isLeaf()) {			    	
			    	checkUpdates.setEnabled(false);
			    }
			    else {
			    	checkUpdates.setEnabled(true);			    	
			    }
			    if(node.toString().startsWith("Add new") || node.toString().equals("Annotations")) {
			    	remove.setEnabled(false);
			    }
			    else {
			    	remove.setEnabled(true);
			    }
			    genometable.clearSelection();
				download.setEnabled(false);				
			}
			catch(Exception ex) {
				ex.printStackTrace();
			}
			}
		});
	    tree.setToggleClickCount(1);
	    tree.setRootVisible(false);		
	    treescroll = new JScrollPane(tree);
		checkGenomes();
		genomeFileText = new JLabel("Select reference fasta-file");
		annotationFileText = new JLabel("Select annotation gff3-file");
		genomeName = new JTextField("Give name of the genome");	
		openRef = new JButton("Browse");	
		openAnno = new JButton("Browse");
		add = new JButton("Add");
		download = new JButton("Download");
		checkEnsembl = new JButton("Ensembl fetch");
		checkEnsembl.setMinimumSize(Main.buttonDimension);
		checkEnsembl.addActionListener(this);
		getLinks = new JButton("Get file links.");
		remove = new JButton("Remove");
		checkUpdates = new JButton("Check updates");
		download.setEnabled(false);
		getLinks.setEnabled(false);
		getLinks.addActionListener(this);
		remove.setEnabled(false);
		download.addActionListener(this);
		remove.addActionListener(this);
		panel.setBackground(Draw.sidecolor);
		checkUpdates.addActionListener(this);
		this.setBackground(Draw.sidecolor);
		frame.getContentPane().setBackground(Draw.sidecolor);
		GridBagConstraints c = new GridBagConstraints();
		
		c.gridx = 0;
		c.gridy = 0;
		c.insets = new Insets(2,4,2,4);		
		c.gridwidth = 2;
		genometable.setSelectionMode(0);
		genometable.setShowGrid(false);		
		remtable.setSelectionMode(0);
		remtable.setShowGrid(false);	
		JScrollPane scroll = new JScrollPane();		
		scroll.getViewport().setBackground(Color.white);
		scroll.getViewport().add(genometable);
		
		
		remscroll.getViewport().setBackground(Color.white);
		
		genometable.addMouseListener(this);
		remtable.addMouseListener(this);
	//	panel.add(welcomeLabel,c);
	//	c.gridy++;
		c.anchor = GridBagConstraints.NORTHWEST;	
		panel.add(new JLabel("Download genome reference and annotation"), c);
		c.gridx++;
		c.anchor = GridBagConstraints.NORTHEAST;	
		panel.add(checkEnsembl, c);
		c.anchor = GridBagConstraints.NORTHWEST;	
		c.weightx = 1;
		c.weighty = 1;
		c.fill = GridBagConstraints.BOTH;
		c.gridx = 0;
		c.gridy++;		
		//c.fill = GridBagConstraints.NONE;
		panel.add(scroll,c);		
		c.gridy++;
		c.fill = GridBagConstraints.NONE;
		panel.add(download,c);
		c.gridx = 1;
		panel.add(sizeError, c);
		c.gridx = 1;
		panel.add(getLinks, c);
		c.gridy++;
		c.gridx = 0;
		c.fill = GridBagConstraints.BOTH;
		panel.add(new JLabel("Add/Remove installed genomes manually"),c);
		c.gridy++;
		panel.add(treescroll,c);
		c.gridy++;
		
		c.fill = GridBagConstraints.NONE;
		c.gridwidth = 1;
		remove.setMinimumSize(Main.buttonDimension);
		panel.add(remove, c);
		c.gridx = 1;
		panel.add(checkUpdates, c);
		checkUpdates.setMinimumSize(Main.buttonDimension);
		checkUpdates.setEnabled(false);
		c.gridwidth = 2;
		c.gridx = 0;
		c.gridy++;
		try {
			if(Main.genomeDir != null) {
				genomedirectory.setText(Main.genomeDir.getCanonicalPath());
			}
			
			genomedirectory.setEditable(false);
			genomedirectory.setBackground(Color.white);
			genomedirectory.setForeground(Color.black);
		} 
		catch (IOException e1) {			
			e1.printStackTrace();
		}
		panel.add(new JLabel("Genome directory:"), c);
		c.gridy++;
		panel.add(genomedirectory, c);
	/*	c.fill = GridBagConstraints.BOTH;
		c.gridy++;
		panel.add(new JLabel("Add genome manually"),c);
		c.gridy++;	
		c.gridwidth = 2;		
		panel.add(new JSeparator(),c);
		c.gridwidth = 1;
		c.gridy++;
		panel.add(genomeFileText, c);
		c.fill = GridBagConstraints.NONE;
		c.gridx = 1;
		panel.add(openRef, c);		
		
		c.gridx = 0;
		openRef.addActionListener(this);
		c.gridy++;
		panel.add(annotationFileText,c);
		c.gridx=1;
		panel.add(openAnno, c);
		c.gridy++;
		
		panel.add(add,c);
				
		openAnno.addActionListener(this);
		add.addActionListener(this);
		add.setEnabled(false);
		*/
		add(panel,BorderLayout.NORTH);
		if(Main.drawCanvas != null) {
			setFonts(Main.menuFont);
		}
	/*	html.append("<a href=http:Homo_sapiens_GRCh37:Ensembl_genes> Homo sapiens GRCh37 with Ensembl</a> or <a href=http:Homo_sapiens_GRCh37:RefSeq_genes>RefSeq</a> gene annotations<br>");
		html.append("<a href=http:Homo_sapiens_GRCh38:Ensembl_genes> Homo sapiens GRCh38 with Ensembl</a> or <a href=http:Homo_sapiens_GRCh38:RefSeq_genes>RefSeq</a> gene annotations<br><br>");
		
		html.append("<a href=http:Mus_musculus_GRCm38:Ensembl_genes> Mus musculus GRCm38 with Ensembl</a> or <a href=http:Mus_musculus_GRCm38:RefSeq_genes>RefSeq</a> gene annotations<br>");
		html.append("<a href=http:Rattus_norvegicus:Ensembl_genes> Rattus norvegicus with Ensembl gene annotations</a><br>");
		html.append("<a href=http:Saccharomyces_cerevisiae:Ensembl_genes> Saccharomyces cerevisiae with Ensembl gene annotation</a><br>");					
		html.append("<a href=http:Ciona_intestinalis:Ensembl_genes> Ciona intestinalis with Ensembl gene annotation</a><br>");
		Object[] row = {"Homo_sapiens_GRCh37"};
		Object[] row = {"Homo_sapiens_GRCh38"};
		
		model.addRow(row);
		/*	genomeName.setPreferredSize(new Dimension(300,20));
		this.add(genomeName);
		this.add(new JSeparator());
		
		this.add(openRef);
		openRef.addActionListener(this);
		this.add(genomeFileText);
		this.add(openAnno);
		openAnno.addActionListener(this);
		this.add(annotationFileText);
		this.add(add);
		add.addActionListener(this);
		if(annotation) {
			openRef.setVisible(false);
			genomeFileText.setVisible(false);
			genomeName.setEditable(false);
		}
		genomeFileText.setEditable(false);
		annotationFileText.setEditable(false);*/
	}
	static void updateEnsemblList() {
		try {
			
		menu = new JPopupMenu();
		area = new JTextArea();
		menuscroll = new JScrollPane();
			
		area.setFont(Main.menuFont);
		
		menu.add(menuscroll);		
		//menu.setPreferredSize(new Dimension(menu.getFontMetrics(Main.menuFont).stringWidth("0000000000000000000000000000000000000000000000000")+Main.defaultFontSize*10, (int)menu.getFontMetrics(Main.menuFont).getHeight()*4));
		menu.setPreferredSize(new Dimension(300,200));
		
			//area.setMaximumSize(new Dimension(300, 600));
		//area.setLineWrap(true);
		//area.setWrapStyleWord(true);
		//area.setPreferredSize(new Dimension(300,200));
		
		area.revalidate();
		menuscroll.getViewport().add(area);
		menu.pack();
		menu.show(AddGenome.treescroll,0, 0);			
	/*	area.addMouseListener(new MouseListener() {

			@Override
			public void mouseClicked(MouseEvent arg0) {
				// TODO Auto-generated method stub
				
			}

			@Override
			public void mouseEntered(MouseEvent arg0) {
				// TODO Auto-generated method stub
				
			}

			@Override
			public void mouseExited(MouseEvent arg0) {
				// TODO Auto-generated method stub
				
			}

			@Override
			public void mousePressed(MouseEvent arg0) {
				StringBuffer buf = new StringBuffer("");
				for(int i= 0; i<(int)(Math.random()*100); i++) {
					buf.append("O");
				}
				AddGenome.area.append(buf.toString() +"\n");
				AddGenome.area.setCaretPosition(AddGenome.area.getText().length());
				AddGenome.area.revalidate();
				
			}

			@Override
			public void mouseReleased(MouseEvent arg0) {
				// TODO Auto-generated method stub
				
			}
			
		});*/
		
				
		FTPClient f = new FTPClient();		
		news = new ArrayList<String[]>();
		area.append("Connecting to Ensembl...\n");
		//System.out.println("Connecting to Ensembl...");
		f.connect("ftp.ensembl.org");	
		f.enterLocalPassiveMode();
		f.login("anonymous", "");		
		//System.out.println("Connected.");
		area.append("Connected.\n");
		
		FTPFile[] files = f.listFiles("pub");
		String releasedir = "";
		String releasenro;
		for(int i = 0 ; i<files.length;i++) {
			
			if(files[i].isDirectory() && files[i].getName().contains("release")) {
				releasedir = files[i].getName();				
			}
		}
		
		
		files = f.listFiles("pub/"+releasedir +"/fasta/");
		releasenro = releasedir.substring(releasedir.indexOf("-")+1);
		area.append("Searching for new genomes");
		for(int i = 0 ; i<files.length;i++) {
			if(files[i].isDirectory()) {
				FTPFile[] fastafiles =f.listFiles("pub/"+releasedir +"/fasta/"+files[i].getName() +"/dna/");
				String[] urls = new String[5];
				for(int j=0; j<fastafiles.length; j++) {
					if(fastafiles[j].getName().contains(".dna.toplevel.")) {
						urls[0] = "ftp://ftp.ensembl.org/pub/" +releasedir +"/fasta/" +files[i].getName() +"/dna/" +fastafiles[j].getName();	
						
						String filePath = "/pub/" +releasedir +"/fasta/" +files[i].getName() +"/dna/" +fastafiles[j].getName();
						f.sendCommand("SIZE", filePath);
						String reply = f.getReplyString().split("\\s+")[1];
						urls[1] = reply;						
						break;
					}						
				}
				if(urls[0] == null) {
					continue;
				}
				FTPFile[] annofiles = f.listFiles("pub/"+releasedir +"/gff3/" +files[i].getName());
				for(int j=0; j<annofiles.length; j++) {
					if(annofiles[j].getName().contains("."+releasenro+".gff3.gz")) {
						urls[2] = "ftp://ftp.ensembl.org/pub/" +releasedir+"/gff3/" +files[i].getName() +"/" +annofiles[j].getName();
						String filePath = "/pub/"  +releasedir+"/gff3/" +files[i].getName() +"/" +annofiles[j].getName();
						f.sendCommand("SIZE", filePath);
						String reply = f.getReplyString().split("\\s+")[1];
						urls[3] = reply;
						break;
					}
				}
				if(urls[2] == null) {
					continue;
				}
				if(files[i].getName().contains("homo_sapiens")) {
					urls[4] = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz";
				}
				else if(files[i].getName().contains("mus_musculus")) {
					urls[4] = "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/cytoBand.txt.gz";
				}
				String name = urls[0].substring(urls[0].lastIndexOf("/")+1, urls[0].indexOf(".dna."));
				//System.out.print(urls[0]+"\t" +urls[1] +"\t" +urls[2] +"\t" +urls[3]);
				if(genomeHash.containsKey(name) || AddGenome.removables.contains(name)) {
					//System.out.println(name +" already in the list.");
					area.append(".");
				}
				else {
					area.append("\nNew genome " +name +" added.\n");
					AddGenome.area.setCaretPosition(AddGenome.area.getText().length());
					AddGenome.area.revalidate();
					//System.out.println("New reference " +name +" found.");
					organisms.add(name);
					news.add(urls);
					
					if(urls[4] != null) {
						//System.out.println(urls[0] +" " + urls[2] +" " +urls[4]);
						URL[] newurls = {new URL(urls[0]),new URL(urls[2]),new URL(urls[4])};
						genomeHash.put(name, newurls);
					}
					else {
						URL[] newurls = {new URL(urls[0]),new URL(urls[2])};
						genomeHash.put(name, newurls);
					}
					Integer[] sizes = {Integer.parseInt(urls[1]), Integer.parseInt(urls[3])};
					sizeHash.put(name, sizes);
					
				}
				/*if(urls[4] != null) {
					System.out.print("\t" +urls[4]);
				}
				System.out.println();
			*/
			}
		
		}
		
		checkGenomes();
		if(news.size() > 0) {
			
			
			try {
					//File file = new File();
					FileWriter fw = new FileWriter(Main.genomeDir.getCanonicalPath() +"/ensembl_fetched.txt");
				    BufferedWriter bw = new BufferedWriter(fw);
				    
				    for(int i=0;i<news.size(); i++) {
						for(int j=0; j<news.get(i).length; j++) {
							if(news.get(i)[j] == null) {
								break;
							}
							 if(j > 0) {
								 bw.write("\t");
							 }
							bw.write(news.get(i)[j]);
						}
						bw.write("\n");
					}
				    bw.close();
				    fw.close();
				
					  
				} catch (IOException e) {
					
				    e.printStackTrace();
				}
			
			}
		}
		catch(Exception e) {
			Main.showError(e.getMessage(), "Error");
			e.printStackTrace();
		}
		
	}
	
	static void makeGenomes() {
		try {
			FileReader freader = null;
			File file = new File(Main.genomeDir.getCanonicalPath() +"/ensembl.txt");
			if(file.exists()) {
				freader = new FileReader(file);
				
				BufferedReader reader = new BufferedReader(freader);			
				String line,name;
				String[] split;
				
				while((line = reader.readLine()) != null) {
					split = line.split("\t");				
					name = split[0].substring(split[0].lastIndexOf("/")+1, split[0].indexOf(".dna."));
					organisms.add(name);
					if(split[0].contains("GRCh38")) {
						split[0] = split[0].replace("toplevel", "primary_assembly");
					}
					if(split.length == 5) {
						URL[] urls = {new URL(split[0]),new URL(split[2]),new URL(split[4])};
						genomeHash.put(name, urls);
					}
					else {
						URL[] urls = {new URL(split[0]),new URL(split[2])};
						genomeHash.put(name, urls);
					}
					Integer[] sizes = {Integer.parseInt(split[1]), Integer.parseInt(split[3])};
					sizeHash.put(name, sizes);
				}
				
				freader.close();
				reader.close();
			}
			file = new File(Main.genomeDir.getCanonicalPath() +"/ensembl_fetched.txt");
			if(file.exists()) {
				
				freader = new FileReader(file);
			
				BufferedReader reader = new BufferedReader(freader);			
				String line,name;
				String[] split;
				while((line = reader.readLine()) != null) {
					split = line.split("\t");				
					name = split[0].substring(split[0].lastIndexOf("/")+1, split[0].indexOf(".dna."));
					organisms.add(name);
					if(split[0].contains("GRCh38")) {
						split[0] = split[0].replace("toplevel", "primary_assembly");
					}
					if(split.length == 5) {
						URL[] urls = {new URL(split[0]),new URL(split[2]),new URL(split[4])};
						genomeHash.put(name, urls);
					}
					else {
						URL[] urls = {new URL(split[0]),new URL(split[2])};
						genomeHash.put(name, urls);
					}
					Integer[] sizes = {Integer.parseInt(split[1]), Integer.parseInt(split[3])};
					sizeHash.put(name, sizes);
				}
				freader.close();
				reader.close();
			}
	/*		URL[] urls = {new URL("ftp://ftp.ensembl.org/pub/grch37/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"),
						  new URL("ftp://ftp.ensembl.org/pub/grch37/update/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.gff3.gz"),
						  new URL("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz")};
			Integer[] sizes = new Integer[2];
			sizes[0] = 900;
			sizes[1] = 40;
			sizeHash.put("Homo_sapiens_GRCh37", sizes);
			genomeHash.put("Homo_sapiens_GRCh37", urls);
			
		/*	URL[] urls2 = {new URL("ftp://ftp.ensembl.org/pub/grch37/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"),
					  new URL("ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.37.3/GFF/ref_GRCh37.p5_top_level.gff3.gz"),
					  new URL("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz") };
			sizes = new Integer[2];
			sizes[0] = 900;
			sizes[1] = 35;
			sizeHash.put("Homo_sapiens_GRCh37:RefSeq_genes", sizes);
			genomeHash.put("Homo_sapiens_GRCh37:RefSeq_genes", urls2);
			*/
		/*	URL[] urls3 = { new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/ciona_intestinalis/dna/Ciona_intestinalis.KH.dna.toplevel.fa.gz"),
						new URL("ftp://ftp.ensembl.org/pub/release-85/gff3/ciona_intestinalis/Ciona_intestinalis.KH.85.gff3.gz")
			};
			sizes = new Integer[2];
			sizes[0] = 34;
			sizes[1] = 5;
			sizeHash.put("Ciona_intestinalis", sizes);
			genomeHash.put("Ciona_intestinalis", urls3);
			
			URL[] urls4 = { new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"),
					new URL("ftp://ftp.ensembl.org/pub/release-85/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.85.gff3.gz")
			};
			sizes = new Integer[2];
			sizes[0] = 4;
			sizes[1] = 1;
			sizeHash.put("Saccharomyces_cerevisiae", sizes);
			genomeHash.put("Saccharomyces_cerevisiae", urls4);
			
			URL[] urls5 = {new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz"),
					  new URL("ftp://ftp.ensembl.org/pub/release-85/gff3/mus_musculus/Mus_musculus.GRCm38.85.gff3.gz"),
					  new URL("http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/cytoBand.txt.gz")};
			
			sizes = new Integer[2];
			sizes[0] = 801;
			sizes[1] = 58;
			sizeHash.put("Mus_musculus_GRCm38", sizes);
			genomeHash.put("Mus_musculus_GRCm38", urls5);
			
		/*	URL[] urls6 = {new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz"),
					  new URL("ftp://ftp.ncbi.nlm.nih.gov/genomes/M_musculus/GFF/ref_GRCm38.p4_top_level.gff3.gz") };
			sizes = new Integer[2];
			sizes[0] = 801;
			sizes[1] = 28;
			sizeHash.put("Mus_musculus_GRCm38:RefSeq_genes", sizes);
			genomeHash.put("Mus_musculus_GRCm38:RefSeq_genes", urls6);
			
			URL[] urls7 = {new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz"),
					  new URL("ftp://ftp.ensembl.org/pub/release-85/gff3/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.85.gff3.gz") };
			sizes = new Integer[2];
			sizes[0] = 812;
			sizes[1] = 15;
			sizeHash.put("Rattus_norvegicus", sizes);
			genomeHash.put("Rattus_norvegicus", urls7);
			
			
			URL[] urls8 = {new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
					  new URL("ftp://ftp.ensembl.org/pub/release-85/gff3/homo_sapiens/Homo_sapiens.GRCh38.85.gff3.gz"),
					  new URL("http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz")};
			sizes = new Integer[2];
			sizes[0] = 860;
			sizes[1] = 38;
			sizeHash.put("Homo_sapiens_GRCh38", sizes);
			genomeHash.put("Homo_sapiens_GRCh38", urls8);
		
			/*URL[] urls9 = {new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
					  new URL("ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/GFF/ref_GRCh38.p7_top_level.gff3.gz"),
					  new URL("http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz")};
			sizes = new Integer[2];
			sizes[0] = 860;
			sizes[1] = 41;
			sizeHash.put("Homo_sapiens_GRCh38:RefSeq_genes", sizes);
			genomeHash.put("Homo_sapiens_GRCh38:RefSeq_genes", urls9);
		*/
		}
		catch(Exception ex) {
			ex.printStackTrace();
		}
	}
	public static class EnsemblFetch extends SwingWorker<String, Object> {

		@Override
		protected String doInBackground() throws Exception {
			ensemblfetch = true;
			updateEnsemblList();
			ensemblfetch = false;
			return null;
		}
		
	}
	public static class OutputRunner extends SwingWorker<String, Object> {
		File genomefile, annotationFile;
		String genomeName, ref;
		URL annofile;
		Boolean createGenome = false, downloadAnnotation = false;
		
		public OutputRunner(String genomeName, File genomefile, File annotationFile) {
			this.genomefile = genomefile;
			this.genomeName = genomeName;
			this.annotationFile = annotationFile;			
		}
		public OutputRunner(String genomeName, File genomefile, File annotationFile, String ref) {
			this.genomefile = genomefile;
			this.genomeName = genomeName;
			this.annotationFile = annotationFile;
			this.ref = ref;
		
		}
		public OutputRunner(URL annotationFile, String ref) {
			this.annofile = annotationFile;
			this.ref = ref;
		}
		
		protected String doInBackground() {
			
			frame.setState(Frame.ICONIFIED);
			try{
				if(Main.drawCanvas != null) {
					Main.drawCanvas.ready("all");
					Main.drawCanvas.loading("Processing files...");
				}
				if(createGenome) {
					
					createGenome(genomeName, genomefile, annotationFile);
					createGenome = false;
				}
				else if(downloadAnnotation) {
					downloadAnnotation(annofile, ref);
				}
				
				
				if(Main.drawCanvas != null) {
					Main.drawCanvas.ready("Processing files...");
				}
			}
			catch(Exception e) {
				frame.setState(Frame.NORMAL);
				e.printStackTrace();
			}
			frame.setState(Frame.NORMAL);
			return "";
		}
		void createGenome(String genomeName, File genomeFile, File annotationFile) {
			
			try {				
				
			File annofile = null;
			if(Main.drawCanvas != null) {
				Main.drawCanvas.splits.get(0).getGenes().clear();
			}
			Boolean ok = false;
			String annotationUrl = "";		
			File fastatest = new File("test");
			
			if(this.ref != null) {	
				String targetDir = "";
				URL[] urlList = AddGenome.genomeHash.get(ref);
				URL fastafile= urlList[0];
				File f = Main.genomeDir;
				if(f.canWrite()) {
					targetDir = Main.genomeDir.getCanonicalPath() +"/" +ref +"/";
				} else {
				  Main.showError("no access", "Error");
				
				}
				
				fastatest = new File(targetDir +FilenameUtils.getName(fastafile.getFile()).substring(0,FilenameUtils.getName(fastafile.getFile()).indexOf(".gz") ));
				
				if(!new File(targetDir).exists()) {
					File fasta = new File(targetDir +FilenameUtils.getName(fastafile.getFile()));
				//	fastatest = new File(targetDir +FilenameUtils.getName(fastafile.getFile()).substring(0,FilenameUtils.getName(fastafile.getFile()).indexOf(".gz") ));
					File target = new File(targetDir);
					target.mkdir();
					if(!fasta.exists() && !fastatest.exists()) {
						if(Main.drawCanvas != null) {
							Main.drawCanvas.loadingtext ="Downloading " +genomeName;
						}
						else {
							System.out.println("Downloading " +genomeName);
						}
						String test = Main.downloadFile(fastafile, targetDir, AddGenome.sizeHash.get(ref)[0]);	
						if(test == null) {
							target.delete();
							return;
						}
					}
				}
				annotationUrl = annotationFile.getName();
				URL gfffile= urlList[1];
				targetDir = Main.genomeDir.getCanonicalPath() +"/" +ref +"/annotation/";
				if(Main.drawCanvas != null) {
					Main.drawCanvas.loadingtext ="Downloading " +gfffile.getFile();
				}
				String filetest = Main.downloadFile(gfffile, targetDir, AddGenome.sizeHash.get(ref)[1]);		
				if(filetest == null) {
					return;
				}
				if(!filetest.equals(FilenameUtils.getName(gfffile.getFile()))) {
					annotationUrl = filetest;
					annotationFile = new File( Main.genomeDir.getCanonicalPath() +"/" +ref +"/annotation/" +filetest +"/" +filetest);
				}				
				if(Main.drawCanvas != null) {
					Main.drawCanvas.loadingtext ="Downloading " +filetest;
				}
				else {
					System.out.println("Downloading " +filetest);
				}
				
				if(urlList.length == 3) {
					URL bandFile = urlList[2];
					targetDir = Main.genomeDir.getCanonicalPath() +"/" +ref +"/";
					Main.downloadFile(bandFile, targetDir, 0);
					File infile = new File(targetDir+FilenameUtils.getName(bandFile.getFile()));
					File outfile = new File(targetDir+"bands.txt");
					MethodLibrary.unzip(infile, outfile);
				}
				
			}
			
			if(ref == null) {
				if(!new File(Main.genomeDir.getCanonicalPath() +"/" +genomeName +"/").exists()) {
					ok = new File(Main.genomeDir.getCanonicalPath() +"/" +genomeName +"/").mkdir();
				}
				else {
					ok = true;
				}
			}
			else {
				ok = true;
			}
			
			if(ok) {
				String genomedir = Main.genomeDir.getCanonicalPath() +"/" +genomeName +"/";
				if(genomeFile != null && !fastatest.exists()) {
					if(Main.drawCanvas != null) {
						Main.drawCanvas.loadingtext = "Unpacking and indexing: " +genomeName;			
					}
					else {
						System.out.println("Unpacking and indexing: " +genomeName);
					}
							
					indexFasta(genomeFile, genomedir);
					Main.addGenomeFile(genomeName);
				}
			
				if(annotationFile != null) {
					if(Main.drawCanvas != null) {
						Main.drawCanvas.loadingtext = "Unpacking and indexing: " +annotationFile.getName();
					}
					else {
						System.out.println("Unpacking and indexing: " +annotationFile.getName());
					}
					
					if(ref == null) {
						if(annotationFile.getName().indexOf(".gff") > -1) {
							ok = new File(genomedir +"annotation/" +annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gff"))+"/").mkdirs();
						}
						else {
							ok = new File(genomedir +"annotation/" +annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gtf"))+"/").mkdirs();
						}
						
					}
					else {
						ok = true;
					}
					
					if(ref != null) {
						annofile = new File(genomedir +"annotation/" +annotationUrl +"/"+annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gff")) +".bed.gz");
						
					}
					else {
						if(annotationFile.getName().indexOf(".gff") > -1) {
							annofile = new File(genomedir +"annotation/"+annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gff")) +"/" +annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gff")) +".bed.gz");
						}
						else {
							annofile = new File(genomedir +"annotation/"+annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gtf")) +"/" +annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gtf")) +".bed.gz");
							
						}
					}
					if(ok) {
					
						if(genomeFile == null) {
							File dir = new File(Main.genomeDir.getCanonicalPath() +"/" +genomeName);
							File[] files = dir.listFiles();
							for(int i = 0 ; i < files.length; i++) {
								if(files[i].getName().endsWith(".fa") || files[i].getName().endsWith(".fasta")) {
									genomeFile = files[i].getCanonicalFile();
									break;
								}
							}
						}
						if(Main.drawCanvas != null) {
							Main.drawCanvas.loadingtext = "Unpacking and indexing: " +annotationFile.getName();			
						}
						else {
							System.out.println("Unpacking and indexing: " +annotationFile.getName());
						}
						parseGFF(annotationFile, annofile, genomeName, genomedir +genomeFile.getName());
						Main.drawCanvas.ready("all");
						Main.addAnnotationFile(genomeName, annofile);
					}
					else {
						
						Main.addAnnotationFile(genomeName, annofile);
					}
					
				}
				if(ref != null) {
					genomeFile.delete();
					annotationFile.delete();
					
				}
				Main.defaultGenome = genomeName;
				if(genomeName.length() > Main.reflength) {
					Main.reflength = genomeName.length();
					Main.refDropdown.setPopupWidth(Main.reflength*Main.letterlength);
				}
				if(annofile != null) {
					Main.defaultAnnotation = annofile.getName();	
				}
				else {
					Main.defaultAnnotation = "";
				}
				if(Main.drawCanvas != null) {
					Main.setChromDrop(genomeName);
					Main.getBands();
					Main.geneDropdown.setSelectedItem(Main.defaultAnnotation);
					Main.geneDropdown.revalidate();
					Main.getExons();
					Main.chromosomeDropdown.setSelectedIndex(0);
					AddGenome.downloading = false;
				}			
				else {
					System.out.println("Genome ready: " +genomeName +" : " +Main.defaultAnnotation);
				}
				DefaultMutableTreeNode newref = null;
				
				if(genomeFile != null) {
					
					DefaultMutableTreeNode last = (DefaultMutableTreeNode)root.getLastChild();
					root.remove(last);
					newref = new DefaultMutableTreeNode(genomeName);
					newref.add(new DefaultMutableTreeNode("Add new annotation..."));
					root.add(newref);
					root.add(last);
					treemodel.reload(root);
				}
				if(annofile != null) {
					DefaultMutableTreeNode parent = null;
					if(newref == null) {
						parent = (DefaultMutableTreeNode)selectedNode.getParent();
						DefaultMutableTreeNode last = (DefaultMutableTreeNode)parent.getLastChild();						
						parent.remove(last);						
						parent.add(new DefaultMutableTreeNode(annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gff"))));
						parent.add(last);
						treemodel.reload(parent);
						
					}
					else {
						
						parent = newref;
						checkGenomes();
					}					
				}
				
				}
				for(int i  =0; i<AddGenome.root.getChildCount(); i++) {
						
					 if(AddGenome.root.getChildAt(i).toString().equals(genomeName)) {
						 AddGenome.tree.setSelectionRow(i);
						 AddGenome.tree.expandRow(i);
						 break;
					 }
				 }
			}
			catch(Exception e) {
				try {					
					if(new File(Main.genomeDir.getCanonicalPath() +"/" +genomeName).listFiles() == null) {
						FileUtils.deleteDirectory(new File(Main.genomeDir.getCanonicalPath() +"/" +genomeName));
					}
					else {
						File[] files = new File(Main.genomeDir.getCanonicalPath() +"/" +genomeName).listFiles();
						int faifound = -1;
						for(int i = 0; i<files.length; i++) {
							if(files[i].getName().endsWith(".fai")) {
								faifound = i;
								break;
							}
						}
						if(faifound > -1) {
							if(files[faifound].length() == 0) {
								FileUtils.deleteDirectory(new File(Main.genomeDir.getCanonicalPath() +"/" +genomeName));
							}
						}
						else {
							FileUtils.deleteDirectory(new File(Main.genomeDir.getCanonicalPath() +"/" +genomeName));
						}
						files = new File(Main.genomeDir.getCanonicalPath() +"/" +genomeName+"/annotation/").listFiles();
						
						if(files != null) {
							for(int i = 0 ; i<files.length; i++) {
								if(files[i].isDirectory()) {
									if(files[i].listFiles() != null) {
										if(files[i].listFiles().length < 2) {
											FileUtils.deleteDirectory(new File(Main.genomeDir.getCanonicalPath() +"/" +genomeName+"/annotation/" +files[i].getName()));
										}
									}
								}
							}
						}
					}
					downloading = false;
				}
				catch(Exception ex) {
					ex.printStackTrace();
				}
				downloading = false;
				e.printStackTrace();
				Main.showError(e.getMessage(), "Error");
			}
		}
		
	}

static String downloadAnnotation(URL fileurl, String genome) {
	File annotationFile = null;
	try {
		
		String annotationUrl = fileurl.getFile();	
		String targetDir = Main.genomeDir.getCanonicalPath() +"/" +genome +"/annotation/";
		File dir = new File(Main.genomeDir.getCanonicalPath() +"/" +genome);
		File[] files = dir.listFiles();
		String genomestring = "";
		for(int i = 0 ; i < files.length; i++) {
			if(files[i].getName().endsWith(".fa") || files[i].getName().endsWith(".fasta")) {
				genomestring = files[i].getCanonicalPath();
				break;
			}
		}
		if(Main.drawCanvas != null) {
			Main.drawCanvas.loadingtext = "Downloading " +fileurl.getHost()+"/"+annotationUrl;
		}
		else {
			System.out.println("Downloading " +annotationUrl);
		}
		String filetest = Main.downloadFile(fileurl, targetDir, AddGenome.sizeHash.get(genome)[1]);	
		annotationFile = new File( Main.genomeDir.getCanonicalPath() +"/" +genome +"/annotation/" +filetest +"/" +filetest);
		String genomedir = Main.genomeDir.getCanonicalPath() +"/" +genome+"/";
		
		File annofile = new File(genomedir +"annotation/" +filetest +"/"+annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gff")) +".bed.gz");
		
		parseGFF(annotationFile, annofile, genome, genomestring);		
		annotationFile.delete();
		Main.addAnnotationFile(genome, annofile);
		checkGenomes();
	}
	catch(Exception e) {
		e.printStackTrace();
		annotationFile.delete();
	}
	
	return "";
}
	
static SAMSequenceDictionary ReadDict(File fastafile) {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		
		try {
			BufferedReader reader = null;
			if(fastafile.getName().endsWith(".gz")) {
				reader = new BufferedReader(new FileReader(fastafile.getCanonicalPath().substring(0,fastafile.getCanonicalPath().lastIndexOf(".gz") ) +".fai"));
			}
			else {
				reader = new BufferedReader(new FileReader(fastafile.getCanonicalPath() +".fai"));
			}
			
		
			String line;
			String[] split;
			
			while((line = reader.readLine())!=null) {				
				split = line.split("\t");
				SAMSequenceRecord record = new SAMSequenceRecord(split[0], Integer.parseInt(split[1]));
				dict.addSequence(record);				
			}
			reader.close();
		}
		catch(Exception e){
			
			e.printStackTrace();
		}
 		return dict;
	}
	
	static void parseGFF(File gffFile, File outfile, String genome, String genomeFile) {		
		
			try {
				
				SAMSequenceDictionary dict = ReadDict(new File(genomeFile));				
				if(gffFile.getName().contains(".gff")) {
					FileRead.readGFF(gffFile, outfile.getCanonicalPath(), dict);
				}
				else {
					FileRead.readGTF(gffFile, outfile.getCanonicalPath(), dict);
				}
				
			}
			catch(Exception e) {
				e.printStackTrace();
			}
			/*	Main.addGenome
				Main.genomehash.get(genome).add(outfile);	
				JMenuItem additem = new JMenuItem(outfile.getName().substring(0,outfile.getName().indexOf(".gff")));
				additem.setName(outfile.getName().substring(0,outfile.getName().indexOf(".gff")));
				additem.addMouseListener(Main.frame);
				Main.addMenu.add(additem);
				
				break;
			*/
			
		
	}
	
	static void indexFasta(File fastafile, String targetDir) {
		try {
			
			ArrayList<String[]> faiList = new ArrayList<String[]>();
			
			BufferedReader reader;
			Boolean zip = false;
			GZIPInputStream gzip = null;
			BufferedWriter fastaWriter = null;
			BufferedWriter indexWriter = null;
			if(fastafile.getName().endsWith(".gz")) {
				zip = true;
				gzip = new GZIPInputStream(new FileInputStream(fastafile));
			    reader = new BufferedReader(new InputStreamReader(gzip));
			    if(!new File (targetDir +fastafile.getName().substring(0, fastafile.getName().indexOf(".gz"))).exists()) {
			    	fastaWriter = new BufferedWriter(new FileWriter(targetDir +fastafile.getName().substring(0, fastafile.getName().indexOf(".gz"))));
			    }
			    if(!new File(targetDir+fastafile.getName().substring(0, fastafile.getName().indexOf(".gz")) +".fai").exists()) {
					indexWriter = new BufferedWriter(new FileWriter(targetDir+fastafile.getName().substring(0, fastafile.getName().indexOf(".gz")) +".fai"));
				}
			}
			else {
				reader = new BufferedReader(new FileReader(fastafile));	
				if(!new File(targetDir+"/" +fastafile.getName()).exists()) {					
					fastaWriter = new BufferedWriter(new FileWriter(targetDir +fastafile.getName()));
				}
				if(!new File(targetDir+"/" +fastafile.getName() +".fai").exists()) {
					indexWriter = new BufferedWriter(new FileWriter(targetDir+fastafile.getName() +".fai"));
				}
				
			}
			String line;
			
			long counter = 0, chromlength = 0, pointer = 0;
			String[] split;
			Main.drawCanvas.loadbarAll = 0;
			
			while((line = reader.readLine()) != null) {
				if(!Main.drawCanvas.loading) {
					break;
				}
				if(fastaWriter != null) {
					fastaWriter.write(line +"\n");
				}
				counter += line.length()+1;
				if(line.startsWith(">")) {
					if(faiList.size() > 0) {
						faiList.get(faiList.size()-1)[1] = ""+chromlength;
					}					
					//Main.drawCanvas.loadbarAll = (int)((counter/(double)filesize)*100);
					//Main.drawCanvas.loadBarSample = Main.drawCanvas.loadbarAll;
					
					String[] row = new String[5];
					faiList.add(row);
					split = line.split("\\s+");
					pointer = counter;
					row[0] = split[0].substring(1);
					
					row[2] = ""+pointer;						
					line = reader.readLine();
					if(fastaWriter != null) {
						fastaWriter.write(line +"\n");
					}
					chromlength = line.length();
					row[3] = ""+line.length();
					row[4] = ""+(line.length()+1);
					
					counter += line.length()+1;						
					
				}
				else {
					chromlength += line.length();
				}					
			}
			faiList.get(faiList.size()-1)[1] = ""+chromlength;
			if(fastaWriter != null) {
				fastaWriter.close();
			}
			reader.close();
		
			if(zip) {
				gzip.close();
			}
			if(indexWriter != null) {
				for(int i = 0 ; i<faiList.size(); i++) {
					for(int j = 0; j<4; j++) {
						indexWriter.write(faiList.get(i)[j] +"\t");
					}
					indexWriter.write(faiList.get(i)[4] +"\n");
				}
				indexWriter.close();
			
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	static void createAndShowGUI() {	
		frame = new JFrame("Add new genome");
		if(Main.drawCanvas != null) {
		 	frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE); 
		 	frame.setVisible(false); 
		}
		else {
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
		 	frame.setVisible(true); 
		}
		 	frame.setResizable(true);    
		 	frame.setTitle("Genome select");
		    JComponent newContentPane = new AddGenome();
		    frame.setContentPane(newContentPane);
		    frame.pack();
		    sizeError.setVisible(false);
	}
	public static void main(String[] args) {
		
		try {
			
			
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		 javax.swing.SwingUtilities.invokeLater(new Runnable() {
	          public void run() {	          
	          		createAndShowGUI();
	          }
	      });	
	}
	public void actionPerformed(ActionEvent event) {
		if(event.getSource() == download) {			
			if(!downloading) {
	    		downloading = true;
	    		downloadGenome(genometable.getValueAt(genometable.getSelectedRow(), 0).toString());	    		
	    		downloading = false;
	    	}	
					}
		else if(event.getSource() == getLinks) {
			URL[] urls = AddGenome.genomeHash.get(genometable.getValueAt(genometable.getSelectedRow(), 0).toString());
			JPopupMenu menu = new JPopupMenu();
			JTextArea area = new JTextArea();
			JScrollPane menuscroll = new JScrollPane();			
			area.setFont(Main.menuFont);			
			menu.add(menuscroll);			
			menu.setPreferredSize(new Dimension(menu.getFontMetrics(Main.menuFont).stringWidth(urls[0].toString())+Main.defaultFontSize*10, (int)menu.getFontMetrics(Main.menuFont).getHeight()*4));
			//area.setMaximumSize(new Dimension(300, 600));
			area.setLineWrap(true);
			area.setWrapStyleWord(true);
			for(int i=0;i<urls.length; i++) {
				area.append(urls[i].toString() +"\n");
			}	
			
			area.setCaretPosition(0);
			area.revalidate();
			menuscroll.getViewport().add(area);
			menu.pack();
			menu.show(this,0, 0);			
			
		}
		else if (event.getSource() == checkEnsembl) {
			if(ensemblfetch ) {
				menu.show(AddGenome.treescroll,0, 0);		
			}
			else {
				EnsemblFetch fetcher = new EnsemblFetch();
				fetcher.execute();
			}
		}
		else if(event.getSource() == checkUpdates) {
			URL testfile = null;
			try {
			// kattoo onko päivityksiä annotaatioon
				String ref = selectedNode.toString();
				if(AddGenome.genomeHash.get(ref) != null) {
					ArrayList<String> testfiles = new ArrayList<String>();
					if(Main.drawCanvas != null) {
						for(int i = 0 ; i<Main.genomehash.get(ref).size();i++) {
							testfiles.add(Main.genomehash.get(ref).get(i).getName().replace(".bed.gz",""));							
						}
					}
					testfile = AddGenome.genomeHash.get(ref)[1];					
					String result = Main.checkFile(testfile, testfiles);
					
					if(result.length() == 0) {						
						Main.showError("You have newest annotation file.", "Note");
					}
					else {
						int n = JOptionPane.showConfirmDialog(Main.drawCanvas, "New annotation file found: " +result +"\nDownload it now?", "Note", JOptionPane.YES_NO_OPTION);
						if(n == JOptionPane.YES_OPTION) {
							URL fileurl = new URL(testfile.getProtocol() +"://" +testfile.getHost() +testfile.getPath().substring(0,testfile.getPath().lastIndexOf("/")+1)+result);
							OutputRunner runner = new OutputRunner(fileurl, ref);
							runner.downloadAnnotation = true;
							runner.execute();						
						}						
					}				
				}
				else {
					Main.showError("This genome is not from Ensembl list, could not check for updates.", "Note", AddGenome.genometable);
				}
			}
			catch(Exception e) {
				Main.showError("Cannot connect to " +testfile.getHost() +".\nTry again later.", "Error");
				e.printStackTrace();
			}				
		}
		else if(event.getSource() == remove) {		
			if(!selectedNode.isLeaf()) {
				String removeref = selectedNode.toString();			
			//	Boolean same = false;
				try {
					if(Main.drawCanvas != null ) {
						if(removeref.equals(Main.refDropdown.getSelectedItem().toString())) {
							Main.referenceFile.close();
					//		same = true;
							if(ChromDraw.exonReader != null) {								
								ChromDraw.exonReader.close();
							}							
						}
					}
					if(Main.genomehash.containsKey(removeref)) {
						for(int i= Main.genomehash.get(removeref).size()-1; i>=0;i--) {
							Main.genomehash.get(removeref).remove(i);
						}
						Main.genomehash.remove(removeref);
						
					}
					if(Main.drawCanvas != null ) {
						Main.refModel.removeElement(removeref);
						Main.refDropdown.removeItem(removeref);
						Main.refDropdown.revalidate();
					}
					
					for(int i = 0 ; i<Main.genome.getItemCount(); i++) {
						if(Main.genome.getItem(i).getName() != null) {
							
							if(Main.genome.getItem(i).getName().equals(removeref)) {
								Main.genome.remove(Main.genome.getItem(i));							
								break;
							}
						}
					}					
									
					FileUtils.deleteDirectory(new File(Main.genomeDir.getCanonicalPath() +"/" +removeref));
					checkGenomes();
					Main.setAnnotationDrop("");
					
					if(Main.genomehash.size() == 0) {
						Main.refDropdown.setSelectedIndex(0);
						Main.setChromDrop("-1");
					}
				}
				catch(Exception e) {
					e.printStackTrace();
					try {
						Main.showError("Could not delete genome folder.\nYou can do it manually by deleting folder " +Main.genomeDir.getCanonicalPath() +"/"+removeref, "Note");
					} catch (IOException e1) {
						
						e1.printStackTrace();
					}
				}
			}
			else  {
				try {
					if(Main.drawCanvas != null ) {						
						if(ChromDraw.exonReader != null) {
							ChromDraw.exonReader.close();
						}					
					}
				
					Main.removeAnnotationFile(selectedNode.getParent().toString(), selectedNode.toString());		

				
					FileUtils.deleteDirectory(new File(Main.genomeDir.getCanonicalPath() +"/" +selectedNode.getParent().toString() +"/annotation/"+selectedNode.toString()));
					
					
				//	root.remove(selectedNode.getParent().getIndex(selectedNode));
				//	root.remove
				//	checkGenomes();
					
				}
				catch(Exception e) {
					e.printStackTrace();
					try {
						Main.showError("Could not delete genome folder.\nYou can do it manually by deleting folder " +Main.genomeDir.getCanonicalPath() +"/"+selectedNode.getParent().toString() +"/annotation/" +selectedNode.toString(), "Note");
					} catch (IOException e1) {
						
						e1.printStackTrace();
					}
				}
				treemodel.removeNodeFromParent(selectedNode);
			}
			
		}
		else if(event.getSource() == add) {
		
			if(genomeFile == null) {
				if(new File(genomeFileText.getText()).exists()) {
					genomeFile = new File(genomeFileText.getText());		
					
				}	
				else {
					genomeFileText.setText("Select reference genome fasta-file.");
					genomeFileText.setForeground(Color.red);
					return;
				}
			}
			
			
			/*if(genomeName.getText().contains("Give name") || genomeName.getText().length() == 0) {
				genomeName.setText("Give name of the genome");
				genomeName.setForeground(Color.red);
				genomeName.revalidate();
				
			}
			else if(!annotation && new File(Main.userDir +"/genomes/"+genomeName.getText().trim().replace("\\s+", "_")).exists()) {
				genomeName.setText("This genome exists already.");
				genomeName.setForeground(Color.red);
				genomeName.revalidate();
			}
			else */
		
			if((genomeFileText.getText().length() == 0 || genomeFileText.getText().startsWith("Select reference"))) {				
				genomeFileText.setText("Select reference genome fasta-file.");
				genomeFileText.setForeground(Color.red);
				genomeFileText.revalidate();
			}
			
			else {	
				
				OutputRunner runner = new OutputRunner(genomeFile.getName().replace(".fasta", "").replace(".gz", ""), genomeFile, annotationFile);
				runner.execute();
			}
			
		}	
		else if(event.getSource() == openRef) {
			try {
				
			 JFileChooser chooser = new JFileChooser(Main.downloadDir);	 
	    	  chooser.setMultiSelectionEnabled(false);
	    	  chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
	    	  chooser.setAcceptAllFileFilterUsed(false);	    	
	    	  MyFilterFasta fastaFilter = new MyFilterFasta();	    
	    	 
	    	  chooser.addChoosableFileFilter(fastaFilter); 	    	
	    	  chooser.setDialogTitle("Select reference fasta-file");
	    	  if(Main.screenSize != null) {
	    		  chooser.setPreferredSize(new Dimension((int)Main.screenSize.getWidth()/3, (int)Main.screenSize.getHeight()/3));
	    	  }
	   	  
	          int returnVal = chooser.showOpenDialog((Component)this.getParent());	         
	          
	         if (returnVal == JFileChooser.APPROVE_OPTION) {	        	 
	        	 genomeFile = chooser.getSelectedFile(); 
	        	 Main.downloadDir = genomeFile.getParent();
	        	 Main.writeToConfig("DownloadDir=" +genomeFile.getParent());
	        	 genomeFileText.setText(genomeFile.getName());
	        	 genomeFileText.revalidate();  	
	        	 frame.pack();
	          }
		 }
		 catch(Exception ex) {
			 ex.printStackTrace();
		 }		
		}
		else if(event.getSource() == openAnno) {
			try {
				
				 JFileChooser chooser = new JFileChooser(Main.downloadDir);	 
		    	  chooser.setMultiSelectionEnabled(false);
		    	  chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
		    	  chooser.setAcceptAllFileFilterUsed(false);	    	
		    	  MyFilterGFF gffFilter = new MyFilterGFF();	    
		    	 
		    	  chooser.addChoosableFileFilter(gffFilter); 	    	
		    	  chooser.setDialogTitle("Select annotation gff3-file");
		    	  if(Main.screenSize != null) {
		    		  chooser.setPreferredSize(new Dimension((int)Main.screenSize.getWidth()/3, (int)Main.screenSize.getHeight()/3));
		    	  }
		          int returnVal = chooser.showOpenDialog((Component)this.getParent());	         
		          
		         if (returnVal == JFileChooser.APPROVE_OPTION) {
		        	 if(genomeFile == null) {
		        		 genomeFile = Main.fastahash.get(Main.hoverGenome);
		        	 }
		        	 annotationFile = chooser.getSelectedFile(); 
		        	 Main.downloadDir = annotationFile.getParent();
		        	 Main.writeToConfig("DownloadDir=" +annotationFile.getParent());
		        	
		        	
		        	 
		        	OutputRunner runner = new OutputRunner(genomeFile.getName().replace(".fasta", "").replace(".gz", ""), genomeFile, annotationFile);
					runner.execute();
		          }
			 }
			 catch(Exception ex) {
				 ex.printStackTrace();
			 }		
		}
	}
	
	static class MyFilterFasta extends javax.swing.filechooser.FileFilter {
		public boolean accept(File file) { 
				if (file.isDirectory()) {
					return true;
			    }	
				if (file.getName().endsWith(".fa")) {
					return true;
				}					
				if (file.getName().endsWith(".fa.gz")) {
					return true;
				}
				if (file.getName().endsWith(".fasta")) {
					return true;
				}
				if (file.getName().endsWith(".fasta.gz")) {
					return true;
				}
				else {
					return false;
				}	
		} 
		
		public String getDescription() { return "*.fa, *.fasta"; }
}
	static class MyFilterGFF extends javax.swing.filechooser.FileFilter {
		public boolean accept(File file) { 
				if (file.isDirectory()) {
					return true;
			    }	
				if (file.getName().endsWith(".gff3")) {
					return true;
				}
				if(file.getName().endsWith(".gff3.gz")) {
					return true;
				}
				if(file.getName().endsWith(".gtf.gz")) {
					return true;
				}
				if(file.getName().endsWith(".gtf")) {
					return true;
				}
				else {
					return false;
				}	
		} 
		
		public String getDescription() { return "*.gff3, *.gff3.gz"; }
}

	@Override
	public void mouseClicked(MouseEvent event) {
		
		
	}


	@Override
	public void mousePressed(MouseEvent e) {
		
		if(e.getSource() == tree) {
		
			if(selectedNode!= null && selectedNode.toString().contains("Add new refe")) {
				try {
					FileDialog fs = new FileDialog(frame, "Select reference fasta-file", FileDialog.LOAD);
			  		  fs.setDirectory(Main.downloadDir);	  		
			  		  fs.setVisible(true);
			  		  String filename = fs.getFile();			
			          fs.setFile("*.fasta;*.fa");
			          fs.setFilenameFilter(new FilenameFilter() {
				  			public boolean accept(File dir, String name) {
						        return name.toLowerCase().contains(".fasta") || name.toLowerCase().contains(".fa");
						     }
						 });
			         if (filename != null) {
			        	 File addfile = new File(fs.getDirectory() +"/" +filename);
			        	
			        	 if(addfile.exists()) {
			        		     	
			        		 genomeFile = addfile; 
				        	 Main.downloadDir = genomeFile.getParent();
				        	 Main.writeToConfig("DownloadDir=" +genomeFile.getParent());
				        	 OutputRunner runner = new OutputRunner(genomeFile.getName().replace(".fasta", "").replace(".fa", "").replace(".gz", ""), genomeFile, null);
				        	 runner.createGenome = true;
				        	 runner.execute();
			        	 }
			        	 else {
			        		 Main.showError("File does not exists.", "Error", frame);
			        	 }
			        	 
			         }
			         if(1==1) {
			        	 return;
			         }
					 JFileChooser chooser = new JFileChooser(Main.downloadDir);	 
			    	  chooser.setMultiSelectionEnabled(false);
			    	  chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
			    	  chooser.setAcceptAllFileFilterUsed(false);	    	
			    	  MyFilterFasta fastaFilter = new MyFilterFasta();	    
			    	 
			    	  chooser.addChoosableFileFilter(fastaFilter); 	    	
			    	  chooser.setDialogTitle("Select reference fasta-file");
			    	  if(Main.screenSize != null) {
			    		  chooser.setPreferredSize(new Dimension((int)Main.screenSize.getWidth()/3, (int)Main.screenSize.getHeight()/3));
			    	  }
			   	  
			          int returnVal = chooser.showOpenDialog((Component)this.getParent());	         
			          
			         if (returnVal == JFileChooser.APPROVE_OPTION) {
			        	 
			        	 genomeFile = chooser.getSelectedFile(); 
			        	 Main.downloadDir = genomeFile.getParent();
			        	 Main.writeToConfig("DownloadDir=" +genomeFile.getParent());
			        	 OutputRunner runner = new OutputRunner(genomeFile.getName().replace(".fasta", "").replace(".gz", ""), genomeFile, null);
			        	 runner.createGenome = true;
			        	 runner.execute();
			        	 
			          }
				 }
				 catch(Exception ex) {
					 ex.printStackTrace();
				 }		
			}
			else if (selectedNode!= null && selectedNode.isLeaf() && selectedNode.toString().contains("Add new anno")) {
				try {
					  FileDialog fs = new FileDialog(frame, "Select annotation gff3/gtf-file", FileDialog.LOAD);
			  		  fs.setDirectory(Main.downloadDir);	  		
			  		  fs.setVisible(true);
			  		  fs.setFile("*.gff3;*.gtf");
			  		  fs.setFilenameFilter(new FilenameFilter() {
				  			public boolean accept(File dir, String name) {
						        return name.toLowerCase().contains(".gff3") || name.toLowerCase().contains(".gtf");
						     }
						 });
			  		  String filename = fs.getFile();			
			       
			         if (filename != null) {
			        	 File addfile = new File(fs.getDirectory() +"/" +filename);
			        	
			        	 if(addfile.exists()) {			        		     	
			        		 annotationFile = addfile; 
			        		 Main.downloadDir = annotationFile.getParent();					        	
				        	 Main.writeToConfig("DownloadDir=" +annotationFile.getParent());		        	 
				        	 OutputRunner runner = new OutputRunner(selectedNode.getParent().toString(), null, annotationFile);
				        	 runner.createGenome = true;
				        	 runner.execute();
			        	 }
			        	 else {
			        		 Main.showError("File does not exists.", "Error", frame);
			        	 }
			        	 
			         }
			         if(1==1) {
			        	 return;
			         }
					 JFileChooser chooser = new JFileChooser(Main.downloadDir);	 
			    	  chooser.setMultiSelectionEnabled(false);
			    	  chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
			    	  chooser.setAcceptAllFileFilterUsed(false);	    	
			    	  MyFilterGFF gffFilter = new MyFilterGFF();	    
			    	
			    	  chooser.addChoosableFileFilter(gffFilter); 	    	
			    	  chooser.setDialogTitle("Select annotation gff3-file");
			    	  if(Main.screenSize != null) {
			    		  chooser.setPreferredSize(new Dimension((int)Main.screenSize.getWidth()/3, (int)Main.screenSize.getHeight()/3));
			    	  }
			          int returnVal = chooser.showOpenDialog((Component)this.getParent());	         
			          
			         if (returnVal == JFileChooser.APPROVE_OPTION) {
			        	
			        	 annotationFile = chooser.getSelectedFile(); 
			        	 Main.downloadDir = annotationFile.getParent();
			        	
			        	 Main.writeToConfig("DownloadDir=" +annotationFile.getParent());		        	 
			        	 OutputRunner runner = new OutputRunner(selectedNode.getParent().toString(), null, annotationFile);
			        	 runner.createGenome = true;
			        	 runner.execute();
			          }
				 }
				 catch(Exception ex) {
					 ex.printStackTrace();
				 }				
			}	      
		}
		if(e.getSource() == genometable) {
			
			if(new File(".").getFreeSpace()/1048576 < sizeHash.get(genometable.getValueAt(genometable.getSelectedRow(), 0))[0]/1048576) {
				sizeError.setVisible(true);
				download.setEnabled(false);
				AddGenome.getLinks.setEnabled(false);
			}
			else {			
				sizeError.setVisible(false);
				download.setEnabled(true);
				AddGenome.getLinks.setEnabled(true);
			}
			tree.clearSelection();
			remove.setEnabled(false);
			checkUpdates.setEnabled(false);
		}
		
	}
	static void downloadGenome(String urls) {
		try {		
			
			URL fastafile= AddGenome.genomeHash.get(urls)[0];
			String targetDir = "";
			//boolean writable = true;
			
			File test = new File(Main.genomeDir.getCanonicalPath() +"/test");
				
			
			
			//File f = new File(Main.genomeDir.getCanonicalPath());			
			
			if(test.mkdir()) {
				 /*if(fastafile.getFile().contains("GRCh38")) {
					
		    		   if((test.getFreeSpace()/1048576) < (60000000000L/1048576)) {
			    		   Main.showError("Sorry, you need more than 60GB of disk space to install GRCh38.\nGRCh38 FASTA file size is ~50GB uncompressed.\nThis drive has " +test.getFreeSpace()/1048576/1000 +"GB.", "Note");
			    		   test.delete();
			    		   return;
			    	   }						    	   					    	  
		    	   }
		    else*/ if((test.getFreeSpace()/1048576) < (5000000000L/1048576)) {
		    		   Main.showError("Sorry, you need more than 5GB of disk space.\nThis drive has " +test.getFreeSpace()/1048576/1000 +"GB.", "Note");
		    		   test.delete();
		    		   return;
		    	   }			
				test.delete();
				
				targetDir = Main.genomeDir.getCanonicalPath() +"/" +urls +"/";
				File fasta = new File(targetDir +FilenameUtils.getName(fastafile.getFile()));				
				URL gfffile= AddGenome.genomeHash.get(urls)[1];
				targetDir = Main.genomeDir.getCanonicalPath() +"/" +urls +"/annotation/" +FilenameUtils.getName(gfffile.getFile()) +"/";
				File gff = new File(targetDir +FilenameUtils.getName(gfffile.getFile()));		
				AddGenome.OutputRunner genomeadd = new AddGenome.OutputRunner(urls, fasta, gff, urls);
				genomeadd.createGenome = true;
				genomeadd.execute();
			} 
			else {				
				try {	  	    
		    		   
			    	  JFileChooser chooser = new JFileChooser();
			    	  chooser.setAcceptAllFileFilterUsed(false);			    	  
			    	  chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			    	  chooser.setDialogTitle("Select a local directory for genome files...");			                 	  
			          File outfile = null;
			          
			          while(true) {
			        	  int returnVal = chooser.showSaveDialog((Component)AddGenome.panel);
			        	  
					      if(returnVal == JFileChooser.APPROVE_OPTION) {   		    	  
						       outfile = chooser.getSelectedFile();
						       File genomedir = new File(outfile.getCanonicalPath() +"/genomes");
						       if(new File(outfile.getCanonicalPath() +"/genomes").mkdir()) {
						    	   if(fastafile.getFile().contains("GRCh38")) {
						    		   if((outfile.getFreeSpace()/1048576) < (200000000000L/1048576)) {
							    		   Main.showError("Please, select local drive with more than 60GB of disk space.\nGRCh38 FASTA file is ~50GB uncompressed.\nThis drive has " +outfile.getFreeSpace()/1048576/1000 +"GB.", "Note");
							    		   genomedir.delete();
							    		   continue;
							    	   }						    	   					    	  
						    	   }
						    	   else if((outfile.getFreeSpace()/1048576) < (5000000000L/1048576)) {
						    		   Main.showError("Please, select local drive with more than 5GB of disk space.\nThis drive has " +outfile.getFreeSpace()/1048576/1000 +"GB.", "Note");
						    		   genomedir.delete();
						    		   continue;
						    	   }						    	   					    	  
						       } 
						       else {						    	  
						    	  Main.showError("No writing permissions for this directory. \nPlease, select new directory for genomes.", "Error");
								  continue;						
						       }
						       
						       /*if (!new File(outfile.getCanonicalPath() +"/genomes").mkdir()) {
						    	   Main.showError("Could not create genome directory in " +outfile.getCanonicalPath(), "Error");
						    	   continue;
						       }*/
						       
						        break;
					      }
					      if(returnVal == JFileChooser.CANCEL_OPTION) {   	
					    	  outfile = null;
					    	  downloading = false;
					    	  break;
					      }
			          }
			          if(outfile != null) {
			        	  Main.genomeDir = new File(outfile.getCanonicalPath() +"/genomes");
			        	  genomedirectory.setText(Main.genomeDir.getCanonicalPath());
					       Main.writeToConfig("genomeDir=" +Main.genomeDir.getCanonicalPath());
					        targetDir = Main.genomeDir.getCanonicalPath() +"/" +urls +"/";
							File fasta = new File(targetDir +FilenameUtils.getName(fastafile.getFile()));				
							URL gfffile= AddGenome.genomeHash.get(urls)[1];
							targetDir = Main.genomeDir.getCanonicalPath() +"/" +urls +"/annotation/" +FilenameUtils.getName(gfffile.getFile()) +"/";
							File gff = new File(targetDir +FilenameUtils.getName(gfffile.getFile()));		
							AddGenome.OutputRunner genomeadd = new AddGenome.OutputRunner(urls, fasta, gff, urls);
							genomeadd.createGenome = true;
							genomeadd.execute();
			          }
			          downloading = false;
				       /* new File(outfile.getCanonicalPath() +"/genomes").mkdir();
					    File fasta = new File(targetDir +FilenameUtils.getName(fastafile.getFile()));							
						
						URL gfffile= AddGenome.genomeHash.get(urls)[1];
						targetDir = AddGenome.userDir +"/genomes/" +urls +"/annotation/" +FilenameUtils.getName(gfffile.getFile()) +"/";
						File gff = new File(targetDir +FilenameUtils.getName(gfffile.getFile()));
				
						AddGenome.OutputRunner genomeadd = new AddGenome.OutputRunner(urls, fasta, gff, urls);
						genomeadd.createGenome = true;
						genomeadd.execute();
				      */
					}
					catch(Exception ex) {
						downloading = false;
						ex.printStackTrace();
					}
			 
			}
		
			
		}
		catch(Exception e) {
			downloading = false;
			e.printStackTrace();
		}
		
	}

	@Override
	public void mouseReleased(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}


	@Override
	public void mouseEntered(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}


	@Override
	public void mouseExited(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
}
