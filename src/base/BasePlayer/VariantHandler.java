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
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.Feature;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;

import java.awt.Color;
import java.awt.Component;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.Serializable;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSlider;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.SwingWorker;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class VariantHandler extends JPanel implements ChangeListener, ActionListener, MouseListener, KeyListener, ComponentListener {
	
	private static final long serialVersionUID = 1L;
	static RangeSlider commonSlider = new RangeSlider(1,1);
	static JSlider qualitySlider = new JSlider(0,60);
	static JSlider gqSlider = new JSlider(0,60);
	static JSlider coverageSlider = new JSlider(1,40);
	static JSlider maxCoverageSlider = new JSlider(1,2000);
	static JSlider callSlider = new JSlider(0,100);
	static JSlider geneSlider = new JSlider(1,1);
	static JLabel geneLabel = new JLabel("At least 1/1 samples share a mutated gene");	
	static JLabel aminoCount = new JLabel("");
	static JLabel callsLabel = new JLabel();
	static JLabel coverageLabel = new JLabel();
	static JLabel maxCoverageLabel = new JLabel();
	static JLabel qualityLabel = new JLabel(), gqLabel = new JLabel(),comparison = new JLabel("Sample comparison");
	static JMenuBar menubar = new JMenuBar();
	static JMenu filters = new JMenu("Variant Filters");
	static JMenuBar aminobar = new JMenuBar();
	static JMenu aminomenu = new JMenu("Variant Annotator");
	static JMenu outputmenu = new JMenu("Save variant table");
	static JMenuItem write = new JMenuItem("Save");
	static JLabel slideLabel = new JLabel("Common variants in 1/1 samples");	
	static JLabel clusterLabel = new JLabel("Window size for variant clusters (use common variant slider)");
	static JTextField clusterBox = new JTextField("0");
	static ButtonGroup outputgroup = new ButtonGroup();
	static JRadioButton tsv = new JRadioButton("TSV");
	static JRadioButton compactTsv = new JRadioButton("Compact TSV");
	static JRadioButton vcf = new JRadioButton("VCF");
	static JCheckBox hidenoncoding = new JCheckBox("Hide non-coding variants");
	static JCheckBox freeze = new JCheckBox("Freeze filters");
	static JCheckBox rscode = new JCheckBox("Hide rs-coded variants");
	static JCheckBox allChroms = new JCheckBox("All chromosomes");
	static JCheckBox allChromsfrom = new JCheckBox("From this chr?");
	static JCheckBox hideSNVs = new JCheckBox("Hide SNVs");
	static JCheckBox hideHomos = new JCheckBox("Hide homozygotes");
	static JCheckBox onlyStats = new JCheckBox("Only stats");
	static JCheckBox hideIndels = new JCheckBox("Hide indels");
	static JCheckBox synonymous = new JCheckBox("Only non-synonymous");
	static JCheckBox nonsense = new JCheckBox("Only truncs");
	static JCheckBox intronic = new JCheckBox("Show intronic");
	static JCheckBox intergenic = new JCheckBox("Show intergenic");
	static JCheckBox utr = new JCheckBox("Show UTR");
	static JCheckBox onlyselected = new JCheckBox("Only selected sample");
	static JCheckBox writetofile = new JCheckBox("Write directly to a file");
	static JFrame frame = new JFrame("Variant handler");    
	static JLabel totalVars = new JLabel("Variant count on screen: 0");
	static JLabel empty = new JLabel("");
	static JButton varcalc = new JButton("Annotate");
	static JButton statcalc = new JButton("Stats");	
	static ArrayList<ControlFile> controlarray = new ArrayList<ControlFile>();
	String userDir;
	static int lastWrittenPos = 0;
	static  ArrayList<String> outputStrings = new ArrayList<String>();
	static MouseWheelListener sliderWheelListener;
	static JTabbedPane tabs = new JTabbedPane();
	static JScrollPane tableScroll = new JScrollPane();
	static JSeparator separator = new JSeparator();
	static AminoTable table;
	static StatsTable stattable;
	static JScrollPane statsScroll = new JScrollPane();
	static ClusterTable clusterTable;	
	static JScrollPane clusterScroll = new JScrollPane();
	static ArrayList<JScrollPane> tablescrolls = new ArrayList<JScrollPane>();
	static ArrayList<BedTable> tables = new ArrayList<BedTable>();
	static NodeSorter nodesorter = new NodeSorter();
	static JPopupMenu menu = new JPopupMenu("Advanced quality control");
	static JButton applyQualities = new JButton("Apply");
	static JButton advQualities = new JButton("More qualities...");
	static OwnVCFCodec vcfCodec= new OwnVCFCodec();
	static String format = "GT:DP:AD:GQ";
	int moveX=0, moveY=0, pressX=0,pressY=0;
	final int buttonHeight = 15, buttonWidth = 40;
	final Dimension buttondimension = new Dimension(buttonWidth, buttonHeight);
	static Color backColor = new Color(230,255,240,230);
	static int clusterSize = 0;
	static JPanel filterpanel = new JPanel() {				
		private static final long serialVersionUID = 1L;

		protected void paintComponent(Graphics g) {
		        super.paintComponent(g);     
		  
		        g.setColor(backColor);
		        g.fillRect(0, 0, this.getWidth(), this.getHeight());        
		       
		    }		
	};
	static JPanel comparepanel = new JPanel() {		
		
		private static final long serialVersionUID = 1L;

		protected void paintComponent(Graphics g) {
		        super.paintComponent(g);     
		        g.setColor(backColor);
		        g.fillRect(0, 0, this.getWidth(), this.getHeight());	        
		        g.setColor(Draw.softColor);
		        g.fillRect(0, 0, this.getWidth(), this.getHeight());	        
		    }		
	};
	
	static JPanel aminopanel = new JPanel(){		
		
		private static final long serialVersionUID = 1L;

		protected void paintComponent(Graphics g) {
		        super.paintComponent(g);     
		        
		        g.setColor(backColor);
		        g.fillRect(0, 0, this.getWidth(), this.getHeight());	        
		       
		    }		
	};

	static int maxCoverage = 1500;

	public VariantHandler() {	
		super(new GridBagLayout());	
		
	 
		GridBagConstraints c = new GridBagConstraints();	
		sliderWheelListener = new MouseWheelListener() {
		    @Override
		    public void mouseWheelMoved(MouseWheelEvent e) {
		        int notches = e.getWheelRotation();
		        JSlider slider = (JSlider)e.getSource();
		        if (notches < 0) {
		          
		        	slider.setValue(slider.getValue() + 1);
		        } else 
		         if(notches > 0) {   
		        	 slider.setValue(slider.getValue() - 1);
		        }
		    }
		};
		c.anchor = GridBagConstraints.NORTHWEST;
		c.fill = GridBagConstraints.BOTH;
		c.insets = new Insets(0,0,0,0);
		c.gridx = 0;
		c.gridy = 0;
		c.gridwidth = 1;		
		userDir = new File(Main.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParent().replace("%20", " ");
		
	
		this.setOpaque(false);
		if(Main.screenSize == null) {
			Main.screenSize = new Dimension(1000,1000);
		}
		stattable = new StatsTable((int)Main.screenSize.width, (int)Main.screenSize.height,statsScroll);
		stattable.setEnabled(true);
		table = new AminoTable((int)Main.screenSize.width, (int)Main.screenSize.height,tableScroll);
		
		table.setEnabled(true);
		clusterTable = new ClusterTable((int)Main.screenSize.width, (int)Main.screenSize.height,clusterScroll);
		clusterTable.setEnabled(true);
		filterpanel.setLayout(new GridLayout(7, 2));
		comparepanel.setLayout(new GridLayout(5, 2));
		aminopanel.setLayout(new GridLayout(1,3));
		
//		tablepanel.setLayout(new GridLayout(1,1));
		createButtons();
		tableScroll.setPreferredSize(new Dimension(500,400));
		
		tableScroll.getViewport().add(table);
		tableScroll.getVerticalScrollBar().addAdjustmentListener(
				
				new AdjustmentListener() {
					
					public void adjustmentValueChanged(AdjustmentEvent event) {
						
						table.repaint();
						
					}										
				}				
		);
		
		tableScroll.getVerticalScrollBar().setUnitIncrement(16);
		
		statsScroll.setPreferredSize(new Dimension(500,400));
		
		statsScroll.getViewport().add(stattable);
		statsScroll.getVerticalScrollBar().addAdjustmentListener(
				
				new AdjustmentListener() {
					
					public void adjustmentValueChanged(AdjustmentEvent event) {
						
						stattable.repaint();
						
					}										
				}				
		);
		
		clusterScroll.getVerticalScrollBar().setUnitIncrement(16);
		clusterScroll.setName("Clusters");
		clusterScroll.getVerticalScrollBar().setUnitIncrement(16);
		
		clusterScroll.setPreferredSize(new Dimension(500,400));
		
		clusterScroll.getViewport().add(clusterTable);
		clusterScroll.getVerticalScrollBar().addAdjustmentListener(
				
				new AdjustmentListener() {
					
					public void adjustmentValueChanged(AdjustmentEvent event) {
						
						clusterTable.repaint();
						
					}										
				}				
		);
		
		clusterScroll.getVerticalScrollBar().setUnitIncrement(16);
		
		add(filterpanel, c);
		
		c.gridy = 1;
		add(comparepanel,c);
		c.gridy = 2;
		add(aminopanel, c);
		c.gridy = 3;
		c.weightx = 1;
		c.weighty = 0.8;
		tabs.addMouseListener(this);
		tabs.setBackground(backColor);
		tabs.add("Genes", tableScroll);	
		tabs.add("Stats", statsScroll);
		//tabs.add("Clusters", clusterScroll);		
		
		add(tabs, c);
		
	}
	
	 protected void paintComponent(Graphics g) {
	        super.paintComponent(g);     
	        
	     //   g.drawImage(image, 0, 0, this.getWidth(), this.getHeight(), null);
	        g.setColor(backColor);
	        g.fillRect(0, 0, this.getWidth(), this.getHeight());	        
	       
	    }
	 
	void createButtons() {
		
		freeze.addChangeListener(
				new ChangeListener() {					
					public void stateChanged(ChangeEvent arg0) {
						
							if(freeze.isSelected()) {
								freezeFilters(true);
							}
							else {
								freezeFilters(false);
							}						
					}
					
				}
		);
		
		
		geneSlider.setPreferredSize(buttondimension);
		geneSlider.setMaximumSize(buttondimension);
		geneSlider.setMinimumSize(buttondimension);
		geneSlider.setValue(1);
		geneSlider.setSnapToTicks(true);
		geneSlider.setMajorTickSpacing(1);
		geneSlider.setMinorTickSpacing(1);
		geneSlider.setOpaque(false);
		geneSlider.addChangeListener(this);
		geneSlider.addMouseWheelListener(sliderWheelListener);
		commonSlider.addMouseWheelListener(sliderWheelListener);
		commonSlider.setPreferredSize(buttondimension);
		commonSlider.setValue(1);
		commonSlider.setUpperValue(1);		
		commonSlider.addChangeListener(this);
		commonSlider.setOpaque(false);
		qualitySlider.addMouseWheelListener(sliderWheelListener);
		qualitySlider.addChangeListener(this);
		qualitySlider.addMouseListener(this);
		qualitySlider.setPreferredSize(buttondimension);
		qualitySlider.setValue(0);
		qualitySlider.setOpaque(false);
		gqSlider.addMouseWheelListener(sliderWheelListener);
		gqSlider.addChangeListener(this);
		gqSlider.addMouseListener(this);
		gqSlider.setPreferredSize(buttondimension);
		gqSlider.setValue(0);
		gqSlider.setOpaque(false);		
		coverageSlider.addMouseWheelListener(sliderWheelListener);
		coverageSlider.addChangeListener(this);
		coverageSlider.addMouseListener(this);
		coverageSlider.setOpaque(false);
		coverageSlider.setValue(4);		
		maxCoverageSlider.addMouseWheelListener(sliderWheelListener);
		maxCoverageSlider.addChangeListener(this);
		maxCoverageSlider.addMouseListener(this);
		maxCoverageSlider.setOpaque(false);
		maxCoverageSlider.setValue(1500);
		callSlider.addMouseWheelListener(sliderWheelListener);
		callSlider.addChangeListener(this);
		callSlider.addMouseListener(this);
		callSlider.setValue(10);
		callSlider.setOpaque(false);
		
		hideSNVs.setOpaque(false);
		hideSNVs.addActionListener(this);
		hideIndels.setOpaque(false);
		hideIndels.addActionListener(this);
		hideHomos.setOpaque(false);
		hideHomos.addActionListener(this);
		rscode.addActionListener(this);
		rscode.setOpaque(false);
		
		hidenoncoding.addActionListener(this);
		hidenoncoding.setOpaque(false);
		frame.getContentPane().setBackground(Color.white);
		frame.setBackground(Color.white);
		frame.addComponentListener(this);
		menubar.setOpaque(true);
		filters.setOpaque(false);
		filters.add(hidenoncoding);
		filters.add(rscode);
		filters.add(hideSNVs);
		filters.add(hideIndels);
		filters.add(hideHomos);
		freeze.setBackground(new Color(100,150,255,50));
		filters.add(freeze);
		menubar.add(filters);
		filterpanel.add(menubar);
		filterpanel.add(totalVars);
		qualityLabel.setToolTipText("Variants below quality threshold will be hidden");
		gqLabel.setToolTipText("Variants below quality threshold will be hidden");
		filterpanel.add(qualityLabel);
		filterpanel.add(qualitySlider);
		filterpanel.add(gqLabel);
		filterpanel.add(gqSlider);
	
		menu.add(new JLabel("Advanced quality filters"));
		menu.addKeyListener(this);
		
		menu.add(applyQualities);
		applyQualities.addActionListener(this);
		advQualities.addActionListener(this);
		filterpanel.add(advQualities);
		filterpanel.add(new JSeparator());
		filterpanel.add(coverageLabel);
		filterpanel.add(coverageSlider);
		filterpanel.add(maxCoverageLabel);
		filterpanel.add(maxCoverageSlider);
		filterpanel.add(callsLabel);
		filterpanel.add(callSlider);		
		comparepanel.add(comparison);
		comparepanel.add(new JLabel(""));
		comparepanel.add(new JSeparator());
		comparepanel.add(new JSeparator());	
		comparepanel.add(geneLabel);
		comparepanel.add(geneSlider);
		comparepanel.add(slideLabel);
		comparepanel.add(commonSlider);
		comparepanel.add(clusterLabel);
		clusterBox.addKeyListener(this);
		comparepanel.add(clusterBox);
		
		varcalc.addActionListener(this);		
		varcalc.setPreferredSize(buttondimension);		
		allChroms.addActionListener(this);
		write.addActionListener(this);
		aminomenu.add(synonymous);
		aminomenu.add(nonsense);
		aminomenu.add(intronic);
		aminomenu.add(intergenic);
		aminomenu.add(utr);
	
		aminomenu.add(allChroms);
		aminomenu.add(onlyselected);
		onlyStats.addActionListener(this);
		aminomenu.add(onlyStats);
		aminomenu.add(writetofile);
		writetofile.addActionListener(this);
		aminomenu.add(varcalc);
	
		aminopanel.add(aminobar);
		aminopanel.add(aminoCount);
		aminobar.add(aminomenu);
		tsv.setSelected(true);
		outputgroup.add(tsv);
		outputgroup.add(compactTsv);
		outputgroup.add(vcf);
		outputmenu.add(tsv);
		outputmenu.add(compactTsv);
		outputmenu.add(vcf);
	//	outputmenu.add(forVEP);
		outputmenu.add(write);
		aminobar.add(outputmenu);
	//	add(table);
	}
	
	private static void createAndShowGUI() {	
		   
		JFrame.setDefaultLookAndFeelDecorated(false);
	   
	
	    frame.setResizable(true);    
	    JComponent newContentPane = new VariantHandler();
	    newContentPane.setOpaque(true); 
	    frame.setContentPane(newContentPane);
	    frame.pack();
	    table.setPreferredSize(new Dimension(tableScroll.getViewport().getWidth(),tableScroll.getViewport().getHeight()));
		table.setMinimumSize(new Dimension(tableScroll.getViewport().getWidth(),tableScroll.getViewport().getHeight()));
		table.resizeTable(tableScroll.getViewport().getWidth());		
	    aminobar.setMinimumSize(new Dimension((int)aminobar.getSize().getWidth(),(int)aminobar.getSize().getHeight()));
	    filters.setMinimumSize(new Dimension((int)filters.getSize().getWidth(),(int)filters.getSize().getHeight()));
	    clusterTable.resizeTable(tableScroll.getViewport().getWidth());
	    if(Main.chromDraw == null) {
	    	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
		    frame.setVisible(true);
	    }
	    else {
		    frame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE); 
		    frame.setVisible(false);
	    }
	    
	}
	
	public static void main(String[] args) {
		
		 javax.swing.SwingUtilities.invokeLater(new Runnable() {
	           public void run() {
	           
	           		createAndShowGUI();
	           	
	           }
	       });		
	}
	@Override
	public void stateChanged(ChangeEvent event) {
		if(event.getSource() == commonSlider){
			slideLabel.setText("Common variants in " +commonSlider.getValue() +"/" +commonSlider.getUpperValue() +" samples");
			Draw.calculateVars = true;
			Draw.updatevars = true;
			
			if((commonSlider.getValue() > 1 || commonSlider.getUpperValue() < Main.varsamples) && (!clusterBox.getText().equals("0") && !clusterBox.getText().equals(""))) {
				clusterSize = Integer.parseInt(clusterBox.getText());
				
				Main.drawCanvas.calcClusters(FileRead.head,1);
				if(tabs.indexOfComponent(clusterScroll) == -1) {
						
					tabs.add(clusterScroll, tabs.indexOfComponent(statsScroll));
				}
			}
			else if(clusterSize != 0) {
				if(clusterBox.getText().equals("0") || clusterBox.getText().equals("")) {
					clusterSize = 0;
				}
				if(tabs.indexOfComponent(clusterScroll) != -1) {
					
					tabs.remove(clusterScroll);
				}
			}
			else {
				if(clusterSize == 0) {
					if(tabs.indexOfComponent(clusterScroll) != -1) {
						
						tabs.remove(clusterScroll);
					}
				}
				
			}
			if(Main.drawCanvas != null) {
				
				Main.drawCanvas.repaint();
				
				if(Main.drawCanvas.splits.get(0).viewLength <=1000000) {
					Main.chromDraw.updateExons = true;
					Main.chromDraw.repaint();
				}
			}
			return;
		}
		else if(event.getSource() == geneSlider){
			geneLabel.setText("At least " +geneSlider.getValue() +"/" +Main.varsamples +" samples share a mutated gene");
		}
		else if(event.getSource() == qualitySlider){
			
			qualityLabel.setText("Min. quality score: " +qualitySlider.getValue());
			Draw.calculateVars = true;
			Draw.updatevars = true;
			if(Main.drawCanvas != null) {
				Main.drawCanvas.repaint();
				if(Main.drawCanvas.splits.get(0).viewLength <=1000000) {
					Main.chromDraw.updateExons = true;
					Main.chromDraw.repaint();
				}
			}
			return;
		}
		else if(event.getSource() == gqSlider) {
			gqLabel.setText("Min. genotype quality score: " +gqSlider.getValue());
			Draw.calculateVars = true;
			Draw.updatevars = true;
			if(Main.drawCanvas != null) {
				Main.drawCanvas.repaint();
				if(Main.drawCanvas.splits.get(0).viewLength <=1000000) {
					Main.chromDraw.updateExons = true;
					Main.chromDraw.repaint();
				}
			}
			return;
		}
		else if(event.getSource() == coverageSlider){
			
			coverageLabel.setText("Min. coverage: " +coverageSlider.getValue());
			Draw.calculateVars = true;
			Draw.updatevars = true;
			if(Main.drawCanvas != null) {
				Main.drawCanvas.repaint();
				if(Main.drawCanvas.splits.get(0).viewLength <=1000000) {
					Main.chromDraw.updateExons = true;
					Main.chromDraw.repaint();
				}
			}
			return;
		}
		else if(event.getSource() == callSlider){
			callsLabel.setText("Min. allelic fraction: " +callSlider.getValue() +"%");
			Draw.calculateVars = true;
			Draw.updatevars = true;
			if(Main.drawCanvas != null) {
				Main.drawCanvas.repaint();
				if(Main.drawCanvas.splits.get(0).viewLength <=1000000) {
					Main.chromDraw.updateExons = true;
					Main.chromDraw.repaint();
				}
			}
			return;
		}
		else if(event.getSource() == maxCoverageSlider){			
			maxCoverageLabel.setText("Max. coverage: " +maxCoverageSlider.getValue());
			maxCoverage = maxCoverageSlider.getValue();
			Draw.calculateVars = true;
			Draw.updatevars = true;
			if(Main.drawCanvas != null) {
				Main.drawCanvas.repaint();
				if(Main.drawCanvas.splits.get(0).viewLength <=1000000) {
					Main.chromDraw.updateExons = true;
					Main.chromDraw.repaint();
				}
			}
			return;
		}
		
	}
	public static class NodeSorter implements Comparator<VarNode> {
		
		public int compare(VarNode o1, VarNode o2) {  
			
				if ( o1.getPosition() <= o2.getPosition()) { 	        	
	        		return -1;  
	        	}
	        	else {
	        		return 1;
	        	}
	                
	        
			
		}
	}
	@Override
	public void actionPerformed(ActionEvent event) {
		
		if(event.getSource() == writetofile) {
			
			if(writetofile.isSelected()) {				
				//tabs.setEnabled(false);
				aminomenu.add(tsv,aminomenu.getItemCount()-1);
				aminomenu.add(compactTsv,aminomenu.getItemCount()-1);
				aminomenu.add(vcf,aminomenu.getItemCount()-1);
				aminomenu.getPopupMenu().pack();
				varcalc.setText("Annotate and write");
				aminomenu.revalidate();
				aminomenu.repaint();
				tabs.revalidate();
			}
			else {				
				//tabs.setEnabled(true);				
				outputmenu.add(vcf,0);
				outputmenu.add(compactTsv,0);
				outputmenu.add(tsv,0);
				outputmenu.revalidate();
				aminomenu.getPopupMenu().pack();
				varcalc.setText("Annotate");
				tabs.revalidate();
			}
		}
		else if(event.getSource() == hidenoncoding) {
		
			Draw.calculateVars = true;
			if(commonSlider.getValue() > 1) {			
			
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}
			Draw.updatevars = true;		
			Main.drawCanvas.repaint();
			
			if(Main.drawCanvas.splits.get(0).viewLength <=1000000) {
				Main.chromDraw.updateExons = true;
				Main.chromDraw.repaint();
			}
		}
		else if(event.getSource() == varcalc) {
			 
			if(writetofile.isSelected()) {
				try {
			  		String outfname = "";
			  		String path =Main.savedir;			 
			  		
		    		//path = new java.io.File(".").getCanonicalPath();		    		    
		    		JFileChooser chooser = new JFileChooser(path);		    
		    		int returnVal = chooser.showSaveDialog((Component)this.getParent());	         
			    	  
		    		
		    		if(returnVal == JFileChooser.APPROVE_OPTION) {  		    	  
		    		outfname = chooser.getSelectedFile().getAbsolutePath();		      	
		    		BufferedWriter output = null;
		    		Main.savedir = chooser.getSelectedFile().getParent();
		        	 Main.writeToConfig("DefaultSaveDir=" +chooser.getSelectedFile().getParent());
		        	 lastWrittenPos = 0;
			    	 if(tsv.isSelected() || compactTsv.isSelected()) {
			    		 if(!outfname.contains(".tsv")) {
					    	   File outfile = new File(outfname +".tsv");
					    	//   statout = new BufferedWriter(new FileWriter(new File (outfname +"_stats.tsv")));
					    	   FileRead.outputName = outfname +".tsv";
					    	   output = new BufferedWriter(new FileWriter(outfile));
					       }
					       else {
					    	   File outfile = new File(outfname);
					    	   FileRead.outputName = outfname;
					    	 //  statout = new BufferedWriter(new FileWriter(new File (outfname.replace(".tsv", "") +"_stats.tsv")));
					    	   output = new BufferedWriter(new FileWriter(outfile));
					       }    
			    		 
			    		 if(output != null) {
			    			String header = createTSVHeader();
			    	    	output.write(header);
			    			VariantHandler.table.clear();
			 				VariantHandler.table.headerHover = 2;
			 				VariantHandler.table.sorter.ascending = true;
			 				VariantHandler.table.createPolygon();
			 				VariantHandler.table.repaint();	
			 				FileRead.output = output;			 				
			 				FileRead calculator = new FileRead();
			 				calculator.varcalc = true;
			 				calculator.execute();						    							    	
			    		 }			    		 
			    	 }
			    	 else {
			    		 if(vcf.isSelected()) {
			    			SAMSequenceDictionary dict = AddGenome.ReadDict(Main.ref);
			    			FileRead.indexCreator = new TabixIndexCreator(dict, TabixFormat.VCF);			    			
			    			FileRead.filepointer = 0;
			    			
			    			
			    			if(!outfname.contains(".vcf")) {
			    				outfname = outfname +".vcf";
						   // 	   File outfile = new File(outfname +".vcf");
						    	//   statout = new BufferedWriter(new FileWriter(new File (outfname +"_stats.tsv")));
						   // 	   FileRead.outputName = outfname +".vcf";
						   // 	   output = new BufferedWriter(new FileWriter(outfile));
						       }
						       else {
						    	   File outfile = new File(outfname);
						    	   FileRead.outputName = outfname;
						    	 //  statout = new BufferedWriter(new FileWriter(new File (outfname.replace(".tsv", "") +"_stats.tsv")));
						  //  	   output = new BufferedWriter(new FileWriter(outfile));
						       }    
				    			
				    		 
				    			
				    			 if(!outfname.endsWith(".gz")) {
				    				 outfname = outfname +".gz";
				    			 }
				    			 FileRead.lastpos = 0;
				    			 FileRead.outputgz = new BlockCompressedOutputStream(outfname);
				    			 FileRead.outFile = new File(outfname);
				    			String header = createVCFHeader();
				    			VCFHeader vcfheader = new VCFHeader();
				    			VCFHeaderLine headerline = new VCFHeaderLine("format","##fileformat=VCFv4.1");
				    			vcfheader.addMetaDataLine(headerline);
				    			vcfCodec.setVCFHeader(vcfheader, VCFHeaderVersion.VCF4_1);
				    			
				    			FileRead.outputgz.write(header.getBytes());
				    			VariantHandler.table.clear();
				 				VariantHandler.table.headerHover = 2;
				 				VariantHandler.table.sorter.ascending = true;
				 				VariantHandler.table.createPolygon();
				 				VariantHandler.table.repaint();	
				 				
				 				//FileRead.output = output;			 				
				 				FileRead calculator = new FileRead();
				 				calculator.varcalc = true;
				 				calculator.execute();					    							    	
				    		 			    		
			    		 }
			    	 }			    	 
		    		}
			     }
			     catch(Exception e) {
			    	 JOptionPane.showMessageDialog(Main.chromDraw, e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
			    	 e.printStackTrace();
			     }
				
			}
			else {
				FileRead.outputgz = null;
				VariantHandler.table.clear();
				VariantHandler.table.headerHover = 2;
				VariantHandler.table.sorter.ascending = true;
				VariantHandler.table.createPolygon();
				VariantHandler.table.repaint();	
				
				FileRead calculator = new FileRead();
				calculator.varcalc = true;
				calculator.execute();	
			}
		}		
		else if(event.getSource() == hideSNVs) {
			Draw.calculateVars = true;
			if(commonSlider.getValue() > 1) {
				
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}
			Draw.updatevars = true;
			Main.drawCanvas.repaint();
			if(Main.drawCanvas.splits.get(0).viewLength <=1000000) {
				Main.chromDraw.updateExons = true;
				Main.chromDraw.repaint();
			}
		}
		else if(event.getSource() == hideIndels) {
			Draw.calculateVars = true;
			Draw.updatevars = true;
			if(commonSlider.getValue() > 1) {
				
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}
			Main.drawCanvas.repaint();
			if(Main.drawCanvas.splits.get(0).viewLength <=1000000) {
				Main.chromDraw.updateExons = true;
				Main.chromDraw.repaint();
			}			
		}
		else if(event.getSource() == hideHomos) {
			Draw.calculateVars = true;
			if(commonSlider.getValue() > 1) {
				
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}
			Draw.updatevars = true;
			Main.drawCanvas.repaint();
			if(Main.drawCanvas.splits.get(0).viewLength <=1000000) {
				Main.chromDraw.updateExons = true;
				Main.chromDraw.repaint();
			}
		}
		else if(event.getSource() == rscode) {
			Draw.calculateVars = true;
			Draw.updatevars = true;
			if(commonSlider.getValue() > 1) {
				
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}
			Main.drawCanvas.repaint();
			if(Main.drawCanvas.splits.get(0).viewLength <=1000000) {
				Main.chromDraw.updateExons = true;
				Main.chromDraw.repaint();
			}			
		}
		else if(event.getSource() == allChroms) {
			if(allChroms.isSelected()) {
				aminomenu.add(allChromsfrom,6);			
				aminomenu.getPopupMenu().pack();
			}
			else {
				aminomenu.remove(allChromsfrom);				
				aminomenu.getPopupMenu().pack();
			}
		}
		
		else if(event.getSource() == write) {
			
			 try {
				
		  		String outfname = "";
		  		String path =Main.savedir;			 
	    		//path = new java.io.File(".").getCanonicalPath();		    		    
	    		JFileChooser chooser = new JFileChooser(path);	
	    		if(tabs.getSelectedComponent().equals(statsScroll)) {
	    			chooser.setDialogTitle("Save variant statistics");
	    		}
	    		else {
	    			chooser.setDialogTitle("Save variants");
	    		}
	    		
	    		int returnVal = chooser.showSaveDialog((Component)this.getParent());	  
	    		
		      
	    		if(returnVal == JFileChooser.APPROVE_OPTION) {  		    	  
	    		outfname = chooser.getSelectedFile().getCanonicalPath();		      	
	    		BufferedWriter output = null;
	    		Main.savedir = chooser.getSelectedFile().getParent();
	        	 Main.writeToConfig("DefaultSaveDir=" +chooser.getSelectedFile().getParent());
		     try {
		    	 File outfile = null;
		    	 if(tsv.isSelected() || compactTsv.isSelected()) {
		    		 if(!outfname.contains(".tsv")) {
				    	  outfile = new File(outfname +".tsv");
				    	  
				    	   output = new BufferedWriter(new FileWriter(outfile));
				       }
				       else {
				    	  outfile = new File(outfname);
				    	 
				    	   output = new BufferedWriter(new FileWriter(outfile));
				       }    
		    		 
		    		 if(output != null) {
		    			 	
					    	//TODO
		    			 	OutputRunner runner = new OutputRunner(output, null,null);
		    			 	runner.execute();
					    	
		    		 }
		    		 
		    	 }
		    	 else if(vcf.isSelected()) {
		    		 BlockCompressedOutputStream outputgz;
		    	
		    		 if(!outfname.contains(".vcf")) {
				    	   outfile = new File(outfname +".vcf.gz");
				    	   outputgz = new BlockCompressedOutputStream(outfile);
				       }
				       else {
				    	   if(!outfname.endsWith(".gz")) {
				    		   outfile = new File(outfname +".gz");
					    	   outputgz = new BlockCompressedOutputStream(outfile);
				    	   }
				    	   else {
					    	   outfile = new File(outfname);
					    	   outputgz = new BlockCompressedOutputStream(outfile);
				    	   }
				       }  
		    		
		    		 if(outputgz != null) {
		    			 	
					    	//TODO
		    			 	OutputRunner runner = new OutputRunner(output, outputgz,outfile);
		    			 	runner.execute();
					    	
		    		 }
		    		 	
			    	
					 
		    	 }
		    	 /*else if(forVEP.isSelected()) {
		    		
				    	   File outfile = new File(outfname);
				    	   output = new BufferedWriter(new FileWriter(outfile));
				    	   
				    	   if(output != null) {
						    	//TODO
						    	 Main.drawCanvas.loading("Writing output...");
						    	 String[] row, position;
						    	 VarNode node;
					//	    	 String rscode, uniprot;
						    	 
						    	 for(int trans = 0; trans < table.transarray.size(); trans++) {
						    		 table.getAminos(table.transarray.get(trans));
						    		 
							    	 for(int s = 0; s<table.aminoarray.size(); s++) {
							    		 row = table.aminoarray.get(s).getRow();
							    		 node = table.aminoarray.get(s).getNode();
							    	//	 rscode = table.getRscode(node, row[5]);
							    		 position = row[2].split(":");
							    		 //getaminosiin jo uniprottikoodi
						//	    		 uniprot = node.getExons().get(0).getTranscript().getUniprot();
							    		 output.write(position[0] +" " +position[1] +" " +position[1] +" " +node.getRefBase() +"/" +row[5] +" 1"); 
							 
							    	 }
						    	 
						    	 }
						    	 output.close();
						    	 Main.drawCanvas.ready("Writing output...");
						     }
				    	   
				    	   
		    	 }   */
		       	  
		     }
		     catch(Exception ex) {
		    	 JOptionPane.showMessageDialog(frame, ex.getMessage());		    	
		     }
	     }
			 
			 }
			 catch(Exception e) {
				 ErrorLog.addError(e.getStackTrace());
				e.printStackTrace(); 
			 }
		}
		else if(event.getSource() == advQualities) {
			menu.show(this, 100, 100);
		}
		else if(event.getSource() == applyQualities) {
			Main.drawCanvas.drawVariables.advQDraw = null;
			for(int i = 1; i<menu.getComponentCount(); i++) {
				if(menu.getComponent(i) instanceof JLabel) {
					JLabel label = (JLabel)menu.getComponent(i);
					JTextField field = (JTextField)menu.getComponent(i+1);
					String format = "<";
					try {		
						Float number = null;
						field.setForeground(Color.black);
						if(field.getText().length() > 0) {
							if(field.getText().trim().startsWith("<")) {
								number = Float.parseFloat(field.getText().substring(1).trim());
							}
							else if(field.getText().trim().startsWith("<=")) {
								format += "=";
								 number = Float.parseFloat(field.getText().substring(2).trim());
							}							
							else if(field.getText().trim().startsWith(">=")) {
								format = ">=";
								number = Float.parseFloat(field.getText().substring(2).trim());
							}
							else if(field.getText().trim().startsWith(">")){
								format = ">";
								number = Float.parseFloat(field.getText().substring(1).trim());
							}
							else {
								number = Float.parseFloat(field.getText().trim());
							}							
							
							if(Main.drawCanvas.drawVariables.advQDraw == null) {
								Main.drawCanvas.drawVariables.advQDraw = new ArrayList<QualEntry>();
							}							
							Main.drawCanvas.drawVariables.advQDraw.add(new QualEntry(label.getText(), number, format));							
						}						
					}
					catch(Exception e) {
						field.setForeground(Color.red);
					}
					
				}
			}
				FileRead.search = true;
			
				Main.drawCanvas.gotoPos(Main.drawCanvas.splits.get(0).chrom, Main.drawCanvas.splits.get(0).start, Main.drawCanvas.splits.get(0).end);
			
			}
		
	
}
static void removeMenuComponents() {
	for(int i = menu.getComponentCount()-2; i>0; i--) {
		if(menu.getComponent(i) instanceof JLabel || menu.getComponent(i) instanceof JTextField) {
			menu.remove(i);			
		}
	}
}
/*	static void writeTranscriptsToVCF(ArrayList<Gene> genes, BufferedWriter output) {
		if(genes.size() == 0) {
			return;
		}
		 ArrayList<VarNode> nodes = new ArrayList<VarNode>();
		 int lastpos = 0;
		 String lastChrom;
		 
		 for(int v=0; v<genes.get(0).varnodes.size(); v++) {		    					
			 if(!nodes.contains(genes.get(0).varnodes.get(v))) {		    						
				 nodes.add(genes.get(0).varnodes.get(v));		    								    						
			 }		    					
		 }	  
		 lastpos = genes.get(0).getEnd();
		 lastChrom = genes.get(0).getChrom();
		 for(int i = 1 ; i< genes.size() ; i++) {
			 if(lastpos < genes.get(i).getStart() || !lastChrom.equals( genes.get(i).getChrom()) || i == genes.size()-1) {
				 Collections.sort(nodes, nodesorter);
				 for(int j = 0 ; j<nodes.size(); j++) {
					
					 writeNodeToFile(nodes.get(j),genes.get(i).getChrom(), output);				    					
				 }
				 nodes.clear();
			 }
			 if(lastpos < genes.get(i).getEnd() ) {
				 lastpos = genes.get(i).getEnd();
			 }
			 if(!lastChrom.equals(genes.get(i).getChrom())) {
				 lastChrom = genes.get(i).getChrom();
			 }
			 for(int v=0; v<genes.get(i).varnodes.size(); v++) {		    					
				 if(!nodes.contains(genes.get(i).varnodes.get(v))) {		    						
					 nodes.add(genes.get(i).varnodes.get(v));		    								    						
				 }		    					
			 }	    				
		 }
	}*/
	String createVCFHeader() {
		StringBuffer headerstring = new StringBuffer("##fileformat=VCFv4.1\n"
				+ "##BasePlayer=<Version: " +Main.version +" output " +new SimpleDateFormat("dd.MM.yyyy HH:mm").format(Calendar.getInstance().getTime())+"\n"
				+ "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
				+ "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">\n"
				+ "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n"
				+ "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
				+ "##reference=" +Main.ref.getName() +"\n");
				
				
				if(VariantHandler.onlyselected.isSelected()) {
					headerstring.append("##sample="+Main.drawCanvas.selectedSample.getName()+"\n");
					headerstring.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+Main.drawCanvas.selectedSample.getName() +"\n");
				}
				else {
					StringBuffer headersamples = new StringBuffer("");
					headerstring.append("##files=");
					for(int i = 0; i<Main.drawCanvas.sampleList.size(); i++) {						
						if((Main.drawCanvas.sampleList.get(i).getTabixFile() != null)) {
							headerstring.append(Main.drawCanvas.sampleList.get(i).getName() +",");
							
						}	
						if(!Main.drawCanvas.sampleList.get(i).multiVCF && (Main.drawCanvas.sampleList.get(i).getTabixFile() != null || Main.drawCanvas.sampleList.get(i).multipart) && !Main.drawCanvas.sampleList.get(i).removed) {
							
							headersamples.append("\t" +Main.drawCanvas.sampleList.get(i).getName());
						}
					}					
					headerstring.deleteCharAt(headerstring.length()-1);
					headerstring.append("\n");
					headerstring.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
					headerstring.append(headersamples +"\n");
					headersamples = new StringBuffer();
				}
	
		return headerstring.toString();
	}
	String createTSVHeader() {
		 StringBuffer headerstring = new StringBuffer("");
   	  headerstring.append("##BasePlayer version: " +Main.version +" output " +new SimpleDateFormat("dd.MM.yyyy HH:mm").format(Calendar.getInstance().getTime())+"\n");
   	  
   	  headerstring.append("##Files: ");
	   	if(VariantHandler.onlyselected.isSelected()) {
			headerstring.append(Main.drawCanvas.selectedSample.getName());
		}
		else {
	   	  for(int i = 0; i< Main.samples; i++) {
	   		if((Main.drawCanvas.sampleList.get(i).getTabixFile() != null || Main.drawCanvas.sampleList.get(i).multipart) && !Main.drawCanvas.sampleList.get(i).removed) {
				
	   			headerstring.append(Main.drawCanvas.sampleList.get(i).getName() +",");
	   		  }   		  
	   	  }
	   	 headerstring.deleteCharAt(headerstring.length()-1);
		}
   	 
   	  headerstring.append("\n");
   	  headerstring.append("##Genome:" +Main.ref.getName() +",Annotation:" +Main.annotationfile +"\n");   	  
   	 
   	  headerstring.append("##Exclude synonymous:" +VariantHandler.synonymous.isSelected() +",Show UTR:" +VariantHandler.utr.isSelected() +",+intronic:" +VariantHandler.intronic.isSelected()+",Only truncs:" +VariantHandler.nonsense.isSelected() +"\n"); 
   	  
   	  if(Control.controlData.fileArray.size() > 0 && Control.controlData.controlsOn) {
   		 
	    		  headerstring.append("##Controls:");
		    	  for(int i = 0; i< Control.controlData.fileArray.size(); i++) {
		    		  if(Control.controlData.fileArray.get(i).controlOn) {
		    			  headerstring.append(Control.controlData.fileArray.get(i).getName() +":" +Control.controlData.fileArray.get(i).alleleFreq +",");
		    		  }		    		  
		    	  }	    		  
		    	  headerstring.deleteCharAt(headerstring.length()-1);
		    	  headerstring.append("\n");
	    	  
   	  }
   	  else {
   		  headerstring.append("##No controls applied\n");
   	  }
   	  if(Main.bedCanvas.bedOn) {
   		  headerstring.append("##Tracks:");
   		  for(int i=0; i<Main.bedCanvas.bedTrack.size(); i++) {
   			  if(Main.bedCanvas.bedTrack.get(i).intersect) {
   				  if(Main.bedCanvas.bedTrack.get(i).limitValue > Double.MIN_VALUE) {
   					  headerstring.append(Main.bedCanvas.bedTrack.get(i).file.getName() +">=" +Main.bedCanvas.bedTrack.get(i).limitValue +",");
   				  }
   				  else {
   					  headerstring.append(Main.bedCanvas.bedTrack.get(i).file.getName() +",");
   				  }
   		      }
   		  }
   		  headerstring.deleteCharAt(headerstring.length()-1);
   		  headerstring.append("\n");
   	  }
   	  if(VariantHandler.commonSlider.getValue() > 1) {
   		  
   		  if(VariantHandler.clusterSize > 0) {
	    	headerstring.append("##Variant clusters in " +VariantHandler.commonSlider.getValue() +"/" +Main.varsamples +" samples within " +VariantHandler.clusterSize +"bp\n");
	      }
   		  else {
   			headerstring.append("##Common variants in " +VariantHandler.commonSlider.getValue() +"/" +Main.varsamples +" samples\n");   			
   		  }   		  
   	  }
   	  if(VariantHandler.geneSlider.getValue() > 1) {
   		  headerstring.append("##At least " +VariantHandler.geneSlider.getValue() +"/" +Main.varsamples +" samples share a mutated gene\n");
   	  }
   	  headerstring.append("##Variant filters:\n");
   	  headerstring.append("##Hide rs-coded variants: " +VariantHandler.rscode.isSelected() +"\n");
   	  headerstring.append("##Hide SNVs: " +VariantHandler.hideSNVs.isSelected() +"\n");
   	  headerstring.append("##Hide indels: " +VariantHandler.hideIndels.isSelected() +"\n");
   	  headerstring.append("##Min. coverage: " +VariantHandler.coverageSlider.getValue() +"\n");
   	  headerstring.append("##Min. allelic/fraction: " +VariantHandler.callSlider.getValue() +"%\n");
   	  headerstring.append("##Min. quality score: " +VariantHandler.qualitySlider.getValue() +"\n");   
   	headerstring.append("##Min. genotype quality score: " +VariantHandler.gqSlider.getValue() +"\n");   
   	  headerstring.append("##Max. coverage: " +VariantHandler.maxCoverageSlider.getValue() +"\n");   	  
   	  StringBuffer controls = new StringBuffer("");	 

	 if(Control.controlData.controlsOn) {
		 controlarray.clear();
		 for(int i = 0;i<Control.controlData.fileArray.size(); i++) {
			
			 if(!Control.controlData.fileArray.get(i).controlOn) {
				 continue;
			 }
			 controls.append("AF: "+Control.controlData.fileArray.get(i).getName() +"\tOR\t");
			 controlarray.add(Control.controlData.fileArray.get(i));    		
		 }
	 }
	 StringBuffer tracks = new StringBuffer("");	 
   	 if(Main.bedCanvas.bedOn) {
   		 for(int i = 0; i<Main.bedCanvas.bedTrack.size(); i++) {
   			 if(Main.bedCanvas.bedTrack.get(i).intersect) {
   				 tracks.append(Main.bedCanvas.bedTrack.get(i).file.getName() +"\t");
   			 }
   		 }
   	 }
   	 String clusters = "";
   	 if(commonSlider.getValue() > 1 && clusterSize > 0) {
   		clusters = "ClusterID\tClusterMutCount\tClusterWidth\tClusterMutFreq\t"; 
   	 }
   	 
	   	if(!tabs.getSelectedComponent().equals(statsScroll)) {
	   		headerstring.append("#Sample\tGene\tMutationCount\tSampleCount\tENSG\tENST\tBioType\tPosition\tStrand\tRegion\tEffect\tBaseChange\tGenotype(calls/coverage)\tquality\trs-code\t" +clusters +controls +tracks+"Description\n");
	   	}
	   	else {
	   		headerstring.append("#Sample\tVariants\tSNVs\tDELs\tINSs\tCoding\tHetero/homo-rate\tTS/TV-rate\tT>A\tT>C\tT>G\tC>A\tC>G\tC>T\tAvg.call/cov\n");
	   	}
   	  return headerstring.toString();
	}
	
	void writeOutput(BufferedWriter output, BlockCompressedOutputStream outputgz, File outFile) {
		 frame.getGlassPane().setVisible(true);
		 table.setEnabled(false);		 
		 frame.getGlassPane().setCursor( Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR)); 
		try {
			
		if(tabs.getSelectedComponent().equals(statsScroll)) {
			 String header = createTSVHeader();	    	 
	    	 output.write(header);
	    	 Sample sample;
	   // 	 output.write("#Sample\tVariants\tSNVs\tDELs\tINSs\tCoding\tHetero/homo-rate\tTS/TV-rate\tT>A\tT>C\tT>G\tC>A\tC>G\tC>T\tAvg.call/cov\n");
	 	  	 for(int i = 0 ; i<VariantHandler.stattable.sampleArray.size(); i++) {
	 	  		 sample = (Sample)VariantHandler.stattable.sampleArray.get(i)[0];
	 	  		 output.write(sample.getName() +"\t");
	 	  		 
	 	  		 for(int j = 1; j<VariantHandler.stattable.headerlengths.length-1; j++) {
	 	  			output.write(VariantHandler.stattable.sampleArray.get(i)[j] +"\t");
	 	  		 }
	 	  		 output.write(VariantHandler.stattable.sampleArray.get(i)[VariantHandler.stattable.headerlengths.length-1] +"\n");
	 	  	 }
	    	output.close();
		}
		else {
			
	    	 if(vcf.isSelected()) {
	    		 String header = createVCFHeader();
	    		 if(outputgz != null) {
	    			SAMSequenceDictionary dict = AddGenome.ReadDict(Main.ref);
		    		FileRead.indexCreator = new TabixIndexCreator(dict, TabixFormat.VCF);			    			
		    		FileRead.filepointer = 0;
    			 	VCFHeader vcfheader = new VCFHeader();
	    			VCFHeaderLine headerline = new VCFHeaderLine("format","##fileformat=VCFv4.1");
	    			vcfheader.addMetaDataLine(headerline);
	    			vcfCodec.setVCFHeader(vcfheader, VCFHeaderVersion.VCF4_1);		    			
	    			outputgz.write(header.getBytes());
	    		 }
	    		 else {
	    			output.write(createVCFHeader());  	
	    		 }
	    		
	    		 writeGeneListToVCF(table.genearray, output, outputgz);
	    		 if(outputgz != null) {
	    				
	    				for(int i = 0 ; i<VariantHandler.outputStrings.size(); i++) {
	    					 outputgz.write(VariantHandler.outputStrings.get(i).getBytes());
	    					 
	    					 Feature vcf = VariantHandler.vcfCodec.decode(VariantHandler.outputStrings.get(i));
	    					
	    					 FileRead.indexCreator.addFeature(vcf, FileRead.filepointer);
	    					 FileRead.filepointer = outputgz.getFilePointer();
	    				 }
	    				VariantHandler.outputStrings.clear();
	    				outputgz.flush();
	    				 
	    			    Index index = FileRead.indexCreator.finalizeIndex(outputgz.getFilePointer());
	    			    
	    			    index.writeBasedOnFeatureFile(outFile);
	    				outputgz.close();
	    				
	    			}
	    	 }
	    	 else {
	    		 output.write(createTSVHeader());  	
	    		 for(int gene = 0; gene < table.genearray.size(); gene++) {
	 	 		 	writeTranscriptToFile(table.genearray.get(gene), output);  	    		
	 	    	 }
	    		 output.close();
	    	 }	    	
	    	
	    	
			
		}
		}
		catch(Exception e) {
			ErrorLog.addError(e.getStackTrace());
			frame.getGlassPane().setCursor( Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
			table.setEnabled(true);
			e.printStackTrace();
			
			JOptionPane.showMessageDialog(Main.chromDraw, e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
		}
		 frame.getGlassPane().setVisible(false);
    	 frame.getGlassPane().setCursor( Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
    	 controlarray.clear();
    	 table.setEnabled(true);
	}
	/*
static void	writeVariantToTSVFile(VarNode node, BufferedWriter output) {
	 try {
		 if(Main.drawCanvas.hideNode(node)) {
			 return;
		 }
		 Entry<String, ArrayList<SampleNode>> entry;
		String genes = "", geneIDs = "", transIDs = "";
		 SampleNode varnode;
		  int mutcount = 0;
    	 String rscode, uniprot, strand, aminochange;
    	 Transcript.Exon exon;
    	 int casefreq = 0;
    	 String clusters;
    	 String genotype = "", biotype = ".";
    	 StringBuffer controls = new StringBuffer(""), tracks = new StringBuffer("");
    	 StringBuffer[] bedarray;
    	 HashMap<ControlFile, SampleNode> temphash = new HashMap<ControlFile, SampleNode>();
		    	    		
    		 if(node.isRscode() != null) {
    			 rscode = node.isRscode();
    		 }
    		 else {
    			 rscode = "N/A";
    		 }
    		 strand = "+";
    		 exon = null;    		
    		 strand = "N/A";    		     		
    		 uniprot = "N/A";							    		
    		 biotype = "N/A";
    		 if(node.getTranscripts() == null || node.getTranscripts().size() < 2) {
    			 System.out.println(Main.chromosomeDropdown.getSelectedItem() +":" +node.getPosition());
    		 }
    		 if(node.getTranscripts().get(0).equals(node.getTranscripts().get(1))) {
    			 if(node.getPosition() < node.getTranscripts().get(0).getStart()) {
    				 genes = "..." +node.getTranscripts().get(0).getGenename();
    				 geneIDs =   "..." +node.getTranscripts().get(0).getENSG();
        			 transIDs = "..." +node.getTranscripts().get(0).getENST() ;
    			 }
    			 else {
    				 genes = node.getTranscripts().get(0).getGenename() +"...";
    				 geneIDs =  node.getTranscripts().get(0).getENSG() +"...";
        			 transIDs = node.getTranscripts().get(0).getENST() +"...";
    			 }
    		 }
    		 else {
    			 genes = node.getTranscripts().get(0).getGenename() +"..." +node.getTranscripts().get(1).getGenename();
    			 geneIDs =  node.getTranscripts().get(0).getENSG() +"..." +node.getTranscripts().get(1).getENSG();
    			 transIDs = node.getTranscripts().get(0).getENST() +"..." +node.getTranscripts().get(1).getENST();
    		 }
    		 if(Main.bedCanvas.bedOn) {	    				 
				 	tracks = new StringBuffer("");	
				 	bedarray = MethodLibrary.makeTrackArray(node,null);
					
					for(int b = 0 ; b<bedarray.length; b++) {
							if(!Main.bedCanvas.bedTrack.get(b).intersect) {
								continue;
							}							
							if(bedarray[b] != null) {
								tracks.append(bedarray[b].toString()+"\t");							
							}							
						}
	    		 }
    		 clusters = "";
	    		
    		 if(commonSlider.getValue() > 1 && clusterSize > 0) {
    			 clusters = node.clusterNode.ID+"\t" +node.clusterNode.nodecount +"\t" +node.clusterNode.width +"\t" +MethodLibrary.round(node.clusterNode.nodecount/(double)node.clusterNode.width, 2) +"\t";
    			
    		 }
    		 for(int var = 0; var < node.vars.size(); var++) {
    			 
    			 entry = node.vars.get(var);   
    			 if(Main.drawCanvas.hideNodeVar(node, entry)) {
						continue;
				 }	
    			 aminochange = "N/A";
    			 mutcount = 0;    				    			
    			 
    			 if(Control.controlData.controlsOn && entry.getValue().size() > 0) {
				
						casefreq = 0;
						controls = new StringBuffer("");
						temphash.clear();
						for(int e = entry.getValue().size()-1; e> -1;e-- ) {
							
							if(entry.getValue().get(e).alleles == null) {
								if(entry.getValue().get(e).isHomozygous()) {
									casefreq+=2;
								}
								else {
									casefreq++;
								}										
							}
							else {
								if(!entry.getValue().get(e).getControlSample().controlOn) {
									continue;
								}
								temphash.put(entry.getValue().get(e).getControlSample(), entry.getValue().get(e));
						
							}									
						}		
					
						for(int i = 0 ; i<controlarray.size(); i++) {
							if(temphash.containsKey(controlarray.get(i))) {
								//AF
								controls.append(MethodLibrary.round(temphash.get(controlarray.get(i)).alleles/(double)temphash.get(controlarray.get(i)).allelenumber,5) +"\t");
								//OR
								controls.append(MethodLibrary.round((casefreq/(double)(Main.varsamples*2-casefreq))/(temphash.get(controlarray.get(i)).alleles/(double)(temphash.get(controlarray.get(i)).allelenumber-temphash.get(controlarray.get(i)).alleles)),5)+" (p=" +MethodLibrary.round(VariantHandler.table.fe.getRightTailedP(casefreq, Main.varsamples*2-casefreq, temphash.get(controlarray.get(i)).alleles, temphash.get(controlarray.get(i)).allelenumber-temphash.get(controlarray.get(i)).alleles) ,12) +")\t");
							}
							else {
								controls.append("N/A\tN/A\t");
							}
						}			
    			 }
    			
	    			 for(int i = 0; i<entry.getValue().size(); i++) {
	    				
	    				 if(Main.drawCanvas.hideVar(entry.getValue().get(i))) {
								continue;
						 }
	    				 
	    				 mutcount++;
	    			 }
    			 
	    			 if(mutcount > 0) {
	    				 
		    			for(int i = 0; i<entry.getValue().size(); i++) {
		    				 varnode = entry.getValue().get(i);
		    				 if(Main.drawCanvas.hideVar(varnode)) {
									continue;
							 }
		    				 if(onlyselected.isSelected()) {
		    					 if(!varnode.getSample().equals(Main.drawCanvas.selectedSample)) {
		    						 continue;
		    					 }
		    					 mutcount = 1;
		    				 }
		    				 if(varnode.isHomozygous()) {
		    					 genotype = "Hom(" +varnode.getCalls() +"/" +varnode.getCoverage() +")";
		    				 }
		    				 else {
		    					 genotype = "Het(" +varnode.getCalls() +"/" +varnode.getCoverage() +")";
		    				 } 				 
		    				
		    			
		    				 try {
		    					 if(output != null) {	    				
			    					 if(entry.getKey().length() > 1) {
			    						output.write(varnode.getSample().getName() +"\t" +
	 		    						genes +"\t" +mutcount +"\t" +"N/A" +"\t" +geneIDs +"\t" +transIDs +"\t" +"N/A" +"\t"+
	 		    						node.getTranscripts().get(0).getChrom() +":" +MethodLibrary.formatNumber((node.getPosition()+1)) +"\t" +strand	+"\t" +"Intergenic"+"\t" +"N/A" +"\t" +entry.getKey() +"\t" +genotype +"\t" +rscode +"\t" +clusters +controls +tracks +"N/A" +"\n");	    
	 		    					
			    					 }
			    					 else {
	    						 		output.write(varnode.getSample().getName() +"\t" +
			    						genes +"\t" +mutcount +"\t" +"N/A" +"\t" +geneIDs +"\t" +transIDs +"\t" +"N/A" +"\t"+
			    						node.getTranscripts().get(0).getChrom() +":" +MethodLibrary.formatNumber((node.getPosition()+1)) +"\t" +strand	+"\t" +"Intergenic"+"\t" +"N/A" +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotype +"\t" +rscode +"\t"  +clusters +controls +tracks +"N/A" +"\n");	    
			    					 }
		    					 }
		    					
		    				 }
		    				 catch(Exception ex) {
		    					 ex.printStackTrace();
		    					 ErrorLog.addError(ex.getStackTrace());
		    				 }
		    			 }
    			 }
    		 }
    		 
	 
    		    	
    	 Main.drawCanvas.loadBarSample = (int)((node.getPosition()/(double)Main.drawCanvas.splits.get(0).chromEnd)*100);     
    	 varnode = null;
		 node = null;
    	 
 	 }
	 catch(Exception exc) {
		 exc.printStackTrace();
		 ErrorLog.addError(exc.getStackTrace());
	 }	
}*/
public class QualEntry implements Serializable {

	private static final long serialVersionUID = 1L;
	String key;
	float value;
	String format;
	public QualEntry(String key, Float value, String format) {
		this.key = key;
		this.value = value;
		this.format = format;
	}
}
static void writeNodeToFile(VarNode node, String chrom, BufferedWriter output, BlockCompressedOutputStream outputgz) {
	try {
	if(vcf.isSelected()) {		
		 
		 StringBuffer info = new StringBuffer("");
		 StringBuffer alts = new StringBuffer("");
		 StringBuffer refs = new StringBuffer("");
		 StringBuffer ADs = new StringBuffer("");
		 StringBuffer sampleinfos = new StringBuffer("");
		 double avgquality = 0;
		 int samplecount = 0;
		 String rscode = ".";
		 String ref = Main.getBase.get(node.getRefBase());
		 int AN= Main.varsamples*2, AC =0, allelenro = 1;
		 StringBuffer ACs = new StringBuffer("");
		 StringBuffer AFs = new StringBuffer("");
		 double AF = 0.0;
		 Entry<String, ArrayList<SampleNode>> entry;
		 String sampleinfo = null;
		
		 int[] coverages = null;
		 HashMap<Short, String> samplehash = new HashMap<Short, String>();
	
		boolean set = false, found = false;
		 
		 for(int var = 0; var < node.vars.size(); var++) {
			 
			 entry = node.vars.get(var);
			 if(Main.drawCanvas.hideNodeVar(node, entry)) {
					continue;
			 }	
			 AF = 0;
			 AC = 0;			
			 
			 if( node.vars.size() == 1) {
				 if(node.indel && entry.getKey().length() > 1) {
					 String[] result = MethodLibrary.makeIndelColumns(chrom, node.getPosition(), ref, entry.getKey());
					 ref = result[0];
					 alts.append(result[1]+",");
					 
				 }
				 else {
					
				 	alts.append(entry.getKey()+",");
				 }
			 }
			 else {
				
				 if(node.indel) {
					 if(!set) {
						 String[] result = MethodLibrary.makeMultiAlt(chrom, node.getPosition(), ref, node);
						 set = true;
						 ref = result[0];
						 alts.append(result[1]);
					 }
					 allelenro = var+1;
					
				 }
				 else {
					 if(refs.length() == 0) {
						 refs.append(ref+",");
						
					 }
					 if(alts.length() > 0) {
						 allelenro++;
					 }
					 alts.append(entry.getKey() +",");
				 }				 
			 }
			//AC++;
			//if(node.vars.size() > 1) {
				coverages = new int[node.vars.size()+1];
				for(int i = 0; i<coverages.length; i++) {
					coverages[i] = 0;
				}
		//	}
			 for(int sample = 0; sample < entry.getValue().size(); sample++) {
				 ADs = new StringBuffer("");
				 if(!Main.drawCanvas.hideVar(entry.getValue().get(sample))) {
					 found = true;
					 if(samplehash.containsKey(entry.getValue().get(sample).getSample().getMultiIndex())) {
						 AC++;
						 coverages[0] = 0;
						 coverages[Integer.parseInt(""+samplehash.get(entry.getValue().get(sample).getSample().getMultiIndex()).charAt(2))] = entry.getValue().get(sample).getCoverage()-entry.getValue().get(sample).getCalls();
						 coverages[allelenro] = entry.getValue().get(sample).getCalls();
						 ADs.append(coverages[0]);
						 for(int i = 1; i<coverages.length; i++) {
							 ADs.append(","+coverages[i]);
						 }
						 sampleinfo = samplehash.get(entry.getValue().get(sample).getSample().getMultiIndex()).charAt(2) +"/" +allelenro +":"+entry.getValue().get(sample).getCoverage()+":" +ADs+":" +entry.getValue().get(sample).getGQ();
						 samplehash.put(entry.getValue().get(sample).getSample().getMultiIndex(),sampleinfo);						
					 }
					 else if(entry.getValue().get(sample).isHomozygous()) {
						 AC+=2;
						 coverages[0] = entry.getValue().get(sample).getCoverage()-entry.getValue().get(sample).getCalls();
						 coverages[allelenro] = entry.getValue().get(sample).getCalls();
						 ADs.append(coverages[0]);
						 for(int i = 1; i<coverages.length; i++) {
							 ADs.append(","+coverages[i]);
						 }
						 sampleinfo = (allelenro) +"/" +allelenro +":"+entry.getValue().get(sample).getCoverage()+":"+ADs+":" +entry.getValue().get(sample).getGQ();
						 
						 samplehash.put(entry.getValue().get(sample).getSample().getMultiIndex(),sampleinfo);		
					 }
					 else {
						 AC++;
						 coverages[0] = entry.getValue().get(sample).getCoverage()-entry.getValue().get(sample).getCalls();
						 coverages[allelenro] = entry.getValue().get(sample).getCalls();
						 ADs.append(coverages[0]);
						 for(int i = 1; i<coverages.length; i++) {
							 ADs.append(","+coverages[i]);
						 }
						
						 sampleinfo = "0/" +allelenro +":"+entry.getValue().get(sample).getCoverage()+":"+ ADs.toString() +":" +entry.getValue().get(sample).getGQ();
						 samplehash.put(entry.getValue().get(sample).getSample().getMultiIndex(),sampleinfo);	
						 
					 }
					 samplecount++;
					 avgquality+= entry.getValue().get(sample).getQuality();
				 }				 
			 }
			 AF = MethodLibrary.round(AC/(double)AN,5);
			 ACs.append(AC +",");
			 AFs.append(AF +",");
		 }		 
		
		 if(!found) {
			 return;
		 }
		 if(node.rscode != null) {
			 rscode = node.rscode;
		 }
		
		 for(Short i = 0 ; i <Main.varsamples; i++) {
			 if(samplehash.containsKey(i)) {
				 sampleinfos.append("\t" +samplehash.get(i));
			 }
			 else {
				 sampleinfos.append("\t" +"0/0");
			 }			
		 }
		 
		
		 alts.deleteCharAt(alts.length()-1);
		 AFs.deleteCharAt(AFs.length()-1);
		 ACs.deleteCharAt(ACs.length()-1);
		 info.append("AN=" +AN +";AC=" +ACs +";AF=" +AFs);
		
		// System.out.println(chrom +"\t" +(node.getPosition()+1) +"\t" +rscode +"\t" +refs +"\t" +alts +"\t" +MethodLibrary.round(avgquality/(double)samplecount,2) +"\t.\t" +info +"\t" +format +sampleinfos );
		
		 if(outputgz != null) {
			
			 String writeline = chrom +"\t" +(node.getPosition()+1) +"\t" +rscode +"\t" +ref +"\t" +alts +"\t" +MethodLibrary.round(avgquality/(double)samplecount,2) +"\tPASS\t" +info +"\t" +format +sampleinfos+"\n";
			 outputgz.write(writeline.getBytes());				 
			 Feature vcf = vcfCodec.decode(writeline);				
			 FileRead.indexCreator.addFeature(vcf, FileRead.filepointer);
			 FileRead.filepointer = outputgz.getFilePointer();
			 /*
			if(FileRead.bigcalc) {
				 outputStrings.add(writeline);
				 if(node.getPosition() - VariantHandler.lastWrittenPos > Settings.windowSize) {
					
					 for(int i = 0 ; i<outputStrings.size(); i++) {
						 outputgz.write(outputStrings.get(i).getBytes());
						 
						 Feature vcf = vcfCodec.decode(outputStrings.get(i));
						
						 FileRead.indexCreator.addFeature(vcf, FileRead.filepointer);
						 FileRead.filepointer = outputgz.getFilePointer();
					 }
					 outputStrings.clear();
					 VariantHandler.lastWrittenPos = node.getPosition();
				 }			
			}
			else {
				 outputgz.write(writeline.getBytes());				 
				 Feature vcf = vcfCodec.decode(writeline);				
				 FileRead.indexCreator.addFeature(vcf, FileRead.filepointer);
				 FileRead.filepointer = outputgz.getFilePointer();
			}*/
		 }
		 else if(output != null) {
			 output.write(chrom +"\t" +(node.getPosition()+1) +"\t" +rscode +"\t" +ref +"\t" +alts +"\t" +MethodLibrary.round(avgquality/(double)samplecount,2) +"\t.\t" +info +"\t" +format +sampleinfos+"\n" );
		 }
		 else {
		//	 System.out.println(chrom +"\t" +(node.getPosition()+1) +"\t" +rscode +"\t" +ref +"\t" +alts +"\t" +MethodLibrary.round(avgquality/(double)samplecount,2) +"\t.\t" +info +"\t" +format +sampleinfos);
		 }
	}
	 }
	 catch(Exception e) {
		 e.printStackTrace();
	 }
}

static void	writeTranscriptToFile(Gene gene, BufferedWriter output) {
		 try {			 
			
			 Entry<String, ArrayList<SampleNode>> entry;
			 VarNode node = null;
			 SampleNode varnode;
			 String[] row;    	
			 StringBuffer[] bedarray;
	    	 String rscode,aminochange, transcripts, exons;	    	
	    	 StringBuffer strand = new StringBuffer("");
	    	 int casefreq = 0;
	    	 String genotype = "", biotype = ".", geneID = "", description = "";
	    	 StringBuffer controls = new StringBuffer(""), tracks = new StringBuffer("");
	    	 String clusters = "";    	
	    	 HashMap<ControlFile, SampleNode> temphash = new HashMap<ControlFile, SampleNode>();
    		 table.getAminos(gene);
    		 StringBuffer samples = new StringBuffer("");
    		 StringBuffer genotypes = new StringBuffer("");
    		 StringBuffer qualities = new StringBuffer("");
    		 
	    	 for(int s = 0; s<table.aminoarray.size(); s++) {
	    		 row = table.aminoarray.get(s).getRow();
	    		 node = table.aminoarray.get(s).getNode();
	    		 strand = new StringBuffer("");
	    		 if(node.isRscode() != null) {
	    			 rscode = node.isRscode();
	    		 }
	    		 else {
	    			 rscode = "N/A";
	    		 }
	    		
	    		 if(row[8].contains("troni")) {
	    			 exons = "Intronic";
	    		 }
	    		 else if(row[8].contains("geni")) {
	    			 exons = "Intergenic";
	    		 }
	    		 else {
	    			 if(row[8].length() == 1) {
		    			 exons = "Exon " +row[8];
		    		 }
		    		 else if(row[8].length() > 1) {
		    			 exons = "Exons " +row[8];
		    		 }
		    		 else {
		    			 exons = "";
		    		 }
	    		 }
	    		 if(gene.intergenic) {
	    			 if(node.getTranscripts().size() == 2) {
	    				 geneID = gene.getID() +";" +node.getTranscripts().get(1).getGene().getID();
	    				 description = gene.getDescription() + ";" +node.getTranscripts().get(1).getGene().getDescription();	    				 	    				 
	    				 strand.append((gene.getStrand()) ? "+;" : "-;");
	    				 strand.append((node.getTranscripts().get(1).getGene().getStrand()) ? "+" : "-");
	    			 }
	    			 else  {
	    				 geneID = gene.getID();
	    				 description = gene.getDescription();
	    			 }
	    		 }
	    		 else {
		    		 if(!gene.getStrand()) {
		    			 strand.append("-");
		    		 }
		    		 else {
		    			 strand.append("+");
		    		 }
		    		 geneID = gene.getID();
		    		 description = gene.getDescription();
	    		 }
	    		 transcripts = row[6];	    		
	    		 biotype = row[7];
	    		
	    		
	    		 clusters = "";
	    		
	    		 if(commonSlider.getValue() > 1 && clusterSize > 0) {	    		
	    			if(node.clusterNode != null) {	    			
	    			 clusters = node.clusterNode.ID+"\t" +node.clusterNode.nodecount +"\t" +node.clusterNode.width +"\t" +MethodLibrary.round(node.clusterNode.nodecount/(double)node.clusterNode.width, 2) +"\t";
	    			}
	    			else {
	    				//System.out.println(node.getPosition());
	    				continue;
	    			}
	    		 }
	    		
	    		else {
	    		 for(int var = 0; var < node.vars.size(); var++) {
	    			 
	    			 entry = node.vars.get(var);
	    			 if(Main.drawCanvas.hideNodeVar(node, entry)) {
							continue;
					 }	
	    			 if(!entry.getKey().equals(row[5])) {
	    				 continue;
	    			 }
	    			
	    			 if(exons.length() > 0) {
	    				 aminochange = row[3];	    				
	    			 }
	    			 else {
	    				 aminochange = "N/A";
	    			 }
	    			
	    			 if(Main.bedCanvas.bedOn) {
	 	    			
	  				 	tracks = new StringBuffer("");	
	  				 	
	  				 	bedarray = MethodLibrary.makeTrackArray(node,entry.getKey());
	 					if(bedarray != null) {
	 						for(int b = 0 ; b<bedarray.length; b++) {
	 							if(!Main.bedCanvas.bedTrack.get(b).intersect) {
	 								continue;
	 							}							
	 							if(bedarray[b] != null) {
	 								tracks.append(bedarray[b].toString()+"\t");							
	 							}
	 							else {
	 								tracks.append("-\t");		
	 							}
	 						}
	 					}
	 					else {
	 						for(int i = 0 ; i<Main.bedCanvas.bedTrack.size(); i++) {
	 							if(Main.bedCanvas.bedTrack.get(i).intersect) {
	 								tracks.append("\t");
	 							}							
	 						}
	 					}
	 	    		 }	    			
	    			 
	    			 if(Control.controlData.controlsOn && entry.getValue().size() > 0) {
					
							casefreq = 0;
							controls = new StringBuffer("");
							temphash.clear();
							for(int e = entry.getValue().size()-1; e> -1;e-- ) {
								
								if(entry.getValue().get(e).alleles == null) {
									if(entry.getValue().get(e).isHomozygous()) {
										casefreq+=2;
									}
									else {
										casefreq++;
									}										
								}
								else {
									if(!entry.getValue().get(e).getControlSample().controlOn) {
										continue;
									}
									temphash.put(entry.getValue().get(e).getControlSample(), entry.getValue().get(e));
							
								}									
							}		
							for(int i = 0 ; i<controlarray.size(); i++) {
								if(temphash.containsKey(controlarray.get(i))) {
									//AF
									controls.append(MethodLibrary.round(temphash.get(controlarray.get(i)).alleles/(double)temphash.get(controlarray.get(i)).allelenumber,5) +"\t");
									//OR
									controls.append(MethodLibrary.round((casefreq/(double)(Main.varsamples*2-casefreq))/(temphash.get(controlarray.get(i)).alleles/(double)(temphash.get(controlarray.get(i)).allelenumber-temphash.get(controlarray.get(i)).alleles)),5)+" (p=" +MethodLibrary.round(VariantHandler.table.fe.getRightTailedP(casefreq, Main.varsamples*2-casefreq, temphash.get(controlarray.get(i)).alleles, temphash.get(controlarray.get(i)).allelenumber-temphash.get(controlarray.get(i)).alleles) ,12) +")\t");
								}
								else {
									controls.append("N/A\tN/A\t");
								}
							}				
	    			 }
    			if(compactTsv.isSelected()) {
    				 samples = new StringBuffer("");
    	    		 genotypes = new StringBuffer("");
    	    		 qualities = new StringBuffer("");
    	    		 
    	    		 for(int i = 0; i<entry.getValue().size(); i++) {
 	    				
	    				 if(Main.drawCanvas.hideVar(entry.getValue().get(i))) {
							continue;
						 }
	    				 
	    				 varnode = entry.getValue().get(i);
	    				
	    				 if(varnode.isHomozygous()) {
	    					 genotypes.append("Hom(" +varnode.getCalls() +"/" +varnode.getCoverage() +");");
	    				 }
	    				 else {
	    					 genotypes.append("Het(" +varnode.getCalls() +"/" +varnode.getCoverage() +");");
	    				 } 		
	    				 samples.append(varnode.getSample().getName() +";");
	    				 qualities.append(varnode.getQuality()+";");
    	    		 }
    	    		 
    	    		genotypes.deleteCharAt(genotypes.length()-1);
    	    		qualities.deleteCharAt(qualities.length()-1);
    	    		samples.deleteCharAt(samples.length()-1);
	    	    		
    	    		 if(output != null) {
    	    			
		    				 if(exons.length() > 0) {
		    					 if(!aminochange.equals("N/A")) {
		    						 output.write(samples +"\t" +
		    						 row[0] +"\t" +row[1] +"\t" +gene.samples.size() +"\t" +geneID +"\t" +transcripts +"\t" +biotype +"\t"+
						  		  	 row[2] +"\t" +strand	+"\t" +exons +"\t" +aminochange +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotypes +"\t" +qualities +"\t"+rscode +"\t" +clusters +controls +tracks +description +"\n");
		    					 }
		    					 else {
		    						 output.write(samples +"\t" +
		    						 row[0] +"\t" +row[1] +"\t" +gene.samples.size() +"\t" +geneID +"\t" +transcripts +"\t" +biotype +"\t"+
		    						 row[2] +"\t" +strand	+"\t" +"UTR" +"\t" +aminochange +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotypes +"\t"+qualities +"\t" +rscode +"\t"+clusters+controls +tracks +description +"\n");
		    					 }
		    				 }
		    				 else {
		    					 if(node.getTranscripts() != null && node.isInGene()) {
		    						 output.write(samples +"\t" +
		    						 row[0] +"\t" +row[1] +"\t" +gene.samples.size() +"\t" +geneID +"\t" +transcripts +"\t" +biotype +"\t"+
		    						 row[2] +"\t" +strand	+"\t" +"Intronic"+"\t" +"N/A" +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotypes +"\t"+qualities +"\t" +rscode+"\t" +clusters+controls +tracks +description +"\n");	    						 
		    					 }
		    					 else {
		    						 output.write(samples +"\t" +
		    						 row[0] +"\t" +row[1] +"\t" +gene.samples.size() +"\t" +geneID +"\t" +transcripts +"\t" +biotype +"\t"+
		    						 row[2] +"\t" +strand	+"\t" +"Intergenic"+"\t" +"N/A" +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotypes +"\t" +qualities +"\t"+rscode+"\t" +clusters +controls +tracks +description +"\n");	    
		    					 }
		    				 }
    	    			 }	    	    		
	    	    		
	    			}
	    			
	    			else if(tsv.isSelected()) {
	    			 for(int i = 0; i<entry.getValue().size(); i++) {
	    				
	    				 if(Main.drawCanvas.hideVar(entry.getValue().get(i))) {
							continue;
						 }
	    				 
	    				 varnode = entry.getValue().get(i);
	    				
	    				 if(varnode.isHomozygous()) {
	    					 genotype = "Hom(" +varnode.getCalls() +"/" +varnode.getCoverage() +")";
	    				 }
	    				 else {
	    					 genotype = "Het(" +varnode.getCalls() +"/" +varnode.getCoverage() +")";
	    				 } 				 	    				
	    				
	    				 try {
	    					 if(output != null) {
			    				 if(exons.length() > 0) {
			    					 if(!aminochange.equals("N/A")) {
			    						 output.write(varnode.getSample().getName() +"\t" +
			    						 row[0] +"\t" +row[1] +"\t" +gene.samples.size() +"\t" +geneID +"\t" +transcripts +"\t" +biotype +"\t"+
							  		  	 row[2] +"\t" +strand	+"\t" +exons +"\t" +aminochange +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotype +"\t" +varnode.getQuality() +"\t"+rscode +"\t" +clusters +controls +tracks +description +"\n");
			    					 }
			    					 else {
			    						 output.write(varnode.getSample().getName() +"\t" +
			    						 row[0] +"\t" +row[1] +"\t" +gene.samples.size() +"\t" +geneID +"\t" +transcripts +"\t" +biotype +"\t"+
			    						 row[2] +"\t" +strand	+"\t" +"UTR" +"\t" +aminochange +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotype +"\t"+varnode.getQuality() +"\t" +rscode +"\t"+clusters+controls +tracks +description +"\n");
			    					 }
			    				 }
			    				 else {
			    					 if(node.getTranscripts() != null && node.isInGene()) {
			    						 output.write(varnode.getSample().getName() +"\t" +
			    						 row[0] +"\t" +row[1] +"\t" +gene.samples.size() +"\t" +geneID +"\t" +transcripts +"\t" +biotype +"\t"+
			    						 row[2] +"\t" +strand	+"\t" +"Intronic"+"\t" +"N/A" +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotype +"\t"+varnode.getQuality() +"\t" +rscode+"\t" +clusters+controls +tracks +description +"\n");	    						 
			    					 }
			    					 else {
			    						 output.write(varnode.getSample().getName() +"\t" +
			    						 row[0] +"\t" +row[1] +"\t" +gene.samples.size() +"\t" +geneID +"\t" +transcripts +"\t" +biotype +"\t"+
			    						 row[2] +"\t" +strand	+"\t" +"Intergenic"+"\t" +"N/A" +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotype +"\t" +varnode.getQuality() +"\t"+rscode+"\t" +clusters +controls +tracks +description +"\n");	    
			    					 }
			    				 }
	    					 }
	    					 else {
	    						/* if(exon != null) {
			    					 if(!aminochange.equals("N/A")) {
				    				 System.out.print(varnode.getSample().getName() +"\t" +
							  		  row[0] +"\t" +row[1] +"\t" +transcript.samples.size() +"\t" +uniprot +"\t" +transcript.getENSG() +"\t" +transcript.getENST() +"\t" +biotype +"\t"+
							  		  "chr" +row[2] +"\t" +strand	+"\t" +"Exon " +exon.getNro()+"\t" +aminochange +"\t" +node.getRefBase() +"->" +entry.getKey() +"\t" +genotype +"\t" +rscode +controls +transcript.getDescription() +"\n");
			    					 }
			    					 else {
			    						 System.out.print(varnode.getSample().getName() +"\t" +
			    						  		  row[0] +"\t" +row[1] +"\t" +transcript.samples.size() +"\t" +uniprot +"\t" +transcript.getENSG() +"\t" +transcript.getENST() +"\t" +biotype +"\t"+
			    						  		"chr" +row[2] +"\t" +strand	+"\t" +"UTR" +"\t" +aminochange +"\t" +node.getRefBase() +"->" +entry.getKey() +"\t" +genotype +"\t" +rscode +controls +transcript.getDescription() +"\n");
			    					 }
			    				 }
			    				 else {
			    					 if(node.getTranscripts().size() > 0) {
			    						 System.out.print(varnode.getSample().getName() +"\t" +
			    						  		  row[0] +"\t" +row[1] +"\t" +transcript.samples.size() +"\t" +uniprot +"\t" +transcript.getENSG() +"\t" +transcript.getENST() +"\t" +biotype +"\t"+
			    						  		"chr" +row[2] +"\t" +strand	+"\t" +"Intronic"+"\t" +"N/A" +"\t" +node.getRefBase() +"->" +entry.getKey() +"\t" +genotype +"\t" +rscode +controls +transcript.getDescription() +"\n");	    						 
			    					 }
			    					 else {
			    						 System.out.print(varnode.getSample().getName() +"\t" +
			    						  		  row[0] +"\t" +row[1] +"\t" +transcript.samples.size() +"\t" +uniprot +"\t" +transcript.getENSG() +"\t" +transcript.getENST() +"\t" +biotype +"\t"+
			    						  		"chr" +row[2] +"\t" +strand	+"\t" +"Intergenic"+"\t" +"N/A" +"\t" +node.getRefBase() +"->" +entry.getKey() +"\t" +genotype +"\t" +rscode  +controls +transcript.getDescription() +"\n");	    
			    					 }
			    				 }*/
	    					 }
	    				 }
	    				 catch(Exception ex) {
	    					 ex.printStackTrace();
	    					 ErrorLog.addError(ex.getStackTrace());
	    				 }
	    			 }	    			
	    		 }    	 
	    		 }
	    		}
	    	 }					    	
	    	 Main.drawCanvas.loadBarSample = (int)((gene.getStart()/(double)Main.drawCanvas.splits.get(0).chromEnd)*100);     
	    	 varnode = null;
			 node = null;
	 	}
	 catch(Exception exc) {
		 JOptionPane.showMessageDialog(Main.chromDraw, exc.getMessage() +" - there were problem writing variants for " +gene.getName(), "Error", JOptionPane.ERROR_MESSAGE);
		 exc.printStackTrace();
		 ErrorLog.addError(exc.getStackTrace());
	 }
		
	}

void writeGeneListToVCF(ArrayList<Gene> genelist, BufferedWriter output, BlockCompressedOutputStream outputgz) {
	boolean found = false;
	VarNode node = VariantHandler.table.genearray.get(0).varnodes.get(0);
	while(node != null) {
		
		if(geneSlider.getValue() > 1) {
			found = false;
			if(node.isInGene()) {
				if(node.getExons() != null) {
					for(int i = 0 ; i<node.getExons().size(); i++) {
						if(node.getExons().get(i).getTranscript().getGene().samples.size() >= geneSlider.getValue()) {							
							found = true;
						}
						else {
							node.getExons().remove(i);
							i--;
						}
					}
				}
				else {
					for(int i = 0 ; i<node.getTranscripts().size(); i++) {
						if(node.getTranscripts().get(i).getGene().samples.size() >= geneSlider.getValue()) {
							found = true;
						}
						else {
							node.getTranscripts().remove(i);
							i--;
						}
					}
				}	
			}
			if(found) {
				writeNodeToFile(node,node.getChrom(), null, outputgz);
			}
		}
		else {
			writeNodeToFile(node,node.getChrom(), null, outputgz);
		}
		 node = node.getNext();
	}
	node = null;
	
}

	public class OutputRunner extends SwingWorker<String, Object> {
		BufferedWriter output;
		BlockCompressedOutputStream outputgz;
		File outFile;
		public OutputRunner(BufferedWriter output, BlockCompressedOutputStream outputgz, File outFile) {
			this.output = output;
			this.outputgz = outputgz;
			this.outFile = outFile;
		}
		
		protected String doInBackground() {
			Main.drawCanvas.loading("Writing output...");
			writeOutput(output, outputgz, outFile);
			Main.drawCanvas.ready("Writing output...");
			return "";
		}
		
	}
	@Override
	public void mouseClicked(MouseEvent arg0) {
		
	}

	@Override
	public void mouseEntered(MouseEvent arg0) {
		
	}

	@Override
	public void mouseExited(MouseEvent arg0) {
		
	}

	@Override
	public void mousePressed(MouseEvent event) {
		if(event.getSource() == tabs) {			
			if(tabs.getSelectedIndex() == 0) {
				VariantHandler.aminoCount.setText(table.variants +" variants");
				outputmenu.setText("Variant output");
				outputmenu.revalidate();
			}
			else if(tabs.getSelectedIndex() == tabs.indexOfComponent(statsScroll)) {
				VariantHandler.aminoCount.setText(stattable.variants +" variants");
				outputmenu.setText("Stats output");
				outputmenu.revalidate();
			}
			else if(tabs.getSelectedIndex() == tabs.indexOfComponent(clusterScroll)) {
				VariantHandler.aminoCount.setText(clusterTable.variants +" variants");
				outputmenu.setText("Variant output");
				outputmenu.revalidate();
			}
			else {
				VariantHandler.aminoCount.setText(tables.get(tabs.getSelectedIndex()-(tabs.indexOfComponent(statsScroll)+1)).variants +" variants");
				outputmenu.setText("Variant output");
				outputmenu.revalidate();
			}
		}		
	}

	@Override
	public void mouseReleased(MouseEvent event) {
		if(event.getSource() == tableScroll.getVerticalScrollBar()) {
			table.repaint();			
			return;
		}
		else if(event.getSource() == coverageSlider) {
			if(coverageSlider.getValue() == coverageSlider.getMaximum()) {
				coverageSlider.setMaximum(coverageSlider.getMaximum()*2);
			}
			if(commonSlider.getValue() > 1) {
				
				
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}
		}
		else if(event.getSource() == qualitySlider) {
			if(qualitySlider.getValue() == qualitySlider.getMaximum()) {
				qualitySlider.setMaximum(qualitySlider.getMaximum()*2);
			}
			if(commonSlider.getValue() > 1) {
				
				
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}
		}
		else if(event.getSource() == gqSlider) {
			if(gqSlider.getValue() == gqSlider.getMaximum()) {
				gqSlider.setMaximum(gqSlider.getMaximum()*2);
			}
			if(commonSlider.getValue() > 1) {
				
				
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}
		}
		else if(event.getSource() == maxCoverageSlider) {
			if(maxCoverageSlider.getValue() == maxCoverageSlider.getMaximum()) {
				maxCoverageSlider.setMaximum(maxCoverageSlider.getMaximum()*2);
			}
			if(commonSlider.getValue() > 1) {
			
				
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}
		}
		else if(event.getSource() == callSlider) {
			if(commonSlider.getValue() > 1) {
				
				
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}
		}
	
	}
	void freezeFilters(boolean value) {
		
		if(value) {
			for(int i =2;i<filterpanel.getComponentCount(); i++) {
				if(filterpanel.getComponent(i).equals(freeze)) {
					continue;
				}
				filterpanel.getComponent(i).setEnabled(false);
				filterpanel.getComponent(i).revalidate();
			}
		}
		else {
			for(int i =0;i<filterpanel.getComponentCount(); i++) {
				filterpanel.getComponent(i).setEnabled(true);
				filterpanel.getComponent(i).revalidate();
			}
		}
		
	}
	
	/*public class Writer extends SwingWorker< String, Object > {
		
		String type;
		
		protected String doInBackground() throws Exception {
			
			
			
			return null;
		}
		
		public Writer(String type) {
			this.type = type;
		}
	}*/
	@Override
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void keyPressed(KeyEvent e) {
		
		if(e.getKeyCode() == KeyEvent.VK_ENTER) {
			
			if(e.getSource() == clusterBox) {
				if(commonSlider.getValue() > 1) {
					clusterSize = Integer.parseInt(clusterBox.getText());
				
					Main.drawCanvas.calcClusters(FileRead.head,1);
					if(tabs.indexOfComponent(clusterScroll) == -1) {						
						tabs.add(clusterScroll, tabs.indexOfComponent(statsScroll));
					}
				}
				if(clusterSize == 0) {
					if(tabs.indexOfComponent(clusterScroll) != -1) {						
						tabs.remove(clusterScroll);
					}
				}
			}
		}
		
	}

	@Override
	public void keyReleased(KeyEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void componentResized(ComponentEvent e) {
		// TODO Auto-generated method stub
		table.resizeTable();
		for(int i = 0 ; i<tables.size(); i++) {
			tables.get(i).resizeTable();
		}
		stattable.resizeTable();
		clusterTable.resizeTable();
	}

	@Override
	public void componentMoved(ComponentEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void componentShown(ComponentEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void componentHidden(ComponentEvent e) {
		// TODO Auto-generated method stub
		
	}

	
	
	
	
	
}
