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
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
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
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;

import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
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
import javax.swing.UIManager;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class VariantHandler extends JPanel implements ChangeListener, ActionListener, MouseListener, KeyListener, ComponentListener {
	
	private static final long serialVersionUID = 1L;
	static RangeSlider commonSlider = new RangeSlider(1,1);
	static RangeSlider callSlider = new RangeSlider(0,100);
	static RangeSlider callSliderIndel = new RangeSlider(0,100);
	static JSlider qualitySlider = new JSlider(0,60);
	static JSlider gqSlider = new JSlider(0,60);
	static JSlider coverageSlider = new JSlider(1,40);
	static JSlider maxCoverageSlider = new JSlider(1,2000);
//	static JSlider callSlider = new JSlider(0,100);
	static JSlider geneSlider = new JSlider(1,1);
	static JLabel geneLabel = new JLabel("At least 1/1 samples share a mutated gene");	
	static JLabel aminoCount = new JLabel("");
	static JLabel callsLabel = new JLabel();
	static JLabel coverageLabel = new JLabel();
	static JLabel maxCoverageLabel = new JLabel();
	static JLabel qualityLabel = new JLabel(), gqLabel = new JLabel(),comparison = new JLabel("Sample comparison");
	static JLabel callsLabelIndel = new JLabel();
	static JLabel coverageLabelIndel = new JLabel();
	static JLabel maxCoverageLabelIndel = new JLabel();
	static JLabel qualityLabelIndel = new JLabel(), gqLabelIndel = new JLabel();
	static JLabel slideLabel = new JLabel("Common variants in 1/1 samples");	
	static JLabel clusterLabel = new JLabel("Window size for variant clusters");
	
	
	static JSlider qualitySliderIndel = new JSlider(0,60);
	static JSlider gqSliderIndel = new JSlider(0,60);
	static JSlider coverageSliderIndel = new JSlider(1,40);
	static JSlider maxCoverageSliderIndel = new JSlider(1,2000);
	//static JSlider callSliderIndel = new JSlider(0,100);		
	
	
//	static JMenuBar menubar = new JMenuBar();
	static JMenu filters;
	
	static JMenuBar aminobar;
	static JMenu aminomenu;
	static JMenu outputmenu;
	static JMenuItem write;
	
	static JTextField clusterBox;
	static ButtonGroup outputgroup;
	static JRadioButton tsv;
	static JRadioButton compactTsv;
	static JRadioButton vcf;
	static JRadioButton oncodrive;
	static JCheckBox hidenoncoding;
	static JCheckBox freeze;
	static JCheckBox rscode;
	static JCheckBox allChroms;
	static JCheckBox allChromsfrom;
	static JCheckBox onlyAutosomes;
	static JCheckBox hideSNVs;
	static JCheckBox hideHomos;
	static JCheckBox onlyStats;
	static JCheckBox hideIndels;
	static JCheckBox synonymous;
	static JCheckBox nonsense;
	static JCheckBox intronic;
	static JCheckBox intergenic;
	static JCheckBox utr;
	static JCheckBox onlyselected;
	static JCheckBox writetofile;
	static JCheckBox indelFilters;
	static JLabel SNVFilters;
	static JFrame frame;    
	static JLabel totalVars;
	static JLabel totalVarsIndel;
	static JLabel empty;	
	static JButton varcalc;
	static JButton statcalc;	
	static ArrayList<ControlFile> controlarray = new ArrayList<ControlFile>();
	String userDir;
	static HashMap<String, Integer> variantSettings;
	static int lastWrittenPos = 0;
	static  ArrayList<String> outputStrings = new ArrayList<String>();
	static MouseWheelListener sliderWheelListener;
	static JTabbedPane tabs;
	static JScrollPane tableScroll = new JScrollPane();
	static JSeparator separator;
	static AminoTable table;
	static StatsTable stattable;
	static JScrollPane statsScroll = new JScrollPane();
	static ClusterTable clusterTable;	
	static JScrollPane clusterScroll = new JScrollPane();
	static ArrayList<JScrollPane> tablescrolls = new ArrayList<JScrollPane>();
	static ArrayList<BedTable> tables = new ArrayList<BedTable>();
	static NodeSorter nodesorter = new NodeSorter();
	static JPopupMenu menu;
	static JPopupMenu menuIndel;
	static JScrollPane menuScroll;
	static JScrollPane menuScrollIndel;
	static JPanel menuPanel;
	static JPanel menuPanelIndel;
	static JButton applyQualities;
	static JButton applyQualitiesIndel;
	static JButton advQualities;
	static JButton advQualitiesIndel;
	static OwnVCFCodec vcfCodec= new OwnVCFCodec();
	static String format = "GT:DP:AD:GQ";
	int moveX=0, moveY=0, pressX=0,pressY=0;
	static JLabel adder, adder2, adder3, adder4, adder5;
	//final int buttonHeight = 15, buttonWidth = 40;
	int buttonHeight = (int)(11);
	int buttonWidth = 12*6;
	static float maxSlideValue = 0;
	Dimension buttondimension = new Dimension(buttonWidth, buttonHeight);		
	//final Dimension buttondimension = new Dimension(Main.buttonWidth, Main.buttonHeight);
	static Color backColor = new Color(228,228,218,255);
	static Color frameColor = new Color(188,188,178,255);
	static int clusterSize = 0;
	static JTabbedPane filterPanes = new JTabbedPane();
	static JPanel filterpanel;
	static JPanel filterpanelIndel;
	static JPanel hidepanel;
	static JPanel comparepanel;	
	static JPanel aminopanel;

	public VariantHandler() {	
		super(new GridBagLayout());		 
		
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
		createButtons();
		GridBagConstraints c = new GridBagConstraints();	
		c.anchor = GridBagConstraints.NORTHWEST;
		
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
	
		
		tabs.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);
		filterPanes.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);
	//	filterpanelIndel.setLayout(new GridLayout(7, 2));
	//	hidepanel.setLayout(new GridBagLayout());
	//	comparepanel.setLayout(new GridLayout(4, 2));
		aminopanel.setLayout(new GridLayout(1,3));
		
		//tableScroll.setPreferredSize(new Dimension(500,400));
		
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
		filterPanes.addKeyListener(this);
//		add(menubar,c);
		c.gridy = 1;
		filterPanes.add("SNVs",filterpanel);
		filterPanes.add("Indels",filterpanelIndel);
		filterPanes.add("Hide", hidepanel);
		
		filterPanes.add("Compare", comparepanel);
	
		add(filterPanes, c);		
		c.gridy = 2;
		//add(comparepanel,c);
	//	c.gridy = 3;
		c.fill = GridBagConstraints.HORIZONTAL;
		add(aminopanel, c);
		c.gridy = 3;
		c.weightx = 1;
		c.weighty = 0.8;
		tabs.addMouseListener(this);
		tabs.setBackground(backColor);
		tabs.add("Genes", tableScroll);	
		tabs.add("Stats", statsScroll);
		
		//tabs.add("Clusters", clusterScroll);		
		c.fill = GridBagConstraints.BOTH;
		
		add(tabs, c);
		if(Main.menuFont == null) {
			Main.menuFont = new Font("SansSerif", Font.PLAIN, 12);
		}
		menuPanel.add(new JLabel("Hard filters"));
		menuPanelIndel.add(new JLabel("Hard filters"));
		setFonts(Main.menuFont);
		
	}
	
	 protected void paintComponent(Graphics g) {
	        super.paintComponent(g);     
	        
	     //   g.drawImage(image, 0, 0, this.getWidth(), this.getHeight(), null);
	        g.setColor(frameColor);
	        g.fillRect(0, 0, this.getWidth(), this.getHeight());	        
	       
	    }
	 
	void createButtons() {
		freeze = new JCheckBox("Freeze filters");
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
		 commonSlider = new RangeSlider(1,1);
	     callSlider = new RangeSlider(0,100);
		 callSliderIndel = new RangeSlider(0,100);
		 qualitySlider = new JSlider(0,60);
		 gqSlider = new JSlider(0,60);
		 coverageSlider = new JSlider(1,40);
		 maxCoverageSlider = new JSlider(1,2000);
		 geneSlider = new JSlider(1,1);
		 geneLabel = new JLabel("At least 1/1 samples share a mutated gene");	
		 aminoCount = new JLabel("");
		 callsLabel = new JLabel();
		 coverageLabel = new JLabel();
		 maxCoverageLabel = new JLabel();
		 qualityLabel = new JLabel();
		 gqLabel = new JLabel();
		 comparison = new JLabel("Sample comparison");
		 callsLabelIndel = new JLabel();
		 coverageLabelIndel = new JLabel();
		 maxCoverageLabelIndel = new JLabel();
		 qualityLabelIndel = new JLabel();
		 gqLabelIndel = new JLabel();
		 slideLabel = new JLabel("Common variants in 1/1 samples");	
		 clusterLabel = new JLabel("Window size for variant clusters");
		
		
		 qualitySliderIndel = new JSlider(0,60);
		 gqSliderIndel = new JSlider(0,60);
		 coverageSliderIndel = new JSlider(1,40);
	
		 filters = new JMenu("Variant Filters");
		
		 aminobar = new JMenuBar();
		 aminomenu = new JMenu("Options");
		 outputmenu = new JMenu("Variant output");
		 write = new JMenuItem("Save");
		
		 clusterBox = new JTextField("0");
		 outputgroup = new ButtonGroup();
		  tsv = new JRadioButton("TSV");
		  compactTsv = new JRadioButton("Compact TSV");
		  vcf = new JRadioButton("VCF");
		  oncodrive = new JRadioButton("Oncodrive");
		hidenoncoding = new JCheckBox("Hide non-coding variants");
		
		rscode = new JCheckBox("Hide rs-coded variants");
		allChroms = new JCheckBox("All chromosomes");
		allChromsfrom = new JCheckBox("From this chr?");
		onlyAutosomes = new JCheckBox("Only autosomes?");
		hideSNVs = new JCheckBox("Hide SNVs");
		hideHomos = new JCheckBox("Hide homozygotes");
		onlyStats = new JCheckBox("Only stats");
		hideIndels = new JCheckBox("Hide indels");
		synonymous = new JCheckBox("Only non-synonymous");
		nonsense = new JCheckBox("Only truncs");
		intronic = new JCheckBox("Show intronic");
		intergenic = new JCheckBox("Show intergenic");
		utr = new JCheckBox("Show UTR");
		onlyselected = new JCheckBox("Only selected sample");
		writetofile = new JCheckBox("Write directly to a file");
		indelFilters = new JCheckBox("Use indel specific filters");
		SNVFilters = new JLabel("SNV & indel filters");
		SNVFilters.setName("header");
		indelFilters.setName("header");
		comparison.setName("header");
		totalVars = new JLabel("Variant count on screen: 0");
		totalVarsIndel = new JLabel("Variant count on screen: 0");
		empty = new JLabel("");	
		varcalc = new JButton("Annotate");
		statcalc = new JButton("Stats");	 
		tabs = new JTabbedPane();
		
		separator = new JSeparator();	
		 menu = new JPopupMenu("Advanced quality control");
		 menuIndel = new JPopupMenu("Advanced quality control");		 
		 menuPanel = new JPanel(new GridBagLayout());
		 menuPanelIndel = new JPanel(new GridBagLayout());
		 applyQualities = new JButton("Apply");
		 applyQualitiesIndel = new JButton("Apply");
		 advQualities = new JButton("Hard filters");
		 advQualitiesIndel = new JButton("Hard filters");
		
		 filterPanes = new JTabbedPane();
			filterpanel = new JPanel(new GridBagLayout()) {				
				private static final long serialVersionUID = 1L;

				protected void paintComponent(Graphics g) {
				        super.paintComponent(g);     
				  
				        g.setColor(backColor);
				        g.fillRect(0, 0, this.getWidth(), this.getHeight());        
				       
				    }		
			};
			filterpanelIndel = new JPanel(new GridBagLayout()) {				
				private static final long serialVersionUID = 1L;

				protected void paintComponent(Graphics g) {
				        super.paintComponent(g);     
				  
				        g.setColor(backColor);
				        g.fillRect(0, 0, this.getWidth(), this.getHeight());        
				       
				    }		
			};
			hidepanel = new JPanel(new GridBagLayout()) {				
				private static final long serialVersionUID = 1L;

				protected void paintComponent(Graphics g) {
				        super.paintComponent(g);     
				  
				        g.setColor(backColor);
				        g.fillRect(0, 0, this.getWidth(), this.getHeight());        
				       
				    }		
			};
			comparepanel = new JPanel(new GridBagLayout()) {		
				
				private static final long serialVersionUID = 1L;

				protected void paintComponent(Graphics g) {
				        super.paintComponent(g);     
				         
				        g.setColor(backColor);
				        g.fillRect(0, 0, this.getWidth(), this.getHeight());	        
				    }		
			};
			
			aminopanel = new JPanel() {		
				
				private static final long serialVersionUID = 1L;

				protected void paintComponent(Graphics g) {
				        super.paintComponent(g);     
				        
				        g.setColor(backColor);
				        g.fillRect(0, 0, this.getWidth(), this.getHeight());	        
				       
				    }		
			};
		
		geneSlider.setValue(1);
		geneSlider.setSnapToTicks(true);
		geneSlider.setMajorTickSpacing(1);
		geneSlider.setMinorTickSpacing(1);
		geneSlider.setOpaque(false);
		geneSlider.addChangeListener(this);
		geneSlider.addMouseWheelListener(sliderWheelListener);
		commonSlider.addMouseWheelListener(sliderWheelListener);
	//	commonSlider.setPreferredSize(buttondimension);
		commonSlider.setValue(1);
		commonSlider.setUpperValue(1);		
		commonSlider.addChangeListener(this);
		commonSlider.setOpaque(false);
		qualitySlider.addMouseWheelListener(sliderWheelListener);
		qualitySlider.addChangeListener(this);
		qualitySlider.addMouseListener(this);
	//	qualitySlider.setPreferredSize(buttondimension);
		qualitySlider.setValue(0);
		qualitySlider.setOpaque(false);
		qualitySliderIndel.addMouseWheelListener(sliderWheelListener);
		qualitySliderIndel.addChangeListener(this);
		qualitySliderIndel.addMouseListener(this);
	//	qualitySliderIndel.setPreferredSize(buttondimension);
		qualitySliderIndel.setValue(0);
		qualitySliderIndel.setOpaque(false);
		gqSlider.addMouseWheelListener(sliderWheelListener);
		gqSlider.addChangeListener(this);
		gqSlider.addMouseListener(this);
//		gqSlider.setPreferredSize(buttondimension);		
		 
		gqSlider.setValue(0);
		gqSlider.setOpaque(false);		
		gqSliderIndel.addMouseWheelListener(sliderWheelListener);
		gqSliderIndel.addChangeListener(this);
		gqSliderIndel.addMouseListener(this);
	//	gqSliderIndel.setPreferredSize(buttondimension);
		gqSliderIndel.setValue(0);
		gqSliderIndel.setOpaque(false);		
		coverageSlider.addMouseWheelListener(sliderWheelListener);
		coverageSlider.addChangeListener(this);
		coverageSlider.addMouseListener(this);
		coverageSlider.setOpaque(false);
		coverageSlider.setValue(4);		
		coverageSliderIndel.addMouseWheelListener(sliderWheelListener);
		coverageSliderIndel.addChangeListener(this);
		coverageSliderIndel.addMouseListener(this);
		coverageSliderIndel.setOpaque(false);
		coverageSliderIndel.setValue(4);		
		maxCoverageSlider.addMouseWheelListener(sliderWheelListener);
		maxCoverageSlider.addChangeListener(this);
		maxCoverageSlider.addMouseListener(this);
		maxCoverageSlider.setOpaque(false);
		maxCoverageSlider.setValue(1500);
		maxCoverageSliderIndel.addMouseWheelListener(sliderWheelListener);
		maxCoverageSliderIndel.addChangeListener(this);
		maxCoverageSliderIndel.addMouseListener(this);
		maxCoverageSliderIndel.setOpaque(false);
		maxCoverageSliderIndel.setValue(1500);
		callSlider.addMouseWheelListener(sliderWheelListener);
		callSlider.addChangeListener(this);
		callSlider.addMouseListener(this);
		callSlider.setValue(10);
		callSlider.setUpperValue(100);
		callSlider.setOpaque(false);
		callSliderIndel.addMouseWheelListener(sliderWheelListener);
		callSliderIndel.addChangeListener(this);
		callSliderIndel.addMouseListener(this);
		callSliderIndel.setValue(10);
		callSliderIndel.setUpperValue(100);
		callSliderIndel.setOpaque(false);
		indelFilters.setOpaque(false);
		indelFilters.addActionListener(this);
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
		menuScroll = new JScrollPane(menuPanel);	
	//	menuScroll.setMinimumSize(new Dimension(200, 100));
		menuPanel.setBackground(Color.white);	
		
		menu.addKeyListener(this);		
		menuScroll.addKeyListener(this);
	//	menuScroll.setMaximumSize(new Dimension(250,500));
		menu.add(menuScroll);	
		menu.add(applyQualities);
		
		menuScrollIndel = new JScrollPane(menuPanelIndel);	
//		menuScrollIndel.setMinimumSize(new Dimension(200, 20));
		menuPanelIndel.setBackground(Color.white);	
		menuIndel.add(menuScrollIndel);	
		menuIndel.addKeyListener(this);		
		menuScrollIndel.addKeyListener(this);
//		menuScrollIndel.setMaximumSize(new Dimension(250,500));
		menuIndel.add(applyQualitiesIndel);
		
		applyQualitiesIndel.addActionListener(this);
		applyQualities.addActionListener(this);
		advQualities.addActionListener(this);
		advQualitiesIndel.addActionListener(this);
		qualityLabel.setToolTipText("Variants below quality threshold will be hidden");		
		gqLabel.setToolTipText("Variants below quality threshold will be hidden");
		GridBagConstraints c = new GridBagConstraints();	
		c.anchor = GridBagConstraints.WEST;		
		c.fill = GridBagConstraints.HORIZONTAL;
		c.insets = new Insets(2,4,0,4);
		c.gridx = 0;
		c.gridy = 0;
		c.weightx = 0;
		c.weighty = 0;
		c.gridwidth = 1;		
		filterpanel.add(SNVFilters,c);
		qualityLabel.setToolTipText("Click for hard filters (advanced)");
		qualityLabel.addMouseListener(this);
		c.gridx = 1;
		adder2 = new JLabel("__________________________________");
		adder2.setForeground(Draw.sidecolor);
		filterpanel.add(adder2,c);
		filterpanel.add(totalVars,c);		
		c.gridy++;
		c.gridx = 0;
		filterpanel.add(qualityLabel,c);
		c.gridx = 1;
		filterpanel.add(qualitySlider,c);
		c.gridy++;
		c.gridx = 0;
		filterpanel.add(gqLabel,c);
		c.gridx = 1;
		filterpanel.add(gqSlider,c);
		//filterpanel.add(advQualities);
		//filterpanel.add(new JSeparator());
		c.gridy++;
		c.gridx = 0;
		filterpanel.add(coverageLabel,c);
		c.gridx = 1;
		filterpanel.add(coverageSlider,c);
		c.gridy++;
		c.gridx = 0;
		filterpanel.add(maxCoverageLabel,c);
		c.gridx = 1;
		filterpanel.add(maxCoverageSlider,c);
		c.gridy++;
		c.gridx = 0;		
		adder = new JLabel("_______________________________________");
		adder.setForeground(Draw.sidecolor);
		filterpanel.add(adder,c);
		filterpanel.add(callsLabel,c);
		c.gridx = 1;
		filterpanel.add(callSlider,c);	
		c.gridy++;
		c.weightx = 1;
		c.weighty = 1;
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.gridx = 0;
		filterpanel.add(new JLabel(),c);	
		
		c.gridx = 0;
		c.gridy = 0;
		c.gridwidth = 1;
		c.weightx = 0;
		c.weighty = 0;
		c.insets = new Insets(0,0,0,0);
		filterpanelIndel.add(indelFilters,c);
		qualityLabelIndel.setToolTipText("Click for hard filters (advanced)");
		qualityLabelIndel.addMouseListener(this);
		c.insets = new Insets(2,4,0,4);
		c.gridx = 1;
		adder4 = new JLabel("__________________________________");
		adder4.setForeground(Draw.sidecolor);
		filterpanelIndel.add(adder4,c);
		filterpanelIndel.add(totalVarsIndel,c);		
		c.gridy++;
		c.gridx = 0;
		filterpanelIndel.add(qualityLabelIndel,c);
		c.gridx = 1;
		filterpanelIndel.add(qualitySliderIndel,c);
		c.gridy++;
		c.gridx = 0;
		filterpanelIndel.add(gqLabelIndel,c);
		c.gridx = 1;
		filterpanelIndel.add(gqSliderIndel,c);
		c.gridy++;
//		c.gridx = 0;
//		filterpanelIndel.add(advQualitiesIndel,c);
//		c.gridx = 1;
//		filterpanelIndel.add(new JSeparator(),c);
//		c.gridy++;
		c.gridx = 0;
		filterpanelIndel.add(coverageLabelIndel,c);
		c.gridx = 1;
		filterpanelIndel.add(coverageSliderIndel,c);
		c.gridy++;
		c.gridx = 0;
		filterpanelIndel.add(maxCoverageLabelIndel,c);
		c.gridx = 1;
		filterpanelIndel.add(maxCoverageSliderIndel,c);
		c.gridy++;
		c.gridx = 0;
		adder3 = new JLabel("__________________________________");
		adder3.setForeground(Draw.sidecolor);
		filterpanelIndel.add(adder3,c);
		filterpanelIndel.add(callsLabelIndel,c);
		c.gridx = 1;
		filterpanelIndel.add(callSliderIndel,c);	
		c.gridy++;
		c.weightx = 1;
		c.weighty = 1;
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.gridx = 0;
		filterpanelIndel.add(new JLabel(),c);			
		
		
		c.insets = new Insets(0,4,0,4);
		
	//	gb.anchor = GridBagConstraints.WEST;
		c.gridx = 0;
		c.gridy = 0;
		c.gridwidth = 1;		
		hidepanel.add(hidenoncoding,c);
		c.gridy++;
		hidepanel.add(rscode,c);
		c.gridy++;
		hidepanel.add(hideSNVs,c);
		c.gridy++;
		hidepanel.add(hideIndels,c);
		c.gridy++;
		hidepanel.add(hideHomos,c);
		c.gridy++;
		freeze.setBackground(new Color(170,220,255));
		hidepanel.add(freeze,c);
		freezeIndels(true);
		c.gridy++;
		c.weightx = 1;
		c.weighty = 1;
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.gridx = 0;
		hidepanel.add(new JLabel(),c);			
		
		c.gridx = 0;
		c.gridy = 0;
		c.gridwidth = 1;
		c.insets = new Insets(2,4,0,4);
		c.weightx = 0;
		c.weighty = 0;
		comparepanel.add(comparison,c);
		c.gridy++;
	//	comparepanel.add(new JLabel(""));
		
		comparepanel.add(geneLabel,c);
		adder5 = new JLabel("__________________________________________");
		adder5.setForeground(Draw.sidecolor);
		comparepanel.add(adder5,c);
		c.gridx = 1;
		comparepanel.add(geneSlider,c);
		c.gridx = 0;
		c.gridy++;
		comparepanel.add(slideLabel,c);
		c.gridx = 1;
		comparepanel.add(commonSlider,c);
		c.gridx = 0;
		c.gridy++;
		comparepanel.add(clusterLabel,c);
		c.gridx = 1;
		clusterBox.addKeyListener(this);
		comparepanel.add(clusterBox,c);
		c.gridy++;
		c.weightx = 1;
		c.weighty = 1;
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.gridx = 0;
		comparepanel.add(new JLabel(),c);		
		varcalc.addActionListener(this);		
	//	varcalc.setPreferredSize(buttondimension);		
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
		if(Main.settingsIcon == null) {
			URL imgUrl = getClass().getResource("icons/settings.png");
			Main.settingsIcon = new ImageIcon(imgUrl);
			imgUrl = getClass().getResource("icons/save.gif");
			Main.save = new ImageIcon(imgUrl);
		}
		aminomenu.setIcon(Main.settingsIcon);
		writetofile.addActionListener(this);
		writetofile.setIcon(Main.save);
		outputmenu.setIcon(Main.save);
		aminopanel.add(aminobar);
		aminopanel.add(aminoCount);
		aminobar.add(varcalc);
		aminobar.add(aminomenu);
		tsv.setSelected(true);
		outputgroup.add(tsv);
		outputgroup.add(compactTsv);
		outputgroup.add(vcf);
		outputgroup.add(oncodrive);
		outputmenu.add(tsv);
		outputmenu.add(compactTsv);
		outputmenu.add(vcf);
		outputmenu.add(oncodrive);
	//	outputmenu.add(forVEP);
		write.setIcon(Main.save);
		outputmenu.add(write);
		aminobar.add(outputmenu);
		setValues();
		//setFonts();
	//	add(table);
	}
	
	static void setValues() {
		if(variantSettings == null) {
			variantSettings = new HashMap<String, Integer>();
			variantSettings.put("SNVquality", 0);
			variantSettings.put("SNVgq", 0);
			variantSettings.put("mincoverage", 4);
			variantSettings.put("maxcoverage", 1500);
			variantSettings.put("fraction", 10);
			variantSettings.put("upperFraction", 100);
			variantSettings.put("Indelquality", 0);
			variantSettings.put("Indelgq", 0);
			variantSettings.put("Indelmincoverage", 4);
			variantSettings.put("Indelmaxcoverage", 1500);
			variantSettings.put("Indelfraction", 10);
			variantSettings.put("IndelUpperFraction", 100);
			variantSettings.put("IndelFilters", 0);
			variantSettings.put("geneSlider", 1);
			variantSettings.put("commonSliderMin", 1);
			variantSettings.put("commonSliderMax", 1);
			variantSettings.put("clusterSize", 0);
			variantSettings.put("hideNoncoding", 0);
			variantSettings.put("hideRscoded", 0);
			variantSettings.put("hideSNVs", 0);
			variantSettings.put("hideIndels", 0);
			variantSettings.put("hideHomozygotes", 0);
			variantSettings.put("freeze", 0);
			variantSettings.put("nonSyn", 0);
			variantSettings.put("truncs", 0);
			variantSettings.put("intronic", 0);
			variantSettings.put("intergenic", 0);
			variantSettings.put("utr", 0);
			variantSettings.put("allChroms", 0);
			variantSettings.put("onlySel", 0);
			variantSettings.put("onlyStats", 0);
			variantSettings.put("writeFile", 0);			
			variantSettings.put("tsv", 1);
			variantSettings.put("compact", 0);
			variantSettings.put("vcf", 0);
			variantSettings.put("oncodrive", 0);			
		}
		
		VariantHandler.qualitySlider.setValue(variantSettings.get("SNVquality"));
		VariantHandler.gqSlider.setValue(variantSettings.get("SNVgq"));
		VariantHandler.coverageSlider.setValue(variantSettings.get("mincoverage"));
		VariantHandler.maxCoverageSlider.setValue(variantSettings.get("maxcoverage"));
		VariantHandler.callSlider.setValue(variantSettings.get("fraction"));
	
		VariantHandler.callSlider.setUpperValue(variantSettings.get("upperFraction"));
		
		VariantHandler.qualitySliderIndel.setValue(variantSettings.get("Indelquality"));
		VariantHandler.gqSliderIndel.setValue(variantSettings.get("Indelgq"));
		VariantHandler.coverageSliderIndel.setValue(variantSettings.get("Indelmincoverage"));
		VariantHandler.maxCoverageSliderIndel.setValue(variantSettings.get("Indelmaxcoverage"));
		VariantHandler.callSliderIndel.setValue(variantSettings.get("Indelfraction"));
		
		VariantHandler.callSliderIndel.setUpperValue(variantSettings.get("IndelUpperFraction"));
		
		VariantHandler.indelFilters.setSelected(variantSettings.get("IndelFilters") == 1);
		checkIndelFilters();
		VariantHandler.maxSlideValue = VariantHandler.maxCoverageSlider.getValue();
		
		if(Main.varsamples > 1) {
			VariantHandler.commonSlider.setMinimum(1);
			VariantHandler.commonSlider.setMaximum(Main.varsamples);
			VariantHandler.geneSlider.setMinimum(1);
			VariantHandler.geneSlider.setMaximum(Main.varsamples);
		}
		else {
			VariantHandler.commonSlider.setMinimum(1);
			VariantHandler.commonSlider.setMaximum(1);
			VariantHandler.geneSlider.setMinimum(1);
			VariantHandler.geneSlider.setMaximum(1);
		}
		VariantHandler.geneSlider.setValue(variantSettings.get("geneSlider"));
		VariantHandler.commonSlider.setValue(variantSettings.get("commonSliderMin"));
		VariantHandler.commonSlider.setUpperValue(variantSettings.get("commonSliderMax"));
		VariantHandler.clusterSize = variantSettings.get("clusterSize");
		VariantHandler.clusterBox.setText(""+VariantHandler.clusterSize);
		VariantHandler.hidenoncoding.setSelected(variantSettings.get("hideNoncoding") == 1);
		VariantHandler.rscode.setSelected(variantSettings.get("hideRscoded") == 1);
		VariantHandler.hideSNVs.setSelected(variantSettings.get("hideSNVs") == 1);
		VariantHandler.hideIndels.setSelected(variantSettings.get("hideIndels") == 1);
		VariantHandler.hideHomos.setSelected(variantSettings.get("hideHomozygotes") == 1);
		VariantHandler.freeze.setSelected(variantSettings.get("freeze") == 1);
		
		VariantHandler.synonymous.setSelected(variantSettings.get("nonSyn") == 1);
		VariantHandler.nonsense.setSelected(variantSettings.get("truncs") == 1);
		VariantHandler.intronic.setSelected(variantSettings.get("intronic") == 1);
		VariantHandler.intergenic.setSelected(variantSettings.get("intergenic") == 1);
		VariantHandler.utr.setSelected(variantSettings.get("utr") == 1);
		VariantHandler.allChroms.setSelected(variantSettings.get("allChroms") == 1);
		VariantHandler.onlyselected.setSelected(variantSettings.get("onlySel") == 1);
		VariantHandler.onlyStats.setSelected(variantSettings.get("onlyStats") == 1);
		VariantHandler.writetofile.setSelected(variantSettings.get("writeFile") == 1);
		checkWriteFiles();
		
		VariantHandler.tsv.setSelected(variantSettings.get("tsv") == 1);
		VariantHandler.compactTsv.setSelected(variantSettings.get("compact") == 1);
		VariantHandler.vcf.setSelected(variantSettings.get("vcf") == 1);
		VariantHandler.oncodrive.setSelected(variantSettings.get("oncodrive") == 1);
		
	}
	static void saveValues() {
		variantSettings.put("SNVquality", VariantHandler.qualitySlider.getValue());
		variantSettings.put("SNVgq", VariantHandler.gqSlider.getValue());
		variantSettings.put("mincoverage", VariantHandler.coverageSlider.getValue());
		variantSettings.put("maxcoverage", VariantHandler.maxCoverageSlider.getValue());
		variantSettings.put("fraction", VariantHandler.callSlider.getValue());
		
		variantSettings.put("upperFraction", VariantHandler.callSlider.getUpperValue());
		
		variantSettings.put("Indelquality", VariantHandler.qualitySliderIndel.getValue());
		variantSettings.put("Indelgq", VariantHandler.gqSliderIndel.getValue());
		variantSettings.put("Indelmincoverage", VariantHandler.coverageSliderIndel.getValue());
		variantSettings.put("Indelmaxcoverage", VariantHandler.maxCoverageSliderIndel.getValue());
		variantSettings.put("Indelfraction", VariantHandler.callSliderIndel.getValue());
		
		variantSettings.put("IndelUpperFraction", VariantHandler.callSliderIndel.getUpperValue());
		
		variantSettings.put("IndelFilters", indelFilters.isSelected() ? 1 : 0);
		
		variantSettings.put("geneSlider", VariantHandler.geneSlider.getValue());
		variantSettings.put("commonSliderMin",VariantHandler.commonSlider.getValue());
		variantSettings.put("commonSliderMax",VariantHandler.commonSlider.getUpperValue());
		variantSettings.put("clusterSize", clusterSize);
		variantSettings.put("hideNoncoding", VariantHandler.hidenoncoding.isSelected() ? 1 : 0);
		variantSettings.put("hideRscoded", VariantHandler.rscode.isSelected() ? 1 : 0);
		variantSettings.put("hideSNVs", VariantHandler.hideSNVs.isSelected() ? 1 : 0);
		variantSettings.put("hideIndels", VariantHandler.hideIndels.isSelected() ? 1 : 0);
		variantSettings.put("hideHomozygotes", VariantHandler.hideHomos.isSelected() ? 1 : 0);
		variantSettings.put("freeze", VariantHandler.freeze.isSelected() ? 1 : 0);
		
		variantSettings.put("nonSyn", VariantHandler.synonymous.isSelected() ? 1 : 0);
		variantSettings.put("truncs", VariantHandler.nonsense.isSelected() ? 1 : 0);
		variantSettings.put("intronic", VariantHandler.intronic.isSelected() ? 1 : 0);
		variantSettings.put("intergenic", VariantHandler.intergenic.isSelected() ? 1 : 0);
		variantSettings.put("utr", VariantHandler.utr.isSelected() ? 1 : 0);
		variantSettings.put("allChroms", VariantHandler.allChroms.isSelected() ? 1 : 0);
		variantSettings.put("onlySel", VariantHandler.onlyselected.isSelected() ? 1 : 0);
		variantSettings.put("onlyStats", VariantHandler.onlyStats.isSelected() ? 1 : 0);
		variantSettings.put("writeFile", VariantHandler.writetofile.isSelected() ? 1 : 0);
	
		variantSettings.put("tsv", VariantHandler.tsv.isSelected() ? 1 : 0);
		variantSettings.put("compact", VariantHandler.compactTsv.isSelected() ? 1 : 0);
		variantSettings.put("vcf", VariantHandler.vcf.isSelected() ? 1 : 0);
		variantSettings.put("oncodrive", VariantHandler.oncodrive.isSelected() ? 1 : 0);
		
	}
	private static void createAndShowGUI() {	
		try {
			
			frame = new JFrame("Variant Manager");
		    frame.setResizable(true);    
		    JComponent newContentPane = new VariantHandler();
		    newContentPane.setOpaque(true); 
		    frame.setContentPane(newContentPane);
		    frame.pack();
		    VariantHandler.filterPanes.setMinimumSize(filterPanes.getSize());
		    table.setPreferredSize(new Dimension(tableScroll.getViewport().getWidth(),tableScroll.getViewport().getHeight()));
			table.setMinimumSize(new Dimension(tableScroll.getViewport().getWidth(),tableScroll.getViewport().getHeight()));
			table.resizeTable(tableScroll.getViewport().getWidth());		
		    aminobar.setMinimumSize(new Dimension((int)aminobar.getSize().getWidth(),(int)aminobar.getSize().getHeight()));
		    frame.setAlwaysOnTop(true);
		    filters.setMinimumSize(new Dimension((int)filters.getSize().getWidth(),(int)filters.getSize().getHeight()));
		    clusterTable.resizeTable(tableScroll.getViewport().getWidth());
		    if(Main.chromDraw == null) {
		    	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
			    frame.setVisible(true);
			   String add =  "##INFO=<ID=ABHom,Number=1,Type=Float,Description=\"Allele Balance for homs (A/(A+O))\">";
			   addMenuComponents(add);
			   add = "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">";
			   addMenuComponents(add);
			   
		    }
		    else {
			    frame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE); 
			    frame.setVisible(false);
		    }
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	   
	}
	
	public static void main(String[] args) {
		try {
			//UIManager.setLookAndFeel("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel"); 
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
			System.setProperty("sun.java2d.d3d", "false"); 	
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
			geneLabel.setText("At least " +geneSlider.getValue() +"/" +geneSlider.getMaximum() +" samples share a mutated gene");
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
			
			callsLabel.setText("Min/max allelic fraction: " +callSlider.getValue() +"% - " +callSlider.getUpperValue() +"%");
			if(Settings.selectedVarDraw == 1) {
				maxSlideValue = callSlider.getUpperValue()/(float)100;
			}
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
			//maxCoverage = maxCoverageSlider.getValue();
			if(Settings.selectedVarDraw == 0) {
				maxSlideValue = maxCoverageSlider.getValue();
			}
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
		else if(event.getSource() == qualitySliderIndel){
			
			qualityLabelIndel.setText("Min. quality score: " +qualitySliderIndel.getValue());
			
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
		else if(event.getSource() == gqSliderIndel) {
			gqLabelIndel.setText("Min. genotype quality score: " +gqSliderIndel.getValue());
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
		else if(event.getSource() == coverageSliderIndel){
			
			coverageLabelIndel.setText("Min. coverage: " +coverageSliderIndel.getValue());
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
		else if(event.getSource() == callSliderIndel){
			callsLabelIndel.setText("Min/max allelic fraction: " +callSliderIndel.getValue() +"% - " +callSliderIndel.getUpperValue() +"%");
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
		else if(event.getSource() == maxCoverageSliderIndel){			
			maxCoverageLabelIndel.setText("Max. coverage: " +maxCoverageSliderIndel.getValue());
	//		maxCoverageIndel = maxCoverageSliderIndel.getValue();
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
	static void checkIndelFilters() {
		if(indelFilters.isSelected()) {			
			freezeIndels(false);
			SNVFilters.setText("SNV filters");			
		}
		else {
			freezeIndels(true);
			SNVFilters.setText("SNV & indel filters");			
		}
		if(Main.drawCanvas != null) {
			Draw.updatevars = true;		
			Main.drawCanvas.repaint();
		}
	}
	static void checkWriteFiles() {
		if(writetofile.isSelected()) {
			writetofile.setBackground(Color.white);
			//tabs.setEnabled(false);
			aminomenu.add(tsv,aminomenu.getItemCount()-1);
			aminomenu.add(compactTsv,aminomenu.getItemCount()-1);
			aminomenu.add(vcf,aminomenu.getItemCount()-1);
			aminomenu.add(oncodrive,aminomenu.getItemCount()-1);
			aminomenu.getPopupMenu().pack();
			//varcalc.setText("Annotate and write");
			aminomenu.revalidate();
			aminomenu.repaint();
			tabs.revalidate();
		}
		else {				
			writetofile.setBackground(Color.lightGray);
			//tabs.setEnabled(true);	
			outputmenu.add(oncodrive,0);
			outputmenu.add(vcf,0);
			outputmenu.add(compactTsv,0);
			outputmenu.add(tsv,0);
			outputmenu.revalidate();
			aminomenu.getPopupMenu().pack();
		//	varcalc.setText("Annotate");
			tabs.revalidate();
		}
	}
	@Override
	public void actionPerformed(ActionEvent event) {
		
		if(event.getSource() == writetofile) {
			checkWriteFiles();	
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
		else if(event.getSource() == indelFilters) {
			
			checkIndelFilters();
		}
		
		else if(event.getSource() == varcalc) {
			/*if(VariantHandler.onlyStats.isSelected()) {
				try {
			  		String outfname = "";
			  		String path =Main.savedir;			 		    		   		    
		    		JFileChooser chooser = new JFileChooser(path);		    
		    		int returnVal = chooser.showSaveDialog((Component)this.getParent());		    		
		    		
		    		if(returnVal == JFileChooser.APPROVE_OPTION) {  		    	  
			    		 outfname = chooser.getSelectedFile().getAbsolutePath();		      	
			    		 BufferedWriter output = null, sigOutput =null;
			    		 Main.savedir = chooser.getSelectedFile().getParent();
			        	 Main.writeToConfig("DefaultSaveDir=" +chooser.getSelectedFile().getParent());
			        	 lastWrittenPos = 0;
				    	 if(tsv.isSelected() || compactTsv.isSelected() || oncodrive.isSelected()) {
				    		   if(!outfname.contains(".tsv")) {
						    	   File outfile = new File(outfname +".tsv");						    	
						    	   FileRead.outputName = outfname +".tsv";
						    	   output = new BufferedWriter(new FileWriter(outfile));
						       }
						       else {
						    	   File outfile = new File(outfname);
						    	   FileRead.outputName = outfname;						    	
						    	   output = new BufferedWriter(new FileWriter(outfile));
						       }    
				    		 
								 if(output != null) {
									sigOutput = new BufferedWriter(new FileWriter(outfname +"_signatures.tsv"));
									String header = createTSVHeader();
									output.write(header);
									VariantHandler.table.clear();
									VariantHandler.table.headerHover = 2;
									VariantHandler.table.sorter.ascending = true;
									VariantHandler.table.createPolygon();
									VariantHandler.table.repaint();	
									FileRead.output = output;	
									FileRead.sigOutput = sigOutput;
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
			else */if(writetofile.isSelected()) {
				try {
			  		String outfname = "";
			  		String path =Main.savedir;			 		    		   		    
		    		JFileChooser chooser = new JFileChooser(path);		    
		    		int returnVal = chooser.showSaveDialog((Component)this.getParent());		    		
		    		if(returnVal == JFileChooser.APPROVE_OPTION) {  		    	  
		    		outfname = chooser.getSelectedFile().getAbsolutePath();		      	
		    		BufferedWriter output = null;
		    		Main.savedir = chooser.getSelectedFile().getParent();
		        	 Main.writeToConfig("DefaultSaveDir=" +chooser.getSelectedFile().getParent());
		        	 lastWrittenPos = 0;
			    	 if(tsv.isSelected() || compactTsv.isSelected() || oncodrive.isSelected()) {
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
						       }
						       else {
						       	   FileRead.outputName = outfname;						    	
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
				aminomenu.add(onlyAutosomes,7);
				aminomenu.getPopupMenu().pack();
			}
			else {
				aminomenu.remove(allChromsfrom);		
				aminomenu.remove(onlyAutosomes);	
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
		    	 if(tsv.isSelected() || compactTsv.isSelected() || oncodrive.isSelected()) {
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
			
			//setFonts();
			
			menu.show(this, 100, 100);
		}
		else if(event.getSource() == advQualitiesIndel) {
			menuIndel.show(this, 100, 100);
		}
		else if(event.getSource() == applyQualities) {
				Main.drawCanvas.drawVariables.advQDraw = null;
			if(VariantHandler.indelFilters.isSelected()) {	
				for(int i = 0; i<menuPanel.getComponentCount(); i++) {
					if(menuPanel.getComponent(i) instanceof JLabel) {
						JLabel label = (JLabel)menuPanel.getComponent(i);
						if(menuPanel.getComponent(i+1) instanceof JTextField) {
							JTextField field = (JTextField)menuPanel.getComponent(i+1);
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
						else {
							
							JCheckBox field = (JCheckBox)menuPanel.getComponent(i+1);
							
													
							if(field.isSelected()) {
								if(Main.drawCanvas.drawVariables.advQDraw == null) {
									Main.drawCanvas.drawVariables.advQDraw = new ArrayList<QualEntry>();
								}	
								Main.drawCanvas.drawVariables.advQDraw.add(new QualEntry(field.getText(), 1F, ""));				
							}
							
						}
					}
				}
			}
			else {
				for(int i = 0; i<menuPanel.getComponentCount(); i++) {
					if(menuPanel.getComponent(i) instanceof JLabel) {
						JLabel label = (JLabel)menuPanel.getComponent(i);
						if(menuPanel.getComponent(i+1) instanceof JTextField) {
							JTextField field = (JTextField)menuPanel.getComponent(i+1);
							
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
									JTextField indelfield = (JTextField)menuPanelIndel.getComponent(i+1);
									indelfield.setText(format +number);
								}						
							}
							catch(Exception e) {
								field.setForeground(Color.red);
							}
						}
						else {
							if(menuPanel.getComponent(i+1) instanceof JCheckBox) {
								JCheckBox field = (JCheckBox)menuPanel.getComponent(i+1);
								
								if(Main.drawCanvas.drawVariables.advQDraw == null) {
									Main.drawCanvas.drawVariables.advQDraw = new ArrayList<QualEntry>();
								}							
								if(field.isSelected()) {
									Main.drawCanvas.drawVariables.advQDraw.add(new QualEntry(field.getText(), 1F, ""));	
									JCheckBox fieldIndel = (JCheckBox)menuPanelIndel.getComponent(i+1);
									fieldIndel.setSelected(true);
								}
							}
						}
					}					
				}
			}
			if(FileRead.head.getNext() != null) {
				FileRead.search = true;
				Main.drawCanvas.forcereload = true;
				Main.drawCanvas.gotoPos(Main.drawCanvas.splits.get(0).chrom, Main.drawCanvas.splits.get(0).start, Main.drawCanvas.splits.get(0).end);
				Main.drawCanvas.forcereload = false;
			}
		}
		else if(event.getSource() == applyQualitiesIndel) {
			Main.drawCanvas.drawVariables.advQDrawIndel = null;
			for(int i = 0; i<menuPanelIndel.getComponentCount(); i++) {
				if(menuPanelIndel.getComponent(i) instanceof JLabel) {
					JLabel label = (JLabel)menuPanelIndel.getComponent(i);
					if(menuPanel.getComponent(i+1) instanceof JTextField) {
						JTextField field = (JTextField)menuPanelIndel.getComponent(i+1);
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
								
								if(Main.drawCanvas.drawVariables.advQDrawIndel == null) {
									Main.drawCanvas.drawVariables.advQDrawIndel = new ArrayList<QualEntry>();
								}							
								Main.drawCanvas.drawVariables.advQDrawIndel.add(new QualEntry(label.getText(), number, format));							
							}						
						}
						catch(Exception e) {
							field.setForeground(Color.red);
						}
						
					}
					else {
						JCheckBox field = (JCheckBox)menuPanelIndel.getComponent(i+1);
						
						if(Main.drawCanvas.drawVariables.advQDrawIndel == null) {
							Main.drawCanvas.drawVariables.advQDrawIndel = new ArrayList<QualEntry>();
						}							
						if(field.isSelected()) {
							Main.drawCanvas.drawVariables.advQDrawIndel.add(new QualEntry(field.getText(), 1F, ""));	
							
						}
					}
				}
				
			}
			if(FileRead.head.getNext() != null) {
				FileRead.search = true;
				Main.drawCanvas.forcereload = true;
				Main.drawCanvas.gotoPos(Main.drawCanvas.splits.get(0).chrom, Main.drawCanvas.splits.get(0).start, Main.drawCanvas.splits.get(0).end);
				Main.drawCanvas.forcereload = false;
			}
		}
	
}
static void removeMenuComponents() {
	
	
		menuPanel.removeAll();
		menuPanelIndel.removeAll();
		menuPanel.add(new JLabel("Hard filters"));
		menuPanelIndel.add(new JLabel("Hard filters"));
	/*for(int i = menu.getComponentCount()-2; i>0; i--) {
		if(menuScroll.getComponent(i) instanceof JLabel || menu.getComponent(i) instanceof JTextField) {
			menu.remove(i);			
		}
	}*/
}
static void addMenuComponents(String line) {
	if(line.startsWith("##FILTER")) {
		String key = line.substring(line.indexOf("<ID=")+4, line.indexOf(","));
		if(Main.drawCanvas != null && Main.drawCanvas.advQualities == null) {
			Main.drawCanvas.advQualities = new HashMap<String, Float>();
			Main.drawCanvas.advQualitiesIndel = new HashMap<String, Float>();
		}
		if(Main.drawCanvas != null) {
			if(!Main.drawCanvas.advQualities.containsKey(key)) {
				Main.drawCanvas.advQualities.put(key, 0F);
				String description = line.substring(line.indexOf("ion=\"")+5);
				description = description.substring(0, description.indexOf("\""));
				JCheckBox addCheck = new JCheckBox(key);
				JCheckBox addCheckIndel = new JCheckBox(key);
				addCheck.setToolTipText(description);
				
				GridBagConstraints constr = new GridBagConstraints();			
				constr.fill = GridBagConstraints.HORIZONTAL;		
				constr.anchor = GridBagConstraints.NORTHWEST;
				constr.gridx = 0;
				addCheck.setOpaque(false);
				constr.gridy = menuPanel.getComponentCount()/2;
				VariantHandler.menuPanel.add(new JLabel(),constr);	
				VariantHandler.menuPanelIndel.add(new JLabel(),constr);	
				
				constr.gridx = 1;
				VariantHandler.menuPanel.add(addCheck,constr);		
				VariantHandler.menuPanelIndel.add(addCheckIndel,constr);		
				if(Main.drawCanvas.drawVariables.advQDrawIndel != null && Main.drawCanvas.drawVariables.advQDrawIndel.size() > 0) {
					if(Main.drawCanvas.drawVariables.advQDraw != null && Main.drawCanvas.drawVariables.advQDraw.size() > 0) {
						
						for(int i = 0 ;i<Main.drawCanvas.drawVariables.advQDraw.size(); i++) {
							if(Main.drawCanvas.drawVariables.advQDraw.get(i).key.equals(key)) {
								addCheck.setSelected(Main.drawCanvas.drawVariables.advQDraw.get(i).value == 1);
								
								break;
							}
						}
					}
					if(Main.drawCanvas.drawVariables.advQDrawIndel != null && Main.drawCanvas.drawVariables.advQDrawIndel.size() > 0) {
						
						for(int i = 0 ;i<Main.drawCanvas.drawVariables.advQDrawIndel.size(); i++) {
							if(Main.drawCanvas.drawVariables.advQDrawIndel.get(i).key.equals(key)) {
								addCheckIndel.setSelected(Main.drawCanvas.drawVariables.advQDrawIndel.get(i).value == 1);
								break;
							}
						}
					}
				}
				else {
					if(Main.drawCanvas.drawVariables.advQDraw != null && Main.drawCanvas.drawVariables.advQDraw.size() > 0) {
						
						for(int i = 0 ;i<Main.drawCanvas.drawVariables.advQDraw.size(); i++) {
							if(Main.drawCanvas.drawVariables.advQDraw.get(i).key.equals(key)) {
								addCheck.setSelected(Main.drawCanvas.drawVariables.advQDraw.get(i).value == 1);
								
								addCheckIndel.setSelected(Main.drawCanvas.drawVariables.advQDraw.get(i).value == 1);
								
								break;
							}
						}
					}
				}
			}
		}
	}
	else {
		
		if(line.indexOf("Number=1")<0) {
			return;
		}
		String key = line.substring(line.indexOf("<ID=")+4, line.indexOf(","));
		if(Main.drawCanvas != null) {
			if(Main.drawCanvas.advQualities == null) {
				Main.drawCanvas.advQualities = new HashMap<String, Float>();
				Main.drawCanvas.advQualitiesIndel = new HashMap<String, Float>();
			}
			if(!Main.drawCanvas.advQualities.containsKey(key)) {
				GridBagConstraints constr = new GridBagConstraints();	
				
				constr.fill = GridBagConstraints.HORIZONTAL;
				
				constr.anchor = GridBagConstraints.NORTHWEST;
				constr.gridx = 0;
				constr.gridy = menuPanel.getComponentCount()/2;
			
				String[] split = line.substring(line.indexOf("<")+1).replace(">", "").split(",");
				StringBuffer tooltip = new StringBuffer("<html>");
				for(int i = 0;i<split.length; i++) {
					
					tooltip.append(split[i] +"<br>");
				}
				tooltip.append("</html>");
				
				Main.drawCanvas.advQualities.put(key, 0F);
				
				JLabel addLabel = new JLabel(key);
				JLabel addLabelIndel = new JLabel(key);
				addLabel.setToolTipText(tooltip.toString());
			
				VariantHandler.menuPanel.add(addLabel,constr);
				VariantHandler.menuPanelIndel.add(addLabelIndel,constr);
				
				JTextField field = new JTextField("<");
				JTextField fieldIndel = new JTextField("<");
				if(Main.drawCanvas.drawVariables.advQDrawIndel != null && Main.drawCanvas.drawVariables.advQDrawIndel.size() > 0) {
					if(Main.drawCanvas.drawVariables.advQDraw != null && Main.drawCanvas.drawVariables.advQDraw.size() > 0) {
						
						for(int i = 0 ;i<Main.drawCanvas.drawVariables.advQDraw.size(); i++) {
							if(Main.drawCanvas.drawVariables.advQDraw.get(i).key.equals(key)) {
								field.setText(Main.drawCanvas.drawVariables.advQDraw.get(i).format +" " +Main.drawCanvas.drawVariables.advQDraw.get(i).value);	
								break;
							}
						}
					}
					if(Main.drawCanvas.drawVariables.advQDrawIndel != null && Main.drawCanvas.drawVariables.advQDrawIndel.size() > 0) {
						
						for(int i = 0 ;i<Main.drawCanvas.drawVariables.advQDrawIndel.size(); i++) {
							if(Main.drawCanvas.drawVariables.advQDrawIndel.get(i).key.equals(key)) {
								fieldIndel.setText(Main.drawCanvas.drawVariables.advQDrawIndel.get(i).format +" " +Main.drawCanvas.drawVariables.advQDrawIndel.get(i).value);	
								break;
							}
						}
					}
				}
				else {
					if(Main.drawCanvas.drawVariables.advQDraw != null && Main.drawCanvas.drawVariables.advQDraw.size() > 0) {
						
						for(int i = 0 ;i<Main.drawCanvas.drawVariables.advQDraw.size(); i++) {
							if(Main.drawCanvas.drawVariables.advQDraw.get(i).key.equals(key)) {
								field.setText(Main.drawCanvas.drawVariables.advQDraw.get(i).format +" " +Main.drawCanvas.drawVariables.advQDraw.get(i).value);	
								fieldIndel.setText(Main.drawCanvas.drawVariables.advQDraw.get(i).format +" " +Main.drawCanvas.drawVariables.advQDraw.get(i).value);	
								break;
							}
						}
					}
				}
				field.setPreferredSize(new Dimension(100,(int)(Main.defaultFontSize*2)));
				field.setToolTipText(tooltip.toString());
				fieldIndel.setPreferredSize(new Dimension(100,(int)(Main.defaultFontSize*2)));
				fieldIndel.setToolTipText(tooltip.toString());
				constr.gridx = 1;
				VariantHandler.menuPanel.add(field,constr);		
				VariantHandler.menuPanelIndel.add(fieldIndel,constr);		
			}	
		}			
	}
	
	for(int i = 0 ; i<VariantHandler.menuPanel.getComponentCount(); i++) {
		VariantHandler.menuPanel.getComponent(i).setFont(Main.menuFont);
		
		if(VariantHandler.menuPanel.getComponent(i) instanceof JTextField) {			
			VariantHandler.menuPanel.getComponent(i).setPreferredSize(new Dimension(Main.defaultFontSize*6,(int)(Main.defaultFontSize*2)));
		}
	}
	for(int i = 0 ; i<VariantHandler.menuPanelIndel.getComponentCount(); i++) {
		VariantHandler.menuPanelIndel.getComponent(i).setFont(Main.menuFont);
		if(VariantHandler.menuPanelIndel.getComponent(i) instanceof JTextField) {			
			VariantHandler.menuPanelIndel.getComponent(i).setPreferredSize(new Dimension(Main.defaultFontSize*6,(int)(Main.defaultFontSize*2)));
		}
	}
	VariantHandler.menu.pack();
	VariantHandler.menuIndel.pack();
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
		StringBuffer headerstring = new StringBuffer("##fileformat=VCFv4.1"+Main.lineseparator
				+ "##BasePlayer=<Version: " +Main.version +" output " +new SimpleDateFormat("dd.MM.yyyy HH:mm").format(Calendar.getInstance().getTime())+Main.lineseparator
				+ "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"+Main.lineseparator
				+ "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">"+Main.lineseparator
				+ "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">" +Main.lineseparator
				+ "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"+Main.lineseparator
				+ "##reference=" +Main.ref.getName() +Main.lineseparator);
				
				
				if(VariantHandler.onlyselected.isSelected()) {
					headerstring.append("##sample="+Main.drawCanvas.selectedSample.getName()+Main.lineseparator);
					headerstring.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+Main.drawCanvas.selectedSample.getName() +Main.lineseparator);
				}
				else {
					StringBuffer headersamples = new StringBuffer("");
					headerstring.append("##files=");
					for(int i = 0; i<Main.drawCanvas.sampleList.size(); i++) {						
						if(Main.drawCanvas.sampleList.get(i).getTabixFile() != null || Main.drawCanvas.sampleList.get(i).calledvariants || (VariantCaller.inanno.isSelected() && Main.drawCanvas.sampleList.get(i).samFile != null)) {
							headerstring.append(Main.drawCanvas.sampleList.get(i).getName() +",");
							
						}	
						if((VariantCaller.inanno.isSelected() && Main.drawCanvas.sampleList.get(i).samFile != null) || Main.drawCanvas.sampleList.get(i).calledvariants || !Main.drawCanvas.sampleList.get(i).multiVCF && (Main.drawCanvas.sampleList.get(i).getTabixFile() != null || Main.drawCanvas.sampleList.get(i).multipart) && !Main.drawCanvas.sampleList.get(i).removed) {
							
							headersamples.append("\t" +Main.drawCanvas.sampleList.get(i).getName());
						}
					}					
					headerstring.deleteCharAt(headerstring.length()-1);
					headerstring.append(Main.lineseparator);
					headerstring.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
					headerstring.append(headersamples +Main.lineseparator);
					headersamples = new StringBuffer();					
				}
	
		return headerstring.toString();
	}
	String createTSVHeader() {
		 StringBuffer headerstring = new StringBuffer("");
   	  headerstring.append("##BasePlayer version: " +Main.version +" output " +new SimpleDateFormat("dd.MM.yyyy HH:mm").format(Calendar.getInstance().getTime())+Main.lineseparator);
   	  
   	  headerstring.append("##Files: ");
	   	if(VariantHandler.onlyselected.isSelected()) {
			headerstring.append(Main.drawCanvas.selectedSample.getName());
		}
		else {
	   	  for(int i = 0; i< Main.samples; i++) {
	   		if(Main.drawCanvas.sampleList.get(i).calledvariants || (Main.drawCanvas.sampleList.get(i).getTabixFile() != null || Main.drawCanvas.sampleList.get(i).multipart) && !Main.drawCanvas.sampleList.get(i).removed) {
				
	   			headerstring.append(Main.drawCanvas.sampleList.get(i).getName() +",");
	   		  }   		  
	   	  }
	   	 headerstring.deleteCharAt(headerstring.length()-1);
		}
   	 
   	  headerstring.append(Main.lineseparator);
   	  headerstring.append("##Genome:" +Main.ref.getName() +",Annotation:" +Main.annotationfile +Main.lineseparator);   	  
   	 
   	  headerstring.append("##" +VariantHandler.synonymous.getText() +":" +VariantHandler.synonymous.isSelected() +"," +VariantHandler.utr.getText() +":" +VariantHandler.utr.isSelected() +"," +VariantHandler.intronic.getText() +":"+VariantHandler.intronic.isSelected()+","+VariantHandler.nonsense.getText() +":" +VariantHandler.nonsense.isSelected() +"," +VariantHandler.intergenic.getText() +":" +VariantHandler.intergenic.isSelected()+Main.lineseparator); 
   	  
   	  if(Control.controlData.fileArray.size() > 0 && Control.controlData.controlsOn) {
   		 
	    		  headerstring.append("##Controls:");
		    	  for(int i = 0; i< Control.controlData.fileArray.size(); i++) {
		    		  if(Control.controlData.fileArray.get(i).controlOn) {
		    			  headerstring.append(Control.controlData.fileArray.get(i).getName() +":" +Control.controlData.fileArray.get(i).alleleFreq);
		    			  if(Control.controlData.fileArray.get(i).remOverlaps.isSelected()) {
		    				  headerstring.append(",Overlap_indels");
		    			  }
		    			  headerstring.append(",");
		    		  }		    		  
		    	  }	    		  
		    	  headerstring.deleteCharAt(headerstring.length()-1);
		    	  headerstring.append(Main.lineseparator);
	    	  
   	  }
   	  else {
   		  headerstring.append("##No controls applied"+Main.lineseparator);
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
   		  headerstring.append(Main.lineseparator);
   	  }
   	  if(VariantHandler.commonSlider.getValue() > 1) {
   		  
   		  if(VariantHandler.clusterSize > 0) {
	    	headerstring.append("##Variant clusters in " +VariantHandler.commonSlider.getValue() +"/" +Main.varsamples +" samples within " +VariantHandler.clusterSize +"bp"+Main.lineseparator);
	      }
   		  else {
   			headerstring.append("##Common variants in " +VariantHandler.commonSlider.getValue() +"/" +Main.varsamples +" samples"+Main.lineseparator);   			
   		  }   		  
   	  }
   	  if(VariantHandler.geneSlider.getValue() > 1) {
   		  headerstring.append("##At least " +VariantHandler.geneSlider.getValue() +"/" +Main.varsamples +" samples share a mutated gene"+Main.lineseparator);
   	  }
   	  headerstring.append("##Variant filters:"+Main.lineseparator);
   	  headerstring.append("##Hide rs-coded variants: " +VariantHandler.rscode.isSelected() +Main.lineseparator);
   	  headerstring.append("##Hide SNVs: " +VariantHandler.hideSNVs.isSelected() +Main.lineseparator);
   	  headerstring.append("##Hide indels: " +VariantHandler.hideIndels.isSelected() +Main.lineseparator);
   	  if(indelFilters.isSelected()) {
	   	  headerstring.append("##Min. coverage SNVs: " +VariantHandler.coverageSlider.getValue()+", Indels:  "  +VariantHandler.coverageSliderIndel.getValue() +Main.lineseparator);
	   	  headerstring.append("##Min. allelic/fraction SNVs: " +VariantHandler.callSlider.getValue()+"%, Indels:  " +VariantHandler.callSliderIndel.getValue() +"%"+Main.lineseparator);
	   	  headerstring.append("##Min. quality score SNVs: " +VariantHandler.qualitySlider.getValue()+", Indels:  "+VariantHandler.qualitySliderIndel.getValue()  +Main.lineseparator);   
	      headerstring.append("##Min. genotype quality score SNVs: " +VariantHandler.gqSlider.getValue()+", Indels:  " +VariantHandler.gqSliderIndel.getValue() +Main.lineseparator);   
	   	  headerstring.append("##Max. coverage SNVs: " +VariantHandler.maxCoverageSlider.getValue() +", Indels:  "+VariantHandler.maxCoverageSliderIndel.getValue() +Main.lineseparator);   
   	  }
   	  else {
   		  headerstring.append("##Min. coverage: " +VariantHandler.coverageSlider.getValue() +Main.lineseparator);
	   	  headerstring.append("##Min. allelic/fraction: " +VariantHandler.callSlider.getValue() +"%"+Main.lineseparator);
	   	  headerstring.append("##Min. quality score: " +VariantHandler.qualitySlider.getValue() +Main.lineseparator);   
	      headerstring.append("##Min. genotype quality score: " +VariantHandler.gqSlider.getValue() +Main.lineseparator);   
	   	  headerstring.append("##Max. coverage: " +VariantHandler.maxCoverageSlider.getValue() +Main.lineseparator);   
   	  }
   	  
   	  
   	  if(Main.drawCanvas.drawVariables.advQDraw != null) {
   		  for(int i = 0 ; i<Main.drawCanvas.drawVariables.advQDraw.size(); i++) {
   			headerstring.append("##" +Main.drawCanvas.drawVariables.advQDraw.get(i).key +": " +Main.drawCanvas.drawVariables.advQDraw.get(i).format +Main.drawCanvas.drawVariables.advQDraw.get(i).value +Main.lineseparator);
   		  }
   	  }
   	 if(Main.drawCanvas.drawVariables.advQDrawIndel != null) {
  		  for(int i = 0 ; i<Main.drawCanvas.drawVariables.advQDrawIndel.size(); i++) {
  			headerstring.append("##Indel: " +Main.drawCanvas.drawVariables.advQDrawIndel.get(i).key +": " +Main.drawCanvas.drawVariables.advQDrawIndel.get(i).format +Main.drawCanvas.drawVariables.advQDrawIndel.get(i).value +Main.lineseparator);
  		  }
  	  }
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
	   		if(oncodrive.isSelected()) {
	   			headerstring.append("#CHROM\tPOS\tREF\tALT\tSAMPLE"+Main.lineseparator);
	   		 }
	   		else {
	   			headerstring.append("#Sample\tGene\tMutationCount\tSampleCount\tENSG\tENST\tBioType\tPosition\tStrand\tRegion\tEffect\tBaseChange\tGenotype(calls/coverage)\tQuality\tGQ\trs-code\t" +clusters +controls +tracks+"Description"+Main.lineseparator);
	   		   
	   		}
	   	}
	   	else {
	   		if(onlyStats.isSelected()) {
	   			headerstring.append("#Sample\tVariants\tSNVs\tDELs\tINSs\tCoding\tHetero/homo-rate\tTS/TV-rate\tT>A\tT>C\tT>G\tC>A\tC>G\tC>T\tAvg.call/cov\tSynonymous\tNonsynonymous\tMissense\tSplice-site\tNonsense\tFrameShift\tInframe"+Main.lineseparator);
	   		}
	   		else {
	   	 		headerstring.append("#Sample\tVariants\tSNVs\tDELs\tINSs\tCoding\tHetero/homo-rate\tTS/TV-rate\tT>A\tT>C\tT>G\tC>A\tC>G\tC>T\tAvg.call/cov"+Main.lineseparator);
	   	 	
	   		}
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
	   // 	 output.write("#Sample\tVariants\tSNVs\tDELs\tINSs\tCoding\tHetero/homo-rate\tTS/TV-rate\tT>A\tT>C\tT>G\tC>A\tC>G\tC>T\tAvg.call/cov");
	 	  	 for(int i = 0 ; i<VariantHandler.stattable.sampleArray.size(); i++) {
	 	  		 sample = (Sample)VariantHandler.stattable.sampleArray.get(i)[0];
	 	  		 output.write(sample.getName());
	 	  		 
	 	  		 for(int j = 1; j<VariantHandler.stattable.headerlengths.length; j++) {
	 	  			output.write("\t" +VariantHandler.stattable.sampleArray.get(i)[j]);
	 	  		 }
	 	  		 if(onlyStats.isSelected()) {
	 	  			 output.write("\t" +sample.syn +"\t" +sample.nonsyn +"\t" +sample.missense +"\t" +sample.splice +"\t" +sample.nonsense +"\t" +sample.fshift +"\t" +sample.inframe);
	 	  		 }
	 	  		output.write(Main.lineseparator);
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
	 		    						node.getTranscripts().get(0).getChrom() +":" +MethodLibrary.formatNumber((node.getPosition()+1)) +"\t" +strand	+"\t" +"Intergenic"+"\t" +"N/A" +"\t" +entry.getKey() +"\t" +genotype +"\t" +rscode +"\t" +clusters +controls +tracks +"N/A" +Main.lineseparator);	    
	 		    					
			    					 }
			    					 else {
	    						 		output.write(varnode.getSample().getName() +"\t" +
			    						genes +"\t" +mutcount +"\t" +"N/A" +"\t" +geneIDs +"\t" +transIDs +"\t" +"N/A" +"\t"+
			    						node.getTranscripts().get(0).getChrom() +":" +MethodLibrary.formatNumber((node.getPosition()+1)) +"\t" +strand	+"\t" +"Intergenic"+"\t" +"N/A" +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotype +"\t" +rscode +"\t"  +clusters +controls +tracks +"N/A" +Main.lineseparator);	    
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
				 if(!Main.drawCanvas.hideVar(entry.getValue().get(sample), entry.getKey().length() > 1)) {
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
						 sampleinfo = samplehash.get(entry.getValue().get(sample).getSample().getMultiIndex()).charAt(2) +"/" +allelenro +":"+entry.getValue().get(sample).getCoverage()+":" +ADs+":" +entry.getValue().get(sample).getGQString();
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
						 sampleinfo = (allelenro) +"/" +allelenro +":"+entry.getValue().get(sample).getCoverage()+":"+ADs+":" +entry.getValue().get(sample).getGQString();
						 
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
						
						 sampleinfo = "0/" +allelenro +":"+entry.getValue().get(sample).getCoverage()+":"+ ADs.toString() +":" +entry.getValue().get(sample).getGQString();
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
			
			 String writeline = chrom +"\t" +(node.getPosition()+1) +"\t" +rscode +"\t" +ref +"\t" +alts +"\t" +MethodLibrary.round(avgquality/(double)samplecount,2) +"\tPASS\t" +info +"\t" +format +sampleinfos+Main.lineseparator;
			 
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
			
			 output.write(chrom +"\t" +(node.getPosition()+1) +"\t" +rscode +"\t" +ref +"\t" +alts +"\t" +MethodLibrary.round(avgquality/(double)samplecount,2) +"\t.\t" +info +"\t" +format +sampleinfos+Main.lineseparator );
		 }
		 else {
		//	 System.out.println(chrom +"\t" +(node.getPosition()+1) +"\t" +rscode +"\t" +ref +"\t" +alts +"\t" +MethodLibrary.round(avgquality/(double)samplecount,2) +"\t.\t" +info +"\t" +format +sampleinfos);
		 }
	}
	else if(oncodrive.isSelected()) {
		Entry<String, ArrayList<SampleNode>> entry;
			String[] result = null;
			 for(int var = 0; var < node.vars.size(); var++) {
				 
				 entry = node.vars.get(var);
				 if(Main.drawCanvas.hideNodeVar(node, entry)) {
						continue;
				 }	
			if(node.indel && entry.getKey().length() > 1) {
				 result = MethodLibrary.makeIndelColumns(node.getChrom(), node.getPosition(), Main.getBase.get(node.getRefBase()), entry.getKey());
					
			}
			 for(int i = 0; i<entry.getValue().size(); i++) {
 				
				 if(Main.drawCanvas.hideVar(entry.getValue().get(i), entry.getKey().length() > 1)) {
					continue;
				 }
				 
				// varnode = entry.getValue().get(i);
				
				  				
				 try {
					 if(output != null) {
						 
						 
						 
						 if(node.indel && entry.getKey().length() > 1) {
								
							 output.write(node.getChrom() +"\t" +(node.getPosition()) +"\t" +result[0] +"\t" +result[1] +"\t" +entry.getValue().get(i).getSample().getName() +Main.lineseparator);
							 
						 }
						 else {
							
							 output.write(node.getChrom() +"\t" +(node.getPosition()+1) +"\t" +Main.getBase.get(node.getRefBase()) +"\t" +entry.getKey() +"\t" +entry.getValue().get(i).getSample().getName() +Main.lineseparator);
	    						
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
    		 StringBuffer GQualities = new StringBuffer("");
    		
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
	    				//continue;
	    			}
	    		 }
	    		
	    		
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
    	    		 GQualities = new StringBuffer("");
    	    		 for(int i = 0; i<entry.getValue().size(); i++) {
 	    				
	    				 if(Main.drawCanvas.hideVar(entry.getValue().get(i), entry.getKey().length() > 1)) {
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
	    				 GQualities.append(varnode.getGQString() +";");
    	    		 }
    	    		 
    	    		genotypes.deleteCharAt(genotypes.length()-1);
    	    		qualities.deleteCharAt(qualities.length()-1);
    	    		samples.deleteCharAt(samples.length()-1);
    	    		GQualities.deleteCharAt(GQualities.length()-1);
    	    		 if(output != null) {
    	    			
		    				 if(exons.length() > 0) {
		    					 if(!aminochange.equals("N/A")) {
		    						 output.write(samples +"\t" +
		    						 row[0] +"\t" +row[1] +"\t" +gene.samples.size() +"\t" +geneID +"\t" +transcripts +"\t" +biotype +"\t"+
						  		  	 row[2] +"\t" +strand	+"\t" +exons +"\t" +aminochange +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotypes +"\t" +qualities +"\t" +GQualities +"\t"+rscode +"\t" +clusters +controls +tracks +description +Main.lineseparator);
		    					 }
		    					 else {
		    						 output.write(samples +"\t" +
		    						 row[0] +"\t" +row[1] +"\t" +gene.samples.size() +"\t" +geneID +"\t" +transcripts +"\t" +biotype +"\t"+
		    						 row[2] +"\t" +strand	+"\t" +"UTR" +"\t" +aminochange +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotypes +"\t"+qualities +"\t" +GQualities +"\t" +rscode +"\t"+clusters+controls +tracks +description +Main.lineseparator);
		    					 }
		    				 }
		    				 else {
		    					 if(node.getTranscripts() != null && node.isInGene()) {
		    						 output.write(samples +"\t" +
		    						 row[0] +"\t" +row[1] +"\t" +gene.samples.size() +"\t" +geneID +"\t" +transcripts +"\t" +biotype +"\t"+
		    						 row[2] +"\t" +strand	+"\t" +"Intronic"+"\t" +"N/A" +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotypes +"\t"+qualities +"\t" +GQualities +"\t" +rscode+"\t" +clusters+controls +tracks +description +Main.lineseparator);	    						 
		    					 }
		    					 else {
		    						 output.write(samples +"\t" +
		    						 row[0] +"\t" +row[1] +"\t" +gene.samples.size() +"\t" +geneID +"\t" +transcripts +"\t" +biotype +"\t"+
		    						 row[2] +"\t" +strand	+"\t" +"Intergenic"+"\t" +"N/A" +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotypes +"\t" +qualities +"\t" +GQualities +"\t"+rscode+"\t" +clusters +controls +tracks +description +Main.lineseparator);	    
		    					 }
		    				 }
    	    			 }	    	    		
	    	    		
	    			}
	    			
	    			else if(tsv.isSelected()) {
	    			 for(int i = 0; i<entry.getValue().size(); i++) {
	    				
	    				 if(Main.drawCanvas.hideVar(entry.getValue().get(i), entry.getKey().length() > 1)) {
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
							  		  	 row[2] +"\t" +strand	+"\t" +exons +"\t" +aminochange +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotype +"\t" +varnode.getQuality() +"\t" +varnode.getGQString() +"\t"+rscode +"\t" +clusters +controls +tracks +description +Main.lineseparator);
			    					 }
			    					 else {
			    						 output.write(varnode.getSample().getName() +"\t" +
			    						 row[0] +"\t" +row[1] +"\t" +gene.samples.size() +"\t" +geneID +"\t" +transcripts +"\t" +biotype +"\t"+
			    						 row[2] +"\t" +strand	+"\t" +"UTR" +"\t" +aminochange +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotype +"\t"+varnode.getQuality() +"\t" +varnode.getGQString() +"\t" +rscode +"\t"+clusters+controls +tracks +description +Main.lineseparator);
			    					 }
			    				 }
			    				 else {
			    					 if(node.getTranscripts() != null && node.isInGene()) {
			    						 output.write(varnode.getSample().getName() +"\t" +
			    						 row[0] +"\t" +row[1] +"\t" +gene.samples.size() +"\t" +geneID +"\t" +transcripts +"\t" +biotype +"\t"+
			    						 row[2] +"\t" +strand	+"\t" +"Intronic"+"\t" +"N/A" +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotype +"\t"+varnode.getQuality() +"\t" +varnode.getGQString() +"\t" +rscode+"\t" +clusters+controls +tracks +description +Main.lineseparator);	    						 
			    					 }
			    					 else {
			    						 output.write(varnode.getSample().getName() +"\t" +
			    						 row[0] +"\t" +row[1] +"\t" +gene.samples.size() +"\t" +geneID +"\t" +transcripts +"\t" +biotype +"\t"+
			    						 row[2] +"\t" +strand	+"\t" +"Intergenic"+"\t" +"N/A" +"\t" +Main.getBase.get(node.getRefBase()) +"->" +entry.getKey() +"\t" +genotype +"\t" +varnode.getQuality() +"\t" +varnode.getGQString() +"\t"+rscode+"\t" +clusters +controls +tracks +description +Main.lineseparator);	    
			    					 }
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
	VarNode node = null;
	if(VariantHandler.table.genearray.size() > 0 && VariantHandler.table.genearray.get(0).varnodes != null && VariantHandler.table.genearray.get(0).varnodes.size() > 0) {
		node = VariantHandler.table.genearray.get(0).varnodes.get(0);
	}
	
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
	public void mouseClicked(MouseEvent event) {
		
		if(event.getSource() == qualityLabel) {
			
			//setFonts();
		
			menu.show(this, 100, 100);
		}
		else if(event.getSource() == qualityLabelIndel) {
			
			menuIndel.show(this, 100, 100);
			
			
		}
	//	filterPanes.setBackground(Main.panel.getBackground());
		
	}

	@Override
	public void mouseEntered(MouseEvent event) {
		if(event.getSource() == qualityLabel) {
			qualityLabel.setForeground(Color.white);
			
		}
		else if(event.getSource() == qualityLabelIndel) {
			qualityLabelIndel.setForeground(Color.white);
		}
	}

	@Override
	public void mouseExited(MouseEvent event) {
		if(event.getSource() == qualityLabel) {
			qualityLabel.setForeground(Color.black);
			
		}
		else if(event.getSource() == qualityLabelIndel) {
			qualityLabelIndel.setForeground(Color.black);
		}
	}

	@Override
	public void mousePressed(MouseEvent event) {
		
		if(event.getSource() == tabs) {			
			if(tabs.getSelectedIndex() == 0) {
				VariantHandler.aminoCount.setText(table.variants +" variants");
				outputmenu.setText("Variant output");
				vcf.setVisible(true);
				VariantHandler.compactTsv.setVisible(true);
				VariantHandler.oncodrive.setVisible(true);
				outputmenu.revalidate();
			}
			else if(tabs.getSelectedIndex() == tabs.indexOfComponent(statsScroll)) {
				if(stattable.bufImage.getWidth() == 1) {
					
					stattable.bufImage = MethodLibrary.toCompatibleImage(new BufferedImage(Main.screenSize.width,Main.screenSize.height, BufferedImage.TYPE_INT_ARGB));	
					stattable.buf = (Graphics2D)stattable.bufImage.getGraphics();
				}
				VariantHandler.aminoCount.setText(stattable.variants +" variants");
				outputmenu.setText("Stats output");
				outputmenu.revalidate();
				vcf.setVisible(false);
				VariantHandler.compactTsv.setVisible(false);
				VariantHandler.oncodrive.setVisible(false);
			}
			else if(tabs.getSelectedIndex() == tabs.indexOfComponent(clusterScroll)) {
				if(clusterTable.bufImage.getWidth() == 1) {
					clusterTable.bufImage = MethodLibrary.toCompatibleImage(new BufferedImage(Main.screenSize.width,Main.screenSize.height, BufferedImage.TYPE_INT_ARGB));	
					clusterTable.buf = (Graphics2D)clusterTable.bufImage.getGraphics();
				}
				VariantHandler.aminoCount.setText(clusterTable.variants +" variants");
				outputmenu.setText("Variant output");
				vcf.setVisible(true);
				VariantHandler.compactTsv.setVisible(true);
				VariantHandler.oncodrive.setVisible(true);
				outputmenu.revalidate();
			}
			else {
				VariantHandler.aminoCount.setText(tables.get(tabs.getSelectedIndex()-(tabs.indexOfComponent(statsScroll)+1)).variants +" variants");
				outputmenu.setText("Variant output");
				vcf.setVisible(true);
				VariantHandler.compactTsv.setVisible(true);
				VariantHandler.oncodrive.setVisible(true);
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
		else if(event.getSource() == coverageSliderIndel) {
			if(coverageSliderIndel.getValue() == coverageSliderIndel.getMaximum()) {
				coverageSliderIndel.setMaximum(coverageSliderIndel.getMaximum()*2);
			}
			if(commonSlider.getValue() > 1) {
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}
		}
		else if(event.getSource() == qualitySliderIndel) {
			if(qualitySliderIndel.getValue() == qualitySliderIndel.getMaximum()) {
				qualitySliderIndel.setMaximum(qualitySliderIndel.getMaximum()*2);
			}
			if(commonSlider.getValue() > 1) {
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}
		}
		else if(event.getSource() == gqSliderIndel) {
			if(gqSliderIndel.getValue() == gqSliderIndel.getMaximum()) {
				gqSliderIndel.setMaximum(gqSliderIndel.getMaximum()*2);
			}
			if(commonSlider.getValue() > 1) {
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}
		}
		else if(event.getSource() == maxCoverageSliderIndel) {
			if(maxCoverageSliderIndel.getValue() == maxCoverageSliderIndel.getMaximum()) {
				maxCoverageSliderIndel.setMaximum(maxCoverageSliderIndel.getMaximum()*2);
			}
			if(commonSlider.getValue() > 1) {
			
				
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}
		}
		else if(event.getSource() == callSliderIndel) {
			if(commonSlider.getValue() > 1) {			
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}
		}
	}
	void freezeFilters(boolean value) {
		
		if(value) {
			for(int i =2;i<filterpanel.getComponentCount(); i++) {
				if(filterpanel.getComponent(i).equals(freeze) ) {
					if(filterpanelIndel.getComponent(i) instanceof JLabel) {
						JLabel label = (JLabel)filterpanelIndel.getComponent(i);
						if(label.getText().contains("___")) {
							continue;
						}
					}
					continue;
				}
				filterpanel.getComponent(i).setEnabled(false);
				filterpanel.getComponent(i).revalidate();
				
			}
			for(int i =2;i<filterpanelIndel.getComponentCount(); i++) {
				if(filterpanelIndel.getComponent(i).equals(freeze)) {
					if(filterpanelIndel.getComponent(i) instanceof JLabel) {
						JLabel label = (JLabel)filterpanelIndel.getComponent(i);
						if(label.getText().contains("___")) {
							continue;
						}
					}
					continue;
				}
				filterpanelIndel.getComponent(i).setEnabled(false);
				filterpanelIndel.getComponent(i).revalidate();
			}
		}
		else {
			for(int i =0;i<filterpanel.getComponentCount(); i++) {
				filterpanel.getComponent(i).setEnabled(true);
				filterpanel.getComponent(i).revalidate();
				
			}
			if(indelFilters.isSelected()) {
				for(int i =0;i<filterpanelIndel.getComponentCount(); i++) {
					filterpanelIndel.getComponent(i).setEnabled(true);
					filterpanelIndel.getComponent(i).revalidate();
				}
			}
		}
		
	}
	
static void setFonts(Font menuFont) {
	try {
		
		for(int i = 0 ; i<VariantHandler.filterpanel.getComponentCount(); i++) {	 
			if(VariantHandler.filterpanel.getComponent(i).getName() != null) {
				VariantHandler.filterpanel.getComponent(i).setFont(Main.menuFontBold);	 
			}
			else {
				VariantHandler.filterpanel.getComponent(i).setFont(menuFont);	    	
			}
		}
		
	    for(int i = 0 ; i<VariantHandler.filterpanelIndel.getComponentCount(); i++) {	
	    	if(VariantHandler.filterpanelIndel.getComponent(i).getName() != null) {
	    		VariantHandler.filterpanelIndel.getComponent(i).setFont(Main.menuFontBold);	    	
	    	}
	    	else {
	    		VariantHandler.filterpanelIndel.getComponent(i).setFont(menuFont);	    	
	    	}
	    	
	    }
	    for(int i = 0 ; i<VariantHandler.comparepanel.getComponentCount(); i++) {	
	    	if(VariantHandler.comparepanel.getComponent(i).getName() != null) {
	    		VariantHandler.comparepanel.getComponent(i).setFont(Main.menuFontBold);	    	
	    	}
	    	else {
	    		VariantHandler.comparepanel.getComponent(i).setFont(menuFont);	   
	    	}
	    }
	    for(int i = 0 ; i<VariantHandler.aminopanel.getComponentCount(); i++) {	   
	    	VariantHandler.aminopanel.getComponent(i).setFont(menuFont);	    	
	    }
	    for(int i = 0 ; i<VariantHandler.hidepanel.getComponentCount(); i++) {	   
	    	VariantHandler.hidepanel.getComponent(i).setFont(menuFont);	    	
	    }
	    VariantHandler.filterPanes.setFont(menuFont);
	    VariantHandler.filters.setFont(menuFont);
	   
	    for(int i = 0 ; i<VariantHandler.filters.getPopupMenu().getComponentCount(); i++) {	   
	    	VariantHandler.filters.getPopupMenu().getComponent(i).setFont(menuFont);    		
	    }
	    VariantHandler.aminomenu.setFont(menuFont);
	    VariantHandler.outputmenu.setFont(menuFont);
	    for(int i = 0; i<VariantHandler.aminomenu.getPopupMenu().getComponentCount(); i++) {
	    	VariantHandler.aminomenu.getPopupMenu().getComponent(i).setFont(menuFont);
	    }
	    for(int i = 0; i<VariantHandler.outputmenu.getPopupMenu().getComponentCount(); i++) {
	    	VariantHandler.outputmenu.getPopupMenu().getComponent(i).setFont(menuFont);
	    }   
	    varcalc.setFont(menuFont);
	    VariantHandler.allChromsfrom.setFont(menuFont.deriveFont(Font.ITALIC));
	    VariantHandler.onlyAutosomes.setFont(menuFont.deriveFont(Font.ITALIC));
	    VariantHandler.tabs.setFont(menuFont);
	    VariantHandler.tabs.revalidate();
	    
	    if(VariantHandler.table != null) {
	    	
	    	VariantHandler.table.buf.setFont(menuFont);
	    
		    VariantHandler.table.rowHeight = menuFont.getSize() +5;
		    VariantHandler.table.fm = VariantHandler.table.buf.getFontMetrics();
		  
		    VariantHandler.stattable.buf.setFont(menuFont);
		    VariantHandler.stattable.rowHeight = menuFont.getSize() +5;
		    VariantHandler.clusterTable.buf.setFont(menuFont);
		    VariantHandler.clusterTable.rowHeight = menuFont.getSize() +5;
		    for(int i = 0 ; i<VariantHandler.tables.size(); i++) {
		    	VariantHandler.tables.get(i).buf.setFont(menuFont);
		    	VariantHandler.tables.get(i).rowHeight = menuFont.getSize() +5;
		    }
	    }
	    
	    for(int i = 0 ; i<VariantHandler.menuPanel.getComponentCount(); i++) {
			VariantHandler.menuPanel.getComponent(i).setFont(Main.menuFont);
			
			if(VariantHandler.menuPanel.getComponent(i) instanceof JTextField) {
				
				VariantHandler.menuPanel.getComponent(i).setPreferredSize(new Dimension(Main.defaultFontSize*6,(int)(Main.defaultFontSize*1.5)));
			}
		}
		for(int i = 0 ; i<VariantHandler.menuPanelIndel.getComponentCount(); i++) {
			VariantHandler.menuPanelIndel.getComponent(i).setFont(Main.menuFont);
				if(VariantHandler.menuPanelIndel.getComponent(i) instanceof JTextField) {
				
				VariantHandler.menuPanelIndel.getComponent(i).setPreferredSize(new Dimension(Main.defaultFontSize*6,(int)(Main.defaultFontSize*1.5)));
			}
		}
		VariantHandler.menu.setFont(menuFont);
		VariantHandler.menuIndel.setFont(menuFont);
		for(int i = 0 ; i<VariantHandler.menu.getComponentCount(); i++) {
			VariantHandler.menu.getComponent(i).setFont(Main.menuFont);
			
		}
		for(int i = 0 ; i<VariantHandler.menuIndel.getComponentCount(); i++) {
			
			VariantHandler.menuIndel.getComponent(i).setFont(Main.menuFont);
			
		}
		VariantHandler.menu.pack();
		VariantHandler.menuIndel.pack();
	    VariantHandler.frame.pack();
	    aminobar.setMinimumSize(new Dimension((int)aminobar.getSize().getWidth(),(int)aminobar.getSize().getHeight()));
	   
	    
	}
	catch(Exception e) {
		e.printStackTrace();
	}
   
}
static void freezeIndels(boolean value) {		
	if(value) {
		for(int i =2;i<filterpanelIndel.getComponentCount(); i++) {
			if(filterpanelIndel.getComponent(i).equals(freeze)) {
				if(filterpanelIndel.getComponent(i) instanceof JLabel) {
					JLabel label = (JLabel)filterpanelIndel.getComponent(i);
					if(label.getText().contains("___")) {
						continue;
					}
				}
				continue;
			}
			filterpanelIndel.getComponent(i).setEnabled(false);
			filterpanelIndel.getComponent(i).revalidate();
		}
	}
	else {
		for(int i =0;i<filterpanelIndel.getComponentCount(); i++) {
			filterpanelIndel.getComponent(i).setEnabled(true);
			filterpanelIndel.getComponent(i).revalidate();
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
		if(e.getKeyCode() == KeyEvent.VK_0) {
			
		//	setFonts();
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
