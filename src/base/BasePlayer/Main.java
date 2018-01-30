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

import java.awt.Component;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GraphicsDevice;
import java.awt.GraphicsEnvironment;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.Insets;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.tribble.readers.TabixReader;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectStreamClass;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.lang.reflect.Field;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;



import javax.swing.BorderFactory;
import javax.swing.DefaultComboBoxModel;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JEditorPane;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSplitPane;
import javax.swing.JTextField;
import javax.swing.SwingWorker;
import javax.swing.ToolTipManager;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;
import javax.swing.filechooser.FileSystemView;
import javax.swing.plaf.basic.BasicSplitPaneDivider;
import javax.swing.plaf.basic.BasicSplitPaneUI;

import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;

import java.security.SecureRandom;
import java.security.cert.CertificateException;
import java.security.cert.X509Certificate;

import javax.net.ssl.HostnameVerifier;
import javax.net.ssl.HttpsURLConnection;
import javax.net.ssl.KeyManager;
import javax.net.ssl.SSLContext;
import javax.net.ssl.SSLSession;
import javax.net.ssl.TrustManager;
import javax.net.ssl.X509TrustManager;
	public class Main extends JPanel implements ActionListener, ChangeListener, ComponentListener, MouseListener, KeyListener, MouseMotionListener {
	private static final long serialVersionUID = 1L;		
    static JFrame frame;  
	//UI
    
    static ImageIcon save, open, settingsIcon;
    static String version = "1.0.0";
    static int sidebarWidth = 200;
    static String[] argsit = {}, args;
    static int defaultFontSize;
    static GraphicsDevice gd;   
    static boolean loaddraw = false;
    static HashMap<Byte, String> getBase = new HashMap<Byte,String>();
    static int width = Toolkit.getDefaultToolkit().getScreenSize().width;
    static int height= Toolkit.getDefaultToolkit().getScreenSize().height;
    static ArrayList<JLabel> labels = new ArrayList<JLabel>();
    static int loadTextWidth = 200;
    static Dimension screenSize;
    static int drawWidth, drawHeight, chromHeight = 150, bedHeight = 200;
    static GradientPaint gradient = new GradientPaint(drawWidth/2-loadTextWidth/2,0,Color.red,drawWidth/2+loadTextWidth,height,Color.green,true);
    private String[] searchList;
    static Short samples = 0, varsamples = 0, readsamples = 0;
    static boolean readingControls = false, readingbeds = false;
    public static String gerp;
    static int searchStart=-1, searchEnd=-1;
    static int textlength = 0, annolength = 0, reflength = 0;
    static int coverages = 0, alleles = 1, qualities = 2, gqs = 3;
    static int letterlength = 0 ;
	//Labels	
    static JTextField chromlabel = new JTextField(" Chrom ");
    static HashMap<String, String> factorNames = new HashMap<String, String>();
    static boolean configChanged = false;
    static int trackdivider = 0;
    static java.util.List<String> chromnamevector = Collections.synchronizedList(new ArrayList<String>());
    static Image A, C, G, T;
    static HashMap<Byte, Integer> baseMap = new HashMap<Byte, Integer>();
    static Dimension chromDimensions, bedDimensions, drawDimensions, buttonDimension;
    static int buttonHeight = 20, buttonWidth = 60;
    static FileRead fileReader = new FileRead();
    static String refchrom = "", selectedGenome;
    static int selectedChrom = 0;
    static Loader loading; 
    static Hashtable<String, Long[]> chromIndex = new Hashtable<String, Long[]>();
    static Hashtable<String, ArrayList<File>> genomehash = new Hashtable<String, ArrayList<File>>();
    static Hashtable<String, File> fastahash = new Hashtable<String, File>();
    static String[] chromnames;
    static String path;
    static boolean cancelhover = false, cancel = false;
    static JLabel backImage;
    static BedCanvas bedCanvas;
    static VarMaster varmaster;
    public static ChromDraw chromDraw;
    static Draw drawCanvas;
    static ControlCanvas controlDraw;
    static java.util.List<String> undoList = Collections.synchronizedList(new ArrayList<String>());
    static int undoPointer = 0;
    static Hashtable<String, Integer> mutTypes = new Hashtable<String, Integer>();
    static Hashtable<String, String> bases;
    static HashMap<Byte, Double> background = new HashMap<Byte, Double>();
	//Buttons etc.
    static JMenuItem variantCaller;
    static String defaultGenome = "";
    static String downloadDir = "";
    static JSplitPane splitPane, trackPane, varpane, drawpane;
    public static boolean nothread = false, noreadthread = false;	
    static Hashtable<String, String[]> searchTable = new Hashtable<String, String[]>();
    static Hashtable<String, String> geneIDMap = new Hashtable<String, String>();
    static JMenuItem addGenome = new JMenuItem("Add new genome...");	   
    static JButton zoomout;
    static JButton dosomething = new JButton("Do stuff!");
    static JButton back, forward;
  //  static double[][] snow = new double[200][4];
    static String[] snowtable = {"A", "C", "G", "T" };
    static String userDir;
    static File genomeDir;
    static boolean clicked = false, clickedAnno = false;
    VarNode prevTemp;
    static boolean rightclick = false;
    static JButton setbut;
    static JPanel chrompan;
	static JLabel memLabel = new JLabel(""), chromLabel = new JLabel("");		
    static JSplitPane upPanel;
	static JPanel panel;
	static JMenuBar menubar;
	static JMenu filemenu;
	static JMenu toolmenu;
	static JMenu help;
	static JMenu about;
	static JButton manage;
	static JMenuItem average;
	static JMenuItem settings;
	static JMenuItem update;
	static JMenuItem errorlog;
	static JLabel helpLabel = new JLabel("This is pre-release version of BasePlayer\nContact: help@baseplayer.fi\n\nUniversity of Helsinki");
	static JMenu genome;
	static JMenu addURL;
	static JTextField urlField;
	static boolean updatelauncher = false;
	static JMenuItem opensamples, exit;
	static JMenuItem addtracks;
	static JMenuItem addcontrols;
	static JMenuItem pleiadesButton;
	static JMenuItem saveProject;
	static JMenuItem saveProjectAs;
	static JMenuItem openProject;
	static JMenuItem clear;
	static JMenuItem clearMemory;
	static JEditorPane area;
	static RandomAccessFile referenceFile;    
	static Font menuFont, menuFontBold;
	static String[] empty = {""};	
	static DefaultComboBoxModel<String> chromModel = new DefaultComboBoxModel<String>(empty);  
	static SteppedComboBox chromosomeDropdown = new SteppedComboBox(chromModel);		
	static SteppedComboBox refDropdown;	
	static SteppedComboBox geneDropdown;	
	static File[] genomes;
	//UI	
    static MouseListener thisMainListener;    
    static DefaultComboBoxModel<String> refModel;
    static DefaultComboBoxModel<String> geneModel;
    static Hashtable<String, int[][]> SELEXhash = new Hashtable<String, int[][]>();
	static JScrollPane drawScroll;
	static JScrollPane chromScroll;
	static JScrollPane bedScroll;
	static JScrollPane controlScroll;
	static String hoverGenome = "", hoverAnnotation = "";
	static Image glass;	 
    static JTextField positionField = new JTextField();
	static JTextField widthLabel = new JTextField();
	static JTextField searchField = new JTextField("Search by position or gene") {
	       
		private static final long serialVersionUID = 1L;

		protected void paintComponent(Graphics g) {
            super.paintComponent(g);
           if(glass != null) {
            g.drawImage(glass, 4, (Main.searchField.getHeight()-(Main.defaultFontSize+4))/2, Main.defaultFontSize+4,Main.defaultFontSize+4, this);
           }
           }
    };
	ActionListener annoDropActionListener = new ActionListener() {			    	

		public void actionPerformed(ActionEvent actionEvent) {			    	
    	
    	  if (actionEvent.getActionCommand() == "comboBoxChanged") {	    
    		  if(clickedAnno) {
    			  if(Main.geneDropdown.getSelectedIndex() < Main.geneDropdown.getItemCount()-1) {
    				  changeAnnotation(geneDropdown.getSelectedIndex());
    				  geneDropdown.setToolTipText(geneDropdown.getSelectedItem().toString());
    			  }
    			  else {
    				  clickedAnno = false;    
    				 
     				 Main.geneDropdown.setSelectedItem(Main.defaultAnnotation);    				 
     				 Main.geneDropdown.revalidate();
     				 clickedAnno = true;
    				 for(int i=0; i<AddGenome.root.getChildCount(); i++) {    					
    					 if(AddGenome.root.getChildAt(i).toString().equals(refDropdown.getSelectedItem().toString())) {
    						 AddGenome.tree.setSelectionRow(i);
    						 AddGenome.tree.expandRow(i);
    						 break;
    					 }
    				 } 
    				 if(AddGenome.frame == null) {
   					  AddGenome.createAndShowGUI();
   				  	}
    				  AddGenome.annotation = true;	   
    				  AddGenome.remove.setEnabled(false);
    				  AddGenome.download.setEnabled(false);
    				  AddGenome.frame.setVisible(true);
    				  AddGenome.frame.setLocation(frame.getLocationOnScreen().x+frame.getWidth()/2 - AddGenome.frame.getWidth()/2, frame.getLocationOnScreen().y+frame.getHeight()/6);
    				  AddGenome.frame.setState(JFrame.NORMAL);
    				 				  
    				 
    				  /*if(geneDropdown.getSelectedItem() != null) {
    					  geneDropdown.setToolTipText(geneDropdown.getSelectedItem().toString());
    				  } */   				  
    			  }
    		  }
    	  }
		}
	};
	
	ActionListener refDropActionListener = new ActionListener() {		    	

				public void actionPerformed(ActionEvent actionEvent) {			    	
		    	
		    	  if (actionEvent.getActionCommand() == "comboBoxChanged") {	    		
		    		  if(clicked) {
		    			  if(!rightclick) {
			    			  if(Main.refDropdown.getSelectedIndex() < Main.refDropdown.getItemCount()-1) {
			    				  changeRef(Main.refDropdown.getSelectedItem().toString());
			    			  }
			    			  else {			    
			    				  clicked = false;
			    				  Main.refDropdown.setSelectedItem(Main.defaultGenome);
			    				  Main.refDropdown.revalidate();
			    				  clicked = true;
			    				  if(AddGenome.frame == null) {
			    					  AddGenome.createAndShowGUI();
			    				  }
			    				  AddGenome.remove.setEnabled(false);
			    				  AddGenome.download.setEnabled(false);
			    				  AddGenome.annotation = false;			    				
			    				  AddGenome.frame.setVisible(true);
			    				  AddGenome.frame.setLocation(frame.getLocationOnScreen().x+frame.getWidth()/2 - AddGenome.frame.getWidth()/2, frame.getLocationOnScreen().y+frame.getHeight()/6);
			    				  AddGenome.frame.setState(JFrame.NORMAL);
			    				 
			    			  }
		    			  }
		    			  else {
		    				
		    			  }
		    		  }
		    	  }
				}
	};
    	  
	ActionListener ChromoDropActionListener = new ActionListener() {
		
    	

		public void actionPerformed(ActionEvent actionEvent) {
    		
    	  try {
    	
    	  if (actionEvent.getActionCommand() == "comboBoxChanged") {	    		 
    		  if(chromosomeDropdown.getItemCount() == 1) {
    			  return;
    		  }
    		  
    		  
    		  drawCanvas.splits.get(0).getReadBuffer().setComposite( drawCanvas.composite);					
  			  drawCanvas.splits.get(0).getReadBuffer().fillRect(0,0, drawCanvas.splits.get(0).getReadImage().getWidth(),Main.drawScroll.getViewport().getHeight());	
  			  drawCanvas.splits.get(0).getReadBuffer().setComposite(drawCanvas.splits.get(0).getBackupr());		
  			  drawCanvas.rbuf.setComposite( drawCanvas.composite);				
  			  drawCanvas.rbuf.fillRect(0,0, drawCanvas.getWidth(),Main.drawScroll.getViewport().getHeight());	
  			  drawCanvas.rbuf.setComposite(drawCanvas.backupr);		
    		  drawCanvas.clearReads();
    		  selectedChrom = chromosomeDropdown.getSelectedIndex();    
    		  if(chromosomeDropdown.getSelectedItem() == null) {
    			  drawCanvas.chrom = "";
    			  drawCanvas.splits.get(0).chromEnd = 100;
    			  drawCanvas.splits.get(0).viewLength = 100;
    			  drawCanvas.splits.get(0).start = 0;
    		  }
    		  else {
    			 
    			  drawCanvas.chrom = chromosomeDropdown.getSelectedItem().toString();
    			  drawCanvas.splits.get(0).chromEnd = chromIndex.get(Main.refchrom + Main.chromosomeDropdown.getSelectedItem().toString())[1].intValue();
    			  chromLabel.setText("Chromosome " +drawCanvas.chrom);
        		  chromDraw.cytoImage = null;	    		 
        		  drawCanvas.splits.get(0).setCytoImage(null);
        		  drawCanvas.splits.get(0).chrom = drawCanvas.chrom;
        		  drawCanvas.splits.get(0).transStart = 0;
        		  drawCanvas.splits.get(0).nullRef();
        		 
        		  FileRead filereader = new FileRead();
        		  filereader.chrom = Main.chromosomeDropdown.getSelectedItem().toString();	    		  
    		 
    		 
    		
    		  
    		  if(!FileRead.search) {
    			 
    			  FileRead.searchStart = 0;
    			  FileRead.searchEnd = drawCanvas.splits.get(0).chromEnd;
    			  drawCanvas.setStartEnd(0.0,(double)drawCanvas.splits.get(0).chromEnd);
    		  }
    		  try {
	    		  for (int i = 0; i<Control.controlData.fileArray.size(); i++) {
	    			  Control.controlData.fileArray.get(i).controlled = false;
	    		  }
	    		  }
	    		  catch(Exception ex) {
	    			  Control.controlData.fileArray = Collections.synchronizedList(new ArrayList<ControlFile>());
	    		  }
	    		  if(nothread) {		  
	    			 
	    			  filereader.changeChrom(filereader.chrom);		    			 
	    			  nothread = false;
	    		  }
	    		  else {
	    			 
	    			  filereader.changeChrom = true;
	 	        	  filereader.execute();	 	        	
	    		  }
    	  		}   	  
    		  } 			
    	
		  chromosomeDropdown.revalidate();
		  chromosomeDropdown.repaint();	    
			}
			catch(Exception e) {
			
				e.printStackTrace();
			}		    			
      	}
    };

//		private String searchString;
	static Dimension fieldDimension;
	private File[] annofiles;

	private int keyCode;

	private String[] returnlist;

	private JTextField chooserTextField;

	private Image iconImage;

	private boolean pleiades = false;

	static String lineseparator;

	static String savedir = "";

	static String controlDir ="";
	static String trackDir = "";
	static String projectDir = "";
	static BasicSplitPaneDivider splitPaneDivider;

	static BasicSplitPaneDivider varPaneDivider;

	static BasicSplitPaneDivider trackPaneDivider;
	
	static BasicSplitPaneDivider chromPaneDivider;
	static String annotationfile = "";

	static String searchChrom = "";
	static File ref;

	protected static int canceltextwidth;  
	static JPanel glassPane = new JPanel()
	    {
	      
		private static final long serialVersionUID = 1L;
		private String updateAvail = "Update available (goto File -> Update)";	
		
		
		public void paintComponent(Graphics g) {
			paintComponent((Graphics2D)g);
			
		}
		public void paintComponent(Graphics2D g)
	       {
			g.setRenderingHints(Draw.rh);
				
			
			if(drawCanvas.loading) {
				
			 if(!loaddraw && drawCanvas.loadingtext.length() > 0 && !drawCanvas.loadingtext.equals("note")) {					 
				 loaddraw = true;
				  g.setFont(Draw.loadingFont);
		         	
		          if(Main.loadTextWidth != (int)(g.getFontMetrics().getStringBounds(drawCanvas.loadingtext, g).getWidth())) {
		        	  Main.loadTextWidth = (int)(g.getFontMetrics().getStringBounds(drawCanvas.loadingtext, g).getWidth());
		        	  gradient = new GradientPaint(drawScroll.getWidth()/2-loadTextWidth/2-drawScroll.getWidth()/6,0,Color.red,drawScroll.getWidth()/2+loadTextWidth,0,Color.green,true);			        	 
		        	  Main.canceltextwidth = (int)(g.getFontMetrics().getStringBounds("Cancel", g).getWidth());		        	  
		          }
		         
		          g.setColor(Draw.intronColor);	
		          g.fillRoundRect(drawScroll.getWidth()/2-loadTextWidth/2-9, frame.getHeight()*1/3-Draw.loadingFont.getSize()-4, loadTextWidth+18, Draw.loadingFont.getSize()*3+12,20,20);	
		          g.setColor(Draw.sidecolor);	
		            
		          g.fillRoundRect(drawScroll.getWidth()/2-loadTextWidth/2-5, frame.getHeight()*1/3-Draw.loadingFont.getSize(), loadTextWidth+10, Draw.loadingFont.getSize()*3+4,20,20);					 					  				  
		          g.drawRoundRect(drawScroll.getWidth()/2-loadTextWidth/2-8, frame.getHeight()*1/3-Draw.loadingFont.getSize()-3, loadTextWidth+15, Draw.loadingFont.getSize()*3+9,20,20);	
				     
		          g.setPaint(gradient);					  
				  g.fillRoundRect(drawScroll.getWidth()/2-loadTextWidth/2, frame.getHeight()*1/3+Main.defaultFontSize/2, (int)(loadTextWidth*(drawCanvas.loadBarSample/100.0)), Main.defaultFontSize, Main.defaultFontSize/2,Main.defaultFontSize/2);
				  g.fillRoundRect(drawScroll.getWidth()/2-loadTextWidth/2, frame.getHeight()*1/3+Main.defaultFontSize*2, (int)(loadTextWidth*(drawCanvas.loadbarAll/100.0)), Main.defaultFontSize, Main.defaultFontSize/2,Main.defaultFontSize/2);
				  g.setColor(Color.black);
				  g.drawRoundRect(drawScroll.getWidth()/2-loadTextWidth/2, frame.getHeight()*1/3+Main.defaultFontSize/2, (int)(loadTextWidth*(drawCanvas.loadBarSample/100.0)), Main.defaultFontSize, 4,4);						 
				  g.drawRoundRect(drawScroll.getWidth()/2-loadTextWidth/2, frame.getHeight()*1/3+Main.defaultFontSize*2, (int)(loadTextWidth*(drawCanvas.loadbarAll/100.0)), Main.defaultFontSize, 4,4);
				  					 
				  g.drawString(drawCanvas.loadingtext, drawScroll.getWidth()/2-loadTextWidth/2, frame.getHeight()*1/3);
				  g.drawRect(drawScroll.getWidth()/2-Main.canceltextwidth/2-Main.defaultFontSize/2, frame.getHeight()*1/3+Draw.loadingFont.getSize()*3-Main.defaultFontSize/4, Main.canceltextwidth+Main.defaultFontSize, Draw.loadingFont.getSize()+Main.defaultFontSize/2);					 
				  g.setColor(Color.lightGray);
				  g.fillRect(drawScroll.getWidth()/2-Main.canceltextwidth/2-Main.defaultFontSize/2, frame.getHeight()*1/3+Draw.loadingFont.getSize()*3-Main.defaultFontSize/4, Main.canceltextwidth+Main.defaultFontSize, Draw.loadingFont.getSize()+Main.defaultFontSize/2);					
				  g.setColor(Color.white);
				  g.drawRect(drawScroll.getWidth()/2-Main.canceltextwidth/2-Main.defaultFontSize/2+1, frame.getHeight()*1/3+Draw.loadingFont.getSize()*3-Main.defaultFontSize/4+1, Main.canceltextwidth+Main.defaultFontSize-2, Draw.loadingFont.getSize()+Main.defaultFontSize/2-2);					 
					 
				  if(cancelhover) {

					if(getCursor().getType() != Cursor.HAND_CURSOR) {
						setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
					}
					 
					 g.setColor(Color.gray);
					  g.fillRect(drawScroll.getWidth()/2-Main.canceltextwidth/2-Main.defaultFontSize/2, frame.getHeight()*1/3+Draw.loadingFont.getSize()*3-Main.defaultFontSize/4, Main.canceltextwidth+Main.defaultFontSize, Draw.loadingFont.getSize()+Main.defaultFontSize/2);					
					 
					  g.setColor(Color.white);
					  g.drawRect(drawScroll.getWidth()/2-Main.canceltextwidth/2-Main.defaultFontSize/2, frame.getHeight()*1/3+Draw.loadingFont.getSize()*3-Main.defaultFontSize/4, Main.canceltextwidth+Main.defaultFontSize, Draw.loadingFont.getSize()+Main.defaultFontSize/2);					 
						
				  }
				  else {
					  if(getCursor().getType() != Cursor.WAIT_CURSOR) {
							setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
						}
					  g.setColor(Color.black);
				  }
				  g.drawString("Cancel", drawScroll.getWidth()/2-Main.canceltextwidth/2, frame.getHeight()*1/3+Draw.loadingFont.getSize()*3+Draw.loadingFont.getSize()-4);
				/*  g.setColor(ChromDraw.backTransparent);
					 for(int i=0;i<snow.length; i++) {
							
							 if(snow[i][3] == 0) {
								 snow[i][2]--;
								 if(snow[i][2] < -5) {
									 snow[i][3] = 1;
								 }
							 }
							 else {
								 snow[i][2]++;
								 if(snow[i][2] > 5) {
									 snow[i][3] = 0;
								 }
							 }
							 if(snow[i][1]/2-Math.abs(snow[i][2]) < 1) {
								 snow[i][0]++;
							 }
							 else {
								 snow[i][0]+=snow[i][1]/2 -Math.abs(snow[i][2]);
								
							 }
						
							if(snow[i][0] > height) {
								snow[i][0] = 1;
								snow[i][1] = (4*Math.random() +1);
							}
						//	g.fillOval(i*10+(int)snow[i][2], (int)snow[i][0], (int)snow[i][1], (int)snow[i][1]);
							g.drawString(snowtable[(int)snow[i][1]-1], (int)(i*(frame.getWidth()/(double)snow.length)+(int)snow[i][2]), (int)snow[i][0]);
					 }*/
				  loaddraw = false;
			 }
			 else {
				/* g.setColor(ChromDraw.backTransparent);
				 for(int i=0;i<snow.length; i++) {
						
						 if(snow[i][3] == 0) {
							 snow[i][2]--;
							 if(snow[i][2] < -3) {
								 snow[i][3] = 1;
							 }
						 }
						 else {
							 snow[i][2]++;
							 if(snow[i][2] > 3) {
								 snow[i][3] = 0;
							 }
						 }
						 if(snow[i][1]/2-Math.abs(snow[i][2]) < 1) {
							 snow[i][0]++;
						 }
						 else {
							 snow[i][0]+=snow[i][1]/2 -Math.abs(snow[i][2]);
							
						 }
					//	snow[i][0]+=snow[i][1]/2;
						if(snow[i][0] > height) {
							snow[i][0] = 1;
							snow[i][1] = (4*Math.random() +1);
						}
			//			g.fillOval(i*10+(int)snow[i][2], (int)snow[i][0], (int)snow[i][1], (int)snow[i][1]);
						g.drawString(snowtable[(int)snow[i][1]-1], (int)(i*(frame.getWidth()/(double)snow.length)+(int)snow[i][2]), (int)snow[i][0]);
				}*/
				 g.setFont(Draw.loadingFont);
				 g.setColor(ChromDraw.backTransparent);
				/* if(chromosomeDropdown.getItemCount() == 1 || Launcher.firstStart) {
					 g.drawString("Welcome", drawScroll.getWidth()/2-20, frame.getHeight()*2/3+70);
				 }
				 else {*/
				//	 g.drawString("Loading done.", drawScroll.getWidth()/2-g.getFontMetrics().stringWidth("Loading done.")/2, frame.getHeight()*2/3+70);
					 if(Main.update.isEnabled()) {
						 
						 g.drawString(updateAvail , drawScroll.getWidth()/2-g.getFontMetrics().stringWidth(updateAvail)/2, frame.getHeight()*2/3+70+g.getFontMetrics(Draw.loadingFont).getHeight()*2);
					 }
					
				 //}
			 }					 
			} 
			else {

				  if(getCursor().getType() != Cursor.DEFAULT_CURSOR) {
						setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
					}
			}
	      }					       
	    };

	private static String chooserText;

	static Double vardivider = 0.0;

	static int annotation = 0;

	public static String defaultAnnotation = "";

	static boolean shift = false;	
	
	
public Main() {	
	
super(new GridBagLayout());	
try {
	
	//UIManager.put("PopupMenu.border", BorderFactory.createMatteBorder(0, 20, 0, 0, new Color(230,230,230)));
	//URL fontUrl = new URL("http://www.webpagepublicity.com/" +
     //       "free-fonts/a/Airacobra%20Condensed.ttf");
//		URL fontUrl = new URL("C:/HY-Data/RKATAINE/WinPython-64bit-3.5.3.1Qt5/python-3.5.3.amd64/share/numdifftools/docs/_build/html/_static/fonts/Inconsolata-Regular.ttf");
//	URL fonturl = this.getClass().getResource("OpenSans-Regular.ttf");
//	menuFont = Font.createFont(Font.TRUETYPE_FONT, new File(fonturl.getFile()));
//	C:\HY-Data\RKATAINE\WinPython-64bit-3.5.3.1Qt5\python-3.5.3.amd64\Lib\site-packages\reportlab\fonts
	Launcher.fromMain = true;
	Launcher.main(args);
	VariantHandler.main(argsit);
	glass = Toolkit.getDefaultToolkit().getImage(getClass().getResource("icons/glass.jpg"));
	ToolTipManager.sharedInstance().setInitialDelay(100);
   // ToolTipManager.sharedInstance().setDismissDelay(2000);
    UIManager.put("ToolTip.background", new Color(255, 255, 214)); 
    UIManager.put( "ToolTip.border", BorderFactory.createCompoundBorder( UIManager.getBorder( "ToolTip.border" ), BorderFactory.createEmptyBorder( 4, 4, 4, 4 ) ) ); 
    lineseparator = System.getProperty("line.separator");
    
	panel = new JPanel(new GridBagLayout());
	//menuFont = menuFont.deriveFont(Font.PLAIN,12);
	Draw.defaultFont = menuFont;
	
	gd = GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice();
	width = gd.getDisplayMode().getWidth();
	height = gd.getDisplayMode().getHeight();
	if(Launcher.fontSize.equals("")) {
		if(width < 1500) {	    	
			defaultFontSize = 11;
			
			buttonHeight = Main.defaultFontSize*2;
			buttonWidth = Main.defaultFontSize*6;
			
		}
		else if(width < 2000) {
			defaultFontSize = 12;
						
			buttonHeight = Main.defaultFontSize*2+4;
			buttonWidth = Main.defaultFontSize*6+4;
			
		}
		else if(width < 3000) {
			defaultFontSize = 15;
			buttonHeight = Main.defaultFontSize*2+4;
			buttonWidth = Main.defaultFontSize*6+4;
		}
		else {
			defaultFontSize = 19;
			buttonHeight = Main.defaultFontSize*2+4;
			buttonWidth = Main.defaultFontSize*6+4;
		}
	}
	else {
		try {
			defaultFontSize = Integer.parseInt(Launcher.fontSize);
		}
		catch(Exception e) {
			defaultFontSize = 12;
		}
	}
	
	menuFont = new Font("SansSerif", Font.PLAIN, Main.defaultFontSize);
	menuFontBold = new Font("SansSerif", Font.BOLD, Main.defaultFontSize);
//	menuFont = new Font("SansSerif", Font.BOLD, Main.defaultFontSize);
}
catch(Exception e) {
	e.printStackTrace();
}
FileSystemView fsv = FileSystemView.getFileSystemView();
File[] paths = File.listRoots();

for(File path:paths){
	if(fsv.getSystemDisplayName(path).contains("merit")) {
		pleiades = true;
	}					    			  
}

screenSize = new Dimension(width, height);

drawHeight = (int)(screenSize.getHeight()*0.6);
sidebarWidth =(int)(screenSize.getWidth()*0.1);
drawWidth = (int)(screenSize.getWidth()-sidebarWidth);
thisMainListener = this;
try {
	htsjdk.samtools.util.Log.setGlobalLogLevel(htsjdk.samtools.util.Log.LogLevel.ERROR);
/*	for(int i=0;i<snow.length; i++) {
		snow[i][0] = (height*Math.random());
		snow[i][1] = (4*Math.random() +1);
		snow[i][2] = (12*Math.random() -6);
		snow[i][3] = (2*Math.random() +1);
	}*/
	frame.addWindowListener(new java.awt.event.WindowAdapter() {
	    @Override
	    public void windowClosing(java.awt.event.WindowEvent windowEvent) {
	        /*if (JOptionPane.showConfirmDialog(frame, "Are you sure to close this window?", "Really Closing?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE) == JOptionPane.YES_OPTION){
	            System.exit(0);
	        }*/
	    	if(configChanged) {
	    		
	    		try {
	    		 BufferedWriter fileWriter = new BufferedWriter(new FileWriter(Launcher.configfile));	
	    		 for(int i = 0 ; i<Launcher.config.size(); i++) {
	    			 fileWriter.write(Launcher.config.get(i) +lineseparator);
	    		 }
	    		 fileWriter.close();
	    		 
	    		}
	    		catch(Exception e) {
	    			e.printStackTrace();
	    		}
	    	}
	    	
	    }
	});
	baseMap.put((byte)'A', 1);
	baseMap.put((byte)'C', 2);
	baseMap.put((byte)'G', 3);
	baseMap.put((byte)'T', 4);
	baseMap.put((byte)'N', 5);
	baseMap.put((byte)'I', 6);
	baseMap.put((byte)'D', 7);
	mutTypes.put("TA", 0);
	mutTypes.put("AT", 0);
	mutTypes.put("TC", 1);
	mutTypes.put("AG", 1);
	mutTypes.put("TG", 2);
	mutTypes.put("AC", 2);
	mutTypes.put("CA", 3);
	mutTypes.put("GT", 3);
	mutTypes.put("CG", 4);
	mutTypes.put("GC", 4);
	mutTypes.put("CT", 5);
	mutTypes.put("GA", 5);
	
	getBase.put((byte)'A', "A");
	getBase.put((byte)'C', "C");
	getBase.put((byte)'G', "G");
	getBase.put((byte)'T', "T");
	getBase.put((byte)'N', "N");
	getBase.put((byte)'a', "A");
	getBase.put((byte)'c', "C");
	getBase.put((byte)'g', "G");
	getBase.put((byte)'t', "T");
	getBase.put((byte)'n', "N");
	java.net.URL imgUrl = getClass().getResource("icons/save.gif");
	save = new ImageIcon(imgUrl);
	imgUrl = getClass().getResource("icons/open.gif");
	open = new ImageIcon(imgUrl);
	imgUrl = getClass().getResource("icons/settings.png");
	settingsIcon = new ImageIcon(imgUrl);
	settings = new JMenuItem("Settings", settingsIcon);
	
    
   
   
 //   Average.frame.setVisible(false);
	
	
	userDir = new File(Main.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParent().replace("%20", " ");
	savedir  = Launcher.defaultSaveDir;
	path = Launcher.defaultDir;		
	gerp = Launcher.gerpfile;
	defaultGenome = Launcher.defaultGenome;
	defaultAnnotation = Launcher.defaultAnnotation;
	
	if(Launcher.genomeDir.equals("")) {
		
		genomeDir = new File(userDir +"/genomes/");
	}
	else {
		if(new File(Launcher.genomeDir).exists()) {
			genomeDir = new File(Launcher.genomeDir);
		}		
		else {
			genomeDir = new File(userDir +"/genomes/");
		}
	}
	File[] genomes = genomeDir.listFiles(new FilenameFilter() {
		 public boolean accept(File dir, String name) {
		        return !name.contains(".txt");
		     }
 		 });
	annotationfile = defaultAnnotation;
	controlDir = Launcher.ctrldir;
	trackDir = Launcher.trackDir;
	projectDir = Launcher.projectDir;
	downloadDir = Launcher.downloadDir;
	chromHeight = (int)(drawHeight * 0.1);
	 drawDimensions = new Dimension(drawWidth,drawHeight-chromHeight);
	 bedDimensions = new Dimension(drawWidth, bedHeight);		   
	 chromDimensions = new Dimension(drawWidth-Main.sidebarWidth-1,drawHeight);
	 drawCanvas = new Draw((int)drawDimensions.getWidth(),(int)drawDimensions.getHeight());
	 controlDraw = new ControlCanvas((int)bedDimensions.getWidth(),(int)bedDimensions.getHeight());
	 iconImage = Toolkit.getDefaultToolkit().getImage(getClass().getResource("icons/icon.png"));
	 frame.setIconImage(iconImage);
/*	if(args.length > 0) {
    	for(int i = 0; i<args.length; i++) {
    		if(args[i].startsWith("-opendir")) {
    			path = args[i].substring(9).replace(" ", "");
    		}
    		else if(args[i].startsWith("-ctrldir")) {
    			Control.path = args[i].substring(9).replace(" ", "");
    		}	    		
    	}	    
    }*/
	
//	BGZIPInputStream in = this.getClass().getResourceAsStream("SELEX_1505_representative_matrices.bedhead.gz");
	 searchField.getDocument().addDocumentListener(new DocumentListener() {
		private String searchstring;		

		public void changedUpdate(DocumentEvent e) {
			if(searchField.getText().contains(";")) {
				searchList = searchField.getText().split(";");
				for(int i = 0 ; i<searchList.length; i++) {
					warn(searchList[i].replace(" ", ""));
				}
			}
			else {
				warn(searchField.getText().replace(" ", ""));
			}
		    
		  }
		  public void removeUpdate(DocumentEvent e) {
			  if(searchField.getText().contains(";")) {
					searchList = searchField.getText().split(";");
					for(int i = 0 ; i<searchList.length; i++) {
						warn(searchList[i].replace(" ", ""));
					}
				}
				else {
					warn(searchField.getText().replace(" ", ""));
				}
		  }
		  public void insertUpdate(DocumentEvent e) {
			  if(searchField.getText().contains(";")) {
					searchList = searchField.getText().split(";");
					for(int i = 0 ; i<searchList.length; i++) {
						warn(searchList[i].replace(" ", ""));
					}
				}
				else {
					warn(searchField.getText().replace(" ", ""));
				}
		  }

		  public void warn(String searchtext) {				 
			
			  	if(searchTable.containsKey(searchtext.toUpperCase())) {
			  		if(searchTable.get(searchtext.toUpperCase())[0].equals(Main.chromosomeDropdown.getSelectedItem())) {
			  			searchChrom = searchTable.get(searchtext.toUpperCase())[0];
			  			searchStart = Integer.parseInt(searchTable.get(searchtext.toUpperCase())[1]);
						searchEnd = Integer.parseInt(searchTable.get(searchtext.toUpperCase())[2]);
			  		}
			  		else {
			  			chromDraw.repaint();
			  			searchStart = -1;
			  			searchEnd = -1;
			  		}
			  		chromDraw.repaint();
			  		 searchField.setForeground(Color.black);
			  	}
			  	else if(searchField.getText().toUpperCase().matches("CHR.{1,2}(?!:)")) {
					
					 if(Main.chromnamevector.contains(searchtext.toUpperCase().substring(3))) {
						 searchField.setForeground(Color.black);
						 
					 }
					 else {
						 chromDraw.repaint();
						 searchField.setForeground(Color.red);
					 }
				}
			  	else if(searchtext.toUpperCase().replace(",","").matches("(CHR)?(.+:)?\\d+(-\\d+)?")) {
					 
					 searchField.setForeground(Color.black);
					 if(searchtext.contains(":")) {
						 searchstring = searchtext.substring(searchtext.indexOf(":")+1).replace(",","");
					 }
					 else {
						 chromDraw.repaint();
						 searchstring = searchtext.replace(",","");
					 }
						 
					 if(!searchstring.contains("-")) {
						 try {
							 searchStart = Integer.parseInt(searchstring);
						 }
						 catch(Exception ex) {
							 
						 }
						 searchEnd = -1;
					 }
					 else {
						 try {
							 searchStart = Integer.parseInt(searchstring.substring(0,searchstring.indexOf("-")));
							 searchEnd = Integer.parseInt(searchstring.substring(searchstring.indexOf("-")+1));
						 }
						 catch(Exception ex) {
							 
						 }
					 }
					chromDraw.repaint();
				
				 }
				 else {
					 chromDraw.repaint();
					 searchField.setForeground(Color.red);
					 searchStart = -1;
					 searchEnd = -1;						
				 }				
		  }
		});
	 
	try {
		
		
		
		
		A=Toolkit.getDefaultToolkit().getImage(getClass().getResource("SELEX/A.png"));
		C=Toolkit.getDefaultToolkit().getImage(getClass().getResource("SELEX/C.png"));
		G=Toolkit.getDefaultToolkit().getImage(getClass().getResource("SELEX/G.png"));
		T=Toolkit.getDefaultToolkit().getImage(getClass().getResource("SELEX/T.png"));
		
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	ErrorLog.main(args);
	
	this.setBackground(Color.black);
	UIManager.put("FileChooser.readOnly", Boolean.TRUE);				
	
			
	panel.setBackground(Draw.sidecolor);
	panel.setBorder(BorderFactory.createLineBorder(Color.white));
	searchField.addKeyListener(this);
	
	frame.addKeyListener(this);
	frame.getContentPane().setBackground(Color.black);
	
	
	glassPane.addMouseListener(this);
	glassPane.addMouseMotionListener(new MouseMotionListener() {
		
		@Override
		public void mouseDragged(MouseEvent arg0) {
			
			
		}

		@Override
		public void mouseMoved(MouseEvent event) {			
			// g.drawRect(drawScroll.getWidth()/2-Main.canceltextwidth/2-Main.defaultFontSize/2, Main.drawScroll.getViewport().getHeight()*2/3+Draw.loadingFont.getSize()*3-Main.defaultFontSize/4, Main.canceltextwidth+Main.defaultFontSize, Draw.loadingFont.getSize()+Main.defaultFontSize/2);					 
				
			if(drawCanvas.loading && event.getX() > drawScroll.getWidth()/2-Main.canceltextwidth/2-Main.defaultFontSize/2 && event.getX()  < drawScroll.getWidth()/2+Main.canceltextwidth/2+Main.defaultFontSize/2 && event.getY() > frame.getHeight()*1/3+Draw.loadingFont.getSize()*3-Main.defaultFontSize/4 && event.getY() < frame.getHeight()*1/3+Draw.loadingFont.getSize()*4+Main.defaultFontSize/2) {
				if(!Main.cancelhover) {
					Main.cancelhover = true;						
			    	Main.glassPane.requestFocus();				    	
				}
			}
			else {
				if(Main.cancelhover) {
					Main.cancelhover = false;
					Main.glassPane.requestFocus(false);
				}
			}				
		}			
	});
	
	
	background.put((byte)'A', 0.3);
	background.put((byte)'C', 0.2);
	background.put((byte)'G', 0.2);
	background.put((byte)'T', 0.3);		
	
	bases = new Hashtable<String, String>();
	bases.put("A", "A");
	bases.put("C", "C");
	bases.put("G", "G");
	bases.put("T", "T");
	bases.put("N", "N");	
	bases.put("delA","delA");	
	bases.put("delC", "delC");	
	bases.put("delG", "delG");	
	bases.put("delT","delT");	
	bases.put("insA", "insA");	
	bases.put("insC", "insC");	
	bases.put("insG", "insG");	
	bases.put("insT", "insT");	

	chromDraw = new ChromDraw(drawWidth,chromHeight);	
	VariantCaller.main(argsit);
	
    
   // frame.setExtendedState(JFrame.MAXIMIZED_BOTH); 
   
  
//	path = "C:/HY-Data/RKATAINE/FilterSomatic/";
//	path = "W:/RikuRator-Samples/Thyroid/11-0176_Kilpi1_5_exome";
//	path = "W:/RikuRator-Samples/Myomas/Tumors/";
//	path = "V:/cg8/projects/CRC/c9_5119_N_LP6005134-DNA_A01/wgspipeline/align";
//	path = "V:/cg8/projects/joint-calling/";
//	    drawCanvas.loading("note");
	try {
		
		File annodir;
		
		File[] annotations;			
		addGenome.addMouseListener(this);			
		genome = new JMenu("Genomes");
		genome.setName("genomeMenu");
		genome.add(addGenome);
		genome.addComponentListener(this);
		File[] fastadir;
		String[] empty = {};
		
		refModel = new DefaultComboBoxModel<String>(empty);
		
		refDropdown = new SteppedComboBox(refModel);
		refDropdown.addMouseListener(this);
		String[] emptygenes = {};
		refDropdown.addActionListener(refDropActionListener);
		
		geneModel = new DefaultComboBoxModel<String>(emptygenes);
		geneDropdown = new SteppedComboBox(geneModel);
		geneDropdown.addMouseListener(this);
	if(genomes != null) {	
		for(int i = 0; i<genomes.length; i++) {
			if(!genomes[i].isDirectory()) {
				continue;
			}
			annodir = new File(genomes[i].getAbsolutePath() +"/annotation/");
			if(genomes[i].isDirectory()) {
				fastadir = genomes[i].listFiles();
				for(int f = 0; f<fastadir.length; f++) {
					if(fastadir[f].isDirectory()) {
						continue;
					}
					if(fastadir[f].getName().contains(".fai")) {
						continue;
					}
					else if(fastadir[f].getName().contains(".fa")) {
						fastahash.put(genomes[i].getName(), fastadir[f]);			
					}
				}
			}
		
			annotations = annodir.listFiles();
			genomehash.put(genomes[i].getName(), new ArrayList<File>());
		
			refModel.addElement(genomes[i].getName());
			if(genomes[i].getName().length() > reflength) {
				reflength = genomes[i].getName().length();
			}
			JMenu addMenu = new JMenu(genomes[i].getName());
			addMenu.addMouseListener(this);
			addMenu.setName(genomes[i].getName());
			JMenuItem addAnnotation = new JMenuItem("Add new annotation file...");
			addAnnotation.addMouseListener(this);
			addAnnotation.setName("add_annotation");
			addMenu.add(addAnnotation);
			JLabel addLabel = new JLabel("  Select annotation: ");
			labels.add(addLabel);
			addMenu.add(addLabel);
			addMenu.add(new JSeparator());
			
			genome.add(addMenu);
			addMenu.addComponentListener(this);
			if(annotations != null) {
				for(int j = 0; j<annotations.length; j++) {
					annofiles = annotations[j].listFiles();
					for(int f = 0; f<annofiles.length; f++) {
						if(annofiles[f].getName().endsWith(".bed.gz")) {
							if(annofiles[f].getName().substring(0,annofiles[f].getName().indexOf(".bed.gz")).length() > annolength) {
								annolength = annofiles[f].getName().length();
							}
							
							genomehash.get(genomes[i].getName()).add(annofiles[f].getAbsoluteFile());	
							JMenuItem additem = new JMenuItem(annofiles[f].getName().substring(0,annofiles[f].getName().indexOf(".bed.gz")));
							additem.setName(annofiles[f].getName().substring(0,annofiles[f].getName().indexOf(".bed.gz")));
							additem.addMouseListener(this);
							addMenu.add(additem);
							additem.addComponentListener(this);
							break;
						}							
					}					
				}	
			}
		}
		refModel.addElement("Add new reference...");
		
		
	}
	
	if(genomes.length == 0) {
		/*if(Launcher.firstStart) {
			Main.writeToConfig("FirstStart=false");
		}*/
		AddGenome.createAndShowGUI();
		AddGenome.frame.setTitle("Add new genome");
		
		AddGenome.remove.setEnabled(false);
		AddGenome.download.setEnabled(false);
		
		AddGenome.frame.setLocation((int)(screenSize.getWidth()/2 - AddGenome.frame.getWidth()/2), (int)(screenSize.getHeight()/6));
		 
	    AddGenome.frame.setState(JFrame.NORMAL);
	    AddGenome.frame.setVisible(true);
	    AddGenome.frame.setAlwaysOnTop(true);
	   /*
		WelcomeScreen.main(args);
		WelcomeScreen.frame.setVisible(true);
		WelcomeScreen.frame.setLocation(frame.getLocationOnScreen().x+frame.getWidth()/2 - WelcomeScreen.frame.getWidth()/2, frame.getLocationOnScreen().y+frame.getHeight()/6);
	   */
		if(genomes.length != 0) {
			if(!genomehash.containsKey(defaultGenome)) {
				
				
				setChromDrop(genomes[0].getName());
				defaultGenome = genomes[0].getName();
				}
				else {
					
					setChromDrop(defaultGenome);
				}
				getBands();			
			    getExons();	
		}
		else {
			setChromDrop(null);
		}
	}
	else {
		if(!genomehash.containsKey(defaultGenome)) {
	
		
			setChromDrop(genomes[0].getName());
			defaultGenome = genomes[0].getName();
		
		}
		else {
			
			setChromDrop(defaultGenome);
		}
		getBands();			
	    getExons();	
	}
	
	if(Launcher.firstStart) {
			
		
		WelcomeScreen.createAndShowGUI();
		WelcomeScreen.frame.setLocation((int)(screenSize.getWidth()/2 - WelcomeScreen.frame.getWidth()/2), (int)(screenSize.getHeight()/6));
		WelcomeScreen.frame.setVisible(true);
	}
	    setMenuBar();
	    setButtons();
	    Settings.main(args);
	   
	  // Settings.main(args);
		frame.requestFocus();	   
	    
		drawCanvas.addKeyListener(this);
		bedCanvas.addKeyListener(this);
		setFonts();
		CheckUpdates check = new CheckUpdates();
		check.execute();
	//	Main.drawCanvas.loading("test");
		Main.drawCanvas.splits.get(0).setCytoImage(Main.chromDraw.createBands(Main.drawCanvas.splits.get(0)));		
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	
}
catch(Exception ex) {
	ex.printStackTrace();
	Main.showError(ex.getMessage(), "Error");
	
}

}
static void getExons() {
	try {
		String s;
		String[] exonSplit;
		Boolean found = false;
		if(ChromDraw.exonReader != null) {
			ChromDraw.exonReader.close();
		}
		for(int i = 0; i<genomehash.get(defaultGenome).size(); i++) {
			if(genomehash.get(defaultGenome).get(i).getName().equals(defaultAnnotation)) {
				ChromDraw.exonReader = new TabixReader(genomehash.get(defaultGenome).get(i).getCanonicalPath());
				
				annotationfile = genomehash.get(defaultGenome).get(i).getName();
				Main.annotation = i;
				found = true;
				break;
			}
			
		}
		if(!found) {
			
			if(genomehash.get(defaultGenome).size() > 0) {
				ChromDraw.exonReader = new TabixReader(genomehash.get(defaultGenome).get(0).getCanonicalPath());
				annotationfile = genomehash.get(defaultGenome).get(0).getName();
				Main.annotation = 0;
			}
			else {
				annotationfile = "";
			}
		}
		if(ChromDraw.exonReader == null) {
			return;
		}
		searchTable.clear();
		while((s = ChromDraw.exonReader.readLine()) != null) {
			if(s.startsWith("#")) {
				continue;
			}
			exonSplit = s.split("\t");
			
			if(!searchTable.containsKey(exonSplit[3].toUpperCase())) {		
				
				String[] adder = {exonSplit[0], exonSplit[1], exonSplit[2]};
				searchTable.put(exonSplit[3].toUpperCase(), adder);				
				if(exonSplit[6].contains(":")) {
					geneIDMap.put(exonSplit[6].split(":")[1].toUpperCase(), exonSplit[3].toUpperCase());
				}
				else {
					geneIDMap.put(exonSplit[6].toUpperCase(), exonSplit[3].toUpperCase());
				}
				
				
			}	
			else {
				try {
					if(Integer.parseInt(searchTable.get(exonSplit[3].toUpperCase())[1]) > Integer.parseInt(exonSplit[1])) {
						searchTable.get(exonSplit[3].toUpperCase())[1] = exonSplit[1];
					}
					if(Integer.parseInt(searchTable.get(exonSplit[3].toUpperCase())[2]) < Integer.parseInt(exonSplit[2])) {
						searchTable.get(exonSplit[3].toUpperCase())[2] = exonSplit[2];
					}
				}
				catch(Exception e) {
					
					System.out.println("error: " +searchTable.get(exonSplit[3].toUpperCase())[1]);
				}
			}
		}	
		
		ChromDraw.exonReader.close();
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}
void setMenuBar() {
	//filemenu.addMouseListener(this);
	
	//toolmenu.addMouseListener(this);
	filemenu = new JMenu("File");
	toolmenu = new JMenu("Tools");
	help = new JMenu("Help");
	about = new JMenu("About");
	menubar = new JMenuBar();
	//help.addMouseListener(this);
	exit = new JMenuItem("Exit");
//	opensamples = new JMenuItem("Add samples");
	zoomout = new JButton("Zoom out");
	back = new JButton("<<");
	forward = new JButton(">>");
	manage = new JButton("Variant Manager");
	opensamples = new JMenuItem("Add samples", open);
	average = new JMenuItem("Coverage calculator");
	update = new JMenuItem("Update");
	update.setVisible(false);
	errorlog = new JMenuItem("View log");
	helpLabel = new JLabel("This is pre-release version of BasePlayer\nContact: help@baseplayer.fi\n\nUniversity of Helsinki");
	
	addURL = new JMenu("Add from URL");
	urlField = new JTextField("Enter URL");
	addtracks = new JMenuItem("Add tracks");
	addcontrols = new JMenuItem("Add controls");
	pleiadesButton = new JMenuItem("PLEIADES");
	saveProject = new JMenuItem("Save project");
	saveProjectAs = new JMenuItem("Save project as...");
	openProject = new JMenuItem("Open project");
	clear = new JMenuItem("Clear data");
	clearMemory = new JMenuItem("Clean memory");
//	welcome = new JMenuItem("Welcome screen");
	filemenu.add(opensamples);
	variantCaller = new JMenuItem("Variant Caller");
	addtracks = new JMenuItem("Add tracks", open);
	filemenu.add(addtracks);
	addcontrols = new JMenuItem("Add Controls", open);
	filemenu.add(addcontrols);
	if(pleiades) {
		pleiadesButton.setPreferredSize(buttonDimension);
		pleiadesButton.addActionListener(this);
		
		filemenu.add(pleiadesButton);
	}
	
	filemenu.add(new JSeparator());
	openProject = new JMenuItem("Open project", open);
	filemenu.add(openProject);
	saveProject = new JMenuItem("Save project", save);
	filemenu.add(saveProject);
	saveProjectAs = new JMenuItem("Save project as...", save);
	filemenu.add(saveProjectAs);		
	filemenu.add(new JSeparator());
	filemenu.add(genome);
	filemenu.add(update);
	filemenu.add(clear);
	filemenu.add(new JSeparator());
	filemenu.add(exit);
	exit.addActionListener(this);
	menubar.add(filemenu);
	manage.addActionListener(this);
	manage.addMouseListener(this);
	update.addActionListener(this);
	average.addActionListener(this);
	average.setEnabled(false);
	average.setToolTipText("No bam/cram files opened");
	toolmenu.add(average);
	toolmenu.add(variantCaller);
	variantCaller.setToolTipText("No bam/cram files opened");
	variantCaller.addActionListener(this);
	variantCaller.setEnabled(false);
	settings.addActionListener(this);
	clearMemory.addActionListener(this);
	errorlog.addActionListener(this);
	toolmenu.add(clearMemory);
	toolmenu.add(errorlog);
	toolmenu.add(new JSeparator());
	toolmenu.add(settings);
	menubar.add(toolmenu);
	menubar.add(manage);
//	help.addActionListener(this);
	
	
	
//	welcome.addActionListener(this);
//	JMenuItem infotable = new JMenuItem(" <html> Line1 <br/> Line2 <br/> Line3 </html> ");
	
//	JLabel aboutText = new JLabel();
	//aboutText.setText (" <html> Line1 <br/> Line2 <br/> Line3 </html> ");
	//infotable.add(aboutText);
	//about.add(helpLabel);
	area = new JEditorPane();
	

//	area.setEditorKit(javax.swing.JEditorPane.createEditorKitForContentType("text/<b>html</b>"));
	String infotext = "<html><h2>BasePlayer</h2>This is a pre-release version of BasePlayer (<a href=https://baseplayer.fi>https://baseplayer.fi</a>)<br/> Author: Riku Katainen <br/> University of Helsinki<br/>"
					+"Tumor Genomics Group (<a href=http://research.med.helsinki.fi/gsb/aaltonen/>http://research.med.helsinki.fi/gsb/aaltonen/</a>) <br/> " 
					+"Contact: help@baseplayer.fi <br/>"
					+"See the BasePlayer manuscript draft at <a href=https://doi.org/10.1101/126482>bioRxiv</a> <br/><br/>"
					+"Supported filetype for variants is VCF and VCF.gz (index file will be created if missing)<br/> "
					+"Supported filetypes for reads are BAM and CRAM. Index files required (.bai or .crai). <br/> "
					+"Supported filetypes for additional tracks are BED(.gz), GFF.gz, BedGraph, BigWig, BigBed.<br/> (tabix index required for bgzipped files). <br/><br/> "
											
					+"For optimal usage, you should have vcf.gz and bam -files for each sample. <br/> "
					+"e.g. in case you have a sample named as sample1, name all files similarly and <br/>"
					+"place in the same folder:<br/>"
					+"sample1.vcf.gz<br/>"
					+"sample1.vcf.gz.tbi<br/>"
					+"sample1.bam<br/>"
					+"sample1.bam.bai<br/><br/>"
					+"When you open sample1.vcf.gz, sample1.bam is recognized and opened<br/>" 
					+"on the same track.<br/><br/>"
					+"Instructional videos can be viewed at our <a href=https://www.youtube.com/channel/UCywq-T7W0YPzACyB4LT7Q3g> Youtube channel</a>"; 
	area = new JEditorPane();
	
	area.setEditable(false);
	area.setEditorKit(JEditorPane.createEditorKitForContentType("text/html"));
	
	area.setText(infotext);
	area.setFont(Main.menuFont);
	area.addHyperlinkListener(new HyperlinkListener() {
		public void hyperlinkUpdate(HyperlinkEvent hyperlinkEvent) {
		    HyperlinkEvent.EventType type = hyperlinkEvent.getEventType();
		    final URL url = hyperlinkEvent.getURL();
		    if (type == HyperlinkEvent.EventType.ACTIVATED) {
		    	Main.gotoURL(url.toString());
		    }
		  }
	});
	
	
	about.add(area);
	about.addMouseListener(this);
	help.add(about);
//	help.add(welcome);
	
	menubar.add(help);
	JLabel emptylab = new JLabel("  ");
	emptylab.setEnabled(false);		
	emptylab.setOpaque(false);
	//	emptylab.setBackground(new Color(250,250,250));
	menubar.add(emptylab);
	
	/*JMenuItem empty = new JMenuItem();
	empty.setEnabled(false);		
	empty.setBackground(new Color(250,250,250));
	menubar.add(empty);*/
//	empty.setMaximumSize(empty.getPreferredSize());
	
	chromosomeDropdown.setBorder(BorderFactory.createMatteBorder(1, 0, 1, 1, Color.lightGray));
	chromosomeDropdown.setBorder(BorderFactory.createCompoundBorder(
			chromosomeDropdown.getBorder(), 
	        BorderFactory.createEmptyBorder(0, 0, 0, 0)));
	
	 chromlabel.setToolTipText("Current chromosome");
	 chromlabel.setFocusable(false);
	 chromlabel.addMouseListener(this);
	 chromlabel.setBackground(Color.white);
    chromlabel.setEditable(false);
    chromlabel.setBorder(BorderFactory.createMatteBorder(1, 1, 1, 0, Color.lightGray));
    chromlabel.setBorder(BorderFactory.createCompoundBorder(
    		chromlabel.getBorder(), 
	        BorderFactory.createEmptyBorder(0, 0, 0, 0)));
	menubar.add(chromlabel);	
	chromosomeDropdown.setBackground(Color.white);
	chromosomeDropdown.setToolTipText("Current chromosome");
    menubar.add(chromosomeDropdown);
    JLabel empty3 = new JLabel("  ");
	empty3.setEnabled(false);
	empty3.setOpaque(false);
	menubar.add(empty3);
	
	menubar.add(back);
	menubar.add(searchField);
	searchField.setForeground(Color.gray);
	
	searchField.setBorder(BorderFactory.createCompoundBorder(
			searchField.getBorder(), 
	        BorderFactory.createEmptyBorder(0, 0, 0, 0)));
	
	searchField.addMouseListener(this);
    menubar.add(back);
    menubar.add(searchField);
    searchField.setForeground(Color.gray);
  
    back.addMouseListener(this);
    back.setToolTipText("Back");
    forward.addMouseListener(this);
    forward.setToolTipText("Forward");
    back.setEnabled(false);
    forward.setEnabled(false);
    
    searchField.addMouseListener(this);

    menubar.add(back);
    menubar.add(searchField);
    searchField.setForeground(Color.gray);
    back.addMouseListener(this);
    forward.addMouseListener(this);
    back.setEnabled(false);
    forward.setEnabled(false);
    forward.setMargin(new Insets(0,2,0,2));
    back.setMargin(new Insets(0,2,0,2));
    menubar.add(forward);
    JLabel empty4 = new JLabel("  ");
    empty4.setOpaque(false);
	empty4.setEnabled(false);		
	menubar.add(empty4);
    menubar.add(zoomout);
    JLabel empty5 = new JLabel("  ");
	empty5.setEnabled(false);		
	empty5.setOpaque(false);
	menubar.add(empty5);	   
    positionField.setEditable(false);
    positionField.setBackground(new Color(250,250,250));
   
    positionField.setMargin(new Insets(0,2,0,0));
    positionField.setBorder(BorderFactory.createCompoundBorder(
    		widthLabel.getBorder(), 
	        BorderFactory.createEmptyBorder(0, 0, 0, 0)));
    menubar.add(positionField);	   
    widthLabel.setEditable(false);
    widthLabel.setBackground(new Color(250,250,250));
    widthLabel.setMargin(new Insets(0,2,0,0));
    widthLabel.setBorder(BorderFactory.createCompoundBorder(
    		widthLabel.getBorder(), 
	        BorderFactory.createEmptyBorder(0, 0, 0, 0)));  
    JLabel empty6 = new JLabel("  ");
    empty6.setEnabled(false);		
    empty6.setOpaque(false);
	menubar.add(empty6);		
    menubar.add(widthLabel);
    JLabel empty7 = new JLabel("  ");
    empty7.setOpaque(false);
    empty7.setEnabled(false);		
	menubar.add(empty7);	
}
void setButtons() {
	try {
	
	
	  //      filemenu.setFont(font);
	GridBagConstraints c = new GridBagConstraints();	
	
	c.anchor = GridBagConstraints.NORTHWEST;
	c.insets = new Insets(1,4,4,2);
	//c.insets = new Insets(5, 5, 2, 5);
	c.gridx = 0;
	c.gridy = 0;
	c.gridwidth = 1;  
/*	opensamples.setMargin(new Insets(0, 0, 0,0));
	
	opensamples.setPreferredSize(buttonDimension);
	addtracks.setMargin(new Insets(0, 0, 0, 0));
	addtracks.setPreferredSize(buttonDimension);
	addcontrols.setMargin(new Insets(0, 0, 0, 0));
	addcontrols.setPreferredSize(buttonDimension);
	*/
	
	menubar.setOpaque(true);
	panel.add(menubar, c);
	c.gridx = 1;
	
	setbut = new JButton("", settingsIcon);
	setbut.setToolTipText("Settings");
	
	setbut.setOpaque(false);
	setbut.setContentAreaFilled(false);
	setbut.setBackground(Main.panel.getBackground());
	setbut.addMouseListener(new MouseListener() {

		@Override
		public void mouseClicked(MouseEvent e) {
							
		}

		@Override
		public void mousePressed(MouseEvent e) {
			Settings.frame.setLocation(frame.getLocationOnScreen().x+frame.getWidth()/2 - Settings.frame.getWidth()/2, frame.getLocationOnScreen().y+frame.getHeight()/6);			
			Settings.frame.setState(JFrame.NORMAL);
			Settings.frame.setVisible(true);
		}

		@Override
		public void mouseReleased(MouseEvent e) {
			
		}

		@Override
		public void mouseEntered(MouseEvent e) {
			setbut.setOpaque(true);
			Main.setbut.setBackground(Color.white);
			Main.setbut.revalidate();
		}

		@Override
		public void mouseExited(MouseEvent e) {
			setbut.setOpaque(false);	
			Main.setbut.revalidate();
		}
		
	});
	setbut.setBorder(null);
	c.insets = new Insets(0,2,0,0);
	menubar.add(setbut, c);
//	c.gridx = 1;
//	c.gridx = 2;		
//	c.gridx = 3;		
	
	
//		zoomout.setMargin(new Insets(0, 0, 0, 0));
//	panel.add(zoomout, c);		
	
	
/*	JMenuItem empty2 = new JMenuItem("");
	empty2.setEnabled(false);		
	menubar.add(empty2);
  
    c.gridx = 4;*/
  //  chromosomeDropdown.setPreferredSize(buttonDimension);
   
   // panel.add(chromosomeDropdown, c);
//	    c.gridx = 5;
   
    
   
    
  //  back.setMargin(new Insets(0, 0, 0, 0));
   // forward.setMargin(new Insets(0, 0, 0, 0));
   // back.setPreferredSize(new Dimension(back.getFontMetrics(back.getFont()).stringWidth("<")+10,(int)fieldDimension.getHeight()));
//    forward.setPreferredSize(new Dimension(forward.getFontMetrics(forward.getFont()).stringWidth("<")+10,(int)fieldDimension.getHeight()));
    //back.setMinimumSize(new Dimension(back.getFontMetrics(back.getFont()).stringWidth("<")+10,(int)fieldDimension.getHeight()));
//    forward.setMinimumSize(new Dimension(forward.getFontMetrics(forward.getFont()).stringWidth("<")+10,(int)fieldDimension.getHeight()));
//    panel.add(back, c);
//	    c.gridx = 7;
 //   panel.add(forward, c);
    
    
    chromosomeDropdown.setMaximumRowCount(25);	   
	chromosomeDropdown.setEnabled(true);
	chromosomeDropdown.addActionListener(ChromoDropActionListener);	
	chromosomeDropdown.addMouseListener(this);
	
	c.gridwidth = 10;  
    c.gridx = 0;
    c.gridy = 1;
    bedScroll = new JScrollPane(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
    bedScroll.getViewport().setPreferredSize(bedDimensions);
    
    drawScroll = new JScrollPane(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);	  	   
    drawScroll.getViewport().setPreferredSize(drawDimensions);	   	   
    
    chromScroll = new JScrollPane(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
    chromScroll.getViewport().setPreferredSize(chromDimensions);
   
	drawScroll.getVerticalScrollBar().setAutoscrolls(false);
//	    chromScroll.getViewport().setPreferredSize(new Dimension(drawWidth,chromHeight-20));
   
   
//    drawScroll.setBorder(BorderFactory.createEmptyBorder());
    //chromScroll.setBorder(BorderFactory.createLoweredBevelBorder());
//    drawScroll.setBorder(BorderFactory.createLoweredBevelBorder());
//    bedScroll.setBorder(BorderFactory.createLoweredBevelBorder());
    
 //   chromScroll.setBorder(BorderFactory.createEmptyBorder());
    
//	    bedScroll.setBorder(BorderFactory.createEmptyBorder());
    
    controlScroll = new JScrollPane(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);	  	   
    controlScroll.getViewport().setPreferredSize(bedDimensions);	 
    controlScroll.getViewport().add(controlDraw);
    controlDraw.setVisible(false);
    controlScroll.setVisible(false);
    
    chromScroll.setBorder(BorderFactory.createEmptyBorder());
    drawScroll.setBorder(BorderFactory.createEmptyBorder());
    bedScroll.setBorder(BorderFactory.createLoweredBevelBorder());
    controlScroll.setBorder(BorderFactory.createLoweredBevelBorder());
   
//	    controlScroll.setBorder(BorderFactory.createEmptyBorder());
    
   /* chromScroll.getVerticalScrollBar().setUI(new BasicScrollBarUI()  {   
    @Override
    protected void configureScrollBarColors(){
        this.thumbColor = Color.black;
    }

    });
    
    drawScroll.getVerticalScrollBar().setUI(new BasicScrollBarUI()  {   
        @Override
        protected void configureScrollBarColors(){
            this.thumbColor = Color.black;
        }

	    });
	    */
    addSplit(chromosomeDropdown.getItemAt(0));
    
    //drawCanvas.splits.add(new SplitClass());
    //drawCanvas.splits.get(0).chrom = Main.chromosomeDropdown.getSelectedItem().toString();
    //drawCanvas.splits.get(0).transcripts = chromDraw.transcripts;
    
//    addSplit("10");
 //   addSplit("20");
   
  //  drawCanvas.splits.get(1).cytoImage;
    
/*    drawCanvas.splits.add(new SplitClass());	
    drawCanvas.splits.get(2).chrom ="20";
    drawCanvas.splits.get(2).transcripts = fileReader.getExons(drawCanvas.splits.get(2).chrom);
    drawCanvas.splits.get(2).chromEnd = chromIndex.get(drawCanvas.splits.get(2).chrom)[1].intValue();
    drawCanvas.splits.get(2).start = 1;
    drawCanvas.splits.get(2).end =  drawCanvas.splits.get(2).chromEnd;
    drawCanvas.splits.get(2).viewLength = drawCanvas.splits.get(2).end-drawCanvas.splits.get(2).start;
  */ // drawCanvas.splits.get(2).cytoImage = chromDraw.cytoImage;
    
   // drawCanvas.resizeCanvas(drawWidth,drawHeight);
    
//    drawCanvas.splits.get(1).pixel = (drawCanvas.getDrawWidth())/(double)(drawCanvas.splits.get(1).viewLength);
 //   drawCanvas.splits.get(2).pixel = (drawCanvas.getDrawWidth())/(double)(drawCanvas.splits.get(2).viewLength);
    
   // chromDraw.setPreferredSize(chromDimensions);	
 //   chromScroll.addMouseListener(this);
  //  chromDraw.addMouseListener(this);
   
   
    chromScroll.getViewport().add(chromDraw);	    
    drawScroll.getViewport().add(drawCanvas);	
    drawScroll.addMouseListener(this);
//	    drawCanvas.addMouseListener(this);
    bedCanvas = new BedCanvas(drawWidth, 200);
   
    bedScroll.getViewport().add(bedCanvas);
   
    frame.setExtendedState(frame.getExtendedState() | JFrame.MAXIMIZED_BOTH);
    
	c.weightx = 1.0;
	c.weighty = 1.0;
	c.fill = GridBagConstraints.BOTH;
//		c.gridy = 5;
	
	trackPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, bedScroll, controlScroll);
	trackPane.setUI(new BasicSplitPaneUI() {
        public BasicSplitPaneDivider createDefaultDivider() {
        return new BasicSplitPaneDivider(this) {
           
			private static final long serialVersionUID = 1L;

			public void setBorder(Border b) {
            }

            @Override
                public void paint(Graphics g) {
                g.setColor(Color.lightGray);
                g.fillRect(0, 0, getSize().width, getSize().height);
                    super.paint(g);
                }
        };
        }
    });
	trackPane.setBorder(null);
	
	trackPane.setDividerSize(0);
	trackPane.setPreferredSize(drawDimensions);
	trackPane.setResizeWeight(0.0);
	trackPane.setContinuousLayout(true);
	trackPane.setVisible(false);
/*	varmaster = new VarMaster((int)bedDimensions.getWidth(),(int)bedDimensions.getHeight());
	drawpane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, varmaster, drawScroll);
	drawpane.setDividerSize(2);
	drawpane.setUI(new BasicSplitPaneUI() {
        public BasicSplitPaneDivider createDefaultDivider() {
        return new BasicSplitPaneDivider(this) {
           
			private static final long serialVersionUID = 1L;

			public void setBorder(Border b) {
            }

            @Override
                public void paint(Graphics g) {
                g.setColor(Color.black);
                g.fillRect(0, 0, getSize().width, getSize().height);
                    super.paint(g);
                }
        };
        }
    });*/
	varpane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, trackPane, drawScroll);
	varpane.setUI(new BasicSplitPaneUI() {
        public BasicSplitPaneDivider createDefaultDivider() {
        return new BasicSplitPaneDivider(this) {
           
			private static final long serialVersionUID = 1L;

			public void setBorder(Border b) {
            }

            @Override
                public void paint(Graphics g) {
                g.setColor(Color.lightGray);
                g.fillRect(0, 0, getSize().width, getSize().height);
                    super.paint(g);
                }
        };
        }
    });
	varpane.setBorder(null);
	
	varpane.setDividerSize(0);
	varpane.setPreferredSize(drawDimensions);
	varpane.setResizeWeight(0.0);
	varpane.setContinuousLayout(true);
//	varpane.addMouseListener(this);
	bedScroll.setVisible(false);
	
//	trackPane.addMouseListener(this);
	
	
	controlScroll.setVisible(false);
	
	//pan.setPreferredSize(new Dimension(Main.sidebarWidth, Main.chromDimensions.height));
//	upPanel.setPreferredSize(new Dimension(Main.chromDimensions.width +Main.sidebarWidth, Main.chromDimensions.height));
//	upPanel.add(pan);
//	upPanel.add(chromScroll);
	chrompan = new JPanel() {
		private static final long serialVersionUID = 1L;
		
		protected void paintComponent(Graphics g) {
		       
			 super.paintComponent(g);     
		        g.setColor(Draw.sidecolor);
		        g.fillRect(0, 0, this.getWidth(),  this.getHeight());	        
		       
		        g.setColor(Color.gray);
		    	g.fillRect(0, 0, 3,this.getHeight());
		    	g.setColor(Color.lightGray);
		    	g.fillRect(2, 0, 2,this.getHeight());
		    	
		    }		
		
	};		
	chrompan.setLayout(new GridBagLayout());
	GridBagConstraints gb = new GridBagConstraints();	
	gb.anchor = GridBagConstraints.NORTHWEST;
	
	gb.insets = new Insets(2,10,2,2);
	gb.gridx = 0;
	gb.gridy = 0;
	gb.gridwidth = 1;	
	gb.fill = GridBagConstraints.HORIZONTAL;
	refDropdown.setBackground(Color.white);
	refDropdown.setBorder(BorderFactory.createMatteBorder(1, 1, 1, 1, Color.lightGray));
	refDropdown.setBorder(BorderFactory.createCompoundBorder(refDropdown.getBorder(), BorderFactory.createEmptyBorder(0, 0, 0, 0)));
	geneDropdown.setBackground(Color.white);
	geneDropdown.setBorder(BorderFactory.createMatteBorder(1, 1, 1, 1, Color.lightGray));
	geneDropdown.setBorder(BorderFactory.createCompoundBorder(geneDropdown.getBorder(), BorderFactory.createEmptyBorder(0, 0, 0, 0)));
	geneDropdown.addActionListener(annoDropActionListener);
	JLabel refLabel = new JLabel("Reference genome:");
	
	JLabel geneLabel = new JLabel("Gene annotation:");
	chromLabel.setName("header");
	chrompan.add(chromLabel,gb);
	gb.gridy++;
	chrompan.add(new JSeparator(),gb);
	gb.gridy++;
	gb.insets = new Insets(0,10,0,2);
	chrompan.add(refLabel,gb);
	gb.gridy++;
	chrompan.add(refDropdown,gb);
	gb.gridy++;
	chrompan.add(geneLabel,gb);
	gb.gridy++;
	chrompan.add(geneDropdown,gb);
	gb.gridy++;
	gb.insets = new Insets(20,10,0,2);
	JLabel memory = new JLabel("Memory usage:");
	memory.setName("header");
	chrompan.add(memory,gb);
	gb.insets = new Insets(0,10,0,2);
	gb.gridy++;
	chrompan.add(memLabel, gb);
	gb.weightx = 1;
	gb.weighty = 1;
	gb.gridwidth = GridBagConstraints.REMAINDER;
	chrompan.add(new JLabel(),gb);
	chrompan.setMinimumSize(new Dimension(1,1));
	chrompan.addComponentListener(this);
	upPanel = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, chrompan, chromScroll);
	drawScroll.addComponentListener(this);
	//upPanel.setBorder(BorderFactory.createLoweredBevelBorder());
	upPanel.setBorder(BorderFactory.createMatteBorder(0, 0,1, 0, Color.white));
	upPanel.setDividerLocation(Main.sidebarWidth-2);
	chrompan.setBackground(Draw.sidecolor);
	BasicSplitPaneUI chromPaneUI = (BasicSplitPaneUI) upPanel.getUI();
	chromPaneDivider = chromPaneUI.getDivider();
	chromPaneDivider.addMouseListener(this);
	
	upPanel.setDividerSize(3);
	
	upPanel.setUI(new BasicSplitPaneUI() {
        public BasicSplitPaneDivider createDefaultDivider() {
        return new BasicSplitPaneDivider(this) {
           
			private static final long serialVersionUID = 1L;

			public void setBorder(Border b) {
            }

            @Override
                public void paint(Graphics g) {
                g.setColor(Color.lightGray);
                g.fillRect(0, 0, getSize().width, getSize().height);
               
                    super.paint(g);
                }
        };
        }
    });
	
	
	
	
	splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, upPanel, varpane);
	
	splitPane.setUI(new BasicSplitPaneUI() {
        public BasicSplitPaneDivider createDefaultDivider() {
        return new BasicSplitPaneDivider(this) {
           
			private static final long serialVersionUID = 1L;

			public void setBorder(Border b) {
            }

            @Override
                public void paint(Graphics g) {
                g.setColor(Color.lightGray);
                g.fillRect(0, 0, getSize().width, getSize().height);
               
                    super.paint(g);
                }
        };
        }
    });
   // splitPane.setBorder(null);
	//splitPane.setResizeWeight(0.5);
	
    BasicSplitPaneUI basicSplitPaneUI = (BasicSplitPaneUI) splitPane.getUI();
    splitPaneDivider = basicSplitPaneUI.getDivider();
//    splitPaneDivider.addMouseListener(this);
    basicSplitPaneUI = (BasicSplitPaneUI) trackPane.getUI();
    trackPaneDivider = basicSplitPaneUI.getDivider();
//    trackPaneDivider.addMouseListener(this);
    BasicSplitPaneUI splitPaneUI = (BasicSplitPaneUI) varpane.getUI();
    varPaneDivider = splitPaneUI.getDivider();
 //   varPaneDivider.addMouseListener(this);
	splitPane.setDividerSize(3);
	
/*	BasicSplitPaneDivider divider = (BasicSplitPaneDivider) splitPane.getComponent(2);
	divider.setBackground(Color.black);
	*/
	splitPane.setPreferredSize(drawDimensions);
	
	
	splitPane.setContinuousLayout(true);

    panel.add(splitPane, c);
    add(panel, c);
   
    opensamples.addActionListener(this);
    addtracks.addActionListener(this);
    addcontrols.addActionListener(this);
   
    openProject.addActionListener(this);
    saveProject.addActionListener(this);
    saveProjectAs.addActionListener(this);
    dosomething.addActionListener(this);
    clear.addActionListener(this);
//    drawScroll.getVerticalScrollBar().addMouseMotionListener(this);
    drawScroll.getVerticalScrollBar().addMouseListener(this);
    drawScroll.getVerticalScrollBar().addMouseMotionListener(this);	 
    drawScroll.getVerticalScrollBar().addMouseWheelListener(new MouseWheelListener() {
		@Override
		public void mouseWheelMoved(MouseWheelEvent e) {
			Draw.setGlasspane(true);
			if ( e.getWheelRotation() < 0 ) {
				if(drawCanvas.drawVariables.visiblestart > 0) {
					drawCanvas.drawVariables.visiblestart--; 
					
				}
				Draw.setScrollbar((int)(drawCanvas.drawVariables.visiblestart*drawCanvas.drawVariables.sampleHeight));
			}
			else {
				if(drawCanvas.drawVariables.visiblestart+drawCanvas.drawVariables.visiblesamples < Main.samples) {
					drawCanvas.drawVariables.visiblestart++; 
				
				}
				Draw.setScrollbar((int)(drawCanvas.drawVariables.visiblestart*drawCanvas.drawVariables.sampleHeight));
			}
			
		}
    	
    }
    		
    		
    );
    drawScroll.getVerticalScrollBar().addAdjustmentListener(new AdjustmentListener() {
    	
    	@Override
		public void adjustmentValueChanged(AdjustmentEvent event) {
    		
			//System.out.println(drawCanvas.drawVariables.visiblestart +" " +(short)(Main.drawScroll.getVerticalScrollBar().getValue()/drawCanvas.drawVariables.sampleHeight));
    		if(drawCanvas.drawVariables.visiblestart != (short)(Main.drawScroll.getVerticalScrollBar().getValue()/drawCanvas.drawVariables.sampleHeight)) {
    			if(!drawCanvas.sidebar) {
    				drawCanvas.drawVariables.visiblestart = (short)(Main.drawScroll.getVerticalScrollBar().getValue()/drawCanvas.drawVariables.sampleHeight);					
    			}
				if(drawCanvas.splits.size() > 1) {
					for(int i = 0; i<drawCanvas.splits.size(); i++) {
						drawCanvas.splits.get(i).updateReads = true;					
					}
				}
				else {
					Draw.updateReads = true;
					Draw.updatevars = true;
				}			
			
				Draw.updatevars = true;
				Main.drawCanvas.repaint();
			} 	
    		
		}   
    	
    });
   /* chromScroll.getVerticalScrollBar().addAdjustmentListener(new AdjustmentListener() {
    	
    	@Override
		public void adjustmentValueChanged(AdjustmentEvent event) {
			
				if(!drawCanvas.scrolldrag) {
					System.out.println("jou");
				}
			
		}    	
    	
    });
  */
    zoomout.addActionListener(this);
     
    FileRead.head = new VarNode(0,  (byte)0,"N", (short)0, (short)0, false,(float)0,(float)0,null,null, new Sample("",(short)1,null), null, null);     
    drawCanvas.current = FileRead.head;
 //   upPanel.addPropertyChangeListener(this);
  //  splitPane.addPropertyChangeListener(this);
   // trackPane.addPropertyChangeListener(this);
   // varpane.addPropertyChangeListener(this);	
    
    frame.addComponentListener(this);
    frame.addMouseListener(this);	    	   
    frame.setGlassPane(glassPane);
    /*if(chromosomeDropdown.getItemCount() > 0) {
    	chromosomeDropdown.setSelectedIndex(0);
    	drawCanvas.chrom = chromosomeDropdown.getItemAt(0);
    }*/
    glassPane.setOpaque(false);
    glassPane.setVisible(false);
   
    positionField.setText("chr1:1-" +MethodLibrary.formatNumber(drawCanvas.splits.get(0).chromEnd));
    positionField.setToolTipText("Current chromosomal region");
    widthLabel.setText(MethodLibrary.formatNumber(drawCanvas.splits.get(0).chromEnd) +"bp");
    widthLabel.setToolTipText("Current region width in base pairs");
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}

static void addSplit(String chrom) {
	
	    SplitClass split = new SplitClass();
	    drawCanvas.selectedSplit = split;
	    split.chrom = chrom;
	    if(drawCanvas.splits.size() > 0 && drawCanvas.selectedSplit.chrom.equals(chrom)) {
	    	split.setGenes(drawCanvas.selectedSplit.getGenes());
	    }
	    else {
	    	
	    	split.setGenes(fileReader.getExons(chrom));
	    }
	    
	    if(Main.samples > 0) {
	    	for(int i = 0; i<drawCanvas.sampleList.size(); i++) {
	    		Reads newReads = new Reads();
	    		newReads.sample = drawCanvas.sampleList.get(i);
	    		drawCanvas.sampleList.get(i).getreadHash().put(split, newReads);
	    	}
	    }
	    try {
	    	if(chromIndex.size() != 0) {
		    	if(chrom == null) {
		    		split.chromEnd =100;
		    	}
		    	else {
		    		split.chromEnd = chromIndex.get(Main.refchrom +chrom)[1].intValue();
		    	}
	    	}
	    	else {
	    		split.chromEnd = 100;
	    	}
	    }
	    catch(Exception e) {
	    	System.out.println(chrom);
	    	e.printStackTrace();
	    }
	    split.start = 0;
	    split.end =  split.chromEnd;
	    split.viewLength =split.end-split.start;
	    drawCanvas.splits.add(split);	
	    Main.drawCanvas.resizeCanvas(Main.drawCanvas.getWidth(), Main.drawCanvas.getHeight());
	  	split.getExonImageBuffer().setFont(Draw.defaultFont);
	  	split.getReadBuffer().setFont(Draw.defaultFont);
	  	split.getSelectbuf().setFont(Draw.defaultFont);
	  
	
		for(int i= 0; i<drawCanvas.splits.size(); i++) {
			drawCanvas.splits.get(i).setCytoImage(null);
			chromDraw.drawCyto(drawCanvas.splits.get(i));
			chromDraw.updateExons = true;
			chromDraw.repaint();
		}
		
		
}



static class MyFilterVCF extends javax.swing.filechooser.FileFilter {
	public boolean accept(File file) { 
			if (file.isDirectory()) {
				return true;
		    }	
			if (file.getName().endsWith(".vcf.gz")) {
				return true;
			}					
			else if(file.getName().endsWith(".vcf")) {
				return true;
			}
			else {
				return false;
			}	
	} 
	
	public String getDescription() { return "*.vcf, *.vcf.gz"; }
}
	static class MyFilterSES extends javax.swing.filechooser.FileFilter {
		public boolean accept(File file) { 
				if (file.isDirectory()) {
					return true;
			    }	
				if (file.getName().contains(".ses")) {
				return true;
			}					
			
			else {
				return false;
			}	
	} 
	
	public String getDescription() { return "*.ses"; }
}
	static class MyFilterBAM extends javax.swing.filechooser.FileFilter {
		public boolean accept(File file) { 
				if (file.isDirectory()) {
					return true;
			    }	
				if (file.getName().endsWith(".bam")) {
				return true;
			}					
			if (file.getName().endsWith(".cram")) {
				return true;
			}
			else {
				return false;
			}	
	} 
	
	public String getDescription() { return "*.bam, *.cram"; }
}
	static class MyFilterBED extends javax.swing.filechooser.FileFilter {
		public boolean accept(File file) { 
				if (file.isDirectory()) {
					return true;
			    }	
				if (file.getName().toLowerCase().endsWith(".bed.gz") || file.getName().toLowerCase().endsWith(".bed")) {
				return true;
			}					
			if (file.getName().toLowerCase().endsWith(".bedgraph.gz")) {
				return true;
			}	
			if(file.getName().toLowerCase().endsWith(".gff.gz")) {
				return true;
			}
			if(file.getName().toLowerCase().endsWith(".bigwig")) {
				return true;
			}
			if(file.getName().toLowerCase().endsWith(".bw")) {
				return true;
			}
			if(file.getName().toLowerCase().endsWith(".bigbed") ||file.getName().toLowerCase().endsWith(".bb")) {
				return true;
			}
			if(file.getName().toLowerCase().endsWith(".tsv.gz") ) {
				return true;
			}
			else {
				return false;
			}	
	} 
	
	public String getDescription() { return "*.bed, *.gff.gz, *.bedgraph.gz, *.bigWig, *.bigBed, *.tsv.gz"; }
}
	static class MyFilterTXT extends javax.swing.filechooser.FileFilter {
		public boolean accept(File file) { 
				if (file.isDirectory()) {
					return true;
			    }	
				if (file.getName().endsWith(".txt")) {
				return true;
			}				
			
			else {
				return false;
			}	
	} 
	
	public String getDescription() { return "Gene list *.txt"; }
}
	
	static class MyFilterCRAM extends javax.swing.filechooser.FileFilter {
		public boolean accept(File file) { 
				if (file.isDirectory()) {
					return true;
			    }	
				if (file.getName().endsWith(".cram")) {
				return true;
			}					
			
			else {
				return false;
			}	
	} 
	
	public String getDescription() { return "*.cram"; }
}
	static class MyFilterALL extends javax.swing.filechooser.FileFilter {
		public boolean accept(File file) { 
				
				return true;
		} 
		
		public String getDescription() { return "All files"; }
}
public static void zoomout() {		
	drawCanvas.setStartEnd(0.0,(double)drawCanvas.splits.get(0).chromEnd);
	
	if(samples > 0) {			
		if(drawCanvas.splits.get(0).chromEnd > Settings.settings.get("readDrawDistance")) {			
			Main.drawCanvas.clearReads();			
		}
		drawCanvas.removeSplits();
		chromDraw.varnode = null;
		chromDraw.vardraw = null;
		
		VariantHandler.table.hoverNode = null;
		VariantHandler.table.selectedNode = null;	
		Draw.updatevars = true;	
		drawCanvas.splits.get(0).getReadBuffer().setComposite( drawCanvas.composite);					
		drawCanvas.splits.get(0).getReadBuffer().fillRect(0,0, drawCanvas.splits.get(0).getReadImage().getWidth(),Main.drawScroll.getViewport().getHeight());	
		drawCanvas.splits.get(0).getReadBuffer().setComposite(drawCanvas.splits.get(0).getBackupr());		
		drawCanvas.rbuf.setComposite( drawCanvas.composite);				
		drawCanvas.rbuf.fillRect(0,0, drawCanvas.getWidth(),Main.drawScroll.getViewport().getHeight());	
		drawCanvas.rbuf.setComposite(drawCanvas.backupr);			
	}	
	
	if(Main.bedCanvas.bedTrack.size() > 0) {
		for(int i = 0 ; i<Main.bedCanvas.bedTrack.size(); i++) {			
			if(Main.bedCanvas.bedTrack.get(i).graph) {		
				Main.bedCanvas.bedTrack.get(i).setCurrent(Main.bedCanvas.bedTrack.get(i).getHead());
				Main.bedCanvas.calcScale(Main.bedCanvas.bedTrack.get(i));				
			}
			Main.bedCanvas.getMoreBeds(Main.bedCanvas.bedTrack.get(i));			
		}		
	}
	bedCanvas.repaint();
	Main.chromDraw.updateExons = true;
	drawCanvas.repaint();
	Main.chromDraw.repaint();

}

boolean checkGenome() {
	if(chromosomeDropdown.getItemAt(0) == null) {
		Main.showError("Add reference genome first.", "Note");
		if(AddGenome.frame == null) {
			  AddGenome.createAndShowGUI();
		 }
		AddGenome.frame.setVisible(true);
	    AddGenome.frame.setLocation(frame.getLocationOnScreen().x+frame.getWidth()/2 - AddGenome.frame.getWidth()/2, frame.getLocationOnScreen().y+frame.getHeight()/6);
	    AddGenome.frame.setState(JFrame.NORMAL);
		return false;
	}
	return true;
}
public void actionPerformed(ActionEvent e) {
	//Logo.frame.setVisible(false);
	if(e.getSource() == pleiadesButton) {	
		gotoURL("http://kaptah.local.lab.helsinki.fi/pleiades/");
	
	}
	else if(e.getSource() == manage) {			
		
		if(VariantHandler.frame == null) {
			VariantHandler.main(argsit);
		}
		VariantHandler.frame.setLocation(frame.getLocationOnScreen().x+frame.getWidth()/2 - VariantHandler.frame.getWidth()/2, frame.getLocationOnScreen().y+frame.getHeight()/6);			
		
		VariantHandler.frame.setState(JFrame.NORMAL);
		VariantHandler.frame.setVisible(true);
		Draw.calculateVars = true;
		Draw.updatevars = true;
		drawCanvas.repaint();
		
	}
	else if(e.getSource() == variantCaller) {
		
		//FileRead.checkSamples();
		
		if(VariantCaller.frame == null) { 
			
			VariantCaller.main(argsit);
			
		}
		VariantCaller.frame.setLocation(frame.getLocationOnScreen().x+frame.getWidth()/2 - VariantCaller.frame.getWidth()/2, frame.getLocationOnScreen().y+frame.getHeight()/6);			
		VariantCaller.frame.setState(JFrame.NORMAL);
		VariantCaller.frame.setVisible(true);
	}
	else if(e.getSource() == average) {
		if(Average.frame == null) {
		 Average.createAndShowGUI();
		}
		Average.setSamples();
		
		Average.frame.setLocation(frame.getLocationOnScreen().x+frame.getWidth()/2 - Average.frame.getWidth()/2, frame.getLocationOnScreen().y+frame.getHeight()/6);			
		
		Average.frame.setState(JFrame.NORMAL);
		Average.frame.setVisible(true);
	}
	else if(e.getSource() == errorlog) {

		ErrorLog.frame.setLocation(frame.getLocationOnScreen().x+frame.getWidth()/2 - ErrorLog.frame.getWidth()/2, frame.getLocationOnScreen().y+frame.getHeight()/6);			
		
		//		VariantHandler.frame.setAlwaysOnTop(true);	
		ErrorLog.frame.setState(JFrame.NORMAL);
		ErrorLog.frame.setVisible(true);
		
	}
/*	else if(e.getSource() == help) {
		JOptionPane.showMessageDialog(Main.chromDraw, "This is pre-release version of BasePlayer\nContact: help@baseplayer.fi\nUniversity of Helsinki", "Help", JOptionPane.INFORMATION_MESSAGE);
	}*/
	else if(e.getSource() == settings) {
	
		Settings.frame.setLocation(frame.getLocationOnScreen().x+frame.getWidth()/2 - Settings.frame.getWidth()/2, frame.getLocationOnScreen().y+frame.getHeight()/6);			
		
		Settings.frame.setState(JFrame.NORMAL);
		Settings.frame.setVisible(true);
	}
	else if(e.getSource() == update) {
		try {
			Updater update = new Updater();
			update.execute();
		
		}
		catch(Exception ex) {
			ex.printStackTrace();
		}
		
	}
	else if(e.getSource() == clearMemory) {
		FileRead.nullifyVarNodes();
		
		
		//FileRead.removeNonListVariants();f
		System.gc();
		chromDraw.repaint();
	}
	else if (e.getSource() == zoomout) {		
		
		zoomout();			
	}
	else if (e.getSource() == dosomething) {
		
		BedNode currentbed = bedCanvas.bedTrack.get(0).getHead().getNext();
		VarNode currentvar = FileRead.head.getNext();
		
		while(currentbed != null) {
			
			while(currentvar != null && currentvar.getPosition() < currentbed.getPosition()) {					
				currentvar = currentvar.getNext();
			}
			while(currentbed != null && currentvar.getPosition() > currentbed.getPosition()+currentbed.getLength()) {					
				currentbed = currentbed.getNext();				
			}
			
			if(currentvar != null && currentvar.getPosition() >= currentbed.getPosition() && currentvar.getPosition() < currentbed.getPosition()+currentbed.getLength()) {
				currentvar.setBedhit(true);					
				currentvar = currentvar.getNext();
			}				
			if(currentvar == null) {
				break;
			}
			currentbed = currentbed.getNext();				
		}			
		
	}
	else if (e.getSource() == clear) {
		clearData();		
	}
	else if(e.getSource() == exit) {
		
		System.exit(0);
	}	
	else if (e.getSource() == opensamples) {
		 try {			
			 if(!checkGenome()) return;
			 if(VariantHandler.frame != null) {
			  VariantHandler.frame.setState(Frame.ICONIFIED);
			 }
			 
			  JFileChooser chooser = new JFileChooser(path);	 
			  getText(chooser.getComponents());
	    	  chooser.setMultiSelectionEnabled(true);
	    	  chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
	    	  chooser.setAcceptAllFileFilterUsed(false);
	    	  MyFilterBAM bamFilter = new MyFilterBAM();
	    	  MyFilterVCF vcfFilter = new MyFilterVCF();
	    	  chooser.addChoosableFileFilter(vcfFilter); 	
	    	  chooser.addChoosableFileFilter(bamFilter); 
	    	  chooser.setDialogTitle("Add samples");
	    	  chooser.setPreferredSize(new Dimension((int)screenSize.getWidth()/3, (int)screenSize.getHeight()/3));
	          int returnVal = chooser.showOpenDialog((Component)this.getParent());	         
	          
	         if (returnVal == JFileChooser.APPROVE_OPTION) {
	    
	        	  File vcffiles[] = chooser.getSelectedFiles(); 
	        	  if(!vcffiles[0].exists()) {
	        		 
	        		  if(Main.chooserText.contains("`")) {
	        			  Main.chooserText.replace("`", "?");		        			  
	        		  }
	        		  if(Main.chooserText.contains(" ")) {
	        			  Main.chooserText.replace(" ", "%20");
	        		  }
	        		  if(Main.chooserText.contains("pleiades")) {
	        			 try {
	        			  URL url = new URL(Main.chooserText);
	        			  System.out.println(Main.chooserText);
	        			  HttpURLConnection httpConn = (HttpURLConnection)url.openConnection();
	        			  httpConn.connect();
	        			  
	        	//		  String fileURL = url.getPath();
	        			 
	        				int responseCode = httpConn.getResponseCode();
	        		//		int BUFFER_SIZE = 4096;
	        				// always check HTTP response code first
	        				
	        				if (responseCode == HttpsURLConnection.HTTP_OK) {
	        				
	        		//			String fileName = "";
	        					
	        				//	String disposition = httpConn.getHeaderField("Content-Disposition");
	        			//		String contentType = httpConn.getContentType();
	        			//		int contentLength = httpConn.getContentLength();
	        					/*
	        					if (disposition != null) {
	        						// extracts file name from header field
	        						int index = disposition.indexOf("filename=");
	        						if (index > 0) {
	        							fileName = disposition.substring(index + 10,
	        									disposition.length() - 1);
	        						}
	        					} else {
	        						// extracts file name from URL
	        						fileName = fileURL.substring(fileURL.lastIndexOf("/") + 1,
	        						fileURL.length());
	        					}
	        					*/
	        			/*		long downloaded = 0;
				      			int bytesRead = -1, counter=0;
				      		*/	
				      			
				      			String loading = drawCanvas.loadingtext;
				      			InputStream inputStream = httpConn.getInputStream();
				      			BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
				      			Main.drawCanvas.loadingtext = loading + " 0MB";
				      			String line;
				      			StringBuffer buffer = new StringBuffer("");
				      			while((line=reader.readLine()) != null) {		
				      				
				      				buffer.append(line);
				      				
				      			}
				      			inputStream.close();
				      			reader.close();
				      			String split2;
				      			String[] split = buffer.toString().split("request");
				      			String location;
				      			ArrayList<File> array = new ArrayList<File>();
				      			File[] paths;
				    			FileSystemView fsv = FileSystemView.getFileSystemView();
				    			paths = File.listRoots();
				    			String loc = "/mnt";
				    			
				    			for(File path:paths){
				    				if(fsv.getSystemDisplayName(path).contains("merit")) {
				    					loc = path.getCanonicalPath();
				    				}					    			  
				    			}
				    			
				      			for(int i = 0; i<split.length; i++) {
				      				
				      				if(!split[i].contains("lastLocation")) {		
				      					
				      					continue;
				      				}
				      				
				      				split2 = split[i].split("\"lastLocation\":\"")[1];
				      				location = split2.substring(0,split2.indexOf("\"}"));
				      				String filename = location.substring(location.lastIndexOf("/"))+".snps_indels.vcf.gz";
				      				
				      				location = location.replace("/mnt", loc) +"/wgspipeline/align/" +filename;
				      				
				      				if(!new File(location).exists()) {
				      					
				      					location = split2.substring(0,split2.indexOf("\"}"));
					      				filename = location.substring(location.lastIndexOf("/"))+".snps_indels.hc.vcf.gz";
				      					location = location.replace("/mnt", loc) +"/wgspipeline/" +filename;
				      					if(!new File(location).exists()) {
					      					System.out.println(location);
					      				}
				      					else {
				      						array.add(new File(location));
				      					}
				      				}
				      				else {
				      					
				      					array.add(new File(location));
				      				}				      				
				      			}				      			
				      			
				      			File[] files = new File[array.size()];
				      			for(int i = 0; i<files.length; i++) {
				      				
				      				files[i] = array.get(i);
				      			}				      			
				      			 FileRead filereader = new FileRead(files);
				      			 filereader.start = (int)drawCanvas.selectedSplit.start;
				        		 filereader.end = (int)drawCanvas.selectedSplit.end;
				        		 filereader.readVCF = true;
				        		 filereader.execute();        			
	        				}
		        		  }
		        		  catch(Exception ex) {
		        			  ex.printStackTrace();
		        		  }
	        		  	}
	        		  	
	        		  return;
	        	  }
	        	  path = vcffiles[0].getParent();
	        	  writeToConfig("DefaultDir=" +path);
	        	  FileRead filereader = new FileRead(vcffiles);
	        	 
	        	  if(chooser.getFileFilter().equals(chooser.getChoosableFileFilters()[0])) {
	        		  filereader.start = (int)drawCanvas.selectedSplit.start;
	        		  filereader.end = (int)drawCanvas.selectedSplit.end;
	        		  filereader.readVCF = true;
	        		  filereader.execute();
	        	  }
	        	  else if(chooser.getFileFilter().equals(chooser.getChoosableFileFilters()[1])){
	        		  filereader.readBAM = true;
	        		  filereader.execute();
	        	  }
	          }
		 }
		 catch(Exception ex) {
			 ex.printStackTrace();
		 }
		
	}
	else if (e.getSource() == addcontrols) {
		if(!checkGenome()) return;
		if(VariantHandler.frame != null) {
		 VariantHandler.frame.setState(Frame.ICONIFIED);
		}
		 JFileChooser chooser = new JFileChooser(controlDir);	 
 //   	  JFileChooser chooser = new JFileChooser(path);	    	  
    	  chooser.setMultiSelectionEnabled(true);	    	  
    	  chooser.setAcceptAllFileFilterUsed(false);	    	 
    	  Control.MyFilter vcfFilter = new Control.MyFilter();
    	  chooser.addChoosableFileFilter(vcfFilter);		    	
    
    	  chooser.setDialogTitle("Add controls");
    	  chooser.setPreferredSize(new Dimension((int)screenSize.getWidth()/3, (int)screenSize.getHeight()/3));
          int returnVal = chooser.showOpenDialog((Component)this.getParent());	       	         
    	  
          if (returnVal == JFileChooser.APPROVE_OPTION) {
        	  
        	  File vcffiles[] = chooser.getSelectedFiles(); 
        	  controlDir = vcffiles[0].getParent();
        	  writeToConfig("DefaultControlDir=" +controlDir);
        	  
        	  if(chooser.getFileFilter().equals(chooser.getChoosableFileFilters()[0])) {
        		 Control.addFiles(vcffiles);
        	  }        	
          }
	}		
	else if (e.getSource() == addtracks) {
		if(!checkGenome()) return;
		if(VariantHandler.frame != null) {
		 VariantHandler.frame.setState(Frame.ICONIFIED);
		}
		 JFileChooser chooser = new JFileChooser(Main.trackDir);	 
		 getText(chooser.getComponents());
    	  chooser.setMultiSelectionEnabled(true);	    	  
    	  chooser.setAcceptAllFileFilterUsed(false);	    	 
    	  MyFilterBED bedFilter = new MyFilterBED();	
    	  MyFilterTXT txtFilter = new MyFilterTXT();
    	  chooser.addChoosableFileFilter(bedFilter); 
    	  chooser.addChoosableFileFilter(txtFilter); 
    	  chooser.setDialogTitle("Add tracks");
    	  chooser.setPreferredSize(new Dimension((int)screenSize.getWidth()/3, (int)screenSize.getHeight()/3));
          int returnVal = chooser.showOpenDialog((Component)this.getParent());	      
         
          try {
        
       /*  if(chooser.getSelectedFile() != null && chooser.getSelectedFile().getCanonicalPath().contains("http")) {
        	 try {
        		 URL url = new URL("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam");
        	//	 FileRead.samFileReader = SamReaderFactory.makeDefault().open(url);
        		 htsjdk.samtools.SamInputResource resource = htsjdk.samtools.SamInputResource.of(url);
        		 URL indexurl = new URL("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.bai");
        		 resource.index(indexurl);
        		 htsjdk.samtools.SamReader reader = SamReaderFactory.make().open(resource);
        		 
        		 Iterator<SAMRecord> bamIterator = reader.queryOverlapping("11", 5000000, 6000000);	
        		
        		 reader.close();
        	//	 SamReaderFactory.make().
        	//	 htsjdk.samtools.SamInputResource.of(url)
	 //       	 URL url = new URL("http://dx2-tutoh-1.ltdk.helsinki.fi/SELEX_Ensembl_2015/hits/SELEX_Ensembl_dom_set_2015_humanGRCh37_9_or_max1M.gff.gz");
	       
	        	
	   //    	 SeekableStream stream = SeekableStreamFactory.getInstance().getStreamFor(url);
	    //   	 System.out.println(url.getHost() +url.getPath());
	    //    TabixReader tabixReader = new TabixReader(url.getHost() +url.getPath(), "http://dx2-tutoh-1.ltdk.helsinki.fi/SELEX_Ensembl_2015/hits/SELEX_Ensembl_dom_set_2015_humanGRCh37_9_or_max1M.gff.gz.tbi", stream);   	
	 			
	      
	        	 
	 		//	TabixReader.Iterator bedIterator = null;
	 			  try {
	 				 
	 			//	 bedIterator = tabixReader.query(Main.chromosomeDropdown.getSelectedItem().toString() +":" +5000000 +"-" +6000000);
	 			  }
	 			  catch(Exception ex) {
	 				  ex.printStackTrace();
	 			  }
	
        	 }
        	 catch(Exception ex) {
        		 ex.printStackTrace();
        	 }
        	 // File vcffiles[] = chooser.getSelectedFiles(); 
        	 
        	//  FileRead filereader = new FileRead(vcffiles);	        	 
        	//filereader.readBED(vcffiles);
        	     
          }
    	 
         else */
        	  if (returnVal == JFileChooser.APPROVE_OPTION) {
        	
        	  File bedfiles[] = chooser.getSelectedFiles(); 
        	  if(bedfiles[0].exists()) {
        		
        		  trackDir = bedfiles[0].getParent();	 	
        		  writeToConfig("DefaultTrackDir=" +trackDir);
	        	  FileRead filereader = new FileRead(bedfiles);		        	
	        	  filereader.readBED(bedfiles);
	        	  
        	  }
        	  else { 		  
        		 
        		  
        		  if(Main.chooserText.length() > 5 && Main.chooserText.endsWith(".bed.gz")  ||  Main.chooserText.endsWith(".gff.gz") || Main.chooserText.endsWith(".bedgraph.gz")) {	        			  
        			 
	        		  if(Main.chooserText.startsWith("http://") || Main.chooserText.startsWith("ftp://")) {
	        			 
	        			  URL url = new URL(Main.chooserText);		      		       	
	      		      	  SeekableStream stream = SeekableStreamFactory.getInstance().getStreamFor(url);
	      		      	  TabixReader tabixReader  = null;
	      		      	  String index = null;
	      		      	  
	      		      	  try {
	      		      		  tabixReader = new TabixReader(Main.chooserText, Main.chooserText +".tbi", stream);
	      		      		  index = Main.chooserText +".tbi";
	      		      	  }
	      		      	  catch(Exception ex) {
	      		      		  try {
	      		      			  tabixReader = new TabixReader(Main.chooserText, Main.chooserText.substring(0,Main.chooserText.indexOf(".gz")) +".tbi", stream);   
	      		      			index = Main.chooserText.substring(0,Main.chooserText.indexOf(".gz")) +".tbi";
	      		      		  }
	      		      		  catch(Exception exc) {
	      		      			  exc.printStackTrace();
	      		      		  }
	      		      	  }
	        			  if(tabixReader != null && index != null) {
	        				  
	        				  FileRead filereader = new FileRead(bedfiles);
	        				  filereader.readBED(Main.chooserText, index);
	        				  
	        				  tabixReader.close();
	        			  }
	        		  }
        		  }
        		  else {
        			  if(Main.chooserText.contains("://")) {
		        			
	        	//		  URL url = new URL(Main.chooserText);	
	        			
	      		      //	  SeekableStream stream = SeekableStreamFactory.getInstance().getStreamFor(url);
	      		      	  
	      		      //	  BBFileReader bbreader = null;
	      		      	
	      		      	  try {
	      		      //		  bbreader = new BBFileReader(Main.chooserText, stream);
	      		      		
	      		      		
	      		      	  }
	      		      	  catch(Exception ex) {
	      		      		 ex.printStackTrace();
	      		      	  }
	      		      	 
	      		      
	        			//  if(bbreader != null) {
	        				  
	        				  FileRead filereader = new FileRead(bedfiles);
	        				  filereader.readBED(Main.chooserText, "nan");
	        				  
	        				  
	        				
	        		//	  }
	      		      //	  stream.close();
	        		  }
        		  }
        	  }
          }
        
          }
        	 catch(Exception ex) {
        		 ex.printStackTrace();
        	 }
	}
	else if (e.getSource() == openProject) {
		if(!checkGenome()) return;
		if(VariantHandler.frame != null) {
		 VariantHandler.frame.setState(Frame.ICONIFIED);
		}
		 openProject(); 
	}
	else if (e.getSource() == saveProjectAs) {
		if(VariantHandler.frame != null) {
		 VariantHandler.frame.setState(Frame.ICONIFIED);
		}
		try {	  	    
	    		   
    	  JFileChooser chooser = new JFileChooser();
    	  chooser.setAcceptAllFileFilterUsed(false);
    	  MyFilterSES sesFilter = new MyFilterSES();	    	  
    	  chooser.addChoosableFileFilter(sesFilter);
    	  chooser.setDialogTitle("Save project as...");
          int returnVal = chooser.showSaveDialog((Component)this.getParent());	           	  
	       
	      if(returnVal == JFileChooser.APPROVE_OPTION) {   		    	  
		       File outfile = chooser.getSelectedFile();	
		       if(!outfile.getAbsolutePath().endsWith(".ses")) {
		    	   outfile = new File(outfile.getAbsolutePath() +".ses");
		       }
		       
		       Serializer ser = new Serializer();
			   ser.serialize(outfile);
	      }
		}
		catch(Exception ex) {
			ex.printStackTrace();
		}
		
	}
	else if (e.getSource() == saveProject) {
		if(drawCanvas.drawVariables.projectName.equals("Untitled")) {
				saveProjectAs.doClick();
			}
			else {
				
				Serializer ser = new Serializer();
				ser.serialize(drawCanvas.drawVariables.projectFile);
			}
		}
/*		else if(e.getSource() == welcome) {
			WelcomeScreen.main(args);
			WelcomeScreen.frame.setVisible(true);
			WelcomeScreen.frame.setLocation(frame.getLocationOnScreen().x+frame.getWidth()/2 - WelcomeScreen.frame.getWidth()/2, frame.getLocationOnScreen().y+frame.getHeight()/6);
			
		}*/
}
	
	public void getText(Component[] comp)  {
		
	    for(int x = 0; x < comp.length; x++)   {
	      if(comp[x] instanceof JPanel) {
	    	  getText(((JPanel)comp[x]).getComponents());
	    	 
	      }
	      else if(comp[x] instanceof JTextField)
	      {
	     
	        chooserTextField = ((JTextField)comp[x]);
	        chooserTextField.addKeyListener(new KeyListener() {

			@Override
			public void keyTyped(KeyEvent e) {
			
				Main.chooserText = chooserTextField.getText().trim();
				if(Main.chooserText.contains("?")) {
				chooserTextField.setText(Main.chooserText.replace("?", "`"));						
				}
			}

			@Override
			public void keyPressed(KeyEvent e) {					
				
			}
			@Override
			public void keyReleased(KeyEvent e) {			
				
			} 					        	
        });
       
      }	     
    }  
}

static void clearData() {
	try {
	FileRead.checkSamples();
	Main.drawCanvas.bam = false;
	undoList.clear();	
	undoPointer = 0;
	bedCanvas.bedOn = false;
	back.setEnabled(false);
	forward.setEnabled(false);
	if(drawCanvas.clusterNodes != null) {
		drawCanvas.clusterNodes.clear();
	}
	Main.variantCaller.setEnabled(false);
	Main.average.setEnabled(false);
	Average.outVector.clear();
	zoomout();
	FileRead.nullifyVarNodes();
	samples = 0;			
	varsamples = 0;
	readsamples = 0;
	FileRead.head.putNext(null);
	drawCanvas.current = FileRead.head;
	Main.drawCanvas.sampleList.clear();
	Main.drawCanvas.selectedSample = null;
	Main.drawCanvas.selectedSampleIndex = -1;
	
	Draw.updatevars = true;
	chromDraw.varnode = null;
	chromDraw.vardraw = null;
	drawCanvas.currentDraw = null;
	
	if(drawCanvas.mismatches != null) {
		drawCanvas.mismatches.clear();
	}
	
	if(drawCanvas.splits.size() > 1) {
		drawCanvas.removeSplits();
	}
	
	for(int i = 0 ; i<bedCanvas.bedTrack.size(); i++) {
		bedCanvas.bedTrack.get(i).getHead().putNext(null);
		bedCanvas.bedTrack.get(i).setCurrent(null);
		bedCanvas.bedTrack.get(i).setDrawNode(null);
		MethodLibrary.removeHeaderColumns(bedCanvas.bedTrack.get(i));
		FileRead.removeTable(bedCanvas.bedTrack.get(i));
		
		if(bedCanvas.bedTrack.get(i).getTable() != null && bedCanvas.bedTrack.get(i).getTable().bedarray != null) {
			bedCanvas.bedTrack.get(i).getTable().bedarray.clear();
			bedCanvas.bedTrack.get(i).getTable().hoverNode = null;
			bedCanvas.bedTrack.get(i).getTable().selectedNode = null;
		}
	}
	
	SplitClass split = Main.drawCanvas.splits.get(0);
	
	if(split.getGenes() != null) {
		for(int g = 0 ; g<split.getGenes().size(); g++) {
			split.getGenes().get(g).mutations = 0;
			split.getGenes().get(g).missense = 0;
			split.getGenes().get(g).nonsense = 0;
			split.getGenes().get(g).synonymous = 0;
			split.getGenes().get(g).intronic = 0;
			split.getGenes().get(g).utr = 0;
			split.getGenes().get(g).samples.clear();
			split.getGenes().get(g).varnodes.clear();
			split.getGenes().get(g).transcriptString = new StringBuffer();
			/*for(int i = 0; i<split.getGenes().get(g).getTranscripts().size(); i++) {
			split.getGenes().get(g).getTranscripts().get(i).mutations = 0;
			split.getGenes().get(g).getTranscripts().get(i).missense = 0;
			split.getGenes().get(g).getTranscripts().get(i).nonsense = 0;
			split.getGenes().get(g).getTranscripts().get(i).synonymous = 0;
			split.getGenes().get(g).getTranscripts().get(i).intronic = 0;
			split.getGenes().get(g).getTranscripts().get(i).utr = 0;
			split.getGenes().get(g).getTranscripts().get(i).samples.clear();
			split.getGenes().get(g).getTranscripts().get(i).varnodes.clear();
		}*/
	}
}
split = null;
Main.drawCanvas.splits.get(0).resetSplits();
bedCanvas.bedTrack.clear();
Main.bedCanvas.trackDivider.clear();
Main.controlDraw.trackDivider.clear();
if(Main.drawCanvas.advQualities != null) {
	Main.drawCanvas.advQualities.clear();
	if(Main.drawCanvas.drawVariables.advQDraw != null) {
		Main.drawCanvas.drawVariables.advQDraw.clear();
	}
	if(Main.drawCanvas.drawVariables.advQDrawIndel != null) {
		Main.drawCanvas.drawVariables.advQDrawIndel.clear();
	}
	
}
VariantHandler.removeMenuComponents();
bedCanvas.track = null;
bedCanvas.infoNode = null;
bedCanvas.preInfoNode = null;
bedScroll.setVisible(false);
controlScroll.setVisible(false);
varpane.setDividerSize(0);
trackPane.setVisible(false);
varpane.setResizeWeight(0.0);   
trackPane.setDividerSize(0);
if(VariantHandler.table.vararray !=null) {
	VariantHandler.table.vararray.clear();		
}

if(VariantHandler.tables.size() > 0) {
	for(int i = 0; i<VariantHandler.tables.size(); i++) {
		if(VariantHandler.tables.get(i).vararray !=null) {
			VariantHandler.tables.get(i).vararray.clear();				
		}
		VariantHandler.tables.get(i).clear();
	}
}
VariantHandler.stattable.clear();
VariantHandler.commonSlider.setMaximum(1);
VariantHandler.commonSlider.setValue(1);
VariantHandler.commonSlider.setUpperValue(1);	
VariantHandler.geneSlider.setMaximum(1);
VariantHandler.geneSlider.setValue(1);

VariantHandler.table.clear();
VariantHandler.table.repaint();
Draw.setScrollbar(0);
drawCanvas.drawVariables.visiblestart = 0;
drawCanvas.drawVariables.visiblesamples = 1;
drawCanvas.sidebar = false;
drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), drawScroll.getViewport().getHeight());
drawCanvas.repaint();
//	Control.controlData.data.clear();¨
if(VariantHandler.table.geneheader.size() > VariantHandler.table.geneheaderlength) {
	while(VariantHandler.table.geneheader.size() > VariantHandler.table.geneheaderlength) {
		VariantHandler.table.geneheader.remove(VariantHandler.table.geneheader.size()-1);
	}
	VariantHandler.table.repaint();
}
for(int i = 0 ; i<Control.controlData.fileArray.size(); i++) {
	MethodLibrary.removeHeaderColumns(Control.controlData.fileArray.get(i));
}
Control.controlData.fileArray.clear();
//	Control.controlData.sampleArray.clear();
Control.controlData.total = 0;
Control.controlData.sampleCount = 0;
//	Control.useCheck.setSelected(false);
drawCanvas.drawVariables.projectName = "Untitled";
drawCanvas.drawVariables.projectFile = null;
frame.setTitle("BasePlayer - Untitled Project");
//	System.gc();
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}
	public static class runner extends SwingWorker<String, Object> {		
		
		protected String doInBackground() {			
			return "";
	}
	
}
private static void createAndShowGUI() {
	try {
		
	/*	if(System.getProperty("os.name").toLowerCase().contains("nux")) {
			
		}
		else {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		}
	*/
	//	SystemUtils.IS_OS_WINDOWS
	//	UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	
	JFrame.setDefaultLookAndFeelDecorated(false);
	frame = new JFrame("BasePlayer");
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
   
    frame.setResizable(true);
    
    JComponent newContentPane = new Main();
    newContentPane.setOpaque(true); 
    
    frame.setContentPane(newContentPane);
    frame.pack();
    frame.setVisible(true);
   
    
  
     
}
public static void main(String[] args) {
	try {
		UIManager.put("Slider.paintValue", false);
	//	UIManager.setLookAndFeel("com.sun.java.swing.plaf.windows.WindowsLookAndFeel"); 
		UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		System.setProperty("sun.java2d.d3d", "false"); 	
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	Main.args = args;
	 javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {	            	
            		//Logo.main(argsit);
            		createAndShowGUI();	            	
            }
        });			
}

static void updatePositions(double start, double end) {
		Main.widthLabel.setBackground(new Color(250,250,250));
		positionField.setText("chr" +drawCanvas.selectedSplit.chrom +":" +MethodLibrary.formatNumber((int)start+1) +"-" +MethodLibrary.formatNumber((int)end));
		Main.widthLabel.setText("" +MethodLibrary.formatNumber((int)(end-start)) +"bp");
	
}
@Override
public void stateChanged(ChangeEvent event) {
	
	
}

public void componentResized(ComponentEvent e) {
	if(e.getSource() == drawScroll) {
		if(Main.sidebarWidth != Main.upPanel.getDividerLocation()+2) {				
			if(samples*drawCanvas.drawVariables.sampleHeight < drawScroll.getViewport().getHeight()) {						
				drawDimensions.setSize(drawScroll.getViewport().getWidth(),drawScroll.getViewport().getSize().height );	
				drawCanvas.setPreferredSize(drawDimensions);					
				drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), drawScroll.getViewport().getHeight());	
				
				if(drawCanvas.splits.size() > 0) {
					for(int i = 0 ; i<drawCanvas.splits.size(); i++) {
						drawCanvas.splits.get(i).updateReads = true;						
					}
				}					
			}		
			Draw.updatevars = true;
			drawCanvas.repaint();				
		}
		Main.sidebarWidth = upPanel.getDividerLocation()+4;
		chromDimensions.setSize(drawScroll.getViewport().getWidth()-upPanel.getDividerLocation(), splitPane.getDividerLocation());
		chromDraw.setPreferredSize(chromDimensions);
		chromDraw.updateExons = true;
		chromDraw.repaint();
		if(Main.samples == 0) {
			drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), drawScroll.getViewport().getHeight());	
		}
		else if(Main.drawCanvas.drawVariables.visiblesamples == Main.samples) {
			
			Main.drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), Main.drawCanvas.getHeight());
			
			if(Main.drawCanvas.drawVariables.sampleHeight < Main.drawScroll.getViewport().getHeight()) {
				
				Main.drawCanvas.drawVariables.visiblesamples = (short)((Main.drawScroll.getViewport().getHeight() /(double)Main.drawCanvas.drawVariables.sampleHeight)+0.5);
				
			}
		}
		
		return;
	}
	else if(e.getSource() == chrompan) {
		if(Main.sidebarWidth != Main.upPanel.getDividerLocation()+2) {
			if(samples*drawCanvas.drawVariables.sampleHeight < drawScroll.getViewport().getHeight()) {		
				
				drawDimensions.setSize(drawScroll.getViewport().getWidth(),drawScroll.getViewport().getSize().height );	
				drawCanvas.setPreferredSize(drawDimensions);
				
				drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), drawScroll.getViewport().getHeight());	
				if(drawCanvas.splits.size() > 0) {
					for(int i = 0 ; i<drawCanvas.splits.size(); i++) {
						drawCanvas.splits.get(i).updateReads = true;
					
					}
				}
				
			}		
			Main.bedCanvas.repaint();
			Main.controlDraw.repaint();
			Draw.updatevars = true;
			drawCanvas.repaint();		
			
		}
		
		Main.sidebarWidth = upPanel.getDividerLocation()+4;
	//	drawCanvas.drawWidth = drawScroll.getViewport().getWidth()-Main.sidebarWidth;
		chromDimensions.setSize(drawScroll.getViewport().getWidth()-upPanel.getDividerLocation(), splitPane.getDividerLocation());
		chromDraw.setPreferredSize(chromDimensions);
		chromDraw.updateExons = true;
		chromDraw.repaint();
		if(Main.drawCanvas.drawVariables.visiblesamples == Main.samples) {
			Main.drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), Main.drawCanvas.getHeight());
			if(Main.drawCanvas.drawVariables.sampleHeight < Main.drawScroll.getViewport().getHeight()) {
				Main.drawCanvas.drawVariables.visiblesamples = (short)((Main.drawScroll.getViewport().getHeight() /(double)Main.drawCanvas.drawVariables.sampleHeight)+0.5);
				//TODO
			}
		}
		else {
			Main.drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), Main.drawCanvas.getHeight());
		}
		
		return;
		
	}
	
	if(e.getComponent().getName() != null && e.getComponent().getName().contains("frame0")) {
	
		if(drawScroll.getViewport().getWidth() > 0) {		
			drawDimensions.setSize(drawScroll.getViewport().getWidth(),drawScroll.getViewport().getSize().height );	
			drawCanvas.setPreferredSize(drawDimensions);
			
			chromDimensions.setSize(drawScroll.getViewport().getWidth()-Main.sidebarWidth-1, splitPane.getDividerLocation());
			
			chromDraw.setPreferredSize(chromDimensions);
			
			drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), drawScroll.getViewport().getHeight());
			
			if(drawCanvas.splits.size() > 0) {
				for(int i = 0 ; i<drawCanvas.splits.size(); i++) {
					drawCanvas.splits.get(i).updateReads = true;				
				}
			}
			drawCanvas.repaint();
			return;
		}			 
    }       
}

	static void cancel() {
	  cancel = true;
	  
	  if(Draw.variantcalculator) {
		  FileRead.cancelvarcount = true;
		  Main.drawCanvas.ready("all");
	  }
	  else if(drawCanvas.loadingtext.contains("Loading variants")) {
		   Main.drawCanvas.ready("all");		
		   drawCanvas.current = null;
		   drawCanvas.currentDraw = null;
		   chromDraw.vardraw = null;
		   chromDraw.varnode = null;
		   drawCanvas.variantsStart = 0;
		   drawCanvas.variantsEnd = 1;
		   FileRead.head.putNext(null);
		   Draw.updatevars = true;
		   FileRead.cancelvarcount = true;
		   FileRead.cancelfileread = true;		  
	  }
	  else if(drawCanvas.loadingtext.contains("Processing variants")) {
		   Main.drawCanvas.ready("Processing variants...");			 
		   FileRead.cancelvarcount = true;
		   drawCanvas.current = null;
		   drawCanvas.currentDraw = null;
		   chromDraw.vardraw = null;
		   chromDraw.varnode = null;
		   drawCanvas.variantsStart = 0;
		   drawCanvas.variantsEnd = 1;
		   FileRead.head.putNext(null);
		   Draw.updatevars = true;
		   FileRead.cancelvarcount = true;
		   FileRead.cancelfileread = true;
		   Main.drawCanvas.ready("all");
	  }
	  else if(drawCanvas.loadingtext.contains("reads")) {
		  FileRead.cancelreadread = true;			   
		  Iterator<Map.Entry<SplitClass, Reads>> it;
		  Map.Entry<SplitClass, Reads> pair;
	      Reads reads;
	      for(int i = 0; i<drawCanvas.sampleList.size(); i++) {
			   if(drawCanvas.sampleList.get(i).getreadHash() != null) {				
					it = drawCanvas.sampleList.get(i).getreadHash().entrySet().iterator();
				    while (it.hasNext()) {
				    	pair = (Map.Entry<SplitClass, Reads>)it.next();
				        reads = pair.getValue();
				        reads.loading = false;
				        reads.setLastRead(null);
				    }						  
			   }
	      }
		  FileRead.cancelfileread = true;
		  Main.drawCanvas.clearReads();		  
		  Main.drawCanvas.ready("all");			  
		 
		  Draw.updateReads = true;
		  Draw.updateCoverages = true;
		  for(int i = 0; i<drawCanvas.splits.size(); i++) {
			  drawCanvas.splits.get(i).updateReads = true;				 
		  }			  
		  
	  }
	  else if(drawCanvas.loadingtext.contains("Processing variants")) {
		   FileRead.cancelfileread = true;
		   Main.drawCanvas.ready("all");
	  }
	  else if(drawCanvas.loadingtext.contains("BED")) { 
		   FileRead.cancelfileread = true;
		   Main.drawCanvas.ready("all");
	  }
	  else if(drawCanvas.loadingtext.contains("samples")) {
		  FileRead.cancelfileread = true;
		   Main.drawCanvas.ready("all");
	  }
	  else if(drawCanvas.loadingtext.contains("Updating")) {
		  FileRead.cancelfileread = true;
		   Main.drawCanvas.ready("Updating BasePlayer... (downloading BasePlayer.jar from http://baseplayer.fi/update/");
	  }
	  else if(drawCanvas.loadingtext.contains("Annotating")) {
		  /* FileRead.cancelfileread = true;
		   BedCanvas.annoTrack.intersect = false;
		   BedCanvas.annoTrack.used = false;
		   
		   Main.bedCanvas.removeBedhits(BedCanvas.annoTrack);		  
		   BedCanvas.annoTrack = null;
		   
		   Main.drawCanvas.ready(Main.drawCanvas.loadingtext);
		   Main.drawCanvas.ready("Annotating variants");
		   Main.drawCanvas.ready("Loading variants...");
		   */
		   
		   FileRead.cancelvarcount = true;
		   drawCanvas.current = null;
		   drawCanvas.currentDraw = null;
		   chromDraw.vardraw = null;
		   chromDraw.varnode = null;
		   drawCanvas.variantsStart = 0;
		   drawCanvas.variantsEnd = 1;
		   FileRead.head.putNext(null);
		   Draw.updatevars = true;
		   FileRead.cancelvarcount = true;
		   FileRead.cancelfileread = true;
		   Main.drawCanvas.ready("all");
	  }
	  else {
		  Main.drawCanvas.ready("all");
	  }
	  
	
//	   Main.drawCanvas.readyQueue.clear();
	  Main.bedCanvas.annotator = false;
	  FileRead.readFiles = false;
	   cancel = false;
	   drawCanvas.repaint();
	   chromDraw.repaint();
	   FileRead.cancelfileread = false;
}
static void cancel(int nro) {
	  cancel = true;
	  Main.drawCanvas.ready("all");
	   Draw.updateReads = false;
	   for(int i = 0; i<drawCanvas.splits.size(); i++) {
		   drawCanvas.splits.get(i).updateReads = false;
		  
	   }
	   
	   Iterator<Map.Entry<SplitClass, Reads>> it;
		 Map.Entry<SplitClass, Reads> pair;
		 Reads reads;
	   for(int i = 0; i<drawCanvas.sampleList.size(); i++) {
		   if(drawCanvas.sampleList.get(i).getreadHash() != null) {				
				it = drawCanvas.sampleList.get(i).getreadHash().entrySet().iterator();
			    while (it.hasNext()) {
			    	pair = (Map.Entry<SplitClass, Reads>)it.next();
			        reads = pair.getValue();
			        reads.loading = false;
			    }
				  
		   }
	   }
	   
		   for(int i = 0 ; i<Main.bedCanvas.bedTrack.size(); i++) {
			   Main.bedCanvas.removeBeds(Main.bedCanvas.bedTrack.get(i));
		   }
		   drawCanvas.current = null;
		   drawCanvas.currentDraw = null;
		   chromDraw.vardraw = null;
		   chromDraw.varnode = null;
		   Draw.updatevars = true;
		   FileRead.cancelvarcount = true;
		   FileRead.cancelfileread = true;
		   FileRead.cancelreadread = true;
		   Main.drawCanvas.variantsStart = 0;
			Main.drawCanvas.variantsEnd = 1;
		   cancel = false;
  
	}
	
	public static void writeToConfig(String attribute) {
		Main.configChanged = true;
		Boolean found = false;
		
		for(int i = 0 ; i<Launcher.config.size(); i++) {
			
			if(Launcher.config.get(i).contains(attribute.subSequence(0, attribute.indexOf("=")))) {
			found = true;
			
			Launcher.config.set(i, attribute);
			break;
		}
	}
	if(!found) {
		Launcher.config.add(attribute);
	}	
	
}
@Override
public void mouseClicked(MouseEvent event) {
	
   if(cancelhover && drawCanvas.loading) {
	 cancel();
   }
  if(event.getComponent().getName() != null && genomehash.containsKey(event.getComponent().getName())) {
	 
	  changeRef(event.getComponent().getName());
	  
   }
  
	/*if(event.getSource() == drawScroll.getVerticalScrollBar()) {
		
		drawCanvas.drawVariables.visiblestart = (short)(drawScroll.getVerticalScrollBar().getValue()/drawCanvas.sampleHeight);
		
		drawCanvas.drawVariables.visibleend = (short)((drawCanvas.drawVariables.visiblestart + drawScroll.getHeight()/drawCanvas.sampleHeight));
	
		drawCanvas.scrollbar = true;
		Draw.updatevars = true;
		
		drawCanvas.repaint();
		drawCanvas.scrollbar = false;
	}*/
	if(event.getSource() == searchField) {
		searchField.setFocusable(true);
		searchField.requestFocus();
		if(event.getClickCount() == 2) {
			
			searchField.setText("");
		}
	}
	else {
		searchField.setFocusable(false);
	}
	if(event.getSource() == back && back.isEnabled()) {
		undoPointer--;
		if(undoPointer < 1) {
			undoPointer = 0;
			back.setEnabled(false);
		}
		forward.setEnabled(true);
		returnlist = parseSearch(undoList.get(undoPointer));
		searchField.setText(undoList.get(undoPointer));
		FileRead.search = true;
		drawCanvas.gotoPos(returnlist[0], Integer.parseInt(returnlist[1]), Integer.parseInt(returnlist[2]));
		
	}
	if(event.getSource() == forward && forward.isEnabled()) {
		undoPointer++;
		if(undoPointer >= undoList.size()-1) {
			undoPointer = undoList.size()-1;
			forward.setEnabled(false);
		}
		back.setEnabled(true);
		returnlist = parseSearch(undoList.get(undoPointer));			
		FileRead.search = true;
		searchField.setText(undoList.get(undoPointer));
		drawCanvas.gotoPos(returnlist[0], Integer.parseInt(returnlist[1]), Integer.parseInt(returnlist[2]));
		
	}
	
}

static void addGenomeFile(String genomeName) {
	Main.genomehash.put(genomeName, new ArrayList<File>());				
	JMenu addMenu = new JMenu(genomeName);
	addMenu.addMouseListener(thisMainListener);
	addMenu.setName(genomeName);	
	if(Main.drawCanvas != null) {
		clicked = false;
		refModel.removeElementAt(refModel.getSize()-1);
		refModel.addElement(genomeName);		
		refModel.addElement("Add new reference...");
		clicked = true;
		genome.add(addMenu);
		genome.revalidate();
	}
}
static void removeAnnotationFile(String genomeName, String annotationFile) {
	
	
	try {
	
	for(int i = 1 ; i<genome.getItemCount(); i++) {			
		if(genome.getItem(i).getName().equals(genomeName)) {				
			JMenu addMenu = (JMenu)genome.getItem(i);	
			for(int j = 0 ; j<addMenu.getItemCount();j++) {
				if(addMenu.getItem(j) == null || addMenu.getItem(j).getText() == null) {
					continue;
				}
				
				if(annotationFile.contains(addMenu.getItem(j).getText())) {
					addMenu.remove(j);
					addMenu.revalidate();
					break;
				}
			}
			break;
		}
	}
	for(int i = 0; i<genomehash.get(genomeName).size();i++) {
		
		if(genomehash.get(genomeName).get(i).getName().contains(annotationFile.replace(".gff3.gz", ""))) {
			genomehash.get(genomeName).remove(i);			
			break;
		}		
	}
	
	Main.defaultAnnotation = "";
	setChromDrop(genomeName);
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}
static void addAnnotationFile(String genomeName, File annotationFile) {
	boolean first = false;
	
	if(genomehash.get(genomeName) == null) {
		genomehash.put(genomeName, new ArrayList<File>());		
	}		
	if(genome.getItemCount() == 0) {
		first = true;
	}
	
	genomehash.get(genomeName).add(annotationFile);	
	JMenuItem additem = new JMenuItem(annotationFile.getName().substring(0,annotationFile.getName().indexOf(".bed.gz")));
	additem.setName(annotationFile.getName().substring(0,annotationFile.getName().indexOf(".bed.gz")));
	additem.addMouseListener(Main.thisMainListener);
	
	for(int i = 1 ; i<genome.getItemCount(); i++) {			
		
		if(genome.getItem(i).getName().equals(genomeName)) {				
			JMenu addMenu = (JMenu)genome.getItem(i);	
			addMenu.setFont(menuFont);
			if(first) {
				JMenuItem addAnnotation = new JMenuItem("Add new annotation file...");
				addAnnotation.setFont(menuFont);
				addAnnotation.addMouseListener(Main.thisMainListener);
				addAnnotation.setName("add_annotation");
				addMenu.add(addAnnotation);
				JLabel addLabel = new JLabel("  Select annotation: ");
				addLabel.setFont(menuFont);
				labels.add(addLabel);
				addMenu.add(addLabel);
				addMenu.add(new JSeparator());
			}
			additem.setFont(menuFont);
			addMenu.add(additem);
			addMenu.revalidate();
			genome.revalidate();
			break;
		}
	}	
	Main.defaultAnnotation = annotationFile.getName();
	setChromDrop(genomeName);
}

public static void putMessage(String message) {
	if(message != null) {
		widthLabel.setBackground(new Color(240,210,160));
		widthLabel.setText(message);
		
    	chromDraw.timer = System.currentTimeMillis();
    	chromDraw.repaint();
	}
	else {
		widthLabel.setBackground(new Color(250,250,250));
		widthLabel.setText(MethodLibrary.formatNumber(drawCanvas.splits.get(0).chromEnd) +"bp");
	}
}
public class CheckUpdates extends SwingWorker<String, Object> {
	private class DefaultTrustManager implements X509TrustManager {

        @Override
        public void checkClientTrusted(X509Certificate[] arg0, String arg1) throws CertificateException {}

        @Override
        public void checkServerTrusted(X509Certificate[] arg0, String arg1) throws CertificateException {}

        @Override
        public X509Certificate[] getAcceptedIssuers() {
            return null;
        }
	}

	protected String doInBackground() {
		try {/*
			if(!new File(userDir+"/genomes/ensembl.txt").exists()) {
				
				Main.downloadFile(new URL("https://www.cs.helsinki.fi/u/rkataine/BasePlayer/update/ensembl.txt"), userDir+"/genomes/", 0);
			
		    }*/
			SSLContext ctx = SSLContext.getInstance("TLS");
	        ctx.init(new KeyManager[0], new TrustManager[] {new DefaultTrustManager()}, new SecureRandom());
	        SSLContext.setDefault(ctx);
		    URL file = new URL("https://baseplayer.fi/update/BasePlayer.jar");
			
		    HttpsURLConnection httpCon = (HttpsURLConnection) file.openConnection();	 
		    //HttpURLConnection httpCon = (HttpURLConnection) file.openConnection();	 
		    httpCon.setHostnameVerifier(new HostnameVerifier() {
	            @Override
	            public boolean verify(String arg0, SSLSession arg1) {
	                return true;
	            }
	        });
		    File homefile = new File(userDir +"/BasePlayer.jar");
		   
		    httpCon.connect();
		   
		    if(httpCon.getLastModified() != homefile.lastModified()) {
		    	//putMessage("Updates available. Please click File->Update to get the most recent version.");
		    	update.setVisible(true);
		    	update.setEnabled(true);
		    	update.setForeground(Color.green);
		    	
		    }
		    else {
		    	//putMessage("BasePlayer is up-to-date.");			    	
		    	update.setEnabled(false);
		    	update.setVisible(false);
		    }
		    httpCon.disconnect();
		    file = new URL("https://baseplayer.fi/update/Launcher.jar");
		    homefile = new File(userDir +"/Launcher.jar");
		    httpCon = (HttpsURLConnection) file.openConnection();
		    httpCon.setHostnameVerifier(new HostnameVerifier() {
	            @Override
	            public boolean verify(String arg0, SSLSession arg1) {
	                return true;
	            }
	        });
		    httpCon.connect();
		    if(!homefile.exists() || httpCon.getLastModified() != homefile.lastModified()) {
		    	updatelauncher = true;					
		    }
		    else {
		    	updatelauncher = false;
		    }
		    
		    httpCon.disconnect();
		    
		   
		    
		}
		catch(Exception e) {
			showError(e.getMessage(), "Error");
			
			update.setEnabled(false);
		}
		try {
			BufferedReader selexReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("SELEX/TFbinding_PFMs.txt")));
			String line, factor;
			String[] split, matrix, values;
			int[][] selexmatrix;
			while((line = selexReader.readLine()) != null) {
				split = line.split("\\t");
				factor = split[0].replace(".pfm", "");				
				matrix = split[1].split(";");				
				values = matrix[0].split(",");
				selexmatrix = new int[4][values.length];
				for(int i = 0; i<4; i++ ) {				
					values = matrix[i].split(",");
					for(int j = 0; j< values.length; j++) {						
						selexmatrix[i][j] =  Integer.parseInt(values[j]);
										
					}
				}	
				factorNames.put(factor, factor);
				Main.SELEXhash.put(factor, selexmatrix);
			}				
			selexReader.close();
			
			String id;
			boolean first = true;
			int linepointer = 0;			
			if(new File(userDir +"/additions/motifs/").exists()) {
				File[] files = new File(userDir +"/additions/motifs/").listFiles();
				for (File file : files) { 
					selexReader = new BufferedReader(new FileReader(file));			
				
					while((line = selexReader.readLine()) != null) {
						try {
						if(line.startsWith(" ")) {
							continue;
						}
						if(line.startsWith(">")) {
							first = true;
							split = line.split("\\s+");
							id = split[0].substring(1);					
							factor = split[1];
							
							factorNames.put(id, factor);
						
							linepointer = 0;
							line = selexReader.readLine();
							selexmatrix = null;
						
							while(line.length() > 1 && !line.startsWith(" ")) {
								
								split = line.substring(line.indexOf("[") +1, line.indexOf("]")).trim().split("\\s+");
								
								if (first) {
									first = false;
									selexmatrix = new int[4][split.length];							
								}
								
								for(int j = 0; j< split.length; j++) {						
									selexmatrix[linepointer][j] =  Integer.parseInt(split[j]);												
								}
								linepointer++;
								line = selexReader.readLine();
							}
							Main.SELEXhash.put(id, selexmatrix);
						}				
						}
						catch(Exception e) {
							
							e.printStackTrace();
						}
					}				
					selexReader.close();			
				}
			}
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		return "";
	}
}
static void showError(String error, String dialogtype) {
	VariantHandler.frame.setAlwaysOnTop(false);
	if(dialogtype.equals("Error")) {
	
		JOptionPane.showMessageDialog(Main.drawCanvas, error, dialogtype, JOptionPane.ERROR_MESSAGE);
	}
	
	else {
		JOptionPane.showMessageDialog(Main.drawCanvas, error, dialogtype, JOptionPane.INFORMATION_MESSAGE);
	}
	VariantHandler.frame.setAlwaysOnTop(true);
	
}
public static void gotoURL(String url) {
	 if( !java.awt.Desktop.isDesktopSupported() ) {

         System.err.println( "Desktop is not supported" );
         System.exit( 1 );
     }
     
     java.awt.Desktop desktop = java.awt.Desktop.getDesktop();

     if( !desktop.isSupported( java.awt.Desktop.Action.BROWSE ) ) {

         System.err.println( "Desktop doesn't support the browse action" );
         System.exit( 1 );
     }

         try {

             java.net.URI uri = new java.net.URI(url);
             desktop.browse( uri );
         }
         catch ( Exception ex ) {

             System.err.println( ex.getMessage() );
             ex.printStackTrace();
         }	  
}



public static String checkFile(URL url, ArrayList<String> others) throws IOException {
	
		URLConnection httpConn = (URLConnection) url.openConnection();
		String fileURL = url.getPath();		
		
		// always check HTTP response code first
		httpConn.connect();
		String fileName = "";
		
		String disposition = httpConn.getHeaderField("Content-Disposition");
		
		if (disposition != null) {
			// extracts file name from header field
			int index = disposition.indexOf("filename=");
			if (index > 0) {
				fileName = disposition.substring(index + 10,
						disposition.length() - 1);
			}
		} else {
			// extracts file name from URL
			fileName = fileURL.substring(fileURL.lastIndexOf("/") + 1);
		}
		InputStream inputStream = null;		
		
		// opens input stream from the HTTP connection
		try {
			inputStream = httpConn.getInputStream();
			
		}
		catch(Exception e) {
			
			if(fileName.endsWith(".gff3.gz")) {
				String urldir = fileURL.substring(0,fileURL.lastIndexOf("/")+1);				
				fileName = getNewFile(url.getHost(), urldir, fileName);
				
				if(!others.contains(fileName.replace(".gff3.gz", ""))) {		
				
					return fileName;
				}
				else {									
					return "";
				}
				
			}
		}
		
		
		if(inputStream == null) {
			return fileName;
		}
		else {
			inputStream.close();			
			if(!others.contains(fileName.replace(".gff3.gz", ""))) {		
				return fileName;
			}
			return "";
		}
}

public static String downloadFile(URL url, String saveDir, int size) throws IOException {
	
	URLConnection httpConn = (URLConnection) url.openConnection();
	String fileURL = url.getPath();
	
	
		int BUFFER_SIZE = 4096;
	// always check HTTP response code first
	httpConn.connect();
		String fileName = "";
		
		String disposition = httpConn.getHeaderField("Content-Disposition");
		
		
		
		if (disposition != null) {
			// extracts file name from header field
			int index = disposition.indexOf("filename=");
			if (index > 0) {
				fileName = disposition.substring(index + 10,
						disposition.length() - 1);
			}
		} else {
			// extracts file name from URL
			fileName = fileURL.substring(fileURL.lastIndexOf("/") + 1);
		}
		InputStream inputStream = null;
		
		
		
		// opens input stream from the HTTP connection
		try {
			inputStream = httpConn.getInputStream();
		}
		catch(Exception e) {
			if(fileName.endsWith(".gff3.gz")) {
				String urldir = fileURL.substring(0,fileURL.lastIndexOf("/")+1);				
				fileName = getNewFile(url.getHost(), urldir, fileName);
				url = new URL(url.getProtocol()+"://" +url.getHost() +"/" +urldir +"/" +fileName);
				httpConn = (URLConnection) url.openConnection();
				inputStream = httpConn.getInputStream();
				
			}
		}
		
		if(inputStream == null) {
			return "";
		}
		if(fileName.endsWith(".gff3.gz")) {
			saveDir = saveDir += "/" +fileName +"/";
			new File(saveDir).mkdirs();
		}
		String saveFilePath = saveDir + File.separator + fileName;
		
		// opens an output stream to save into file
		FileOutputStream outputStream = new FileOutputStream(saveFilePath);
		long downloaded = 0;
		int bytesRead = -1, counter=0;
		String loading = "";
		byte[] buffer = new byte[BUFFER_SIZE];
		if(Main.drawCanvas != null) {
			loading = drawCanvas.loadingtext;
			Main.drawCanvas.loadingtext = loading + " 0MB";
		}
		while ((bytesRead = inputStream.read(buffer)) != -1) {
			outputStream.write(buffer, 0, bytesRead);
			downloaded +=buffer.length;
			if(Main.drawCanvas != null) {
				counter++;
				if(counter == 100) {
					
					Main.drawCanvas.loadingtext = loading + " " +(downloaded/1048576) +"/~" +(size/1048576) +"MB";
					Main.drawCanvas.loadBarSample = (int)(downloaded/(double)size*100);
					Main.drawCanvas.loadbarAll = (int)(downloaded/(double)size*100);
					if(Main.drawCanvas.loadBarSample > 100) {
						Main.drawCanvas.loadBarSample = 100;
						Main.drawCanvas.loadbarAll = 100;
					}
					counter = 0;
				}
			}
		}
		if(Main.drawCanvas != null) {
			Main.drawCanvas.loadBarSample = 0;
			Main.drawCanvas.loadbarAll = 0;
		}
		outputStream.close();
		inputStream.close();
		return fileName;
}



static String getNewFile(String server, String folder, String oldfile) {
	String minfile = "";
	try {
		
		String filename = oldfile;
		FTPClient f = new FTPClient();		
		f.connect(server);	
		f.enterLocalPassiveMode();
		f.login("anonymous", "");		
		FTPFile[] files = f.listFiles(folder);
		int left = 0, right = filename.length()-1, distance = 0;		
		int mindistance = 100;
		
		for(int i = 0 ; i<files.length; i++) {			
			if(files[i].getName().endsWith(".gff3.gz")) {
				distance = 0; 
				right = Math.min(filename.length(), files[i].getName().length());
				left = 0;
			
				while(left<right) {
					if(filename.charAt(left) != files[i].getName().charAt(left)) {
						distance++;								
					}		
					left++;
				}
				distance += Math.abs(filename.length() - files[i].getName().length());
				if(distance < mindistance) {
					mindistance = distance;
					minfile = files[i].getName();					
				}
			}					
		}		
	}
	catch(Exception ex) {
		ex.printStackTrace();
	}
	return minfile;
}
public class Updater extends SwingWorker<String, Object> {
	private class DefaultTrustManager implements X509TrustManager {

        @Override
        public void checkClientTrusted(X509Certificate[] arg0, String arg1) throws CertificateException {}

        @Override
        public void checkServerTrusted(X509Certificate[] arg0, String arg1) throws CertificateException {}

        @Override
        public X509Certificate[] getAcceptedIssuers() {
            return null;
        }
	}
	
	public void downloadFile(String fileURL, String saveDir) throws IOException {
		
		System.out.println("Updating file from: " +fileURL +"\nSaving to: " +saveDir);
		System.out.println("Opening connection...");
		try {
		URL url = new URL(fileURL);
	
		//HttpsURLConnection httpConn = (HttpsURLConnection) url.openConnection();
		//httpConn.connect();
		SSLContext ctx = SSLContext.getInstance("TLS");
        ctx.init(new KeyManager[0], new TrustManager[] {new DefaultTrustManager()}, new SecureRandom());
        SSLContext.setDefault(ctx);
		HttpsURLConnection httpConn = (HttpsURLConnection) url.openConnection();
		httpConn.setHostnameVerifier(new HostnameVerifier() {
            @Override
            public boolean verify(String arg0, SSLSession arg1) {
                return true;
            }
        });
		httpConn.connect();
		int responseCode = httpConn.getResponseCode();
		int BUFFER_SIZE = 4096;
		System.out.println("Ready.");
		// always check HTTP response code first
		if (responseCode == HttpsURLConnection.HTTP_OK) {
			String fileName = "";
			String disposition = httpConn.getHeaderField("Content-Disposition");
		//	String contentType = httpConn.getContentType();
			int contentLength = httpConn.getContentLength();

			if (disposition != null) {
				// extracts file name from header field
				int index = disposition.indexOf("filename=");
				if (index > 0) {
					fileName = disposition.substring(index + 10,
							disposition.length() - 1);
				}
			} else {
				// extracts file name from URL
				fileName = fileURL.substring(fileURL.lastIndexOf("/") + 1,
						fileURL.length());
			}

			// opens input stream from the HTTP connection
			InputStream inputStream = httpConn.getInputStream();
			String saveFilePath = saveDir + File.separator + fileName +"_temp";
			File testFile = new File(saveDir + File.separator +"_test");
			// opens an output stream to save into file
			if(!testFile.mkdir()) {
				Main.showError("Could not update BasePlayer. No writing permissions in BasePlayer folder.\nStart BasePlayer as adminstrator and press update again.\n"
						+ "Alternatively, give writing permissions to your BasePlayer folder.", "Error");
				testFile.delete();
				return;
			}
			testFile.delete();
			FileOutputStream outputStream = new FileOutputStream(saveFilePath);

			int bytesRead = -1;
			byte[] buffer = new byte[BUFFER_SIZE];
			while ((bytesRead = inputStream.read(buffer)) != -1) {
				outputStream.write(buffer, 0, bytesRead);
			}
			outputStream.close();
			inputStream.close();
			File tempfile = new File(saveFilePath), newfile = new File(saveDir + File.separator + fileName);
			
			if(tempfile.length() == contentLength) {
				
					 
					  InputStream in = new FileInputStream(tempfile);
					  
					  //For Append the file.
					//  OutputStream out = new FileOutputStream(f2,true);

					  //For Overwrite the file.
					  OutputStream out = new FileOutputStream(newfile);

					  byte[] buf = new byte[1024];
					  int len;
					  while ((len = in.read(buf)) > 0){
						  out.write(buf, 0, len);
					  }
					  in.close();
					  out.close();
					 
										  
				newfile.setLastModified(httpConn.getLastModified());
				tempfile.delete();
				if(!fileName.contains("Launcher")) {
					Main.showError("BasePlayer updated. Please restart program to apply changes.", "Note");
				}
				update.setEnabled(false);
			}
			else {
				Main.showError("BasePlayer could not be updated.", "Note");
			}				

			System.out.println("File downloaded");
		} else {
			System.out.println("No file to download. Server replied HTTP code: " + responseCode);
			ErrorLog.addError("No file to download. Server replied HTTP code: " + responseCode);
		}
			httpConn.disconnect();
		}
		catch(Exception ex) {
			ex.printStackTrace();
		}
	}
	
	protected String doInBackground() {
		
		try {
			Main.drawCanvas.loading("Updating BasePlayer... (downloading BasePlayer.jar from https://baseplayer.fi/update/");			    
		  
		    downloadFile("https://baseplayer.fi/update/BasePlayer.jar", userDir);
			if(updatelauncher) {
				downloadFile("https://baseplayer.fi/update/Launcher.jar", userDir);
			}
			/*if(!new File(userDir+"/genomes/ensembl.txt").exists()) {
				downloadFile("https://www.cs.helsinki.fi/u/rkataine/BasePlayer/update/ensembl.txt", userDir+"/genomes/");
		    }*/
			Main.drawCanvas.ready("Updating BasePlayer... (downloading BasePlayer.jar from https://baseplayer.fi/update/");	
			/*if(homefile.length() > 0) {
				File replacefile = new File(userDir +"/BasePlayer.jar");
				FileUtils.copyFile(homefile, replacefile);				
				FileUtils.forceDelete(homefile);
			}
			else {
				
				JOptionPane.showMessageDialog(Main.chromDraw,"BasePlayer couldn't be updated. Please try again.", "Note", JOptionPane.INFORMATION_MESSAGE);
				return "";
			}
		   /*
		    file = new URL("https://www.cs.helsinki.fi/u/rkataine/Rikurator/update/Launcher.jar");
		    homefile = new File(userDir +"/Launcher.jar");
			
			FileUtils.copyURLToFile(file, homefile);
		    */				
				
			
			
		   
		    
		}
		catch(Exception e) {
			 Main.drawCanvas.ready("Updating BasePlayer... (downloading BasePlayer.jar from https://baseplayer.fi/update/");	
			e.printStackTrace();
			ErrorLog.addError(e.getStackTrace());
			Main.showError(e.getMessage(), "Error");
		
			
		}
		return "";
	}
	
}
@Override
public void mouseEntered(MouseEvent event) {
	
	
	if(event.getComponent().getName() != null) {
		
		if(event.getComponent() instanceof JMenu) {
			hoverGenome = event.getComponent().getName();
			hoverAnnotation = "";
			
		}
		else if(event.getComponent() instanceof JMenuItem) {
			
			hoverAnnotation = event.getComponent().getName();
		}
			
		
		
	}
//	if(event.getSource() == manage) {
		
		//Logo.frame.setVisible(false);
		//Syste
		//filemenu.setSelected(false);
		//toolmenu.setSelected(false);
		//help.setSelected(false);
		//filemenu.getPopupMenu().setVisible(false);
		//toolmenu.getPopupMenu().setVisible(false);
		//help.getPopupMenu().setVisible(false);
		//manage.setFocusPainted(true);
		
//	}

	/*if(event.getSource() == filemenu) {
		//Logo.frame.setVisible(false);
		filemenu.doClick();
	}
	else if(event.getSource() == toolmenu) {
		//Logo.frame.setVisible(false);
		toolmenu.doClick();
	}
	else if(event.getSource() == help) {
	//	Logo.frame.setVisible(false);
		help.doClick();
	}*/
	
	
	
}


@Override
public void mouseExited(MouseEvent event) {
	
//	if(!drawCanvas.loading && !drawCanvas.scrollbar && event.getSource() == drawScroll.getVerticalScrollBar()) {
		
	//	Draw.setGlasspane(false);
//	}
	if(event.getSource() == drawCanvas) {
		Main.drawCanvas.sidebar = false;
	
	//	Main.drawCanvas.selectedSampleIndex = -1;
			Main.drawCanvas.repaint();
	
		}
		
	}

	@Override
public void mousePressed(MouseEvent event) {
	//if(Logo.frame.isVisible()) {
//	frame.requestFocus();
//	}
//	Logo.frame.setVisible(false);
if(event.getSource() == refDropdown) {
	switch(event.getModifiers()) {	
		case InputEvent.BUTTON1_MASK: {	
			if(Main.genomehash.size() == 0) {
				if(AddGenome.frame == null) {
				  AddGenome.createAndShowGUI();
			  	}
				AddGenome.frame.setTitle("Add new genome");
				AddGenome.annotation = false;
				AddGenome.remove.setEnabled(false);
				AddGenome.download.setEnabled(false);
				AddGenome.frame.setVisible(true);
				AddGenome.frame.setLocation(frame.getLocationOnScreen().x+frame.getWidth()/2 - AddGenome.frame.getWidth()/2, frame.getLocationOnScreen().y+frame.getHeight()/6);
				 
			   AddGenome.frame.setState(JFrame.NORMAL);
			}
			rightclick = false;
			break;
		}
		case InputEvent.BUTTON3_MASK: {	
			rightclick = true;
			break;
		}
	}
}
else if(event.getSource() == geneDropdown) {
	switch(event.getModifiers()) {	
	case InputEvent.BUTTON1_MASK: {	
		if(Main.genomehash.size() == 0) {
			if(AddGenome.frame == null) {
				AddGenome.createAndShowGUI();
			}
			AddGenome.frame.setTitle("Add new genome");
			AddGenome.annotation = false;
			AddGenome.remove.setEnabled(false);
			AddGenome.download.setEnabled(false);
			AddGenome.frame.setVisible(true);
			AddGenome.frame.setLocation(frame.getLocationOnScreen().x+frame.getWidth()/2 - AddGenome.frame.getWidth()/2, frame.getLocationOnScreen().y+frame.getHeight()/6);
			 
		   AddGenome.frame.setState(JFrame.NORMAL);
		}
		rightclick = false;
		break;
	}
	case InputEvent.BUTTON3_MASK: {	
		rightclick = true;
		break;
	}
}
}
else if(event.getSource() == chromlabel) {
	chromosomeDropdown.showPopup();
	
}

else if(event.getSource() == splitPaneDivider) {
	Main.vardivider = bedCanvas.nodeImage.getHeight()/(double)varPaneDivider.getY();
	
//		Main.bedCanvas.resize = true;
}
else if(event.getSource() == varPaneDivider) {
//		Main.bedCanvas.resize = true;			
	Main.vardivider = bedCanvas.nodeImage.getHeight()/(double)varPaneDivider.getY();	
	
}
else if(event.getSource() == filemenu) {			
	/*if(!filemenu.isSelected()){				
		filemenu.doClick();			
	}
	*/
}
else if(event.getSource() == toolmenu) {
	/*if(!toolmenu.isSelected()){				
		toolmenu.doClick();			
	}*/
	
}
else if(drawCanvas.loadingtext.equals("note")) {
	Main.drawCanvas.loadingtext = "";
	Main.drawCanvas.ready("note");
	
}

else if(event.getSource() == drawScroll.getVerticalScrollBar()) {
	
	 if(Main.glassPane.getCursor().getType() != Cursor.WAIT_CURSOR) {
		 Main.glassPane.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
		}
	Draw.setGlasspane(true);			
	
}
else if(event.getSource() == searchField) {
	searchField.requestFocus();
	searchField.setForeground(Color.black);		
	if(searchField.getText().contains("Search by")) {
		searchField.setText("");
	}
	
}
else if(event.getSource() == addGenome) {
		if(AddGenome.frame == null) {
			AddGenome.createAndShowGUI();
		}
		AddGenome.frame.setTitle("Add new genome");
		AddGenome.annotation = false;
		AddGenome.remove.setEnabled(false);
		AddGenome.download.setEnabled(false);
		AddGenome.frame.setVisible(true);
		AddGenome.frame.setLocation(frame.getLocationOnScreen().x+frame.getWidth()/2 - AddGenome.frame.getWidth()/2, frame.getLocationOnScreen().y+frame.getHeight()/6);
		AddGenome.frame.setState(JFrame.NORMAL);
	
}		   
else if(event.getComponent().getName() != null) {
	if(event.getComponent().getName().equals("frame0")) {
		return;
	}
	
	try {
		if(event.getComponent().getName().equals("add_annotation")) {
			if(AddGenome.frame == null) {
				AddGenome.createAndShowGUI();
			}
			AddGenome.annotation = true;
			AddGenome.frame.setTitle("Add new annotation file for " +Main.selectedGenome);
			AddGenome.remove.setEnabled(false);
			AddGenome.download.setEnabled(false);
			AddGenome.frame.setVisible(true);
			AddGenome.frame.setLocation(frame.getLocationOnScreen().x+frame.getWidth()/2 - AddGenome.frame.getWidth()/2, frame.getLocationOnScreen().y+frame.getHeight()/6);
			AddGenome.genomeName.setText(hoverGenome);				
			return;
		}
			
		if(hoverAnnotation.length() > 0) {
			for(int j = 0; j<genomehash.get(hoverGenome).size(); j++) {
				if(genomehash.get(hoverGenome).get(j).getName().contains(hoverAnnotation)) {
					annotationfile = genomehash.get(hoverGenome).get(j).getName();
					 Main.annotation = j;
					 break;
				}
			}
		
		
		 defaultGenome = hoverGenome;
		
	     setChromDrop(defaultGenome); 		
	     getBands();
	     if(genomehash.get(defaultGenome).size() > 0 && genomehash.get(defaultGenome).get(annotation) != null) {
	    	 
	    	 changeAnnotation(annotation);
	   
	     }
			 drawCanvas.chrom = chromosomeDropdown.getItemAt(0);
			 chromosomeDropdown.setSelectedIndex(0);
		}
		}
		catch(Exception e) {
			e.printStackTrace();
		}				  
	 }	
}

public void changeRef(String dir) {
	   writeToConfig("DefaultGenome=" +dir);
   defaultGenome = dir;
   setChromDrop(Main.refDropdown.getSelectedItem().toString()); 		  
   getBands();
   try {
	   if(genomehash.get(dir).size() > 0) {
		  
		    ChromDraw.exonReader = new TabixReader(genomehash.get(dir).get(0).getCanonicalPath());
			String s;
			String[] exonSplit;
			searchTable.clear();
			while((s = ChromDraw.exonReader.readLine()) != null) {
				exonSplit = s.split("\t");
				if(!searchTable.containsKey(exonSplit[3].toUpperCase())) {		
					
					String[] adder = {exonSplit[0], exonSplit[1], exonSplit[2]};
					searchTable.put(exonSplit[3].toUpperCase(), adder);			
					if(exonSplit[6].contains(":")) {
						geneIDMap.put(exonSplit[6].split(":")[1].toUpperCase(), exonSplit[3].toUpperCase());
						}
						else {
							geneIDMap.put(exonSplit[6].toUpperCase(), exonSplit[3].toUpperCase());
						}
					}	
					else {
						
						if(Integer.parseInt(searchTable.get(exonSplit[3].toUpperCase())[1]) > Integer.parseInt(exonSplit[1])) {
							searchTable.get(exonSplit[3].toUpperCase())[1] = exonSplit[1];
						}
						if(Integer.parseInt(searchTable.get(exonSplit[3].toUpperCase())[2]) < Integer.parseInt(exonSplit[2])) {
							searchTable.get(exonSplit[3].toUpperCase())[2] = exonSplit[2];
						}
					}					
				}
			
		   }
		   }
		   catch(Exception e) {
			   e.printStackTrace();
		   }
	   changeAnnotation(0);
	   chromosomeDropdown.setSelectedIndex(0);
	   drawCanvas.chrom = chromosomeDropdown.getItemAt(0);
  }
	
	
	public void changeAnnotation(int annotation) {
		try {
			if(genomehash.get(defaultGenome).size() < 1) {
				return;
			}
			Main.annotation = annotation;
			Main.defaultAnnotation = genomehash.get(defaultGenome).get(annotation).getName();
			Main.annotationfile = Main.defaultAnnotation;
			 writeToConfig("DefaultGenome=" +defaultGenome);
	     writeToConfig("DefaultGenes=" +annotationfile);
	 	ChromDraw.exonReader = new TabixReader(genomehash.get(defaultGenome).get(annotation).getCanonicalPath());			    
		String s;
		String[] exonSplit;
		searchTable.clear();
		geneIDMap.clear();
		while((s = ChromDraw.exonReader.readLine()) != null) {
		
			exonSplit = s.split("\t");
			if(!searchTable.containsKey(exonSplit[3].toUpperCase())) {		
				
				String[] adder = {exonSplit[0], exonSplit[1], exonSplit[2]};
				searchTable.put(exonSplit[3].toUpperCase(), adder);	
				if(exonSplit[6].contains(":")) {
					geneIDMap.put(exonSplit[6].split(":")[1].toUpperCase(), exonSplit[3].toUpperCase());
				}
				else {
					geneIDMap.put(exonSplit[6].toUpperCase(), exonSplit[3].toUpperCase());
				}
			}	
			else {
				
				if(Integer.parseInt(searchTable.get(exonSplit[3].toUpperCase())[1]) > Integer.parseInt(exonSplit[1])) {
					searchTable.get(exonSplit[3].toUpperCase())[1] = exonSplit[1];
				}
				if(Integer.parseInt(searchTable.get(exonSplit[3].toUpperCase())[2]) < Integer.parseInt(exonSplit[2])) {
					searchTable.get(exonSplit[3].toUpperCase())[2] = exonSplit[2];
				}
			}
		}
		ChromDraw.exonReader.close();
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	  

}

@Override
public void mouseReleased(MouseEvent event) {
	
	if(event.getSource() == splitPaneDivider) {
		Main.bedCanvas.resize = false;
		Main.bedCanvas.repaint();			
	}
	if(event.getSource() == varPaneDivider) {
		
		Main.bedCanvas.resize = false;
		Main.bedCanvas.repaint();
		Draw.updatevars = true;
		Draw.updateReads = true;
		Draw.updateCoverages = true;
		drawCanvas.repaint();
	}
	
		
	if(event.getSource() == drawScroll.getVerticalScrollBar()) {
	
		if(Main.drawScroll.getVerticalScrollBar().getValue() > (drawCanvas.drawVariables.visiblestart*drawCanvas.drawVariables.sampleHeight+drawCanvas.drawVariables.sampleHeight/2.0)) {
			drawCanvas.drawVariables.visiblestart++;
			
			Draw.setScrollbar((int)(drawCanvas.drawVariables.visiblestart*drawCanvas.drawVariables.sampleHeight));
		}
		else {
			Draw.setScrollbar((int)(drawCanvas.drawVariables.visiblestart*drawCanvas.drawVariables.sampleHeight));
		}
		
	}
	
	if(!drawCanvas.loading) {
		Draw.setGlasspane(false);
	}
	drawCanvas.scrolldrag = false;
	
	if(Main.bedCanvas.bedTrack.size() > 0) {
		for(int i = 0 ; i<Main.bedCanvas.bedTrack.size(); i++) {
			
			Main.bedCanvas.getMoreBeds(Main.bedCanvas.bedTrack.get(i));
		}
	}
}
@Override
public void componentHidden(ComponentEvent arg0) {
	
}
@Override
public void componentMoved(ComponentEvent arg0) {
	
	Logo.frame.setLocation(frame.getLocation().x+(int)width/2-300, frame.getLocation().y+(int)height/2-300);
}
@Override
public void componentShown(ComponentEvent arg0) {
	
}
//	@Override
/*	public void propertyChange(PropertyChangeEvent event) {
	
		if(drawScroll.getViewport().getWidth() > 0) {
		//	System.out.println(event.getSource());
		}
			//chromPaneDivider
			/*if(event.getSource() == upPanel) {
				 
				Main.sidebarWidth = upPanel.getDividerLocation()+2;
				chromDimensions.setSize(drawScroll.getViewport().getWidth()-upPanel.getDividerLocation()-1, splitPane.getDividerLocation());
				chromDraw.setPreferredSize(chromDimensions);
			
			//	if(samples > 0) {					
					if(samples*drawCanvas.drawVariables.sampleHeight < drawScroll.getViewport().getHeight()) {						
						drawDimensions.setSize(drawScroll.getViewport().getWidth(),drawScroll.getViewport().getSize().height );	
						drawCanvas.setPreferredSize(drawDimensions);
						
						drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), drawScroll.getViewport().getHeight());	
						if(drawCanvas.splits.size() > 0) {
							for(int i = 0 ; i<drawCanvas.splits.size(); i++) {
								drawCanvas.splits.get(i).updateReads = true;
							
							}
						}
						
					}					
					drawCanvas.repaint();
					
			//	}
				chromDraw.updateExons = true;
				chromDraw.repaint();
			}
			if(event.getSource() == splitPaneDivider) {		
				
				chromDimensions.setSize(drawScroll.getViewport().getWidth()-Main.sidebarWidth-1, splitPane.getDividerLocation());
				chromDraw.setPreferredSize(chromDimensions);
				
				if(samples > 0) {					
					if(samples*drawCanvas.drawVariables.sampleHeight < drawScroll.getViewport().getHeight()) {						
						drawDimensions.setSize(drawScroll.getViewport().getWidth(),drawScroll.getViewport().getSize().height );	
						drawCanvas.setPreferredSize(drawDimensions);
						
						drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), drawScroll.getViewport().getHeight());	
						if(drawCanvas.splits.size() > 0) {
							for(int i = 0 ; i<drawCanvas.splits.size(); i++) {
								drawCanvas.splits.get(i).updateReads = true;
							
							}
						}
						
					}					
					drawCanvas.repaint();
					
				}
				chromDraw.updateExons = true;
				chromDraw.repaint();
			}
		/*	if(event.getSource() == varpane) {			
				
				if(samples > 0) {		
					
					if(samples*drawCanvas.drawVariables.sampleHeight < drawScroll.getViewport().getHeight()) {						
						drawDimensions.setSize(drawScroll.getViewport().getWidth(),drawScroll.getViewport().getSize().height );	
						drawCanvas.setPreferredSize(drawDimensions);
						chromDimensions.setSize(drawScroll.getViewport().getWidth()-Main.sidebarWidth-1,  splitPane.getDividerLocation());
						chromDraw.setPreferredSize(chromDimensions);
						drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), drawScroll.getViewport().getHeight());	
						if(drawCanvas.splits.size() > 0) {
							for(int i = 0 ; i<drawCanvas.splits.size(); i++) {
								drawCanvas.splits.get(i).updateReads = true;
							
							}
						}
						
					}
					
					drawCanvas.repaint();
				}
				
			}	
		}
		*/
//	}

static void setFonts() {		
	if(Settings.bold.isSelected()) {
		menuFont = new Font("SansSerif", Font.BOLD, Main.defaultFontSize);
		menuFontBold = new Font("SansSerif", Font.BOLD, Main.defaultFontSize+1);
		
	}
	else {
		menuFont = new Font("SansSerif", Font.PLAIN, Main.defaultFontSize);
		menuFontBold = new Font("SansSerif", Font.BOLD, Main.defaultFontSize);
	}
	
	
	Draw.defaultFont = menuFont.deriveFont((float)Main.defaultFontSize-1);
	
	for(int i = 0 ; i<menubar.getComponentCount(); i++) {
		menubar.getComponent(i).setFont(Main.menuFont);
	}

	Draw.loadingFont = menuFont.deriveFont((float)(Main.defaultFontSize*1.5));//		
	buttonHeight = (int)(Main.defaultFontSize*1.5);
	buttonWidth = Main.defaultFontSize*6;
	//searchField.setMargin(new Insets(0,Main.defaultFontSize+8, 0, 0));
	searchField.setBorder(BorderFactory.createCompoundBorder(searchField.getBorder(), BorderFactory.createEmptyBorder(0,Main.defaultFontSize+8, 0, 0)));
	buttonDimension = new Dimension(buttonWidth, buttonHeight);		
	
	ChromDraw.seqFont= ChromDraw.seqFont.deriveFont((float)(Main.defaultFontSize+2));
	bedCanvas.buf.setFont(Draw.defaultFont);
	bedCanvas.nodebuf.setFont(Draw.defaultFont);
	bedCanvas.fm = bedCanvas.nodebuf.getFontMetrics();
	for(int i = 0 ; i<bedCanvas.bedTrack.size(); i++) {
		for(int c = 0 ; c<bedCanvas.bedTrack.get(i).getPopup().getComponentCount(); c++) {
			bedCanvas.bedTrack.get(i).getPopup().getComponent(c).setFont(menuFont);
			
		}			
		if(bedCanvas.bedTrack.get(i).getSelector() != null) {
			bedCanvas.bedTrack.get(i).getSelector().setFonts(menuFont);
		}
		bedCanvas.bedTrack.get(i).getLimitField().setPreferredSize(new Dimension(bedCanvas.buf.getFontMetrics().stringWidth("__Value limit__"), Main.defaultFontSize+6));
		bedCanvas.bedTrack.get(i).getLimitField().setMinimumSize(new Dimension(bedCanvas.buf.getFontMetrics().stringWidth("__Value limit__"), Main.defaultFontSize+6));
		
	}
	for(int i = 0 ; i<Control.controlData.fileArray.size(); i++) {
		for(int c = 0 ;c < Control.controlData.fileArray.get(i).getPopupMenu().getComponentCount(); c++) {
			Control.controlData.fileArray.get(i).getPopupMenu().getComponent(c).setFont(menuFont);
			
		}
	}
	Average.setFonts(menuFont);
	menubar.setMargin(new Insets(0,2,0,2));
    filemenu.setMinimumSize(filemenu.getPreferredSize());
    toolmenu.setMinimumSize(toolmenu.getPreferredSize());
    help.setMinimumSize(help.getPreferredSize());
    manage.setPreferredSize(new Dimension(bedCanvas.buf.getFontMetrics().stringWidth("Variant Managerrrrrrrr") +4, buttonHeight));
  
    manage.setMinimumSize(new Dimension(bedCanvas.buf.getFontMetrics().stringWidth("Variant Managerrrrrrrrr") +4, buttonHeight));
    chromlabel.setPreferredSize(new Dimension(bedCanvas.buf.getFontMetrics().stringWidth("..Chrom..") +4, buttonHeight));
    chromlabel.setMinimumSize(new Dimension(bedCanvas.buf.getFontMetrics().stringWidth("..Chrom..") +4, buttonHeight));
	for(int i = 0 ; i < panel.getComponentCount(); i++) {
		panel.getComponent(i).setFont(Main.menuFont);
	}
	for(int i = 0 ; i < filemenu.getItemCount(); i++) {
		if(filemenu.getItem(i) != null) {
			filemenu.getItem(i).setFont(Main.menuFont);
		}			
	}
	for(int i = 0 ; i < toolmenu.getItemCount(); i++) {
		if(toolmenu.getItem(i) != null) {
			toolmenu.getItem(i).setFont(Main.menuFont);
			
		}			
	}
	Main.area.setFont(Main.menuFont);
	for(int i = 0 ; i < help.getItemCount(); i++) {
		if(help.getItem(i) != null) {
			help.getItem(i).setFont(Main.menuFont);
		}			
	}
	for(int i = 0 ; i < genome.getItemCount(); i++) {			
		genome.getItem(i).setFont(Main.menuFont);
		if(genome.getItem(i) instanceof JMenu) {
			JMenu menu = (JMenu)genome.getItem(i);
			for(int j = 0 ; j<menu.getItemCount(); j++) {					
				if(menu.getItem(j) != null) {
					menu.getItem(j).setFont(Main.menuFont);
				}					
			}
		}			
	}
	for(int i = 0 ; i<labels.size(); i++) {
		labels.get(i).setFont(Main.menuFont);
	}
	
	VariantCaller.setFonts(menuFont);
	
	for(int i = 0 ; i<Main.drawCanvas.splits.size(); i++) {
		Main.drawCanvas.splits.get(i).getExonImageBuffer().setFont(Draw.defaultFont);
		Main.drawCanvas.splits.get(i).getReadBuffer().setFont(Draw.defaultFont);
		Main.drawCanvas.splits.get(i).getSelectbuf().setFont(Draw.defaultFont);
	}
	
	
	for(int i = 0 ; i<chrompan.getComponentCount(); i++) {
		if(chrompan.getComponent(i).getName() != null) {
			chrompan.getComponent(i).setFont(menuFontBold);
		}
		else {
			chrompan.getComponent(i).setFont(menuFont);
		}
	}
	if(AddGenome.tree != null) {
		AddGenome.setFonts(menuFont);
	}
	Settings.setFonts(menuFont);
	chromDraw.selectImageBuffer.setFont(Draw.defaultFont);
	chromDraw.chromImageBuffer.setFont(Draw.defaultFont);
	manage.setToolTipText("No variants on screen");
	manage.setMargin(new Insets(0,4,0,4));
	zoomout.setPreferredSize(new Dimension(bedCanvas.buf.getFontMetrics().stringWidth("Zoom outtttttt") +4, buttonHeight));
	zoomout.setMinimumSize(new Dimension(bedCanvas.buf.getFontMetrics().stringWidth("Zoom outtttttt") +4, buttonHeight));
	zoomout.setMargin(new Insets(0,4,0,4));
    fieldDimension = new Dimension(widthLabel.getFontMetrics(widthLabel.getFont()).stringWidth("chrX:000,000,000-000,000,000bp")+4, buttonHeight);
    positionField.setPreferredSize(fieldDimension);
    positionField.setMinimumSize(fieldDimension);
    controlDraw.buf.setFont(Draw.defaultFont);
    controlDraw.nodebuf.setFont(Draw.defaultFont);
    controlDraw.fm = controlDraw.buf.getFontMetrics();
    controlDraw.repaint();	   
    letterlength = chromosomeDropdown.getFontMetrics(chromosomeDropdown.getFont()).stringWidth("E");
	chromosomeDropdown.setPopupWidth(textlength*letterlength+25);
	chromosomeDropdown.revalidate();
	chromosomeDropdown.repaint();	
	chromosomeDropdown.setPreferredSize(new Dimension(Main.defaultFontSize*5,buttonHeight));
    geneDropdown.setPopupWidth(annolength*letterlength);
    refDropdown.setPopupWidth(reflength*letterlength);
    //searchField.setMargin(new Insets(0,buttonHeight+4, 0, 0));
    searchField.setPreferredSize(fieldDimension);
    searchField.setMinimumSize(fieldDimension);
	widthLabel.setPreferredSize(new Dimension(widthLabel.getFontMetrics(widthLabel.getFont()).stringWidth("000,000,000bp (Right click to cancel zoom)  NNNNNNNNNNNNNNNNNNNNNNNN")+10,buttonHeight));
    widthLabel.setMinimumSize(new Dimension(widthLabel.getFontMetrics(widthLabel.getFont()).stringWidth("000,000,000bp") +10,buttonHeight));		  
    back.setFont(menuFont);
	back.setPreferredSize(new Dimension(back.getFontMetrics(back.getFont()).stringWidth(".<<.")+10,buttonDimension.height));	
	forward.setFont(menuFont);
	forward.setPreferredSize(new Dimension(forward.getFontMetrics(forward.getFont()).stringWidth(".>>.")+10,buttonDimension.height));	  
    chromDraw.bounds = chromDraw.chromImageBuffer.getFontMetrics().getStringBounds("K", chromDraw.chromImageBuffer).getWidth();
    chromDraw.cytoHeight = defaultFontSize +10;
    chromDraw.exonDrawY = defaultFontSize*2 +10;
    drawCanvas.sidebuf.setFont(Draw.defaultFont);
    drawCanvas.buf.setFont(Draw.defaultFont);
    drawCanvas.varStringLen = drawCanvas.buf.getFontMetrics().stringWidth(drawCanvas.varloadString);		
    if(VariantHandler.filters != null) {
    	VariantHandler.setFonts(menuFont);
    }
	chromDraw.updateExons = true;
	chromDraw.repaint();
	
    for(int i = 0 ; i <Main.drawCanvas.sampleList.size(); i++ ) {
    	if(Main.drawCanvas.sampleList.get(i).getreadHash() != null) {	    		
    	
	    	for(int j = 0 ; j<Main.drawCanvas.splits.size(); j++) {
	    		if(Main.drawCanvas.sampleList.get(i).getreadHash().get(Main.drawCanvas.splits.get(j))==null) {
	    			continue;
	    		}
	    		double temp = (Main.drawCanvas.sampleList.get(i).getreadHash().get(Main.drawCanvas.splits.get(j)).readHeight+2)/(double)Main.drawCanvas.sampleList.get(i).getreadHash().get(Main.drawCanvas.splits.get(j)).readwheel;
	    		 Main.drawCanvas.sampleList.get(i).getreadHash().get(Main.drawCanvas.splits.get(j)).readfont =new Font("SansSerif", Font.BOLD, defaultFontSize);
	    		 Main.drawCanvas.sampleList.get(i).getreadHash().get(Main.drawCanvas.splits.get(j)).readHeight = defaultFontSize;
	    		 Main.drawCanvas.sampleList.get(i).getreadHash().get(Main.drawCanvas.splits.get(j)).readwheel = (int)((Main.drawCanvas.sampleList.get(i).getreadHash().get(Main.drawCanvas.splits.get(j)).readHeight+2)/(double)temp);
	    		 Draw.updateReads = true;
	    		 Main.drawCanvas.repaint();	    	
	    	}
    	}
    }		
    
    splitPane.setDividerLocation(Main.chrompan.getComponentCount()*(Main.defaultFontSize+6));
    splitPane.revalidate();
}

static void setAnnotationDrop(String ref) {
	if(Main.drawCanvas != null) {
		geneModel.removeAllElements();
		int maxlength = 0, letterLength = chromosomeDropdown.getFontMetrics(chromosomeDropdown.getFont()).stringWidth("E");
		if(genomehash.get(ref) != null) {
			for(int i = 0; i<genomehash.get(ref).size(); i++) {
				if(genomehash.get(ref).get(i).getName().length() > maxlength) {
					maxlength= genomehash.get(ref).get(i).getName().length();
				}
				geneModel.addElement(genomehash.get(ref).get(i).getName());			
			}
		}
		String addAnno = "Add new annotation...";
		if(addAnno.length() > maxlength) {
			maxlength = addAnno.length();
		}
		geneModel.addElement("Add new annotation...");			
		geneDropdown.setPopupWidth(maxlength*letterLength);
	}
}

static void setChromDrop(String dir) {
	try {
		
		if(!new File(genomeDir.getCanonicalPath()+"/" +dir).exists() || dir.length() == 0) {
			
			/*String[] empty = {""};			
			chromModel = new DefaultComboBoxModel<String>(empty);
			chromosomeDropdown = new SteppedComboBox(chromModel);
			*/
			if(chromModel != null) {
				chromModel.removeAllElements();
				chromosomeDropdown.removeAllItems();
				chromosomeDropdown.revalidate();
				chromosomeDropdown.repaint();
				Main.searchTable.clear();
				if(Main.drawCanvas.splits.size() > 0) {
					Main.drawCanvas.splits.get(0).clearGenes();
				}
				Main.chromDraw.updateExons = true;
				Main.chromDraw.repaint();
			}
			else {
				String[] empty = {""};				
				chromModel = new DefaultComboBoxModel<String>(empty);
				chromosomeDropdown = new SteppedComboBox(chromModel);
				chromosomeDropdown.revalidate();
				chromosomeDropdown.repaint();
			}
		}
		else {
			
			File chromindex = null;
			selectedGenome = dir;
			File defdir = new File(genomeDir.getCanonicalPath() +"/" +dir);
			File[] files = defdir.listFiles();
			String chromtemp = "";
			chromnamevector = Collections.synchronizedList(new ArrayList<String>());
	//		boolean faifound = false;
			for(int i = 0; i< files.length; i++) {
				if(files[i].isDirectory()) {
					continue;
				}
				if(files[i].getName().contains(".gz")) {
					continue;
				}
				else if(files[i].getName().contains(".fai")) {
					chromindex = new File(defdir.getCanonicalPath() +"/" +files[i].getName());
		//			faifound = true;
				}
				else if(files[i].getName().contains(".fa")) {
					chromtemp = defdir.getCanonicalPath() +"/" +files[i].getName();					
				}
			}
			if(chromtemp.equals("")) {			
				String[] empty = {""};
				
				chromModel = new DefaultComboBoxModel<String>(empty);
				chromosomeDropdown = new SteppedComboBox(chromModel);
				chromosomeDropdown.revalidate();
				chromosomeDropdown.repaint();
				return;
			}
		
			if(referenceFile != null) {
				referenceFile.close();
			}
		    referenceFile = new RandomAccessFile(chromtemp, "r");
		   
		    ref = new File(chromtemp);
			
			BufferedReader reader = new BufferedReader(new FileReader(chromindex));
			String line;
			String[] split;
			ChromDraw.chromPos = new HashMap<String, Integer>();
			chromIndex = new Hashtable<String, Long[]>();
			textlength = 0;
			while((line = reader.readLine()) != null) {
				
				split = line.split("\t");
				chromnamevector.add(split[0]);
				Long[] add = {Long.parseLong(split[2]), Long.parseLong(split[1]),  Long.parseLong(split[3])};
				if(split[0].equals("MT")) {
					ChromDraw.chromPos.put("M", Integer.parseInt(split[1]));
					chromIndex.put("M", add);
				}
				chromIndex.put(split[0], add);
				if(split[0].length() > textlength) {
					textlength = split[0].length();
				}
				ChromDraw.chromPos.put(split[0], Integer.parseInt(split[1]));
				
			}
			reader.close();
			MethodLibrary.ChromSorter sorter = new MethodLibrary.ChromSorter();
			Collections.sort(chromnamevector, sorter);
			chromnames = new String[chromnamevector.size()];
			
			if(chromnamevector.get(0).contains("chr")) {
				
				refchrom = "chr";
			}
			else {
				refchrom = "";
			}
			for(int i = 0; i< chromnames.length; i++) {
				chromnames[i] = chromnamevector.get(i).replace(refchrom, "");
				
			}
			
			
			chromModel = new DefaultComboBoxModel<String>(chromnames);
			if(chromosomeDropdown == null) {
				chromosomeDropdown = new SteppedComboBox(chromModel);
			
			}
			else {
				
				chromosomeDropdown.setModel(chromModel);
			}
			clicked = false;
			
			refDropdown.setSelectedItem(dir);
			clicked = true;
			clickedAnno = false;
			setAnnotationDrop(dir);
			
			if(defaultAnnotation.length() == 0) {
				geneDropdown.setSelectedIndex(0);
			}
			else {
				geneDropdown.setSelectedItem(defaultAnnotation);
			}
			clickedAnno = true;				
			
			
			int letterlength = chromosomeDropdown.getFontMetrics(chromosomeDropdown.getFont()).stringWidth("E");
			chromosomeDropdown.setPopupWidth(textlength*letterlength+25);		
			chromosomeDropdown.revalidate();
			chromosomeDropdown.repaint();
			refDropdown.setToolTipText(refDropdown.getSelectedItem().toString());
			geneDropdown.setToolTipText(geneDropdown.getSelectedItem().toString());
		}
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}
public static class OpenProject extends SwingWorker<String, Object> {		
	
	File projectfile;
	
	
	public OpenProject(File file) {
		
		this.projectfile = file;
		
	}
	@SuppressWarnings("unchecked")
	protected String doInBackground() {		
		Main.drawCanvas.loading("Loading project...");
	 	Boolean missing = false;
		
		try {
			FileInputStream fin = new FileInputStream(projectfile);
			ObjectInputStream ois = new ObjectInputStream(fin) {
				
										 
				protected ObjectStreamClass readClassDescriptor() throws IOException, ClassNotFoundException {
				    ObjectStreamClass read = super.readClassDescriptor();
				    
				    if (read.getName() != null && !read.getName().contains("[") && !read.getName().contains(".")) {
				    	 try {
					    	Field f = read.getClass().getDeclaredField("name");
				            f.setAccessible(true);
				            f.set(read, "base.BasePlayer." +read.getName());
				    	 }
				    	 catch(Exception e) {
				    		 System.out.println(read.getName());
				    		 e.printStackTrace();
				    		
				    	 }
				    }
				    return read;
				}
			
				
			};
			
			try {
					drawCanvas.sampleList = (ArrayList<Sample>)ois.readObject();			
					drawCanvas.splits = (ArrayList<SplitClass>) ois.readObject();
					
					for(int i = 0; i<drawCanvas.splits.size(); i++) {
						drawCanvas.splits.get(i).resetSplits();
					}
				
					Main.samples = (short)drawCanvas.sampleList.size();
					drawCanvas.drawVariables = (DrawVariables)ois.readObject();		
					
				for(int i = 0 ; i<drawCanvas.sampleList.size(); i++) {
					
					drawCanvas.sampleList.get(i).resetreadHash();
					if(drawCanvas.sampleList.get(i).getTabixFile() != null) {		
						if(!new File(drawCanvas.sampleList.get(i).getTabixFile()).exists()) {
							
							ErrorLog.addError(drawCanvas.sampleList.get(i).getTabixFile() +" not found.");
							missing = true;
							drawCanvas.removeSample(drawCanvas.sampleList.get(i));									
							i--;
							continue;
						}								
						
						if(drawCanvas.sampleList.get(i).vcfchr == null) {
							drawCanvas.sampleList.get(i).vcfchr = "";
						}
						//if(!drawCanvas.sampleList.get(i).multipart) {
							if(drawCanvas.sampleList.get(i).getVCFInput() == null) {
								drawCanvas.sampleList.get(i).setInputStream();
							}
							FileRead.checkMulti(drawCanvas.sampleList.get(i));
							if(!drawCanvas.sampleList.get(i).multipart && !drawCanvas.sampleList.get(i).multiVCF) {
								Main.varsamples++;
							}
						/*}							
						else {
							
						}*/
					//	Main.varsamples++;
						
					}
					else if(drawCanvas.sampleList.get(i).calledvariants) {
						Main.varsamples++;
					}
					
					if(drawCanvas.sampleList.get(i).samFile != null) {
						readsamples++;
						if(drawCanvas.sampleList.get(i).samFile.getName().endsWith("cram")) {
							drawCanvas.sampleList.get(i).readString = "CRAM";
						}
						else {
							drawCanvas.sampleList.get(i).readString = "BAM";
						}
					}
					else {
						drawCanvas.sampleList.get(i).readString = "No BAM/CRAM";
					}
					if(drawCanvas.sampleList.get(i).longestRead == null) {
						
						drawCanvas.sampleList.get(i).longestRead = 0;
					}
				}
				for(int i = 0; i<drawCanvas.splits.size(); i++) {
					drawCanvas.splits.get(i).setDivider(4.0);
				}									
				
				try {
					readingControls = true;
					Control.controlData = (ControlData)ois.readObject();
					if(Control.controlData.fileArray.size() > 0) {
						
						if(Control.controlData.fileArray.get(0) instanceof ControlFile) {
							if(Control.controlData.fileArray.size() > 0) {								
								  Main.trackPane.setVisible(true);							
								  Main.varpane.setDividerSize(3);	 									
								  Main.varpane.setDividerLocation(0.1);
								  Main.controlScroll.setVisible(true);
								  Main.controlDraw.setVisible(true);	
								  varpane.revalidate();
								
							}
						  for(int i = 0 ; i<Control.controlData.fileArray.size(); i++) {								 
							  
							 if(!new File(Control.controlData.fileArray.get(i).tabixfile).exists()) {
								 ErrorLog.addError(Control.controlData.fileArray.get(i).tabixfile +" not found.");
								 controlDraw.removeControl(i);								
								 i--;
								 missing = true;
								 continue;
							 }
							 controlDraw.trackDivider.add(0.0);
							  Control.controlData.fileArray.get(i).setMenu();
				//			  MethodLibrary.addHeaderColumns(Control.controlData.fileArray.get(i));
						
						  }
					}
					else {									
						  Control.controlData.fileArray = Collections.synchronizedList(new ArrayList<ControlFile>());
						 
					}
				}
					
					readingControls = false;		
		
		//		if(ois.available() > 0) {
					readingbeds = true;
				try {
					
					bedCanvas.bedTrack = (ArrayList<BedTrack>)ois.readObject();
					
				}
				catch(EOFException excep) {
					
				}
				try {
					Settings.settings = (HashMap<String,Integer>)ois.readObject();
					Settings.setValues();
				}
				catch(Exception excep) {
					
				}
				try {
					VariantHandler.variantSettings = (HashMap<String,Integer>)ois.readObject();
					
					VariantHandler.setValues();
				}
				catch(Exception excep) {
					
				}
				
				ois.close();
			//	}
				if(bedCanvas.bedTrack != null && bedCanvas.bedTrack.size() > 0) {					
					  
					boolean first = true;
					for(int i = 0 ; i< bedCanvas.bedTrack.size(); i++) {
					
						if(!bedCanvas.bedTrack.get(i).file.exists()) {
							bedCanvas.bedTrack.remove(i);							
							i--;
							continue;
						}
						if(first ) {
							if(Main.trackPane.isVisible()) {
								
								 Main.varpane.setDividerLocation(Main.varpane.getDividerLocation()+100);
								  Main.trackPane.setDividerLocation(Main.varpane.getDividerLocation()/2);
								 if(Main.controlScroll.isVisible()) {
									 Main.trackPane.setDividerSize(3);
								 }
							}
							else {
								Main.trackPane.setVisible(true);
								Main.varpane.setDividerLocation(0.1);
								Main.varpane.setDividerSize(3);	
							}
							 Main.bedScroll.setVisible(true);
							 Main.bedCanvas.setVisible(true);	
							first = false;
						}
						bedCanvas.bedTrack.get(i).setHead();
						bedCanvas.bedTrack.get(i).setColors();
						bedCanvas.bedTrack.get(i).first = true;
						
						bedCanvas.trackDivider.add(0.0);
						bedCanvas.bedTrack.get(i).setBedLevels();
						bedCanvas.bedTrack.get(i).setmenu();
						FileRead.setTable(bedCanvas.bedTrack.get(i));
						if(bedCanvas.bedTrack.get(i).file == null) {
							//SeekableStream stream = SeekableStreamFactory.getInstance().getStreamFor(bedCanvas.bedTrack.get(i).url);									
							//bedCanvas.bedTrack.get(i).setBBfileReader(new BBFileReader(bedCanvas.bedTrack.get(i).url.toString(), stream, bedCanvas.bedTrack.get(i)));
						}
						else if(bedCanvas.bedTrack.get(i).file.getName().toLowerCase().endsWith("bigwig") || bedCanvas.bedTrack.get(i).file.getName().toLowerCase().endsWith("bw")) {									
							bedCanvas.bedTrack.get(i).setZoomlevel(1);									
							bedCanvas.bedTrack.get(i).setBBfileReader(new BBFileReader(bedCanvas.bedTrack.get(i).file.getCanonicalPath(), bedCanvas.bedTrack.get(i)));
							bedCanvas.bedTrack.get(i).getSelectorButton().setVisible(false);
						}
						else if(bedCanvas.bedTrack.get(i).file.getName().toLowerCase().endsWith(".bedgraph") || bedCanvas.bedTrack.get(i).file.getName().toLowerCase().endsWith(".bedgraph.gz")) {
							 bedCanvas.bedTrack.get(i).setSelector();
							bedCanvas.bedTrack.get(i).getSelectorButton().setVisible(true);
					    	  
					    }
						else {							 
							 bedCanvas.bedTrack.get(i).setSelector();
							 bedCanvas.bedTrack.get(i).getSelectorButton().setVisible(true);
						}
						if((bedCanvas.bedTrack.get(i).file != null &&bedCanvas.bedTrack.get(i).file.length() / 1048576 < Settings.settings.get("bigFile")) || bedCanvas.bedTrack.get(i).getZoomlevel() != null) {
							 bedCanvas.bedTrack.get(i).small = true;					    	  	    	 
					    }	
					    else {
					    	 bedCanvas.bedTrack.get(i).small = false;	
					    	 FileRead.setBedTrack(bedCanvas.bedTrack.get(i));
					    }
						if(bedCanvas.bedTrack.get(i).graph) {
							if(bedCanvas.bedTrack.get(i).getCollapseBox() == null) {
								bedCanvas.bedTrack.get(i).setCollapsebox();
							}
							bedCanvas.bedTrack.get(i).getCollapseBox().setText("Auto scale");
						}
					}					
				}
				readingbeds = false;
				 for(int i = 0 ; i<Control.controlData.fileArray.size(); i++) {				 
					 MethodLibrary.addHeaderColumns(Control.controlData.fileArray.get(i));			
				 }
				 if(Average.frame != null) {
					 if(Average.frame.isVisible()) {
					 	Average.setSamples();
				 	}
				 }
			}
			catch(Exception ex) {
				ex.printStackTrace();
				clearData();
			}
			
			
	    	frame.setTitle("BasePlayer - Project: " +drawCanvas.drawVariables.projectName);
	    	FileRead.checkSamples();
	    	
			Main.drawCanvas.drawVariables.visiblesamples = Main.samples;
			Main.drawCanvas.checkSampleZoom();
			drawCanvas.resizeCanvas(Main.drawCanvas.getWidth(), Main.drawCanvas.getHeight());
			
			Draw.setScrollbar(drawCanvas.drawVariables.scrollbarpos);
			
			for(int i= 0; i<drawCanvas.splits.size(); i++) {
			/*	drawCanvas.splits.get(i).setCytoImage(null);
				chromDraw.drawCyto(drawCanvas.splits.get(i));						
				chromDraw.updateExons = true;
				*/
				FileRead.search = true;						
				drawCanvas.gotoPos(drawCanvas.splits.get(i).chrom, drawCanvas.splits.get(i).start, drawCanvas.splits.get(i).end);		
				
				//chromDraw.repaint();
			}
			if(missing) {
				Main.showError("Missing files. Goto Tools->View log.", "Note");
				
			}
		}
			catch(Exception ex) {
				Main.showError("Sorry, your project must be created again.", "Note");
	    		clearData();
	    		Main.drawCanvas.ready("Loading project...");  
	    		ex.printStackTrace();
			}
		}
		catch(Exception ex) {
			ex.printStackTrace();
			Main.drawCanvas.ready("Loading project...");  
			clearData();
		}
		Main.drawCanvas.ready("Loading project...");  
		return "";
	}
}
void openProject() {
	
	 JFileChooser chooser = new JFileChooser(Main.projectDir);	
	
	 chooser.setAcceptAllFileFilterUsed(false);
  	 MyFilterSES sesFilter = new MyFilterSES();	    	  
  	 chooser.addChoosableFileFilter(sesFilter);
  	 chooser.setDialogTitle("Open project");
  	 chooser.setPreferredSize(new Dimension((int)screenSize.getWidth()/3, (int)screenSize.getHeight()/3));
     int returnVal = chooser.showOpenDialog((Component)this.getParent());	         	         
  
     if (returnVal == JFileChooser.APPROVE_OPTION) {
    	 	clearData();
    	 	projectDir = chooser.getSelectedFile().getParent();	 	
			writeToConfig("DefaultProjectDir=" +projectDir);
			OpenProject opener = new OpenProject(chooser.getSelectedFile());
			opener.execute();
     }
}


static void getBands() {
	try {
		File bandfile = new File(genomeDir.getCanonicalPath() +"/" +selectedGenome +"/bands.txt");
		
		if(bandfile.exists()) {
			ChromDraw.bandVector.clear();	
			BufferedReader in = new BufferedReader(new FileReader(bandfile));
			String line;
			while((line = in.readLine()) !=null) {			
				ChromDraw.bandVector.add(line.split("\\s"));						
			}
			in.close();
		}
		else {
			ChromDraw.bandVector.clear();
		}
		
	}
	catch(Exception e) {		
		
		e.printStackTrace();
	}
}
@Override
/*	public boolean dispatchKeyEvent(KeyEvent e) {
		
		int keyCode = e.getKeyCode();
		
		if(e.getID() == KeyEvent.KEY_PRESSED) {
			
			if((e.getModifiers() & KeyEvent.CTRL_MASK) != 0) {
				Main.drawCanvas.ctrlpressed = 100;
				if(keyCode == KeyEvent.VK_S) {
					if(drawCanvas.drawVariables.projectName.equals("Untitled")) {
						saveProjectAs.doClick();
					}
					else {
						Serializer ser = new Serializer();
						ser.serialize(drawCanvas.drawVariables.projectFile);
					}
				}
			
			}
			else if(keyCode == KeyEvent.VK_7) {
				
			}		
			else if(keyCode == KeyEvent.VK_O && ((e.getModifiers() & KeyEvent.CTRL_MASK) != 0)) {
		
			}
			else if(keyCode == KeyEvent.VK_ENTER) {
				
				if(e.getSource() == searchField) {
					drawCanvas.scrollbar = false;
					searchString = searchField.getText();
					if(searchField.getText().toUpperCase().startsWith("S ") && searchField.getText().length() > 2) {
						
						for(int i = 0; i<Main.drawCanvas.sampleList.size(); i++) {
							
							if(Main.drawCanvas.sampleList.get(i).getName().toUpperCase().contains(searchField.getText().toUpperCase().substring(2))) {
								drawCanvas.drawVariables.visiblestart = (short)i;
								drawCanvas.drawVariables.visibleend = (short)i;
								
								drawCanvas.resizeCanvas(this.getWidth(), (int)(Main.samples*drawCanvas.drawVariables.sampleHeight));
								Draw.setScrollbar((int)(i*drawCanvas.drawVariables.sampleHeight));								
								break;
							}
						}
						return true;
					}
					if(searchField.getText().replace(" ", "").toUpperCase().matches("CHR.{1,2}(?!:)")) {
						
						 if(Main.chromnamevector.contains(searchField.getText().replace(" ", "").toUpperCase().substring(3))) {
							 Main.chromosomeDropdown.setSelectedItem(searchField.getText().toUpperCase().substring(3));							
						 }						
						 return true;
					}
					if(searchField.getText().contains(",")) {
						searchString = searchField.getText().replace(" ", "").replace(",", "");
					}
					if(searchString.contains("chr")) {
						searchString = searchString.replace(" ", "").replace("chr", "");
					}
					
					if(searchTable.containsKey(searchString.toUpperCase())) {
						
						FileRead.search = true;
						String[] result = searchTable.get(searchString.toUpperCase());
						drawCanvas.clearReads();
						FileRead.searchStart = Integer.parseInt(result[1]);
						
						FileRead.searchEnd = Integer.parseInt(result[2]);
						searchChrom = result[0];
						
						searchStart = Integer.parseInt(result[1]);
						
						searchEnd = Integer.parseInt(result[2]);
						
						drawCanvas.gotoPos(result[0], Integer.parseInt(result[1]), Integer.parseInt(result[2]));
					}
					else if(searchString.replace(" ", "").matches("\\d+-?\\d+?")) {
						FileRead.search = true;
						if(searchString.contains("-")) {
							drawCanvas.clearReads();
							drawCanvas.gotoPos(Main.drawCanvas.splits.get(0).chrom,Integer.parseInt(searchString.replace(" ", "").split("-")[0]), Integer.parseInt(searchString.replace(" ", "").split("-")[1]));		
						}
						else {
							drawCanvas.clearReads();
							drawCanvas.gotoPos(Main.drawCanvas.splits.get(0).chrom,Integer.parseInt(searchString.replace(" ", ""))-200, Integer.parseInt(searchString.replace(" ", ""))+200);		
						}	
					
					}
					
					else if(searchString.replace(" ", "").matches(".+:\\d+-?\\d+?")) {
						String[] result = searchString.replace(" ", "").split(":");
						FileRead.search = true;
						if(result[1].contains("-")) {
							drawCanvas.clearReads();
							drawCanvas.gotoPos(result[0].replace(" ", ""), Integer.parseInt(result[1].split("-")[0]), Integer.parseInt(result[1].split("-")[1]));
						}
						else {
							drawCanvas.clearReads();
							drawCanvas.gotoPos(result[0].replace(" ", ""), Integer.parseInt(result[1])-200, Integer.parseInt(result[1])+200);
						}
					}
					
				}
			
			}
		}
		
		return false;
	}
	@Override*/
public void mouseDragged(MouseEvent event) {
	
	if(event.getSource() == drawScroll.getVerticalScrollBar()) {
		glassPane.setVisible(true);
		drawCanvas.scrolldrag = true;
/*		drawCanvas.drawVariables.visiblestart = (short)(Main.drawScroll.getVerticalScrollBar().getValue()/drawCanvas.drawVariables.sampleHeight);		
		//	if(drawCanvas.drawVariables.visiblestart + drawCanvas.drawVariables.visiblesamples < Main.samples-1) {
		
			drawCanvas.scrolldrag = true;
			
			if(drawCanvas.splits.size() > 1) {
				for(int i = 0; i<drawCanvas.splits.size(); i++) {
					drawCanvas.splits.get(i).updateReads = true;					
				}
			}
			else {
				Draw.updateReads = true;
				Draw.updatevars = true;
			}
			Draw.updatevars = true;
			
//		}
 */
		Draw.updatevars = true;
		//drawCanvas.repaint();
	// TODO Auto-generated method stub
	}
	
}
@Override
public void mouseMoved(MouseEvent event) {
	
	// TODO Auto-generated method stub
	
}
@Override
public void keyTyped(KeyEvent e) {
	// TODO Auto-generated method stub
	
}
static boolean zoomtopos(String chrom, String pos, String sample) {			
	
	if(sample.length() > 10) {
		sample = sample.toUpperCase().substring(0, 10);
	}
	else {
		sample = sample.toUpperCase();
	}
	boolean found = false;
	for(int i = 0; i<Main.drawCanvas.sampleList.size(); i++) {
	
		if(Main.drawCanvas.sampleList.get(i).getName().toUpperCase().contains(sample)) {
			drawCanvas.drawVariables.visiblestart = (short)i;
			drawCanvas.drawVariables.visiblesamples = (short)1;
			
			drawCanvas.resizeCanvas(Main.drawCanvas.getWidth(), (int)(Main.samples*drawCanvas.drawVariables.sampleHeight));
			Draw.setScrollbar((int)(i*drawCanvas.drawVariables.sampleHeight));								
			found = true;
			break;
		}
	}	
	if(!found) {
		return false;
	}
	Main.nothread = true;
	   Main.noreadthread = true;
	   FileRead.search = true;	
	Main.drawCanvas.gotoPos(chrom, Integer.parseInt(pos)-100, Integer.parseInt(pos)+100);	
	return true;
}

public class Seqfetcher extends SwingWorker<String, Object> {
	File file, outfile;
	BedTrack track;
	
	public Seqfetcher(File file, File outfile) {
		this.file = file;
		this.outfile = outfile;
		Settings.softClips = 1;
	}

	void fetchSeq() {
		BufferedReader reader = null;
		BufferedWriter writer = null;
		try {
			reader = new BufferedReader(new FileReader(file));
			writer = new BufferedWriter(new FileWriter(outfile));
			String line, chrom, position, sample;
			String[] splitter, positionsplit;
			FileRead readreader = new FileRead();
			/*while(!reader.readLine().split("\\s+")[0].equals("19:52486706")) {
									
			}*/
			boolean firsterror = true;
			StringBuffer nonfounds = new StringBuffer("");
			while((line = reader.readLine()) != null) {
				if(!Main.drawCanvas.loading) {
					writer.close();
					break;
				}
				splitter = line.split("\\s+");
				
				positionsplit = splitter[0].split(":");
				chrom = positionsplit[0];
				position = positionsplit[1];
				sample = splitter[1];			
				VariantHandler.hideIndels.setSelected(true);
				VariantHandler.hideSNVs.setSelected(true);		
				
				if(!zoomtopos(chrom, position, sample)) {
					if(firsterror) {
						Main.showError("Sample: " +sample +" not found in opened samples.\nCheck error log for more missing files.", "Error");
						ErrorLog.addError(sample +" not found.");
						nonfounds.append(sample);
						firsterror = false;
					}
					else {
						if(!nonfounds.toString().contains(sample)) {
							nonfounds.append(sample);
							ErrorLog.addError(sample +" not found.");
						}							
					}
					continue;
				}
				
				
				SplitClass split = drawCanvas.splits.get(0);
				Sample readsample = drawCanvas.sampleList.get(drawCanvas.drawVariables.visiblestart);
				ReadNode read;
				
				int centerpos = Integer.parseInt(position);
				
				int minpos = Integer.MAX_VALUE, maxpos = 0;
				ArrayList<Object[]> readlist = new ArrayList<Object[]>();
				readreader.splitIndex = split;
				if(readsample.getreadHash().get(split) == null) {
					readsample.resetreadHash();
				}
				readreader.getReads(chrom, centerpos-100, centerpos+100, readsample.getreadHash().get(split));
			
				if(readsample.getreadHash().get(split) == null) {
					readsample.resetreadHash();
				}
				for(int i = 0; i<readsample.getreadHash().get(split).getReads().size();i++)  {
					
					read = readsample.getreadHash().get(split).getReads().get(i);
					do {
						if(read.getPosition() < centerpos && read.getPosition()+read.getWidth() > centerpos) {
							if(read.getMismatches() != null && read.getMismatches().size() > 10) {
								if(read.getPosition() < minpos) {
									minpos = read.getPosition();
								}
								if(read.getPosition()+read.getWidth() > maxpos) {
									maxpos = read.getPosition()+read.getWidth();
								}
								SAMRecord readsam = Main.fileReader.getRead(chromosomeDropdown.getSelectedItem().toString(),read.getPosition(), read.getPosition()+read.getWidth(),read.getName(), readsample.getreadHash().get(split));
								try {
									Object[] adder = {read.getPosition(), readsam.getReadString()};
									readlist.add(adder);
								}
								catch(Exception e) {
									continue;
								}
								
							}
						}
					
					}
					while ((read = read.getNext()) != null);				
				}
				if(maxpos-minpos < 10) {
					String error =">" +readsample.getName() +"|BP="+chromosomeDropdown.getSelectedItem().toString() +":" +MethodLibrary.formatNumber(centerpos) +" (Breakpoint mismatches not found))\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
					
					writer.write(error+"\n");
					continue;
				}
				int[][] matrix = new int[5][maxpos-minpos];
				
				for(int i = 0 ; i<5;i++) {
					for(int j = 0;j<maxpos-minpos;j++) {
						matrix[i][j] = 0;
					}
				}
				for(int i = 0 ; i<readlist.size();i++) {
					for(int j =0; j<readlist.get(i)[1].toString().length(); j++) {		
						if((int)readlist.get(i)[0]-minpos+j >= matrix[0].length) {
							break;
						}
						matrix[baseMap.get((byte)(readlist.get(i)[1].toString().charAt(j)))-1][(int)readlist.get(i)[0]-minpos+j]++;
					}
				}
				StringBuffer buffer = new StringBuffer("");//readsample.getName() +"\nPosition: " +chromosomeDropdown.getSelectedItem().toString() +":" +MethodLibrary.formatNumber(minpos) +"\nBreak point: " +chromosomeDropdown.getSelectedItem().toString() +":" +MethodLibrary.formatNumber(centerpos) +"\n" );
		//		System.out.print("A [ ");
				read = null;
			/*	buffer.append("A [ ");
				for(int i = 0 ; i<4;i++) {
					for(int j =0; j<maxpos-minpos; j++) {		
						if(j==(centerpos-minpos)+1) {
			//				System.out.print("| ");
							buffer.append("| ");
						}
		//				System.out.print(matrix[i][j]+" ");
						buffer.append(matrix[i][j]+" ");
					}
		//			System.out.println(" ]");
					buffer.append("]\n");
					if(i == 0) {
		//				System.out.print("C [ ");
						buffer.append("C [ ");
					}
					else if(i == 1) {
		//				System.out.print("G [ ");
						buffer.append("G [ ");
					}
					else if(i == 2) {
		//				System.out.print("T [ ");
						buffer.append("T [ ");
					}
				}
				*/
			//	System.out.println(buffer.toString());
			//	buffer.append("\n");
				buffer.append(">" +readsample.getName() +"|BP="+chromosomeDropdown.getSelectedItem().toString() +":" +MethodLibrary.formatNumber(centerpos) +" (LeftPosition=" +chromosomeDropdown.getSelectedItem().toString() +":" +MethodLibrary.formatNumber(minpos) +")\n");
				int max=0,maxindex=0;
				String[] bases = {"A","C","G","T"};
				StringBuffer fasta = new StringBuffer("");
				for(int j =0; j<maxpos-minpos; j++) {
					max = 0;
					maxindex = 0;
					for(int i = 0; i<4; i++) {
						if(matrix[i][j] > max ) {
							max = matrix[i][j];
							maxindex = i;
						}
						else if(max > 0 && matrix[i][j] == max) {
							max = 0;
							maxindex = -1;
						}
					}
					if(maxindex > -1) {
						fasta.append(bases[maxindex]);
					}
					else {
						fasta.append("N");
					}
				}
				buffer.append(fasta.toString());
			//	System.out.println(buffer.toString());
				writer.write(buffer.toString() +"\n");
				
		}
			reader.close();
			writer.close();
			Main.nothread = false;
			Main.noreadthread = false;
			FileRead.search = false;	
			Draw.variantcalculator = false;
			//VariantHandler.hideIndels.setSelected(false);
			//VariantHandler.hideSNVs.setSelected(false);
			chromDraw.updateExons = true;
			chromDraw.repaint();
			
			Main.showError("Fasta file ready!", "Note");
		}
		catch(Exception e) {
			try {
				if(reader != null && writer != null) {
					reader.close();
					writer.close();
				}
			}
			catch(Exception ex) {
				
			}
			Main.nothread = false;
			Main.noreadthread = false;
			FileRead.search = false;	
			chromDraw.updateExons = true;
			chromDraw.repaint();
			e.printStackTrace();
		}
	}
	protected String doInBackground() {
		Main.drawCanvas.loading("Writing fasta");
		fetchSeq();
		Main.drawCanvas.ready("Writing fasta");
		return "";
	}
	
}

	public static void getConsSeq( ) {
		
				
		Sample readsample = Main.drawCanvas.sampleList.get(Main.drawCanvas.drawVariables.visiblestart);
		SplitClass split = Main.drawCanvas.splits.get(0);
		ReadNode read;
		int centerpos = chromDraw.getPosition((int)(Main.drawCanvas.getDrawWidth()/2.0 +split.pixel/2), split);
		int minpos = Integer.MAX_VALUE, maxpos = 0;
		ArrayList<Object[]> readlist = new ArrayList<Object[]>();
		
		for(int i = 0; i<readsample.getreadHash().get(split).getReads().size();i++)  {
			
			read = readsample.getreadHash().get(split).getReads().get(i);
			do {
				if(read.getPosition() < centerpos && read.getPosition()+read.getWidth() > centerpos) {
					if(read.getMismatches() != null && read.getMismatches().size() > 10) {
						if(read.getPosition() < minpos) {
							minpos = read.getPosition();
						}
						if(read.getPosition()+read.getWidth() > maxpos) {
							maxpos = read.getPosition()+read.getWidth();
						}
						SAMRecord readsam = Main.fileReader.getRead(chromosomeDropdown.getSelectedItem().toString(),read.getPosition(), read.getPosition()+read.getWidth(),read.getName(), readsample.getreadHash().get(split));
						
						Object[] adder = {read.getPosition(), readsam.getReadString()};
						readlist.add(adder);
					}
				}			
			}
			while ((read = read.getNext()) != null);				
		}
		
		int[][] matrix = new int[5][maxpos-minpos];
		
		for(int i = 0 ; i<5;i++) {
			for(int j = 0;j<maxpos-minpos;j++) {
				matrix[i][j] = 0;
			}
		}
		for(int i = 0 ; i<readlist.size();i++) {
			for(int j =0; j<readlist.get(i)[1].toString().length(); j++) {		
				
				matrix[baseMap.get((byte)(readlist.get(i)[1].toString().charAt(j)))-1][(int)readlist.get(i)[0]-minpos+j]++;
			}
		}
		StringBuffer buffer = new StringBuffer("");//readsample.getName() +"\nPosition: " +chromosomeDropdown.getSelectedItem().toString() +":" +MethodLibrary.formatNumber(minpos) +"\nBreak point: " +chromosomeDropdown.getSelectedItem().toString() +":" +MethodLibrary.formatNumber(centerpos) +"\n" );

		read = null;
	
		buffer.append(">" +readsample.getName() +"|BP="+chromosomeDropdown.getSelectedItem().toString() +":" +centerpos  +" (LeftPosition=" +chromosomeDropdown.getSelectedItem().toString() +":" +MethodLibrary.formatNumber(minpos) +")\n");
		int max=0,maxindex=0;
		String[] bases = {"A","C","G","T"};
		StringBuffer fasta = new StringBuffer("");
		for(int j =0; j<maxpos-minpos; j++) {
			max = 0;
			maxindex = 0;
			for(int i = 0; i<4; i++) {
				if(matrix[i][j] > max ) {
					max = matrix[i][j];
					maxindex = i;
				}
				else if(max > 0 && matrix[i][j] == max) {
					max = 0;
					maxindex = -1;
				}
			}
			if(maxindex > -1) {
				fasta.append(bases[maxindex]);
			}
			else {
				fasta.append("N");
			}
		}
		buffer.append(fasta.toString());
		
		System.out.println(buffer.toString());
	}
	

@Override

public void keyPressed(KeyEvent e) {
	keyCode = e.getKeyCode();
	
	
	if(!Main.shift && (e.getModifiers() & KeyEvent.SHIFT_MASK) != 0) {
		
		Main.shift = true;
	}
	else if((e.getModifiers() & KeyEvent.CTRL_MASK) != 0) {
		Main.drawCanvas.ctrlpressed = 100;
	
		if(keyCode == KeyEvent.VK_S) {
			
			if(drawCanvas.drawVariables.projectName.equals("Untitled")) {
				saveProjectAs.doClick();
			}
			else {
				Serializer ser = new Serializer();
				ser.serialize(drawCanvas.drawVariables.projectFile);
			}
		}
		if(keyCode == KeyEvent.VK_W) {
			
			System.out.println("\n----\n");
			
			/*
				int[][] array = new int[Main.varsamples][VariantHandler.callSlider.getUpperValue()+1];
				
					for(int i = 0; i<array.length; i++) {
						for(int j = 0; j<array[i].length; j++) {
							array[i][j] = 0;
						}
					}
					/*VarNode node = FileRead.head.getNext();
					int counter = 0;
					while(node != null) {
						if(!drawCanvas.hideNode(node)) {					
							counter++;
							for(int i = 0; i<node.vars.size(); i++) {
								if(!drawCanvas.hideNodeVar(node, node.vars.get(i))) {
									for(int j = 0;j<node.vars.get(i).getValue().size(); j++) {
										if(!drawCanvas.hideVar(node.vars.get(i).getValue().get(j), false)) {
											array[node.vars.get(i).getValue().get(j).getSample().getIndex()][(int)(MethodLibrary.round(node.vars.get(i).getValue().get(j).getAlleleFraction()*100,2))]++;
										}
										
									}
									
									
								}						
							}
						}
				
						node = node.getNext();
					}*/
				int width = Main.drawCanvas.getWidth()-Main.sidebarWidth;
			
				JPopupMenu menu = new JPopupMenu();
			Plotter plotter = new Plotter(width);			
			plotter.setPreferredSize(new Dimension(width,400));
			menu.add(plotter);
			menu.pack();
			menu.show(Main.drawCanvas,Main.sidebarWidth, drawScroll.getVerticalScrollBar().getValue());
		}
		if(keyCode == KeyEvent.VK_PLUS || keyCode == 107) {
			
	//		defaultFontSize++;
			
	//		setFonts();
		}
		if(keyCode == KeyEvent.VK_M || keyCode == KeyEvent.VK_MINUS || keyCode == 109) {
			
	//		defaultFontSize--;
			
	//		setFonts();
		}
	/*	else if(keyCode == KeyEvent.VK_M) {
			VarNode current = FileRead.head.getNext();
			 while(current != null) {
				 
				 if(current.getNext() != null && current.getNext().getPosition()-current.getPosition() > 100000) {
					 System.out.println(current.getPosition() +" " +(current.getNext().getPosition()-current.getPosition()));
				 }
				 current = current.getNext();
			 }
			 current = null;
		}
		*/
	}
	else if(keyCode == KeyEvent.VK_DELETE) {
		if(Main.drawCanvas.selectedSample != null) {
			Main.drawCanvas.removeSample(Main.drawCanvas.selectedSample);
		}
	}
/*	else if(keyCode == KeyEvent.VK_7) {
		
		
	}*/		
	else if(keyCode == KeyEvent.VK_O && ((e.getModifiers() & KeyEvent.CTRL_MASK) != 0)) {
	/*	clearData();
		
		try {
			FileInputStream fin = new FileInputStream("C:/HY-Data/RKATAINE/test.ses");
			ObjectInputStream ois = new ObjectInputStream(fin);
			drawCanvas.sampleList = (ArrayList<Sample>) ois.readObject();
			Main.samples = (short)drawCanvas.sampleList.size();
			
			drawCanvas.splits = (ArrayList<SplitClass>) ois.readObject();
			for(int i = 0; i<drawCanvas.splits.size(); i++) {
				drawCanvas.splits.get(i).resetSplits();
			}
			drawCanvas.drawVariables = (DrawVariables)ois.readObject();
			
			ois.close();
			drawCanvas.resizeCanvas(Main.drawCanvas.getWidth(), Main.drawCanvas.getHeight());
			
			Draw.setScrollbar(drawCanvas.drawVariables.scrollbarpos);
			
			for(int i= 0; i<drawCanvas.splits.size(); i++) {
				drawCanvas.splits.get(i).setCytoImage(null);
				chromDraw.drawCyto(drawCanvas.splits.get(i));
				chromDraw.updateExons = true;
				FileRead.search = true;
				drawCanvas.gotoPos(drawCanvas.splits.get(i).chrom, drawCanvas.splits.get(i).start, drawCanvas.splits.get(i).end);
				chromDraw.repaint();
			}
			for(int i = 0 ;i<drawCanvas.sampleList.size(); i++) {
				for(int r = 0 ;r<drawCanvas.sampleList.size(); r++) {
					if(drawCanvas.sampleList.get(i).getreadHash().size() > 0) {
						drawCanvas.sampleList.get(i).resetreadHash();
					}
				}
			}
		
		}
		catch(Exception ex) {
			ex.printStackTrace();
		}
		*/
	}

	else if(keyCode == KeyEvent.VK_F9) {
		
	
		FileRead.head.putNext(null);
		drawCanvas.variantsStart = 0;
		drawCanvas.variantsEnd = 1;
		Draw.updatevars = true;
		Main.drawCanvas.repaint();

	}
	else if(keyCode == KeyEvent.VK_F11) {
		
	/*	try {
			BBFileReader reader = new BBFileReader( Main.bedCanvas.bedTrack.get(0).file.getCanonicalPath(),  Main.bedCanvas.bedTrack.get(0));
			int zoomlevel = 1;
			for(int i =2;i<reader.getZoomLevels().getZoomHeaderCount();i++) {
				if(reader.getZoomLevels().getZoomLevelHeader(i).getReductionLevel() < (Main.drawCanvas.splits.get(0).viewLength/(Main.drawCanvas.splits.get(0).pixel*Main.drawCanvas.splits.get(0).viewLength))) {
					zoomlevel = i;
				}
				else {
					break;
				}
			}
			
	
			   
		}
		catch(Exception ex) {
			ex.printStackTrace();
		}*/
	}
	else if(keyCode == KeyEvent.VK_F12) {
		/*VarNode next = Main.drawCanvas.current.getNext();
		
		MethodLibrary.makeMultiAlt("2",next.getPosition(), "G", next);
		
			next = null;
		*/
	/*	try {
			URL urli = new URL("ftp://ftp.ensembl.org/pub/grch37/update/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz");
			System.out.println(urli.getProtocol() +" " +urli.getHost() +" " +urli.getPath().substring(0,urli.getPath().lastIndexOf("/")+1));
		}
		catch(Exception ex) {
			ex.printStackTrace();
		}*/
	}
	else if(keyCode == KeyEvent.VK_F8) {
		
		File file = new File(searchField.getText().replaceAll(" ", ""));
		if(!file.exists()) {
			if(Main.drawCanvas.splits.get(0).viewLength < 1000) {
				
				Main.getConsSeq();
			}
		}
		else {
			try {	  	    
	    		   
		    	  JFileChooser chooser = new JFileChooser(file.getPath());
		    	  chooser.setAcceptAllFileFilterUsed(true);			    	  	  
		    	  
		    	  chooser.setDialogTitle("Save FASTA file as...");
		          int returnVal = chooser.showSaveDialog((Component)this.getParent());	           	  
			       
			      if(returnVal == JFileChooser.APPROVE_OPTION) {  
			    	  String suffix = "";
			    	  if(!chooser.getSelectedFile().getName().endsWith(".fa") && !chooser.getSelectedFile().getName().endsWith(".fasta")) {
			    		  suffix = ".fa";
			    	  }
			    	  
				       File outfile = new File(chooser.getSelectedFile().getCanonicalPath() +suffix);	
				     //  File outfile = new File("test.fa");	
				       Main.nothread = true;
					   Main.noreadthread = true;
					   FileRead.search = true;	
					   Draw.variantcalculator = true;
					
					   Seqfetcher fetcher = new Seqfetcher(file, outfile);
					   fetcher.execute();
				       
			      }
				}
				catch(Exception ex) {
					ex.printStackTrace();
				}
			
		}
		
	}
	else if(keyCode == KeyEvent.VK_ENTER) {
		
		if(e.getSource() == searchField) {
			if(searchField.getText().equals("tati") || searchField.getText().equals("third")) {
				ReadNode read;
				HashMap<String, String[]> chrs = new HashMap<String, String[]>();
				HashMap<String, Integer> temp = new HashMap<String, Integer>();
				for(int j = 0; j<drawCanvas.sampleList.get(0).getreadHash().get(drawCanvas.splits.get(0)).getReads().size();j++) {
					read = drawCanvas.sampleList.get(0).getreadHash().get(drawCanvas.splits.get(0)).getReads().get(j);
					
					while(read != null) {	
						
						if (read.SA != null) {
							String[] SAs = read.SA.split(";");
							temp.clear();
							for(int i = 0; i<SAs.length; i++) {
								String[] sa = SAs[i].split(",");
								if(temp.containsKey(sa[0])) {
									continue;
								}
								else {
									temp.put(sa[0], 1);
								}
								if(!chrs.containsKey(sa[0])) {
									String[] add = {"1", sa[1]};
									chrs.put(sa[0], add);
								}
								else {
									String[] add = chrs.get(sa[0]);
									String[] newString = {""+(Integer.parseInt(add[0])+1), add[1] };
									chrs.put(sa[0], newString); 
								}									
							}								
						}
						read = read.getNext();						
					}												
				}
				Iterator<Map.Entry<String, String[]>> it = chrs.entrySet().iterator();
				String result = "Results for splitted reads:\n\n";
			    while (it.hasNext()) {
			        Map.Entry<String, String[]> pair = it.next();
			       
			        result += "chr" +pair.getKey() + ":" +MethodLibrary.formatNumber(Integer.parseInt(pair.getValue()[1]))  +" = " + pair.getValue()[0] +"\n";
			        it.remove(); // avoids a ConcurrentModificationException
			    }
			    ErrorLog.addError(result);
			    ErrorLog.frame.setLocation(frame.getLocationOnScreen().x+10, frame.getLocationOnScreen().y+10);
				
				ErrorLog.frame.setState(JFrame.NORMAL);
				ErrorLog.frame.setVisible(true);
			 //   JOptionPane.showMessageDialog(Main.chromDraw, result, "Tati's results", JOptionPane.INFORMATION_MESSAGE);
				read = null;
				return;
			}
			drawCanvas.scrollbar = false;
		//	searchString = searchField.getText();
				if(searchField.getText().toUpperCase().startsWith("S ") && searchField.getText().length() > 2) {
					
					for(int i = 0; i<Main.drawCanvas.sampleList.size(); i++) {
						
						if(Main.drawCanvas.sampleList.get(i).getName().toUpperCase().contains(searchField.getText().toUpperCase().substring(2))) {
							drawCanvas.drawVariables.visiblestart = (short)i;
							drawCanvas.drawVariables.visiblesamples = (short)1;
							
							drawCanvas.resizeCanvas(this.getWidth(), (int)(Main.samples*drawCanvas.drawVariables.sampleHeight));
							Draw.setScrollbar((int)(i*drawCanvas.drawVariables.sampleHeight));								
							break;
						}
					}
					return;
				}
			//	drawCanvas.clearReads();
				if(searchField.getText().contains(";")) {
					searchList = searchField.getText().split(";");
					
					for(int i = 0 ; i<searchList.length; i++) {
						
						returnlist = parseSearch(searchList[i]);
						
						
						if(returnlist != null) {
							
							if(i == 0) {
								FileRead.search = true;
								drawCanvas.gotoPos(returnlist[0], Integer.parseInt(returnlist[1]), Integer.parseInt(returnlist[2]));
								FileRead.search = false;
							}
							else {
								
								drawCanvas.addSplit(returnlist[0], Integer.parseInt(returnlist[1]), Integer.parseInt(returnlist[2]));
							}
						}
					}
				}
				else {
					try {
					returnlist = parseSearch(searchField.getText());
					if(returnlist != null) {
						if(undoPointer < undoList.size()-1) {
							undoList.add(undoPointer+1,searchField.getText());
							if(undoPointer < undoList.size()-1) {
								for(int i = undoPointer+2 ; i< undoList.size(); i++) {
									undoList.remove(i);
									i--;
								}
							}
						}
						else {
							undoList.add(searchField.getText());
						}
						undoPointer = undoList.size()-1;
						
							forward.setEnabled(false);
							
						if(undoPointer > 0) {
							back.setEnabled(true);
						}
						FileRead.search = true;
						drawCanvas.gotoPos(returnlist[0], Integer.parseInt(returnlist[1]), Integer.parseInt(returnlist[2]));
					}
					}
					catch(Exception ex) {
						ex.printStackTrace();
					}
				}
				
				
			}
	
	}
	
}
@Override
public void keyReleased(KeyEvent e) {
	Main.drawCanvas.ctrlpressed = 5;
	Main.shift = false;
	
}
static Double calcAffiniyChange(VarNode node, String alt, BedNode bednode) {
			if(alt == null || alt.length() > 1) {
				return 0.0;
			}			
			if(bednode.getTrack().selex) {
				
				int index = node.getPosition()-bednode.getPosition();
				
				if(index < 0) {
					
					return 0.0;
				}
				int[][] matrix=null;
				if(bednode.forward) {
					matrix = Main.SELEXhash.get(bednode.id);
				}
				else {
					matrix = MethodLibrary.reverseMatrix(Main.SELEXhash.get(bednode.id));
				}					
				
				double sum;										
				Double value;					
				double mutatedvalue;
				sum = matrix[0][index] + matrix[1][index] + matrix[2][index] + matrix[3][index];
			//	System.out.println(matrix[0][index] +" " +matrix[1][index] +" " +matrix[2][index] +" " +matrix[3][index]);
				value = matrix[baseMap.get(node.getRefBase())-1][index]/(double)sum;
				
				mutatedvalue = matrix[baseMap.get((byte)alt.charAt(0))-1][index]/(double)sum;
				
				value = value*Math.log(value/background.get(node.getRefBase()));
				if(mutatedvalue != 0) {
					mutatedvalue = mutatedvalue*Math.log(mutatedvalue/background.get((byte)alt.charAt(0)));
				}
				
				
				
				return mutatedvalue-value;
				//System.out.println(mutatedvalue-value);
									
				}
				else {
					return 0.0;
				}
		
		}
String[] parseSearch(String searchstring) {
		if(searchstring.replace(" ", "").toUpperCase().matches("CHR.{1,2}(?!:)")) {
		
		 if(Main.chromnamevector.contains(searchstring.replace(" ", "").toUpperCase().substring(3))) {
			 Main.chromosomeDropdown.setSelectedItem(searchstring.toUpperCase().substring(3));							
		 }						
		return null;
	}
	if(searchstring.contains(",")) {
		searchstring = searchstring.replace(" ", "").replace(",", "");
	}
	if(searchstring.contains("chr")) {
		searchstring = searchstring.replace(" ", "").replace("chr", "");
	}
	
	if(searchTable.containsKey(searchstring.replace(" ", "").toUpperCase()) || geneIDMap.containsKey(searchstring.replace(" ", "").toUpperCase())) {
		
		FileRead.search = true;
		String[] result  = {};
		if(searchTable.containsKey(searchstring.replace(" ", "").toUpperCase())) {
			result = searchTable.get(searchstring.replace(" ", "").toUpperCase());
		}			
		else {
			result = searchTable.get(geneIDMap.get(searchstring.replace(" ", "").toUpperCase()));
		}
	//	drawCanvas.clearReads();
		FileRead.searchStart = Integer.parseInt(result[1]);
		
		FileRead.searchEnd = Integer.parseInt(result[2]);
		searchChrom = result[0];
		
		searchStart = Integer.parseInt(result[1]);
		
		searchEnd = Integer.parseInt(result[2]);
		
	//	drawCanvas.gotoPos(result[0], Integer.parseInt(result[1]), Integer.parseInt(result[2]));
		String[] returnstring = {result[0],result[1],result[2]};
		return returnstring;
	}
	else if(searchstring.replace(" ", "").matches("\\d+-?\\d+?")) {
		
		if(searchstring.contains("-")) {
		//	drawCanvas.clearReads();
		//	drawCanvas.gotoPos(Main.drawCanvas.splits.get(0).chrom,Integer.parseInt(searchString.replace(" ", "").split("-")[0]), Integer.parseInt(searchString.replace(" ", "").split("-")[1]));		
			String[] returnstring = {Main.drawCanvas.splits.get(0).chrom,searchstring.replace(" ", "").split("-")[0],searchstring.replace(" ", "").split("-")[1]};
			return returnstring;
		}
		else {
		//	drawCanvas.clearReads();
			
	//		drawCanvas.gotoPos(Main.drawCanvas.splits.get(0).chrom,Integer.parseInt(searchString.replace(" ", ""))-200, Integer.parseInt(searchString.replace(" ", ""))+200);	
			String[] returnstring = {Main.drawCanvas.splits.get(0).chrom,""+(Integer.parseInt(searchstring.replace(" ", ""))-200),""+(Integer.parseInt(searchstring.replace(" ", ""))+200)};
			return returnstring;
		}	
		
		//chromDraw.updateExons = true;
		//chromDraw.repaint();
	}
	
	else if(searchstring.replace(" ", "").matches(".+:\\d+-?\\d+?")) {
		String[] result = searchstring.replace(" ", "").split(":");
	
		if(result[1].contains("-")) {
	//		drawCanvas.clearReads();
		//	drawCanvas.gotoPos(result[0].replace(" ", ""), Integer.parseInt(result[1].split("-")[0]), Integer.parseInt(result[1].split("-")[1]));
			
			searchChrom = result[0].replace(" ", "");
			String[] returnstring = {result[0].replace(" ", ""),result[1].split("-")[0],result[1].split("-")[1]};
			return returnstring;
		}
		else {
	//		drawCanvas.clearReads();
	//		drawCanvas.gotoPos(result[0].replace(" ", ""), Integer.parseInt(result[1])-200, Integer.parseInt(result[1])+200);
			searchChrom = result[0].replace(" ", "");
			String[] returnstring = {result[0].replace(" ", ""),""+(Integer.parseInt(result[1])-200),""+(Integer.parseInt(result[1])+200)};
			return returnstring;
		}
	}
	return null;
	
}


}





