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
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.tribble.readers.TabixReader;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;
import java.awt.Color;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
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

import javax.swing.BorderFactory;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSplitPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JTextPane;
import javax.swing.SwingWorker;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.plaf.basic.BasicSplitPaneDivider;
import javax.swing.plaf.basic.BasicSplitPaneUI;

	public class Main extends JPanel implements ActionListener, ChangeListener, ComponentListener, MouseListener, PropertyChangeListener, KeyListener, MouseMotionListener {
		private static final long serialVersionUID = 1L;
		
	    static JFrame frame = new JFrame("BasePlayer");  
		//UI
	    static String version = "1.0.0";
	    static int sidebarWidth = 200;
	    static String[] argsit = {}, args;
	    static GraphicsDevice gd = GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice(); 
	    
	    //testqeqwqwr
	    static int width = gd.getDisplayMode().getWidth();
	    static int height = gd.getDisplayMode().getHeight();
	    static int loadTextWidth = 200;
	    static Dimension screenSize = new Dimension(width, height);
	    static int drawWidth = (int)screenSize.getWidth()-200, drawHeight = (int)screenSize.getHeight()-100, chromHeight = 150, bedHeight = 200;
	    static GradientPaint gradient = new GradientPaint(drawWidth/2-loadTextWidth/2,0,Color.red,drawWidth/2+loadTextWidth,height,Color.green,true);
	    private String[] searchList;
	    static Short samples = 0, varsamples = 0;
	    static JLabel widthLabel = new JLabel("");
	    public static String gerp;
	    static int searchStart=-1, searchEnd=-1;
		//Labels	
	    static boolean configChanged = false;
	    static int trackdivider = 0;
	    static java.util.List<String> chromnamevector = Collections.synchronizedList(new ArrayList<String>());
	    static Image A, C, G, T;
	    static HashMap<String, Integer> baseMap = new HashMap<String, Integer>();
	    Dimension chromDimensions, bedDimensions, drawDimensions, buttonDimension;
	    static final int buttonHeight = 20, buttonWidth = 60;
	    static JTextField searchField;
	    static JTextField positionField = new JTextField();
	    static FileRead fileReader = new FileRead();
	    static String refchrom = "", selectedGenome;
	    static int selectedChrom = 0;
	    static Loader loading; 
	    static Hashtable<String, Long[]> chromIndex = new Hashtable<String, Long[]>();
	    static Hashtable<String, ArrayList<File>> genomehash = new Hashtable<String, ArrayList<File>>();
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
		//Buttons etc.
	    static String defaultGenome = "";
	    static JSplitPane splitPane, trackPane, varpane, drawpane;
	    public static boolean nothread;	
	    static Hashtable<String, String[]> searchTable = new Hashtable<String, String[]>();
	    
	    JButton zoomout = new JButton("Zoom out");
	    JButton dosomething = new JButton("Do stuff!");
	    JButton back = new JButton("<"), forward = new JButton(">");
	    static double[][] snow = new double[200][4];
	    static String[] snowtable = {"A", "C", "G", "T" };
	    static String userDir;
	    VarNode prevTemp;
		static JTextArea info = new JTextArea();
//		static JLabel fileLabel = new JLabel("Open VCF-files");
		static JPanel panel = new JPanel(new GridBagLayout());
		static JMenuBar menubar = new JMenuBar();
		JMenu filemenu = new JMenu("File");
		JMenu toolmenu = new JMenu("Tools");
		JMenu help = new JMenu("Help");
		JMenuItem manage = new JMenuItem("Variant handler");
		JMenuItem average = new JMenuItem("Coverage handler");
		JMenuItem settings = new JMenuItem("Settings");
		JMenuItem update = new JMenuItem("Update");
		JMenuItem errorlog = new JMenuItem("View log");
		JLabel helpLabel = new JLabel("This is pre-release version of BasePlayer\nContact: help@baseplayer.fi\n\nUniversity of Helsinki");
		JMenu genome = new JMenu("Change genome");
		static boolean updatelauncher = false;
		static JMenuItem opensamples = new JMenuItem("Open samples");
		static JMenuItem addtracks = new JMenuItem("Add tracks");
		static JMenuItem addcontrols = new JMenuItem("Add controls");
		static JMenuItem saveProject = new JMenuItem("Save project");
		static JMenuItem saveProjectAs = new JMenuItem("Save project as...");
		static JMenuItem openProject = new JMenuItem("Open project");
		JMenuItem clear = new JMenuItem("Clear data");
		JMenuItem clearMemory = new JMenuItem("Clean memory");
		static RandomAccessFile referenceFile;    
	    static JComboBox<String> chromosomeDropdown;		
		//UI		
	    
	    DefaultComboBoxModel<String> chromModel;
	    static Hashtable<String, int[][]> SELEXhash = new Hashtable<String, int[][]>();
		static JScrollPane drawScroll;
		static JScrollPane chromScroll;
		static JScrollPane bedScroll;
		static JScrollPane controlScroll;
		
		ActionListener ChromoDropActionListener = new ActionListener() {
			
	    	public void actionPerformed(ActionEvent actionEvent) {
	    		
	    	  try {
	    		
	    	  if (actionEvent.getActionCommand() == "comboBoxChanged") {
	    		 
	    		  chromosomeDropdown.validate();
	    		  chromosomeDropdown.revalidate();
	    		  chromosomeDropdown.repaint();	    
	    		  drawCanvas.splits.get(0).getReadBuffer().setComposite( drawCanvas.composite);					
	  			  drawCanvas.splits.get(0).getReadBuffer().fillRect(0,0, drawCanvas.splits.get(0).getReadImage().getWidth(),Main.drawScroll.getViewport().getHeight());	
	  			  drawCanvas.splits.get(0).getReadBuffer().setComposite(drawCanvas.splits.get(0).getBackupr());		
	  			  drawCanvas.rbuf.setComposite( drawCanvas.composite);				
	  			  drawCanvas.rbuf.fillRect(0,0, drawCanvas.getWidth(),Main.drawScroll.getViewport().getHeight());	
	  			  drawCanvas.rbuf.setComposite(drawCanvas.backupr);		
	    		  drawCanvas.clearReads();
	    		  selectedChrom = chromosomeDropdown.getSelectedIndex();    		  
	    		  drawCanvas.chrom = chromosomeDropdown.getSelectedItem().toString();
	    		  chromDraw.cytoImage = null;	    		 
	    		  drawCanvas.splits.get(0).setCytoImage(null);
	    		  drawCanvas.splits.get(0).chrom = drawCanvas.chrom;
	    		  drawCanvas.splits.get(0).transStart = 0;
	    		  drawCanvas.splits.get(0).chromEnd = chromIndex.get(Main.refchrom + Main.chromosomeDropdown.getSelectedItem().toString())[1].intValue();
	    		  drawCanvas.splits.get(0).readSeqStart = 0;
	    		  drawCanvas.splits.get(0).readSequence = null;					
	    		  FileRead filereader = new FileRead();
	    		  filereader.chrom = Main.chromosomeDropdown.getSelectedItem().toString();	    		  
	    		  drawCanvas.setStartEnd(1.0,(double)drawCanvas.splits.get(0).chromEnd);
	    		  
	    		  if(!FileRead.search) {
	    			  FileRead.searchStart = 1;
	    			  FileRead.searchEnd = drawCanvas.splits.get(0).chromEnd;
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
	    		 
	        //	  drawCanvas.repaint();
	    		 
	    		 
	    	  }   	  			
    			}
    			catch(Exception e) {
    			
    				e.printStackTrace();
    			}		    			
	      	}
	    };

//		private String searchString;
		private Dimension fieldDimension;
		private File[] annofiles;

		private int keyCode;

		private String[] returnlist;

		private JTextField chooserTextField;

		private Image iconImage;

		static String controlDir ="";
		static String trackDir = "";
		static String projectDir = "";
		static BasicSplitPaneDivider splitPaneDivider;

		static BasicSplitPaneDivider varPaneDivider;

		static BasicSplitPaneDivider trackPaneDivider;

		static String annotationfile = "";

		static String searchChrom = "";
		static File ref;  
		static JPanel glassPane = new JPanel()
		    {
		      
			private static final long serialVersionUID = 1L;	
			
			
			public void paintComponent(Graphics g) {
				paintComponent((Graphics2D)g);
			}
			public void paintComponent(Graphics2D g)
		       {
				
				if(drawCanvas.loading) {
					
				 if(drawCanvas.loadingtext.length() > 0 && !drawCanvas.loadingtext.equals("note")) {
					 
					  
					  g.setFont(Draw.loadingFont);
			          g.setColor(Draw.zoomColor);			
			          if(Main.loadTextWidth != (int)(g.getFontMetrics().getStringBounds(drawCanvas.loadingtext, g).getWidth())) {
			        	  Main.loadTextWidth = (int)(g.getFontMetrics().getStringBounds(drawCanvas.loadingtext, g).getWidth());
			        	  gradient = new GradientPaint(drawScroll.getWidth()/2-loadTextWidth/2-drawScroll.getWidth()/6,0,Color.red,drawScroll.getWidth()/2+loadTextWidth,0,Color.green,true);			        	 
			          }
			         
					  g.fillRect(drawScroll.getWidth()/2-loadTextWidth/2-5, Main.drawScroll.getViewport().getHeight()*2/3-35, loadTextWidth+10, 80);
					 
					  
					  
					  g.setPaint(gradient);
					  
					  g.fillRect(drawScroll.getWidth()/2-loadTextWidth/2, Main.drawScroll.getViewport().getHeight()*2/3+10, (int)(loadTextWidth*(drawCanvas.loadBarSample/100.0)), 10);
					  g.fillRect(drawScroll.getWidth()/2-loadTextWidth/2, Main.drawScroll.getViewport().getHeight()*2/3+25, (int)(loadTextWidth*(drawCanvas.loadbarAll/100.0)), 10);
					  g.setColor(Color.black);
					  g.drawRect(drawScroll.getWidth()/2-loadTextWidth/2, Main.drawScroll.getViewport().getHeight()*2/3+10, (int)(loadTextWidth*(drawCanvas.loadBarSample/100.0)), 10);						 
					  g.drawRect(drawScroll.getWidth()/2-loadTextWidth/2, Main.drawScroll.getViewport().getHeight()*2/3+25, (int)(loadTextWidth*(drawCanvas.loadbarAll/100.0)), 10);
					  g.drawRect(drawScroll.getWidth()/2-loadTextWidth/2-5, Main.drawScroll.getViewport().getHeight()*2/3-35, loadTextWidth+10, 80);						 
					  g.drawString(drawCanvas.loadingtext, drawScroll.getWidth()/2-loadTextWidth/2, Main.drawScroll.getViewport().getHeight()*2/3);
					  g.drawRect(drawScroll.getWidth()/2-40, Main.drawScroll.getViewport().getHeight()*2/3+50, 80, 30);					 
					  g.setColor(Color.lightGray);
					  g.fillRect(drawScroll.getWidth()/2-40, Main.drawScroll.getViewport().getHeight()*2/3+50, 80, 30);
					
					  
					  if(cancelhover) {

						if(getCursor().getType() != Cursor.HAND_CURSOR) {
							setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
						}
						  g.setColor(Color.white);
					  }
					  else {
						  if(getCursor().getType() != Cursor.DEFAULT_CURSOR) {
								setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
							}
						  g.setColor(Color.black);
					  }
					  g.setFont(ChromDraw.seqFont);
					  g.drawString("Cancel", drawScroll.getWidth()/2-20, Main.drawScroll.getViewport().getHeight()*2/3+70);
					  g.setColor(ChromDraw.backTransparent);
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
						 }
				 }
				 else {
					 g.setColor(ChromDraw.backTransparent);
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
					}
				 }					 
				} 
		      }					       
		    };

		private static String chooserText;

		static Double vardivider = 0.0;

		static int annotation = 0;

		public static String defaultAnnotation = "";	
		
		
	public Main() {	
		
	super(new GridBagLayout());	
	try {
		htsjdk.samtools.util.Log.setGlobalLogLevel(htsjdk.samtools.util.Log.LogLevel.ERROR);
		for(int i=0;i<snow.length; i++) {
			snow[i][0] = (height*Math.random());
			snow[i][1] = (4*Math.random() +1);
			snow[i][2] = (12*Math.random() -6);
			snow[i][3] = (2*Math.random() +1);
		}
		frame.addWindowListener(new java.awt.event.WindowAdapter() {
		    @Override
		    public void windowClosing(java.awt.event.WindowEvent windowEvent) {
		        /*if (JOptionPane.showConfirmDialog(frame, "Are you sure to close this window?", "Really Closing?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE) == JOptionPane.YES_OPTION){
		            System.exit(0);
		        }*/
		    	if(configChanged) {
		    		
		    		try {
		    		 BufferedWriter fileWriter = new BufferedWriter(new FileWriter(Launcher.maindir +"/config.txt"));	
		    		 for(int i = 0 ; i<Launcher.config.size(); i++) {
		    			 fileWriter.write(Launcher.config.get(i) +"\n");
		    		 }
		    		 fileWriter.close();
		    		 
		    		}
		    		catch(Exception e) {
		    			e.printStackTrace();
		    		}
		    	}
		    	
		    }
		});
		baseMap.put("A", 1);
		baseMap.put("C", 2);
		baseMap.put("G", 3);
		baseMap.put("T", 4);
		baseMap.put("N", 5);
		baseMap.put("I", 6);
		baseMap.put("D", 7);
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
		final Image glass = Toolkit.getDefaultToolkit().getImage(getClass().getResource("glass.jpg"));
	    searchField = new JTextField("Search by position or gene") {
	        protected void paintComponent(Graphics g) {
	            super.paintComponent(g);
	           
	            g.drawImage(glass, 2, 2, fieldDimension.height-4,fieldDimension.height-4, this);
	        }
	    };
	    Average.main(argsit);
	    Average.frame.setVisible(false);
		Launcher.fromMain = true;
		Launcher.main(args);
		Settings.main(args);
		path = Launcher.defaultDir;		
		gerp = Launcher.gerpfile;
		defaultGenome = Launcher.defaultGenome;
		defaultAnnotation = Launcher.defaultAnnotation;
		
		annotationfile = defaultAnnotation;
		controlDir = Launcher.ctrldir;
		trackDir = Launcher.trackDir;
		projectDir = Launcher.projectDir;
		 drawDimensions = new Dimension(drawWidth,drawHeight-200);
		 drawCanvas = new Draw((int)drawDimensions.getWidth(),(int)drawDimensions.getHeight());
		 iconImage = Toolkit.getDefaultToolkit().getImage(getClass().getResource("icon.png"));
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
			
			BufferedReader selexReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("SELEX/TFbinding_PFMs.txt")));
			String line, factor;
			String[] split, matrix, values;
			int[][] selexmatrix;
			while((line = selexReader.readLine()) != null) {
				split = line.split("\\t");
				factor = split[0];
				
				matrix = split[1].split(";");	
				
				values = matrix[0].split(",");
				selexmatrix = new int[4][values.length];
				for(int i = 0; i<4; i++ ) {				
					values = matrix[i].split(",");
					for(int j = 0; j< values.length; j++) {						
						selexmatrix[i][j] =  Integer.parseInt(values[j]);
										
					}
				}	
				Main.SELEXhash.put(factor, selexmatrix);
			}				
			selexReader.close();
			selexReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("SELEX/oldSELEX.txt")));
			
			while((line = selexReader.readLine()) != null) {
				split = line.split("\\t");
				factor = split[0];
				
				matrix = split[1].split(";");	
				
				values = matrix[0].split(",");
				selexmatrix = new int[4][values.length];
				for(int i = 0; i<4; i++ ) {				
					values = matrix[i].split(",");
					for(int j = 0; j< values.length; j++) {						
						selexmatrix[i][j] =  Integer.parseInt(values[j]);
										
					}
				}	
				Main.SELEXhash.put(factor, selexmatrix);
			}				
			selexReader.close();
			
			
			
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
		
		userDir = new File(Main.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParent().replace("%20", " ");		
		panel.setBackground(Draw.sidecolor);
		searchField.addKeyListener(this);
		
		frame.addKeyListener(this);
		frame.getContentPane().setBackground(Color.black);
		chromosomeDropdown = new JComboBox<String>();
		glassPane.addMouseListener(this);
		glassPane.addMouseMotionListener(new MouseMotionListener() {
			
			@Override
			public void mouseDragged(MouseEvent arg0) {
				
				
			}

			@Override
			public void mouseMoved(MouseEvent event) {			
				
				if(drawCanvas.loading && event.getX() > Main.drawScroll.getWidth()/2-40 && event.getX()  < Main.drawScroll.getWidth()/2+40 && event.getY() > Main.drawScroll.getViewport().getHeight()*2/3 +40 && event.getY() < Main.drawScroll.getViewport().getHeight()*2/3 +70) {
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
		/*
		baseTable.add("A");
		baseTable.add("C");
		baseTable.add("G");
		baseTable.add("T");
		baseTable.add("N");
		baseTable.add("delA");
		baseTable.add("delC");
		baseTable.add("delG");
		baseTable.add("delT");
		baseTable.add("insA");
		baseTable.add("insC");
		baseTable.add("insG");
		baseTable.add("insT");
		*/
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
		
	    VariantHandler.main(argsit);
	    VariantHandler.frame.setVisible(false);
	  
	//	path = "C:/HY-Data/RKATAINE/FilterSomatic/";
	//	path = "W:/RikuRator-Samples/Thyroid/11-0176_Kilpi1_5_exome";
	//	path = "W:/RikuRator-Samples/Myomas/Tumors/";
	//	path = "V:/cg8/projects/CRC/c9_5119_N_LP6005134-DNA_A01/wgspipeline/align";
	//	path = "V:/cg8/projects/joint-calling/";
	    drawCanvas.loading("note");
		try {
			File genomedir = new File(userDir +"/genomes/"), annodir;
			File[] genomes = genomedir.listFiles(), annotations;
		if(genomes != null) {	
			for(int i = 0; i<genomes.length; i++) {
				annodir = new File(genomes[i].getAbsolutePath() +"/annotation/");
				annotations = annodir.listFiles();
				genomehash.put(genomes[i].getName(), new ArrayList<File>());
			//	System.out.println(genomes[i].getName());
				JMenu addMenu = new JMenu(genomes[i].getName());
				addMenu.addMouseListener(this);
				addMenu.setName(genomes[i].getName());
				addMenu.add(new JLabel("  Select annotation: "));
				addMenu.add(new JSeparator());
				genome.add(addMenu);
				if(annotations != null) {
					for(int j = 0; j<annotations.length; j++) {
						annofiles = annotations[j].listFiles();
						for(int f = 0; f<annofiles.length; f++) {
							if(annofiles[f].getName().endsWith(".bed.gz")) {
								genomehash.get(genomes[i].getName()).add(annofiles[f].getAbsoluteFile());	
								JMenuItem additem = new JMenuItem(annofiles[f].getName().substring(0,annofiles[f].getName().indexOf(".bed.gz")));
								additem.setName(annofiles[f].getName().substring(0,annofiles[f].getName().indexOf(".bed.gz")));
								additem.addMouseListener(this);
								addMenu.add(additem);
								break;
							}
							else if (annofiles[f].getName().contains(".gff")) {
								if(!new File(annofiles[f].getCanonicalPath().substring(0,annofiles[f].getCanonicalPath().indexOf(".gff")) +".bed.gz").isFile()) {
								//	drawCanvas.loading("Converting " +annofiles[f].getName());
									FileRead.readGFF(annofiles[f].getCanonicalFile(), annofiles[f].getCanonicalPath().substring(0,annofiles[f].getCanonicalPath().indexOf(".gff")) +".bed.gz");
							//	    drawCanvas.ready("Converting " +annofiles[f].getName());
									genomehash.get(genomes[i].getName()).add(new File(annofiles[f].getCanonicalPath().substring(0,annofiles[f].getCanonicalPath().indexOf(".gff")) +".bed.gz"));	
									JMenuItem additem = new JMenuItem(annofiles[f].getName().substring(0,annofiles[f].getName().indexOf(".gff")));
									additem.setName(annofiles[f].getName().substring(0,annofiles[f].getName().indexOf(".gff")));
									additem.addMouseListener(this);
									addMenu.add(additem);
									
									break;
								}
								else {
									genomehash.get(genomes[i].getName()).add(new File(annofiles[f].getCanonicalPath().substring(0,annofiles[f].getCanonicalPath().indexOf(".gff")) +".bed.gz"));	
									JMenuItem additem = new JMenuItem(annofiles[f].getName().substring(0,annofiles[f].getName().indexOf(".gff")));
									additem.setName(annofiles[f].getName().substring(0,annofiles[f].getName().indexOf(".gff")));
									additem.addMouseListener(this);
									addMenu.add(additem);
									break;
								}
							}
						}					
					}	
				}
			}
			if(!genomehash.containsKey(defaultGenome)) {
				
				setChromDrop(genomes[0].getName());
				defaultGenome = genomes[0].getName();
			}
			else {
				
				setChromDrop(defaultGenome);
			}
		}
			getBands();			
		    getExons();			
				
		    Draw.image=Toolkit.getDefaultToolkit().getImage(getClass().getResource("background.jpg"));
	
			setButtons();
			drawCanvas.addKeyListener(this);
			bedCanvas.addKeyListener(this);
			CheckUpdates check = new CheckUpdates();
			check.execute();
		
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	catch(Exception ex) {
		ex.printStackTrace();
		JOptionPane.showMessageDialog(null, ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
	}
		}
	void getExons() {
		try {
		String s;
		String[] exonSplit;
		Boolean found = false;
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
		
		while((s = ChromDraw.exonReader.readLine()) != null) {
			exonSplit = s.split("\t");
			
			if(exonSplit[9].equals("1")) {		
				
				String[] adder = {exonSplit[0], exonSplit[1], exonSplit[2]};
				searchTable.put(exonSplit[3].toUpperCase(), adder);					
			}						
		}	
		ChromDraw.exonReader.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	void setMenuBar() {
		filemenu.addMouseListener(this);
		toolmenu.addMouseListener(this);
		
		help.addMouseListener(this);
		filemenu.setPreferredSize(buttonDimension);
		toolmenu.setPreferredSize(buttonDimension);
		toolmenu.setMinimumSize(buttonDimension);
		filemenu.setMinimumSize(buttonDimension);
		filemenu.add(opensamples);
		filemenu.add(addtracks);
		filemenu.add(addcontrols);
		filemenu.add(new JSeparator());
		filemenu.add(openProject);
		filemenu.add(saveProject);
		filemenu.add(saveProjectAs);		
		filemenu.add(new JSeparator());
		filemenu.add(genome);
		filemenu.add(update);
		filemenu.add(clear);
		
		menubar.add(filemenu);
		manage.addActionListener(this);
		update.addActionListener(this);
		average.addActionListener(this);
		toolmenu.add(manage);
		toolmenu.add(average);
		toolmenu.add(settings);
		settings.addActionListener(this);
		clearMemory.addActionListener(this);
		errorlog.addActionListener(this);
		toolmenu.add(clearMemory);
		toolmenu.add(errorlog);
		menubar.add(toolmenu);
		help.setPreferredSize(buttonDimension);
		help.setMinimumSize(buttonDimension);
	//	help.addActionListener(this);
		
		JMenu about = new JMenu("About");
	//	JMenuItem infotable = new JMenuItem(" <html> Line1 <br/> Line2 <br/> Line3 </html> ");
		
	//	JLabel aboutText = new JLabel();
		//aboutText.setText (" <html> Line1 <br/> Line2 <br/> Line3 </html> ");
		//infotable.add(aboutText);
		//about.add(helpLabel);
		JTextPane area = new JTextPane();
		area.setContentType("text/html");
		area.setEditable(false);
	//	area.setEditorKit(javax.swing.JEditorPane.createEditorKitForContentType("text/<b>html</b>"));
		String infotext = "<html><h2>BasePlayer</h2>This is pre-release version of BasePlayer<br/> Author: Riku Katainen <br/> University of Helsinki<br/>"
						+"Tumor Genomics Group (<a href=http://research.med.helsinki.fi/gsb/aaltonen/>http://research.med.helsinki.fi/gsb/aaltonen/</a>) <br/> " 
						+"Contact: help@baseplayer.fi <br/><br/>"
						+"Package includes GRCh37 / HG19 reference genome and Ensembl genes release 78<br/><br/>"
						+"For optimal usage, you should have vcf and bam -files for each sample. <br/> "
						+"e.g. in case you have a sample named as sample1, name all files similarly and <br/>"
						+"place in the same folder:<br/>"
						+"sample1.vcf.gz<br/>"
						+"sample1.vcf.gz.tbi<br/>"
						+"sample1.bam<br/>"
						+"sample1.bam.bai<br/><br/>"
						+"When you open sample1.vcf.gz, sample1.bam is recognized and opened<br/>" 
						+"on the same track.<br/><br/>"
						+"Instructional video will be available in Youtube soon..."; 
		
		area.setText(infotext);
		about.add(area);
		about.addMouseListener(this);
		help.add(about);
		
		menubar.add(help);
	}
	void setButtons() {
		try {
		GridBagConstraints c = new GridBagConstraints();	
		buttonDimension = new Dimension(buttonWidth, buttonHeight);
		fieldDimension = new Dimension(300, buttonHeight);
		setMenuBar();
		c.anchor = GridBagConstraints.NORTHWEST;
		
		c.insets = new Insets(5, 5, 2, 5);
		c.gridx = 0;
		c.gridy = 0;
		c.gridwidth = 1;  
		opensamples.setMargin(new Insets(0, 0, 0, 0));
		
		opensamples.setPreferredSize(buttonDimension);
		addtracks.setMargin(new Insets(0, 0, 0, 0));
		addtracks.setPreferredSize(buttonDimension);
		addcontrols.setMargin(new Insets(0, 0, 0, 0));
		addcontrols.setPreferredSize(buttonDimension);
		panel.add(menubar, c);
		c.gridx = 1;
		
	 //   c.gridy = 1;
		c.gridx = 2;
		
		c.gridx = 3;
		
		
		zoomout.setPreferredSize(new Dimension(zoomout.getFontMetrics(zoomout.getFont()).stringWidth("Zoom out")+8,(int)buttonDimension.getHeight() ));
		zoomout.setMinimumSize(new Dimension(zoomout.getFontMetrics(zoomout.getFont()).stringWidth("Zoom out")+8 ,(int)buttonDimension.getHeight() ));
		zoomout.setMargin(new Insets(0, 0, 0, 0));
	//	zoomout.setBackground(Draw.sidecolor);
		panel.add(zoomout, c);		
		
	   // c.gridx = 4;
	   // panel.add(slider, c);
	  //  slideLabel.setPreferredSize(new Dimension(200, 15));
	  //  c.gridx = 5;
	//    panel.add(slideLabel, c);
	    c.gridx = 4;
	    chromosomeDropdown.setPreferredSize(buttonDimension);
	    panel.add(chromosomeDropdown, c);
	    c.gridx = 5;
	   
	    
	 //   searchField.setText("Search by position or gene");
	    searchField.setMargin(new Insets(0,fieldDimension.height+2, 0, 0));
	    searchField.setPreferredSize(fieldDimension);
	//    searchField.setMinimumSize(fieldDimension);
	    searchField.setMinimumSize(new Dimension(100, (int)fieldDimension.getHeight()));
	    searchField.addMouseListener(this);
	    panel.add(searchField, c);
	    searchField.setForeground(Color.gray);
	    c.gridx = 6;
	    back.addMouseListener(this);
	    forward.addMouseListener(this);
	    back.setEnabled(false);
	    forward.setEnabled(false);
	    back.setMargin(new Insets(0, 0, 0, 0));
	    forward.setMargin(new Insets(0, 0, 0, 0));
	    back.setPreferredSize(new Dimension(back.getFontMetrics(back.getFont()).stringWidth("<")+10,(int)fieldDimension.getHeight()));
	    forward.setPreferredSize(new Dimension(forward.getFontMetrics(forward.getFont()).stringWidth("<")+10,(int)fieldDimension.getHeight()));
	    back.setMinimumSize(new Dimension(back.getFontMetrics(back.getFont()).stringWidth("<")+10,(int)fieldDimension.getHeight()));
	    forward.setMinimumSize(new Dimension(forward.getFontMetrics(forward.getFont()).stringWidth("<")+10,(int)fieldDimension.getHeight()));
	    panel.add(back, c);
	    c.gridx = 7;
	    panel.add(forward, c);
	   
	    c.gridx = 8;
	    positionField.setPreferredSize(fieldDimension);
	    positionField.setMinimumSize(new Dimension(100, (int)fieldDimension.getHeight()));
	    positionField.setEditable(false);
	    panel.add(positionField, c);
	    
	    c.gridx = 9;
	    widthLabel.setPreferredSize(new Dimension(widthLabel.getFontMetrics(widthLabel.getFont()).stringWidth("000,000,000bp")+2,(int)fieldDimension.getHeight()));
	    widthLabel.setMinimumSize(new Dimension(widthLabel.getFontMetrics(widthLabel.getFont()).stringWidth("000,000,000bp")+2,(int)fieldDimension.getHeight()));
	    panel.add(widthLabel, c);
	  
	  //  chromosomeDropdown.setPreferredSize(new Dimension(buttonSize, buttonHeight));
	    chromosomeDropdown.setMaximumRowCount(25);	 
	   
		chromosomeDropdown.setEnabled(true);
		chromosomeDropdown.addActionListener(ChromoDropActionListener);	
		chromosomeDropdown.addMouseListener(this);
		
		c.gridwidth = 10;  
	    c.gridx = 0;
	    c.gridy = 1;
	    
	    bedDimensions = new Dimension(drawWidth, bedHeight);
	   
	    chromDimensions = new Dimension(drawWidth,chromHeight);
	    
	    bedScroll = new JScrollPane(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
	    bedScroll.getViewport().setPreferredSize(bedDimensions);
	    
	    drawScroll = new JScrollPane(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);	  	   
	    drawScroll.getViewport().setPreferredSize(drawDimensions);	   	   
	    
	    chromScroll = new JScrollPane(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
	    chromScroll.getViewport().setPreferredSize(chromDimensions);
	   
		drawScroll.getVerticalScrollBar().setAutoscrolls(false);
//	    chromScroll.getViewport().setPreferredSize(new Dimension(drawWidth,chromHeight-20));
	   
	   
	//    drawScroll.setBorder(BorderFactory.createEmptyBorder());
	    chromScroll.setBorder(BorderFactory.createLineBorder(Color.black, 2));
	    drawScroll.setBorder(BorderFactory.createLineBorder(Color.black, 2));
	 //   chromScroll.setBorder(BorderFactory.createEmptyBorder());
	    
	    bedScroll.setBorder(BorderFactory.createEmptyBorder());
	    controlDraw = new ControlCanvas((int)bedDimensions.getWidth(),(int)bedDimensions.getHeight());
	    controlScroll = new JScrollPane(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);	  	   
	    controlScroll.getViewport().setPreferredSize(bedDimensions);	 
	    controlScroll.getViewport().add(controlDraw);
	    controlDraw.setVisible(false);
	    controlScroll.setVisible(false);
	    controlScroll.setBorder(BorderFactory.createEmptyBorder());
	    
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
	    
	    drawCanvas.resizeCanvas(drawWidth,drawHeight);
	   
	//    drawCanvas.splits.get(1).pixel = (drawCanvas.getDrawWidth())/(double)(drawCanvas.splits.get(1).viewLength);
	 //   drawCanvas.splits.get(2).pixel = (drawCanvas.getDrawWidth())/(double)(drawCanvas.splits.get(2).viewLength);
	    
	    chromDraw.setPreferredSize(chromDimensions);	
	    chromScroll.addMouseListener(this);
	    chromDraw.addMouseListener(this);
	    chromScroll.getViewport().add(chromDraw);	    
	    drawScroll.getViewport().add(drawCanvas);	
	    drawScroll.addMouseListener(this);
	    drawCanvas.addMouseListener(this);
	    bedCanvas = new BedCanvas(drawWidth, 200);
	   
	    bedScroll.getViewport().add(bedCanvas);
	   
	//    frame.setExtendedState(frame.getExtendedState() | JFrame.MAXIMIZED_BOTH);
	   
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
	                g.setColor(Color.black);
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
	                g.setColor(Color.black);
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
		varpane.addMouseListener(this);
		bedScroll.setVisible(false);
		
		trackPane.addMouseListener(this);
		
		
		controlScroll.setVisible(false);
		splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, chromScroll, varpane);
		
		splitPane.setUI(new BasicSplitPaneUI() {
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
	    });
	    splitPane.setBorder(null);
	    
	    BasicSplitPaneUI basicSplitPaneUI = (BasicSplitPaneUI) splitPane.getUI();
	    splitPaneDivider = basicSplitPaneUI.getDivider();
	    splitPaneDivider.addMouseListener(this);
	    basicSplitPaneUI = (BasicSplitPaneUI) trackPane.getUI();
	    trackPaneDivider = basicSplitPaneUI.getDivider();
	    trackPaneDivider.addMouseListener(this);
	    BasicSplitPaneUI splitPaneUI = (BasicSplitPaneUI) varpane.getUI();
	    varPaneDivider = splitPaneUI.getDivider();
	    varPaneDivider.addMouseListener(this);
		splitPane.setDividerSize(3);
	/*	BasicSplitPaneDivider divider = (BasicSplitPaneDivider) splitPane.getComponent(2);
		divider.setBackground(Color.black);
		*/
		splitPane.setPreferredSize(drawDimensions);
		
		splitPane.setResizeWeight(0.2);
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
				glassPane.setVisible(true);
				if ( e.getWheelRotation() < 0 ) {
					if(drawCanvas.drawVariables.visiblestart > 0) {
						drawCanvas.drawVariables.visiblestart--; 
					//	drawCanvas.drawVariables.visibleend--;	
						
					}
					Draw.setScrollbar((int)(drawCanvas.drawVariables.visiblestart*drawCanvas.drawVariables.sampleHeight));
				}
				else {
					if(drawCanvas.drawVariables.visiblestart+drawCanvas.drawVariables.visiblesamples < Main.samples) {
						drawCanvas.drawVariables.visiblestart++; 
				//		drawCanvas.drawVariables.visibleend++;
					
					}
					Draw.setScrollbar((int)(drawCanvas.drawVariables.visiblestart*drawCanvas.drawVariables.sampleHeight));
				}
				
			}
	    	
	    }
	    		
	    		
	    );
	  
	/*    chromScroll.getVerticalScrollBar().addAdjustmentListener(new AdjustmentListener() {
	    	
	    	@Override
			public void adjustmentValueChanged(AdjustmentEvent event) {
				
			//		chromDraw.repaint();
				
			}    	
	    	
	    });
	  */  
	    zoomout.addActionListener(this);
	     
	    FileRead.head = new VarNode(0, "N", "N", (short)0, (short)0, false,(short)0,null, null, null, null);     
	    drawCanvas.current = FileRead.head;
	  
	//  splitPane.addComponentListener(this);
	    splitPane.addPropertyChangeListener(this);
	    trackPane.addPropertyChangeListener(this);
	    varpane.addPropertyChangeListener(this);
	  //  drawScroll.getVerticalScrollBar().addAdjustmentListener(this);
//	    frame.addPropertyChangeListener(this);
	    frame.addComponentListener(this);
	    frame.addMouseListener(this);	    	   
	    frame.setGlassPane(glassPane);
	    if(chromosomeDropdown.getItemCount() > 0) {
	    	chromosomeDropdown.setSelectedIndex(0);
	    	drawCanvas.chrom = chromosomeDropdown.getItemAt(0);
	    }
	    glassPane.setOpaque(false);
	    glassPane.setVisible(false);
	   
	    positionField.setText("chr1:1-" +MethodLibrary.formatNumber(drawCanvas.splits.get(0).chromEnd));
	    widthLabel.setText(MethodLibrary.formatNumber(drawCanvas.splits.get(0).chromEnd) +"bp");
	    
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
		    	if(chrom == null) {
		    		split.chromEnd =1;
		    	}
		    	else {
		    	split.chromEnd = chromIndex.get(Main.refchrom +chrom)[1].intValue();
		    	}
		    }
		    catch(Exception e) {
		    	System.out.println(chrom);
		    	e.printStackTrace();
		    }
		    split.start = 1;
		    split.end =  split.chromEnd;
		    split.viewLength =split.end-split.start;
		    drawCanvas.splits.add(split);	
		    Main.drawCanvas.resizeCanvas(Main.drawCanvas.getWidth(), Main.drawCanvas.getHeight());
			
			
			for(int i= 0; i<drawCanvas.splits.size(); i++) {
				drawCanvas.splits.get(i).setCytoImage(null);
				chromDraw.drawCyto(drawCanvas.splits.get(i));
				chromDraw.updateExons = true;
				chromDraw.repaint();
			}
			//Main.drawCanvas.repaint();
	}
	
	
	
	static class MyFilterVCF extends javax.swing.filechooser.FileFilter {
		public boolean accept(File file) { 
				if (file.isDirectory()) {
					return true;
			    }	
				if (file.getName().contains(".vcf.gz") && !file.getName().contains(".tbi")) {
					return true;
				}					
				
				else {
					return false;
				}	
		} 
		
		public String getDescription() { return "*.vcf.gz"; }
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
				if (file.getName().endsWith(".bed.gz")) {
					return true;
				}					
				if (file.getName().endsWith(".bedgraph.gz")) {
					return true;
				}	
				if(file.getName().endsWith(".gff.gz")) {
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
				else {
					return false;
				}	
		} 
		
		public String getDescription() { return "*.bed.gz, *.gff.gz, *.bedgraph.gz, *.bigWig, *.bigBed"; }
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
		
		
		drawCanvas.setStartEnd(1.0,(double)drawCanvas.splits.get(0).chromEnd);
		
		if(samples > 0) {
			
			if(drawCanvas.splits.get(0).chromEnd > Settings.readDrawDistance) {
				
				Main.drawCanvas.clearReads();	
				
				
			}
			drawCanvas.removeSplits();
			chromDraw.varnode = null;
			chromDraw.vardraw = null;
	//		drawCanvas.drawVariables.visiblestart = 0;
	//		drawCanvas.drawVariables.visiblesamples = (short)(Main.samples);
	//		Main.drawCanvas.checkSampleZoom();
		//	drawCanvas.drawVariables.visibleend = (short)(Main.samples-1);
			/*if(Main.samples*drawCanvas.sampleHeight > drawScroll.getViewport().getHeight()) {
				drawCanvas.drawVariables.visibleend = (short)(drawScroll.getViewport().getHeight()/drawCanvas.sampleHeight);			}
			else {
				drawCanvas.drawVariables.visibleend = (short)(samples-1);
			}
			
			*/
			VariantHandler.table.hoverNode = null;
			VariantHandler.table.selectedNode = null;
		//	Main.drawCanvas.get(i).resizeCanvas(drawCanvas.get(i).getWidth(), drawScroll.getViewport().getHeight());
			Draw.updatevars = true;
			
		//	drawCanvas.eraseReads();
	//		drawCanvas.eraseReads = true;
			drawCanvas.splits.get(0).getReadBuffer().setComposite( drawCanvas.composite);					
			drawCanvas.splits.get(0).getReadBuffer().fillRect(0,0, drawCanvas.splits.get(0).getReadImage().getWidth(),Main.drawScroll.getViewport().getHeight());	
			drawCanvas.splits.get(0).getReadBuffer().setComposite(drawCanvas.splits.get(0).getBackupr());		
			drawCanvas.rbuf.setComposite( drawCanvas.composite);				
			drawCanvas.rbuf.fillRect(0,0, drawCanvas.getWidth(),Main.drawScroll.getViewport().getHeight());	
			drawCanvas.rbuf.setComposite(drawCanvas.backupr);		
			/*drawCanvas.rbuf.setComposite( drawCanvas.composite);				
			drawCanvas.rbuf.fillRect(0,0, (int)Main.screenSize.getWidth(),Main.drawScroll.getViewport().getHeight());	
			drawCanvas.rbuf.setComposite(drawCanvas.backupr);		
		*/
		}
		
		
	//	if(drawScroll.getViewport().getHeight() > 0) {
	//		Main.drawCanvas.resizeCanvas(drawCanvas.getWidth(), drawScroll.getViewport().getHeight());
	//	}
			bedCanvas.repaint();
			Main.chromDraw.updateExons = true;
			drawCanvas.repaint();
			
	//		Draw.setScrollbar(0);
			Main.chromDraw.repaint();
		
	}
	public void actionPerformed(ActionEvent e) {
		Logo.frame.setVisible(false);
		
		if(e.getSource() == manage) {			
			
			VariantHandler.frame.setLocation(frame.getLocationOnScreen().x+10, frame.getLocationOnScreen().y+10);
			VariantHandler.frame.setState(JFrame.NORMAL);
			VariantHandler.frame.setVisible(true);
			Draw.calculateVars = true;
			Draw.updatevars = true;
			drawCanvas.repaint();
			
		}
		else if(e.getSource() == average) {
		
		
			Average.setSamples();
			Average.frame.setLocation(frame.getLocationOnScreen().x+10, frame.getLocationOnScreen().y+10);
			Average.frame.setState(JFrame.NORMAL);
			Average.frame.setVisible(true);
		}
		else if(e.getSource() == errorlog) {
			ErrorLog.frame.setLocation(frame.getLocationOnScreen().x+10, frame.getLocationOnScreen().y+10);
			//		VariantHandler.frame.setAlwaysOnTop(true);	
			ErrorLog.frame.setState(JFrame.NORMAL);
			ErrorLog.frame.setVisible(true);
			
		}
		else if(e.getSource() == help) {
			JOptionPane.showMessageDialog(Main.chromDraw, "This is pre-release version of BasePlayer\nContact: help@baseplayer.fi\nUniversity of Helsinki", "Help", JOptionPane.INFORMATION_MESSAGE);
		}
		else if(e.getSource() == settings) {
			Settings.frame.setLocation(frame.getLocationOnScreen().x+10, frame.getLocationOnScreen().y+10);			
			Settings.frame.setState(JFrame.NORMAL);
			Settings.frame.setVisible(true);
		}
		else if(e.getSource() == update) {
			try {
				Updater update = new Updater();
				update.execute();
				/*	}
	       
	        InputStream html = null;

	        html = url.openStream();
	        
	        int c = 0;
	        StringBuffer buffer = new StringBuffer("");

	        while(c != -1) {
	            c = html.read();
	            
	        buffer.append((char)c);
	        }
	        System.out.println(buffer.toString());*/
			}
			catch(Exception ex) {
				ex.printStackTrace();
			}
			
		}
		else if(e.getSource() == clearMemory) {
			
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
					currentvar.setBedhit();					
					currentvar = currentvar.getNext();
				}				
				if(currentvar == null) {
					break;
				}
				currentbed = currentbed.getNext();				
			}
			
			/*
			if(trackPane.getDividerSize() != 0) {			
				trackPane.setDividerSize(0);
				
				trackPane.setResizeWeight(0.0);
				splitPane.revalidate();
			
			}
			else {
				trackPane.setDividerSize(5);				
				trackPane.setResizeWeight(0.2);
				bedCanvas.setSize(drawScroll.getWidth(), 200);
				
			}
			/*
			Transcript transcript;
			VarNode current = FileRead.head.getNext();
			Transcript.Exon exon;
			
			while(current != null) {
				if(current.getExons().size() == 0) {
					current = current.getNext();
					continue;
				}
				current.amino = chromDraw.getAminoChange(current, true);
				//System.out.println(chromDraw.getAminoChange(current, true));
				current = current.getNext();
				
			}
			Draw.updatevars = true;
			drawCanvas.repaint();
					
		
				  Main.drawCanvas.insertList.clear();
				  FileRead filereader = new FileRead();       	  
				  filereader.searchInsSites = true;
				  filereader.execute();
			 */
		}
		else if (e.getSource() == clear) {
			clearData();			
			
		}
		
		else if (e.getSource() == opensamples) {
			 try {					    	
		    	
		//    	  JFileChooser chooser = new JFileChooser(path);
				 JFileChooser chooser = new JFileChooser(path);	 
		    	  chooser.setMultiSelectionEnabled(true);
		    	  chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
		    	  chooser.setAcceptAllFileFilterUsed(false);
		    	  MyFilterBAM bamFilter = new MyFilterBAM();
		    	  MyFilterVCF vcfFilter = new MyFilterVCF();
		    	//  MyFilterBED bedFilter = new MyFilterBED();
		 //   	  MyFilterCRAM cramFilter = new MyFilterCRAM();
		    	//  MyFilterALL allFilter = new MyFilterALL();
		    	  chooser.addChoosableFileFilter(vcfFilter); 	
		    	  chooser.addChoosableFileFilter(bamFilter); 
		    	//  chooser.addChoosableFileFilter(bedFilter); 
		    //	  chooser.addChoosableFileFilter(cramFilter);
		    //	  chooser.addChoosableFileFilter(allFilter);
		    	  chooser.setDialogTitle("Add samples");
		    	  chooser.setPreferredSize(new Dimension((int)screenSize.getWidth()/3, (int)screenSize.getHeight()/3));
		          int returnVal = chooser.showOpenDialog((Component)this.getParent());	         
		          
		         if (returnVal == JFileChooser.APPROVE_OPTION) {
		        	  
		        	  File vcffiles[] = chooser.getSelectedFiles(); 
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
		        	/*  else if(chooser.getFileFilter().equals(chooser.getChoosableFileFilters()[2])){
		        		  filereader.readBED(vcffiles);
		        		  
		        	  }*/
		        	/*  else if(chooser.getFileFilter().equals(chooser.getChoosableFileFilters()[2])) {
		        		 
		        	  }*/
		          }
			 }
			 catch(Exception ex) {
				 ex.printStackTrace();
			 }
			
		}
		else if (e.getSource() == addcontrols) {
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
			
	    	    	  
			 JFileChooser chooser = new JFileChooser(Main.trackDir);	 
			 getText(chooser.getComponents());
	    	  chooser.setMultiSelectionEnabled(true);	    	  
	    	  chooser.setAcceptAllFileFilterUsed(false);	    	 
	    	  MyFilterBED bedFilter = new MyFilterBED();	    	
	    	  chooser.addChoosableFileFilter(bedFilter); 
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
	        		 System.out.println(bamIterator.next().getAlignmentStart());
	        		 reader.close();
	        	//	 SamReaderFactory.make().
	        	//	 htsjdk.samtools.SamInputResource.of(url)
		 //       	 URL url = new URL("http://dx2-tutoh-1.ltdk.helsinki.fi/SELEX_Ensembl_2015/hits/SELEX_Ensembl_dom_set_2015_humanGRCh37_9_or_max1M.gff.gz");
		        //	 System.out.println(url.getPath());
		        	
		   //    	 SeekableStream stream = SeekableStreamFactory.getInstance().getStreamFor(url);
		    //   	 System.out.println(url.getHost() +url.getPath());
		    //    TabixReader tabixReader = new TabixReader(url.getHost() +url.getPath(), "http://dx2-tutoh-1.ltdk.helsinki.fi/SELEX_Ensembl_2015/hits/SELEX_Ensembl_dom_set_2015_humanGRCh37_9_or_max1M.gff.gz.tbi", stream);   	
		 			
		        //	 System.out.println(tabixReader.getSource());
		        	 
		 		//	TabixReader.Iterator bedIterator = null;
		 			  try {
		 				 
		 			//	 bedIterator = tabixReader.query(Main.chromosomeDropdown.getSelectedItem().toString() +":" +5000000 +"-" +6000000);
		 			  }
		 			  catch(Exception ex) {
		 				  ex.printStackTrace();
		 			  }
		   //     	 System.out.println(bedIterator.next());
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
	        	
	        	  File vcffiles[] = chooser.getSelectedFiles(); 
	        	  if(vcffiles[0].exists()) {
	        		  trackDir = vcffiles[0].getParent();	 	
	        		  writeToConfig("DefaultTrackDir=" +trackDir);
		        	  FileRead filereader = new FileRead(vcffiles);
		        	  if(chooser.getFileFilter().equals(chooser.getChoosableFileFilters()[0])) {
		        		  filereader.readBED(vcffiles);
		        	  }	
	        	  }
	        	  else { 		  
	        		 
	        		  
	        		  if(Main.chooserText.length() > 5 && Main.chooserText.endsWith(".bed.gz") || Main.chooserText.endsWith(".gff.gz") || Main.chooserText.endsWith(".bedgraph.gz")) {	        			  
	        			 
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
		        				  
		        				  FileRead filereader = new FileRead(vcffiles);
		        				  filereader.readBED(Main.chooserText, index);
		        				  
		        				  tabixReader.close();
		        			  }
		        		  }
	        		  }
	        		  else {
	        			  if(Main.chooserText.contains("://")) {
			        			
		        			  URL url = new URL(Main.chooserText);	
		        			
		      		      	  SeekableStream stream = SeekableStreamFactory.getInstance().getStreamFor(url);
		      		      	
		      		    //  	  BBFileReader bbreader = null;
		      		      	
		      		      	  try {
		      		      //		  bbreader = new BBFileReader(Main.chooserText, stream);
		      		      		
		      		      	  }
		      		      	  catch(Exception ex) {
		      		      		 ex.printStackTrace();
		      		      	  }
		      		      	 
		      		      	FileRead filereader = new FileRead(vcffiles);
	        				  filereader.readBED(Main.chooserText, "nan");
		        		/*	  if(bbreader != null) {
		        				  System.out.println(Main.chooserText);
		        				  
		        				  
		        				 
		        			  }*/
		      		      	  stream.close();
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
			  
			 JFileChooser chooser = new JFileChooser(Main.projectDir);	    	
			 chooser.setAcceptAllFileFilterUsed(false);
	    	  MyFilterSES sesFilter = new MyFilterSES();	    	  
	    	  chooser.addChoosableFileFilter(sesFilter);
	    	  chooser.setDialogTitle("Open project");
	    	  chooser.setPreferredSize(new Dimension((int)screenSize.getWidth()/3, (int)screenSize.getHeight()/3));
	          int returnVal = chooser.showOpenDialog((Component)this.getParent());	         	         
	    	  
	          if (returnVal == JFileChooser.APPROVE_OPTION) {
	        	  clearData();
					VariantHandler.maxCoverage = 1500;
					projectDir = chooser.getSelectedFile().getParent();	 	
	        		writeToConfig("DefaultProjectDir=" +projectDir);
					try {
						FileInputStream fin = new FileInputStream(chooser.getSelectedFile());
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
							
							Main.samples = (short)drawCanvas.sampleList.size();			
						
							drawCanvas.splits = (ArrayList<SplitClass>) ois.readObject();
							for(int i = 0 ; i<drawCanvas.sampleList.size(); i++) {
								drawCanvas.sampleList.get(i).resetreadHash();
								if(drawCanvas.sampleList.get(i).getTabixFile() != null || drawCanvas.sampleList.get(i).multipart) {		
									if(drawCanvas.sampleList.get(i).vcfchr == null) {
										drawCanvas.sampleList.get(i).vcfchr = "";
									}
									if(!drawCanvas.sampleList.get(i).multipart) {
										if(drawCanvas.sampleList.get(i).getVCFInput() == null) {
											drawCanvas.sampleList.get(i).setInputStream();
										}
									}
									Main.varsamples++;
								}
								if(drawCanvas.sampleList.get(i).longestRead == null) {
									
									drawCanvas.sampleList.get(i).longestRead = 0;
								}
							}
							for(int i = 0; i<drawCanvas.splits.size(); i++) {
								drawCanvas.splits.get(i).setDivider(4.0);
							}
						
							for(int i = 0; i<drawCanvas.splits.size(); i++) {
								drawCanvas.splits.get(i).resetSplits();
							}
							drawCanvas.drawVariables = (DrawVariables)ois.readObject();						
							VariantHandler.commonSlider.setMaximum(Main.varsamples);
							VariantHandler.commonSlider.setValue(1);
							VariantHandler.commonSlider.setUpperValue(Main.varsamples);
							VariantHandler.geneSlider.setMaximum(Main.varsamples);
							VariantHandler.geneSlider.setValue(1);
							try {
								Control.controlData = (ControlData)ois.readObject();
								
								if(Control.controlData.fileArray.size() > 0) {
									if(Control.controlData.fileArray.get(0) instanceof ControlFile) {
									  Main.trackPane.setVisible(true);
								//	  Main.trackPane.setDividerSize(3);
									  Main.varpane.setDividerSize(3);	 									
									  Main.varpane.setDividerLocation(0.1);
									  Main.controlScroll.setVisible(true);
									  Main.controlDraw.setVisible(true);	
									  varpane.revalidate();
									  
									  for(int i = 0 ; i<Control.controlData.fileArray.size(); i++) {
										
										  VariantHandler.table.addRowGeneheader("AF: " +Control.controlData.fileArray.get(i).getName());
										  VariantHandler.table.addRowGeneheader("OR");
										  if(VariantHandler.tables.size() > 0) {
											  for(int t = 0; t < VariantHandler.tables.size(); t++) {
												  VariantHandler.tables.get(t).addRowGeneheader("AF: " +Control.controlData.fileArray.get(i).getName());
												  VariantHandler.tables.get(t).addRowGeneheader("OR");
											  }
										  }
										  VariantHandler.clusterTable.addRowGeneheader("AF: " +Control.controlData.fileArray.get(i).getName());
										  VariantHandler.clusterTable.addRowGeneheader("OR");
									  }
								}
								else {									
									  Control.controlData.fileArray = Collections.synchronizedList(new ArrayList<ControlFile>());
								}
							}
												
					
					//		if(ois.available() > 0) {
							try {
								bedCanvas.bedTrack = (ArrayList<BedTrack>)ois.readObject();
							}
							catch(EOFException excep) {
								
							}
						//	}
							if(bedCanvas.bedTrack != null && bedCanvas.bedTrack.size() > 0) {
								
								if(Main.trackPane.isVisible()) {
									
									 Main.varpane.setDividerLocation(Main.varpane.getDividerLocation()+100);
									  Main.trackPane.setDividerLocation(Main.varpane.getDividerLocation()/2);
									 if(Main.controlScroll.isVisible()) {
										 Main.trackPane.setDividerSize(4);
									 }
								}
								else {
									Main.trackPane.setVisible(true);
									Main.varpane.setDividerLocation(0.1);
									Main.varpane.setDividerSize(3);	
								}
								 Main.bedScroll.setVisible(true);
								 Main.bedCanvas.setVisible(true);	
								  
								 
								for(int i = 0 ; i< bedCanvas.bedTrack.size(); i++) {
									bedCanvas.bedTrack.get(i).setHead();
									//bedCanvas.bedTrack.get(i).setNames();
									bedCanvas.bedTrack.get(i).setColors();
									bedCanvas.bedTrack.get(i).first = true;
									bedCanvas.trackDivider.add(0.0);
									bedCanvas.bedTrack.get(i).setBedLevels();
									FileRead.setTable(bedCanvas.bedTrack.get(i));
								}
								
							}
						}
						catch(Exception ex) {
							ex.printStackTrace();
							clearData();
						}
						
						ois.close();
						
					/*	if(Control.controlData.fileArray.size() > 0) {
							Control.useCheck.setEnabled(true);
						}*/
						frame.setTitle("BasePlayer - Project: " +drawCanvas.drawVariables.projectName);
						
						Main.drawCanvas.drawVariables.visiblesamples = Main.samples;
						Main.drawCanvas.checkSampleZoom();
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
						
					//	for(int i = 0 ;i<drawCanvas.sampleList.size(); i++) {
							//for(int r = 0 ;r<drawCanvas.sampleList.size(); r++) {
							//	if(drawCanvas.sampleList.get(i).getReadClass().size() > 0) {
					//				drawCanvas.sampleList.get(i).resetReadClass();
							//	}
							//}
					//	}
						}
						catch(Exception ex) {
							JOptionPane.showMessageDialog(drawCanvas, "Sorry, your project must be created again.", "Error", JOptionPane.ERROR_MESSAGE);
				    		clearData();
				    		ex.printStackTrace();
						}
					}
					catch(Exception ex) {
						ex.printStackTrace();
						clearData();
					}
	        	  	        	  
	          }
		}
		else if (e.getSource() == saveProjectAs) {
			
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
}
	
	public void getText(Component[] comp)  {
		
	    for(int x = 0; x < comp.length; x++)   {
	      if(comp[x] instanceof JPanel) {
	    	  getText(((JPanel)comp[x]).getComponents());
	    	 
	      }
	      else if(comp[x] instanceof JTextField)
	      {
	     //   ((JTextField)comp[x]).setEditable(false);
	    	// System.out.println();
	     
	        chooserTextField = ((JTextField)comp[x]);
	        chooserTextField.addKeyListener(new KeyListener() {

				@Override
				public void keyTyped(KeyEvent e) {
				
					Main.chooserText = chooserTextField.getText().replace(" ", "");
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

void clearData() {
	undoList.clear();
	undoPointer = 0;
	bedCanvas.bedOn = false;
	back.setEnabled(false);
	forward.setEnabled(false);
	if(drawCanvas.clusterNodes != null) {
		drawCanvas.clusterNodes.clear();
	}
	Average.outVector.clear();
	zoomout();
	Main.bedCanvas.trackDivider.clear();
	samples = 0;			
	varsamples = 0;
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
		FileRead.removeTable(bedCanvas.bedTrack.get(i));	
		if(bedCanvas.bedTrack.get(i).getTable() != null) {
			bedCanvas.bedTrack.get(i).getTable().bedarray.clear();
		}
	}
	SplitClass split = Main.drawCanvas.splits.get(0);
	for(int g = 0 ; g<split.getGenes().size(); g++) {
		for(int i = 0; i<split.getGenes().get(g).getTranscripts().size(); i++) {
			split.getGenes().get(g).getTranscripts().get(i).mutations = 0;
			split.getGenes().get(g).getTranscripts().get(i).missense = 0;
			split.getGenes().get(g).getTranscripts().get(i).nonsense = 0;
			split.getGenes().get(g).getTranscripts().get(i).synonymous = 0;
			split.getGenes().get(g).getTranscripts().get(i).intronic = 0;
			split.getGenes().get(g).getTranscripts().get(i).utr = 0;
			split.getGenes().get(g).getTranscripts().get(i).samples.clear();
			split.getGenes().get(g).getTranscripts().get(i).varnodes.clear();
		}
	}
	split = null;
	bedCanvas.bedTrack.clear();
	bedCanvas.drawNode = null;
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
	public static class runner extends SwingWorker<String, Object> {		
		
		protected String doInBackground() {			
			return "";
		}
		
	}
	private static void createAndShowGUI() {
		try {
		
			if(System.getProperty("os.name").toLowerCase().contains("nux")) {
				 UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
			}
			else {
				UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
			}
		//	SystemUtils.IS_OS_WINDOWS
		//	UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		JFrame.setDefaultLookAndFeelDecorated(false);
	    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	   
	    frame.setResizable(true);
	    
	    JComponent newContentPane = new Main();
	    newContentPane.setOpaque(true); 
	    
	    frame.setContentPane(newContentPane);

	    frame.pack();
	    frame.setVisible(true);
	    
	    if(Main.menubar.getWidth() > 100){
	    	Main.menubar.setMinimumSize(new Dimension(Main.menubar.getWidth(), Main.menubar.getHeight()));
	    }
	    else {
	    	Main.menubar.setMinimumSize(new Dimension(180, 22));
	    }
	}
	public static void main(String[] args) {
		
		

		System.setProperty("sun.java2d.d3d", "false"); 
	
		Main.args = args;
		 javax.swing.SwingUtilities.invokeLater(new Runnable() {
	            public void run() {	     
	            	
	            		createAndShowGUI();	            	
	            }
	        });			
	}
	
	static void updatePositions(double start, double end) {
			positionField.setText("chr" +drawCanvas.selectedSplit.chrom +":" +MethodLibrary.formatNumber((int)start) +"-" +MethodLibrary.formatNumber((int)end));
			Main.widthLabel.setText("" +MethodLibrary.formatNumber((int)(end-start)) +"bp");
		
	}
	@Override
	public void stateChanged(ChangeEvent event) {
		
		
	}

	public void componentResized(ComponentEvent e) {
		
		if(e.getComponent().getName().contains("frame0")) {
		
		if(drawScroll.getViewport().getWidth() > 0) {		
			drawDimensions.setSize(drawScroll.getViewport().getWidth(),drawScroll.getViewport().getSize().height );	
			drawCanvas.setPreferredSize(drawDimensions);
			chromDimensions.setSize(drawScroll.getViewport().getWidth(), chromDimensions.getHeight());
			chromDraw.setPreferredSize(chromDimensions);
			drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), drawScroll.getViewport().getHeight());
			
			if(drawCanvas.splits.size() > 0) {
				for(int i = 0 ; i<drawCanvas.splits.size(); i++) {
					drawCanvas.splits.get(i).updateReads = true;				
				}
			}
			drawCanvas.repaint();
			
		}
			
			/*
			if(samples > 0) {			
				
				if(Main.samples*Main.drawCanvas.sampleHeight < drawScroll.getViewport().getHeight()) {
					Main.drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), drawScroll.getViewport().getHeight());
					Main.drawCanvas.revalidate();
					
				//	drawCanvas.resizeCanvas(drawScroll.getHeight());
//				Main.drawPanel.setPreferredSize(new Dimension(Main.drawCanvas.getWidth(), drawScroll.getHeight()));
				}
				else {
					Main.drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), Main.samples*drawCanvas.sampleHeight);
					
					Main.drawCanvas.revalidate();
					
				//	drawCanvas.resizeCanvas( Main.samples*Main.drawCanvas.sampleHeight);
//					Main.drawPanel.setPreferredSize(new Dimension(Main.drawCanvas.getWidth(), Main.samples*Main.drawCanvas.sampleHeight));
				}
			
			 }
			 else {
				
				 drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), drawScroll.getViewport().getHeight());		
			 }
			*/
			//Main.drawCanvas.setPreferredSize(new Dimension(frame.getWidth(), frame.getHeight()));	
			 
	//		 Main.drawScroll.setPreferredSize(new Dimension(frame.getWidth(), frame.getHeight()));
			 
	     }       
	}
	static void cancel() {
		  cancel = true;
		  if(drawCanvas.loadingtext.contains("Loading variants")) {
			  Main.drawCanvas.ready("Loading variants...");
			  drawCanvas.current = null;
			   drawCanvas.currentDraw = null;
			   chromDraw.vardraw = null;
			   chromDraw.varnode = null;
			   drawCanvas.variantsStart = 0;
			   drawCanvas.variantsEnd = 0;
			   Draw.updatevars = true;
			   FileRead.cancelvarcount = true;
			   FileRead.cancelfileread = true;
		  }
		  else if(drawCanvas.loadingtext.contains("Processing variants")) {
			  Main.drawCanvas.ready("Processing variants...");
			 
			   FileRead.cancelvarcount = true;
			 
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
					    }
						  
				   }
			   }
			   Main.drawCanvas.clearReads();
			   Main.drawCanvas.ready("Loading reads");
			//   Main.drawCanvas.timer = System.currentTimeMillis();
			   Draw.updateReads = true;
			   Draw.updateCoverages = true;
			   for(int i = 0; i<drawCanvas.splits.size(); i++) {
				   drawCanvas.splits.get(i).updateReads = true;
				 
			   }
			  
			   FileRead.cancelfileread = true;
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
			   Main.drawCanvas.ready("Updating BasePlayer...");
		  }
		  else {
			  Main.drawCanvas.ready("all");
		  }
		  
		  
	//	   Main.drawCanvas.readyQueue.clear();
		
		   cancel = false;
		   drawCanvas.repaint();
		   chromDraw.repaint();
	}
	static void cancel(int nro) {
		  cancel = true;
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
		   Main.drawCanvas.ready("all");
		//   Main.drawCanvas.timer = System.currentTimeMillis();
		   drawCanvas.current = null;
		   drawCanvas.currentDraw = null;
		   chromDraw.vardraw = null;
		   chromDraw.varnode = null;
		   Draw.updatevars = true;
		   FileRead.cancelvarcount = true;
		   FileRead.cancelfileread = true;
		   FileRead.cancelreadread = true;
	//	   Main.drawCanvas.readyQueue.clear();
	//	   JOptionPane.showMessageDialog(Main.glassPane, "Try to search more specific region to visualize all data", "Information", JOptionPane.ERROR_MESSAGE);
		  
		   cancel = false;
		//   chromDraw.repaint();
		   
		   
		   
	}
	
	public void writeToConfig(String attribute) {
		Main.configChanged = true;
		Boolean found = false;
		if(Launcher.config.get(0).contains("Riku")) {
			Launcher.config.set(0, "#Config file for BasePlayer");
		}
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
		   System.out.println("Genome changed to: " +event.getComponent().getName());
		   writeToConfig("DefaultGenome=" +event.getComponent().getName());
		   
		   setChromDrop(event.getComponent().getName()); 		  
		   getBands();
		   try {
			   if(genomehash.get(event.getComponent().getName()).size() > 0) {
				    ChromDraw.exonReader = new TabixReader(genomehash.get(event.getComponent().getName()).get(0).getCanonicalPath());
				//	annotationfile = genomehash.get(event.getComponent().getName()).get(0).getName();
					
				//	writeToConfig("DefaultGenes=" +annotationfile);
					String s;
					String[] exonSplit;
					searchTable.clear();
					while((s = ChromDraw.exonReader.readLine()) != null) {
						exonSplit = s.split("\t");
						if(exonSplit[9].equals("1")) {		
							
							String[] adder = {exonSplit[0], exonSplit[1], exonSplit[2]};
							searchTable.put(exonSplit[3].toUpperCase(), adder);					
						}						
					}
				
			   }
			   }
			   catch(Exception e) {
				   e.printStackTrace();
			   }
		   chromosomeDropdown.setSelectedIndex(0);
		   drawCanvas.chrom = chromosomeDropdown.getItemAt(0);
		  
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
	public class CheckUpdates extends SwingWorker<String, Object> {
		

		protected String doInBackground() {
			try {
								
			    URL file = new URL("https://www.cs.helsinki.fi/u/rkataine/BasePlayer/update/BasePlayer.jar");
				
			    HttpURLConnection httpCon = (HttpURLConnection) file.openConnection();	 
			    File homefile = new File(userDir +"/BasePlayer.jar");
			   
			   
			    if(httpCon.getLastModified() != homefile.lastModified()) {
			    	chromDraw.message = "Updates available. Please click File->Update to get the most recent version.";
			    	chromDraw.timer = System.currentTimeMillis();
			    	chromDraw.repaint();
			    	update.setEnabled(true);
					
			    }
			    else {
			    	chromDraw.message = "BasePlayer is up-to-date.";
			    	chromDraw.timer = System.currentTimeMillis();
			    	chromDraw.repaint();
			    	update.setEnabled(false);
			    }
			    httpCon.disconnect();
			    file = new URL("https://www.cs.helsinki.fi/u/rkataine/BasePlayer/update/Launcher.jar");
			    homefile = new File(userDir +"/Launcher.jar");
			    httpCon = (HttpURLConnection) file.openConnection();	 
			    if(!homefile.exists() || httpCon.getLastModified() != homefile.lastModified()) {
			    	updatelauncher = true;					
			    }
			    else {
			    	updatelauncher = false;
			    }
			    httpCon.disconnect();
			    
			}
			catch(Exception e) {
				JOptionPane.showMessageDialog(Main.chromDraw, e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
				update.setEnabled(false);
			}
			return "";
		}
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
	public class Updater extends SwingWorker<String, Object> {
		
		
		public void downloadFile(String fileURL, String saveDir) throws IOException {
			URL url = new URL(fileURL);
			HttpURLConnection httpConn = (HttpURLConnection) url.openConnection();
			
			int responseCode = httpConn.getResponseCode();
			int BUFFER_SIZE = 4096;
			// always check HTTP response code first
			if (responseCode == HttpURLConnection.HTTP_OK) {
				String fileName = "";
				String disposition = httpConn.getHeaderField("Content-Disposition");
				String contentType = httpConn.getContentType();
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

				System.out.println("Content-Type = " + contentType);
				System.out.println("Content-Disposition = " + disposition);
				System.out.println("Content-Length = " + contentLength);
				System.out.println("fileName = " + fileName);

				// opens input stream from the HTTP connection
				InputStream inputStream = httpConn.getInputStream();
				String saveFilePath = saveDir + File.separator + fileName +"_temp";
				
				// opens an output stream to save into file
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
					
					JOptionPane.showMessageDialog(Main.chromDraw,"BasePlayer updated. Please restart program to apply changes.", "Note", JOptionPane.INFORMATION_MESSAGE);
					update.setEnabled(false);
				}
				else {
					JOptionPane.showMessageDialog(Main.chromDraw,"BasePlayer could not be updated.", "Note", JOptionPane.INFORMATION_MESSAGE);
				}				

				System.out.println("File downloaded");
			} else {
				System.out.println("No file to download. Server replied HTTP code: " + responseCode);
				ErrorLog.addError("No file to download. Server replied HTTP code: " + responseCode);
			}
			httpConn.disconnect();
		}
		
		protected String doInBackground() {
			
			try {
				Main.drawCanvas.loading("Updating BasePlayer...");
				   
			    
			  
			    downloadFile("https://www.cs.helsinki.fi/u/rkataine/BasePlayer/update/BasePlayer.jar", userDir);
				if(updatelauncher) {
					downloadFile("https://www.cs.helsinki.fi/u/rkataine/BasePlayer/update/Launcher.jar", userDir);
				}
				Main.drawCanvas.ready("Updating BasePlayer...");	
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
				 Main.drawCanvas.ready("Updating BasePlayer...");	
				e.printStackTrace();
				ErrorLog.addError(e.getStackTrace());
				JOptionPane.showMessageDialog(Main.chromDraw,e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
				
			}
			return "";
		}
		
	}
	@Override
	public void mouseEntered(MouseEvent event) {
		if(event.getSource() == filemenu) {
			filemenu.doClick();
		}
		else if(event.getSource() == toolmenu) {
			toolmenu.doClick();
		}
		else if(event.getSource() == help) {
			help.doClick();
		}
		
	}


	@Override
	public void mouseExited(MouseEvent event) {
		
		if(!drawCanvas.loading && !drawCanvas.scrollbar && event.getSource() == drawScroll.getVerticalScrollBar()) {
			Draw.setGlasspane(false);
		}
		if(event.getSource() == drawCanvas) {
			Main.drawCanvas.sidebar = false;
		
		//	Main.drawCanvas.selectedSampleIndex = -1;
			Main.drawCanvas.repaint();
	
		}
	}


	@Override
	public void mousePressed(MouseEvent event) {
		Logo.frame.setVisible(false);
		if(event.getSource() == splitPaneDivider) {
			Main.vardivider = bedCanvas.nodeImage.getHeight()/(double)varPaneDivider.getY();
			Main.bedCanvas.resize = true;
		}
		if(event.getSource() == varPaneDivider) {
			Main.bedCanvas.resize = true;
			
			Main.vardivider = bedCanvas.nodeImage.getHeight()/(double)varPaneDivider.getY();
			
		}
		if(event.getSource() == filemenu) {			
			if(!filemenu.isSelected()){				
				filemenu.doClick();			
			}
		}
		else if(event.getSource() == toolmenu) {
			if(!toolmenu.isSelected()){				
				toolmenu.doClick();			
			}
		}
		else if(drawCanvas.loadingtext.equals("note")) {
			Main.drawCanvas.loadingtext = "";
			Main.drawCanvas.ready("note");
		}
		else if(event.getSource() == drawScroll.getVerticalScrollBar()) {
			
			drawCanvas.scrollbar = true;
			Draw.setGlasspane(true);			
			if(drawCanvas.drawVariables.visiblestart != (short)(Main.drawScroll.getVerticalScrollBar().getValue()/drawCanvas.drawVariables.sampleHeight)) {
				drawCanvas.drawVariables.visiblestart = (short)(Main.drawScroll.getVerticalScrollBar().getValue()/drawCanvas.drawVariables.sampleHeight);					
				
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
			}
		}
		else if(event.getSource() == searchField) {
			searchField.requestFocus();
			searchField.setForeground(Color.black);		
			if(searchField.getText().contains("Search by")) {
				searchField.setText("");
			}
		}
		else if(event.getComponent().getName() != null) {
				
			   for(int i = 0; i<genome.getItemCount(); i++) {
				   for(int j = 0; j<genomehash.get(genome.getItem(i).getName()).size(); j++) {
					  
					   if(genomehash.get(genome.getItem(i).getName()).get(j).getName().contains(event.getComponent().getName())) {
						 
						     
						   	 System.out.println("Annotation for " +genome.getItem(i).getName() +" changed to " +genomehash.get(genome.getItem(i).getName()).get(j).getName());
						     annotationfile = genomehash.get(genome.getItem(i).getName()).get(j).getName();
						     defaultGenome = genome.getItem(i).getName();
						     
						     Main.annotation = j;
						     writeToConfig("DefaultGenome=" +defaultGenome);
						     writeToConfig("DefaultGenes=" +annotationfile);
						     setChromDrop(genome.getItem(i).getName()); 		
						     getBands();
						     try {
								  
						    	 if(genomehash.get(genome.getItem(i).getName()).get(j).getName().endsWith(".bed.gz")) {
						    		 	
									    ChromDraw.exonReader = new TabixReader(genomehash.get(defaultGenome).get(j).getCanonicalPath());
									    
										String s;
										String[] exonSplit;
									
										while((s = ChromDraw.exonReader.readLine()) != null) {
										
											exonSplit = s.split("\t");
											if(exonSplit[9].equals("1")) {		
												
												String[] adder = {exonSplit[0], exonSplit[1], exonSplit[2]};
												searchTable.put(exonSplit[3].toUpperCase(), adder);					
											}						
										}
										 ChromDraw.exonReader.close();
										 
						    	 }
						    /*	 else if(genomehash.get(genome.getItem(i).getName()).get(j).getName().contains(".gff")) {
						    		
								    FileRead.readGFF(genomehash.get(genome.getItem(i).getName()).get(j), genomehash.get(genome.getItem(i).getName()).get(j).getCanonicalPath().substring(0,genomehash.get(genome.getItem(i).getName()).get(j).getCanonicalPath().indexOf(".gff")) +".bed.gz");
								    	 
								     
						    	 }
								*/   
								   }
								   catch(Exception e) {
									   e.printStackTrace();
								   }
						     drawCanvas.chrom = chromosomeDropdown.getItemAt(0);
							   chromosomeDropdown.setSelectedIndex(0);
							   
						     
					   }
				   }
				  
			   }
			//   System.out.println("Genome changed to: " +event.getComponent().getName());
			    
			   
			  
			  
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
	
	//	if(!drawCanvas.loading && drawCanvas.scrolldrag && event.getSource() == drawScroll.getVerticalScrollBar()) {
		//	drawCanvas.scrollbar = false;
			
		/*	for(int i = drawCanvas.drawVariables.visiblestart; i<drawCanvas.drawVariables.visibleend; i++) {
				drawCanvas.sampleList.get(i).prepixel = 0;
			}*/
		//	Draw.setScrollbar((int)(drawCanvas.drawVariables.visiblestart*drawCanvas.drawVariables.sampleHeight));
			//Draw.updatevars = true;
			//Draw.updateReads = true;
			//Draw.updateCoverages = true;
	//		drawCanvas.repaint();
			
	//	}
		
		if(event.getSource() == drawScroll.getVerticalScrollBar()) {
		/*	for(int i = 0; i < drawCanvas.drawVariables.visiblestart-1; i++) {
				drawCanvas.clearReads(drawCanvas.sampleList.get(i));
			}
			for(int i = drawCanvas.drawVariables.visibleend+1; i < Main.samples; i++) {
				drawCanvas.clearReads(drawCanvas.sampleList.get(i));
			}
			drawCanvas.drawVariables.visiblestart = (short)(Main.drawScroll.getVerticalScrollBar().getValue()/drawCanvas.drawVariables.sampleHeight);					
			drawCanvas.drawVariables.visibleend = (short)(drawCanvas.drawVariables.visiblestart + drawScroll.getViewport().getHeight()/drawCanvas.drawVariables.sampleHeight);	
			
			if(drawCanvas.splits.size() > 1) {
				for(int i = 0; i<drawCanvas.splits.size(); i++) {
					drawCanvas.splits.get(i).updateReads = true;					
				}
			}
			else {
				Draw.updateReads = true;
				Draw.updatevars = true;
			}
			drawCanvas.scrollbar = false;
		//	Draw.updateReads = true;
			drawCanvas.repaint();
		*/
		}
		
		if(!drawCanvas.loading) {
			Draw.setGlasspane(false);
		}
		drawCanvas.scrolldrag = false;
		
		
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
	@Override
	public void propertyChange(PropertyChangeEvent event) {
		
		if(drawScroll.getViewport().getWidth() > 0) {
	
			if(event.getSource() == splitPane) {					
				
				if(samples > 0) {					
					if(samples*drawCanvas.drawVariables.sampleHeight < drawScroll.getViewport().getHeight()) {						
						drawDimensions.setSize(drawScroll.getViewport().getWidth(),drawScroll.getViewport().getSize().height );	
						drawCanvas.setPreferredSize(drawDimensions);
						chromDimensions.setSize(drawScroll.getViewport().getWidth(), chromDimensions.getHeight());
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
				chromDraw.updateExons = true;
				chromDraw.repaint();
			}
			if(event.getSource() == varpane) {			
				
				if(samples > 0) {		
					
					if(samples*drawCanvas.drawVariables.sampleHeight < drawScroll.getViewport().getHeight()) {						
						drawDimensions.setSize(drawScroll.getViewport().getWidth(),drawScroll.getViewport().getSize().height );	
						drawCanvas.setPreferredSize(drawDimensions);
						chromDimensions.setSize(drawScroll.getViewport().getWidth(), chromDimensions.getHeight());
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
				/*
					drawDimensions.setSize(drawScroll.getViewport().getWidth(),drawScroll.getViewport().getHeight() );		
					drawCanvas.setPreferredSize(drawDimensions);			
					drawCanvas.resizeCanvas((int)drawDimensions.getWidth(), (int)drawDimensions.getHeight());						
					bedCanvas.setPreferredSize(new Dimension((int)drawDimensions.getWidth(), bedScroll.getViewport().getHeight() ));
					*/
				/*	if(samples > 0) {
				
					if(samples*drawCanvas.drawVariables.sampleHeight < drawScroll.getViewport().getHeight()) {
					
						drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), drawScroll.getViewport().getHeight());	
					
					}
					Draw.setScrollbar((int)(drawCanvas.drawVariables.visiblestart*drawCanvas.drawVariables.sampleHeight));
				}
				drawCanvas.repaint();
			*/
			}	
		}
		/*
		if(splitPane.getComponent(0).getHeight() > 0) {
			Draw.resize = true;
			
			for(int i = 0; i<splitPane.getComponentCount(); i++) {						
				
				if(splitPane.getComponent(i) == chromScroll) {
					System.out.println(i);
					if(splitPane.getComponent(i).getHeight() > 0) {						
						chromDraw.resizeDraw( drawScroll.getViewport().getWidth(), splitPane.getComponent(i).getHeight());
						
					}
				}			
			}
			
			for(int i = 0; i<trackPane.getComponentCount(); i++) {
				
				if(trackPane.getComponent(i) == drawScroll) {
					
					if(samples > 0) {			
						
						if(Main.samples*Main.drawCanvas.sampleHeight < drawScroll.getViewport().getHeight()) {
							
							Main.drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), drawScroll.getViewport().getHeight());
							Main.drawCanvas.revalidate();
						}
						else {
							
							Main.drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), Main.samples*drawCanvas.sampleHeight);								
							Main.drawCanvas.revalidate();							
							
						}
						}
					 else {		
						
						 	drawCanvas.resizeCanvas(drawScroll.getViewport().getWidth(), drawScroll.getViewport().getHeight());
							Main.drawCanvas.revalidate();
					 }
					}
					 
			//	}
			
			}
			}
			
			Draw.resize = false;
		*/
	}

	void setChromDrop(String dir) {
		try {
			File chromindex = null;
			selectedGenome = dir;
			File defdir = new File(userDir +"/genomes/" +dir);
			File[] files = defdir.listFiles();
			String chromtemp = "";
			chromnamevector = Collections.synchronizedList(new ArrayList<String>());
			boolean faifound = false;
			for(int i = 0; i< files.length; i++) {
				if(files[i].isDirectory()) {
					continue;
				}
				if(files[i].getName().contains(".fai")) {
					chromindex = new File(defdir.getCanonicalPath() +"/" +files[i].getName());
					faifound = true;
				}
				else if(files[i].getName().contains(".fa")) {
					chromtemp = defdir.getCanonicalPath() +"/" +files[i].getName();					
				}
			}
			if(chromtemp.equals("")) {
				
				return;
			}
			if(!faifound) {
				
				//TODO
				
			}
		    referenceFile = new RandomAccessFile(chromtemp, "r");
		    ref = new File(chromtemp);
			
			BufferedReader reader = new BufferedReader(new FileReader(chromindex));
			String line;
			String[] split;
			ChromDraw.chromPos = new HashMap<String, Integer>();
			chromIndex = new Hashtable<String, Long[]>();
			
			while((line = reader.readLine()) != null) {
				split = line.split("\t");
				chromnamevector.add(split[0]);
				Long[] add = {Long.parseLong(split[2]), Long.parseLong(split[1]),  Long.parseLong(split[3])};
				if(split[0].equals("MT")) {
					ChromDraw.chromPos.put("M", Integer.parseInt(split[1]));
					chromIndex.put("M", add);
				}
				chromIndex.put(split[0], add);
				
				ChromDraw.chromPos.put(split[0], Integer.parseInt(split[1]));
				
			}
			
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
			reader.close();		
			
			chromModel = new DefaultComboBoxModel<String>(chromnames);
			chromosomeDropdown.setModel(chromModel);
			chromosomeDropdown.revalidate();
			chromosomeDropdown.repaint();
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	

	void getBands() {
		try {
			File bandfile = new File(userDir +"/genomes/" +selectedGenome +"/bands.txt");
			System.out.println(userDir +"/genomes/" +selectedGenome);
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
			drawCanvas.scrolldrag = true;
			drawCanvas.drawVariables.visiblestart = (short)(Main.drawScroll.getVerticalScrollBar().getValue()/drawCanvas.drawVariables.sampleHeight);					
			if(drawCanvas.splits.size() > 1) {
				for(int i = 0; i<drawCanvas.splits.size(); i++) {
					drawCanvas.splits.get(i).updateReads = true;					
				}
			}
			else {
				Draw.updateReads = true;
				Draw.updatevars = true;
			}
			
		}
		Draw.updatevars = true;
		// TODO Auto-generated method stub
		
	}
	@Override
	public void mouseMoved(MouseEvent event) {
		
		// TODO Auto-generated method stub
		
	}
	@Override
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void keyPressed(KeyEvent e) {
		keyCode = e.getKeyCode();
		try {
		
		}
		catch(Exception ex) {
			ex.printStackTrace();
		}
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
		else if(keyCode == KeyEvent.VK_7) {
			/*SplitClass swap = drawCanvas.splits.get(1);
			
			drawCanvas.splits.remove(1);
			drawCanvas.splits.add(0, swap);
			swap = null;
			drawCanvas.splits.get(0).offset = Main.sidebarWidth;
			drawCanvas.splits.get(1).offset = Main.sidebarWidth + drawCanvas.drawWidth;
			
			Main.drawCanvas.repaint();
			Main.chromDraw.drawCyto(drawCanvas.splits.get(0));
			Main.chromDraw.drawCyto(drawCanvas.splits.get(1));
			Main.chromDraw.repaint();
		*/
		}		
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
		else if(keyCode == KeyEvent.VK_F5) {
			if(!Main.bedCanvas.graph) {				
				Main.bedCanvas.graph = true;
			}
			else {
				Main.bedCanvas.graph = false;
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
		
		if(searchTable.containsKey(searchstring.replace(" ", "").toUpperCase())) {
			
			FileRead.search = true;
			String[] result = searchTable.get(searchstring.replace(" ", "").toUpperCase());
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

	
	
	

