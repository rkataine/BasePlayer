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

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.awt.AlphaComposite;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Font;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.geom.QuadCurve2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Timer;
import java.util.TimerTask;


import javax.swing.JPanel;

import javax.swing.SwingUtilities;

import base.BasePlayer.BedCanvas.bedFeatureFetcher;


public class Draw extends JPanel implements MouseMotionListener, MouseListener {	
	private static final long serialVersionUID = 1L;
	String readString ="";
	static int calc = 0;
	boolean annotationOn = false;
	QuadCurve2D curve = new QuadCurve2D.Double(0,0,40,40,0,0);
	public ArrayList<SplitClass> splits = new ArrayList<SplitClass>();
	MethodLibrary.mateListSorter mateSorter = new MethodLibrary.mateListSorter();
	static boolean calculateVars = true;
	int variantsStart =0, variantsEnd = 0;
	String varloadString = "Click here to load variants";
	String varloadhint = "Note: variant annotation is still possible through variant manager.";
	Rectangle varStartRect = new Rectangle(), varEndRect = new Rectangle();
	int varStringLen = 0;
	int clusterId = 1;
	int readInfoWidth = 200, varCalc = 0, indelCalc = 0,hoverInfoWidth = 200, removesamples =0;	
	
	int width=0, height=0, drawWidth=0;
	SAread saReads = null;
	SampleNode sampleOverLap = null;
	VarNode varOverLap = null;
	final Color chr1color = new Color(153,102,0), 
			chr2color = new Color(102,102,0),
			chr3color = new Color(153,153,30),
			chr4color = new Color(204,0,0),
			chr5color = new Color(255,0,0),
			chr6color = new Color(255,0,204),
			chr7color = new Color(255,204,204),
			chr8color = new Color(255,153,0),
			chr9color = new Color(255,204,0),
			chr10color = new Color(255,255,0),
			chr11color = new Color(204,255,0),
			chr12color = new Color(0,255,0),
			chr13color = new Color(53,128,0),
			chr14color = new Color(0,0,204),
			chr15color = new Color(102,153,255),
			chr16color = new Color(153,204,255),
			chr17color = new Color(0,255,255),
			chr18color = new Color(204,255,255),
			chr19color = new Color(153,0,204),
			chr20color = new Color(204,51,255),
			chr21color = new Color(204,153,255),
			chr22color = new Color(102,102,102),
			chrXcolor = new Color(153,153,153),
			chrYcolor = new Color(204,204,204),
			chrMTcolor = new Color(204,204,153);
	String[] clickedReadInfo = new String[11], hoverMateInfo = new String[11];
	//ArrayList<ReadNode> clickmates = null;
	final static int minzoom = 40;
	final static Color redColor = new Color(226, 37, 24), greenColor = new Color(37,226,24);
	ArrayList<String> readyQueue = new ArrayList<String>();
	Graphics2D g2v, buf, rbuf, sidebuf, /*logobuf,*/ cbuf, triImage;
	BufferedImage backbuffer, varbuffer, readbuffer, sidebuffer, logoImage, coveragebuffer, tribu;
	Composite backupv, backupb, backupr, backups, backupc;
	ArrayList<Sample> sampleList = new ArrayList<Sample>();
	ArrayList<ClusterNode> clusterNodes;
	final static BasicStroke strongStroke = new BasicStroke(4), doubleStroke = new BasicStroke(2), basicStroke = new BasicStroke(1), dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, ChromDraw.dash, 0.0f), doubledashed = new BasicStroke(2.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, ChromDraw.dash, 0.0f);
	Rectangle mouseRect = new Rectangle(), remoBox = new Rectangle(), hoverRect = new Rectangle();
	static final Color zoomColor = new Color(255,150,50,140), softLight = new Color(255,255,255,50), softColor = new Color(0,0,0,50), loadColor = new Color(255,255,200,200), forwardColor = new Color(171,194,171), reverseColor = new Color(194,171,171), reverseText= new Color(255,50,70), forwardText = new Color(255,255,255),reverseTextLow= new Color(150,100,100,200), forwardTextLow = new Color(220,255,240,200), dark = new Color(25,25,25);
	static Color backColor = new Color(90,90,90), curtainColor = new Color(70,70,60);
	Color varColor= Color.red;
	static Color sidecolor;
	static Color backTransparent = new Color(255,255,230,100);
	static boolean updatevars = false;
	public boolean loading = false;
	boolean first = true;
	Image leftTriangle = null;
	Image rightTriangle = null;
	static Font loadingFont = new Font("SansSerif", Font.BOLD, Main.defaultFontSize*2), defaultFont = new Font("SansSerif", Font.BOLD, Main.defaultFontSize), biggerFont = new Font("SansSerif", Font.BOLD, Main.defaultFontSize+2);	
	VarNode current, currentDraw;
	Iterator<Map.Entry<Short, SampleNode>> varIterator;	
	final static double defaultSampleHeight = 15;
	HashMap<Character,Color> baseColors = new HashMap<Character,Color>();
	Rectangle2D.Double rect;	
	String tempreadname = "";
	Boolean drawed = false;
	//static Image image;
	Sample selectedSample = null; 
	//String chrom = "";
	final Composite composite = AlphaComposite.getInstance(AlphaComposite.CLEAR,0.5f);
	int mouseX=0, mouseY=0, pressX=0, pressY=0, dragX=0,tempDrag =0, moveX=0,moveY=0;
	double starttemp, endtemp;	
	boolean mouseDrag = false, moveDrag = false, scrollbar = false, lineZoomer = false;
	int tempSample, selectedSampleIndex = -1, bottom = 0;
	String loadingtext = "";
	int loadbarAll = 0, loadBarSample= 0, varpos;
	static RenderingHints rh = new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
	static RenderingHints rhText = new RenderingHints(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_LCD_HRGB);
	
	boolean zoomDrag =false;
	private Sample sample;	
	private Rectangle drawRect;	
	private int readY;
	int selectedIndex = 0;
	private int readLevel;
	static boolean updateReads = false;
	public static boolean resize = false;
	private static int clustersum;
	private static VarNode tempnode;	
	boolean clusterCalc;
	private ReadNode row;	
	private boolean foundread;	
	private Sample readsample;	
	private Color readColor;
	private ReadNode read;
	private ReadNode selectedRead;	
	ArrayList<Entry<Integer, Byte>> mismatches;
	private CigarElement element;
	private int readpos;
	private Color tempColor;
	private Point p1;
	private Point p2;
	private Point p3;
	private int[] xs = new int[3];
	private int[] ys = new int[3];
	private Polygon triangle;
	private ReadNode clickedRead;
	private Rectangle clickedRect;
	private int readwheel;
	private boolean mouseWheel = false;
	private SampleNode samplenode;	
	private int secondSample = -1;
	boolean sidebar;	
	private Sample sidebarSample;
	private Sample removeSample;
	SplitClass selectedSplit = null;

	private ReadNode selectedMate;
	public long timer = 0;
	Timer navTimer = new Timer();
	TimerTask task = new MyTimerTask();
	private boolean splitremove;
	final static Color intronColor= new Color(170,170,160), interColor = new Color(140,140,140);
	private int prezoom;
	private double pressYrelative;
	boolean sampleZoomer;
	private boolean mousePress;
	private Entry<String, ArrayList<SampleNode>> entry;
//	private boolean foundexon;
	boolean resizeSidebar;	
	public DrawVariables drawVariables = new DrawVariables();
	private int remindex;
	public int ctrlpressed = 5;
	private Entry<Integer, Byte> mismatchentry;
	private Composite backupside;
	private Reads readHash;	
	private int coveragevariant;
	private int coveragevariantheight;
	private int coveragebottom;
	private String overlapText = "";
	private boolean overlap = false;
	private Color overlapColor = Color.black;
	private int prepixel;
	private Rectangle clickedMateRect;	
	private double tempViewLength;
	public boolean scrolldrag;
	private ReadNode hoverMate = null;
	private int readtextWidth;
	private int maxwidth = 0;
	private int mateadd = 0;
	private int previsiblestart = 0;
	
	private ArrayList<splitTuple> splitList = new ArrayList<splitTuple>();
	private int preread;
	private int thisread;
	private int insertion = 0;
	private int yValue;
	private double barHeight;
	private int maxcovtemp;
	private Sample clickedReadSample = null;
	private boolean readScrollDrag = false;
	private int moveYtemp = 0;
	private boolean founddiscordant;
	private int yposition;
	private Rectangle testRect;
	private double tempwidth;
	private boolean getMoreVariants = false;
	private int boxheight;
	HashMap<String, Float> advQualities = null;
	HashMap<String, Float> advQualitiesIndel = null;
	private Rectangle saReadRect = new Rectangle();
	private Color infoBoxColor = new Color(225,225,210,240);
	public boolean bam = false;
	private boolean zoomNote = false;
	private boolean foundvisible;
	public String baseHover = null;
	boolean forcereload = false;
	private int verticalScrollValue;
	
	private int overlapCoverage;
	private boolean readsidebar = false;
	private double heightFrac;
	private double hei;
	private boolean sampleDrag =false;
	private int hoverIndex = -1;
	private int sampleYpos;
	private int bottomYpos;
	private boolean sampleInfo;
	private int samePosCount = 1;
	public boolean intersect = false;
	private boolean bamHover;
	private boolean coverageregion;
	private int clipLength;
	private double sampleheight;
	private int totalheight;
	private int prestart = 0;
	
	private short previsiblesamples;
	
	static int clickedAdd = 0;
	static boolean updateCoverages = false;
	
	public static boolean variantcalculator =false;

	
public Draw(int width, int height) {
	
	super();
	this.width = width;
	this.height = height;
	
	this.drawWidth = (int)((this.width - Main.sidebarWidth)/(double)splits.size());
	addKeyListener(new KeyListener() {
		@Override
		public void keyTyped(KeyEvent e) {
			
		}

		@Override
		public void keyPressed(KeyEvent e) {			
			
			int keyCode = e.getKeyCode();
			if((e.getModifiers() & KeyEvent.CTRL_MASK) != 0) {
				return;
			}
			if(keyCode == KeyEvent.VK_PLUS || keyCode == 107) {
				
				if(selectedSplit.viewLength <= Settings.readDrawDistance && Main.drawCanvas.drawVariables.sampleHeight > 100 && Main.drawCanvas.sampleList.get(selectedIndex).getreadHash() != null) {
					double temp = (Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2)/(double)Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel;
					if(Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight < 20) {
						Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readfont =new Font("SansSerif", Font.BOLD, Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight);
						Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight++;
						Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel = (int)((Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2)/(double)temp);
						updateReads = true;
						repaint();
					}		
				}
				else {
				
					if(Draw.backColor.getAlpha() <= 250) {
				//		System.out.println(Draw.backColor.getAlpha()+5);
						backColor = new Color(backColor.getRed(),backColor.getGreen(),backColor.getBlue(),Draw.backColor.getAlpha()+5);
						
						repaint();
					}
				}
			}
			else if(keyCode == KeyEvent.VK_MINUS || keyCode == 109) {
				
				if(selectedSplit.viewLength <= Settings.readDrawDistance && Main.drawCanvas.drawVariables.sampleHeight > 100 && Main.drawCanvas.sampleList.get(selectedIndex).getreadHash() != null) {
				if(Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight > 1) {
					double temp = (Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2)/(double)Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel;
					Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readfont =new Font("SansSerif", Font.BOLD, Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight);
					
					Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight--;
					Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel = (int)((Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2)/(double)temp);
					updateReads = true;
					repaint();
				}
			}	
			else {
				
				if(Draw.backColor.getAlpha() >= 5) {
				//	System.out.println(Draw.backColor.getAlpha()-5);
					backColor = new Color(backColor.getRed(),backColor.getGreen(),backColor.getBlue(),Draw.backColor.getAlpha()-5);
					repaint();
				}
				
			}
			}
			if(keyCode == KeyEvent.VK_PLUS) {
				if(selectedSplit.viewLength <= Settings.readDrawDistance && Main.drawCanvas.drawVariables.sampleHeight > 100 && Main.drawCanvas.sampleList.get(selectedIndex).getreadHash() != null) {
				/*	double temp = (Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2)/(double)Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel;
					if(Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight < 20) {
						Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readfont =new Font("SansSerif", Font.BOLD, Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight);
						Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight++;
						Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel = (int)((Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2)/(double)temp);
						updateReads = true;
						repaint();
					}	*/	
				}
				else {
				
					if(Draw.backColor.getAlpha() <= 250) {
						
						backColor = new Color(60,60,60,Draw.backColor.getAlpha()+5);
						repaint();
					}
				}
			}
			
		}

		@Override
		public void keyReleased(KeyEvent e) {
			
		}
		
	});
	baseColors.put('A', Color.GREEN);
	baseColors.put('C', Color.BLUE);
	baseColors.put('G', Color.YELLOW);
	baseColors.put('T', Color.RED);
	
      
	//resizeCanvas(width, height);
	backbuffer = MethodLibrary.toCompatibleImage(new BufferedImage((int)Main.screenSize.getWidth(), (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
	buf = (Graphics2D)backbuffer.getGraphics();
	buf.setRenderingHints(rh);
	
	backupb = buf.getComposite();
		
	varbuffer = MethodLibrary.toCompatibleImage(new BufferedImage((int)Main.screenSize.getWidth(), (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
	g2v = (Graphics2D)varbuffer.getGraphics();
	g2v.setRenderingHints(rhText);
	backupv = g2v.getComposite();
	
	readbuffer = MethodLibrary.toCompatibleImage(new BufferedImage((int)Main.screenSize.getWidth(), (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
	rbuf = (Graphics2D)readbuffer.getGraphics();
	backupr = rbuf.getComposite();
	rbuf.setRenderingHints(rh);
	coveragebuffer = MethodLibrary.toCompatibleImage(new BufferedImage((int)Main.screenSize.getWidth(), (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
	cbuf = (Graphics2D)coveragebuffer.getGraphics();
	cbuf.setRenderingHints(rhText);
	backupc = cbuf.getComposite();
	
	sidebuffer = MethodLibrary.toCompatibleImage(new BufferedImage((int)Main.screenSize.getWidth(), (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));
	sidebuf = (Graphics2D)sidebuffer.getGraphics();
	backupside = sidebuf.getComposite();
	sidebuf.setRenderingHints(rhText);
	logoImage = MethodLibrary.toCompatibleImage(new BufferedImage(500, 100, BufferedImage.TYPE_INT_ARGB));
//	logobuf = (Graphics2D)logoImage.getGraphics();
//	logobuf.setRenderingHints(rh);
	addMouseMotionListener(this);
	addMouseListener(this);
	addMouseWheelListener(
			new MouseWheelListener() {				

				public void mouseWheelMoved(MouseWheelEvent mwe) {
					Main.drawCanvas.mouseWheel = true;
					if((sidebar || selectedSplit.viewLength > Settings.readDrawDistance) && Main.samples > 0) {
						
				//		clickedRead = null;
				//		clickmates = null;
						splitList.clear();
						selectedRead = null;					
					
						if ( mwe.getWheelRotation() < 0 ) {
							
							if(drawVariables.visiblestart > 0) {
								drawVariables.visiblestart--; //(short)(Main.drawScroll.getVerticalScrollBar().getValue()/drawVariables.sampleHeight-1);	
							//	drawVariables.visibleend--; // (short)(drawVariables.visiblestart + Main.drawScroll.getViewport().getHeight()/drawVariables.sampleHeight-1);	
								
							}
							setScrollbar((int)(drawVariables.visiblestart*drawVariables.sampleHeight));
						}
						else {
							
					
							if(drawVariables.visiblestart+drawVariables.visiblesamples < Main.samples) {
								drawVariables.visiblestart++; 		
							//	drawVariables.visibleend++; 	
					
							}
							setScrollbar((int)(drawVariables.visiblestart*drawVariables.sampleHeight));
							
						}
						
						
				
					
						if(splits.size() > 1) {
							for(int i = 0; i<splits.size(); i++) {
								splits.get(i).updateReads = true;		
							
							}
							updateCoverages = true;
						}
						else {
							
							Draw.updateReads = true;
							updateCoverages = true;
							Draw.updatevars = true;
							
						}
						
						repaint();
					}
					else {
					if(Main.samples > 0 && selectedSplit.viewLength < Settings.readDrawDistance) {
						if ( mwe.getWheelRotation() < 0 ) {							
							if(Main.drawCanvas.sampleList.get(selectedIndex).getreadHash() != null && Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().size() > 0) {
								if(Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel + (Main.drawCanvas.drawVariables.sampleHeight-(Main.drawCanvas.drawVariables.sampleHeight/selectedSplit.getDivider())) < (Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getReadSize()*(Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2))) {
									Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel+=(Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2)*ctrlpressed;
									selectedSplit.updateReads = true;									
									
								}						
							}						
						}
						else {
							if(Main.drawCanvas.sampleList.get(selectedIndex).getreadHash() != null && Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().size() > 0) {
								Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel -=(Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2)*ctrlpressed;
								selectedSplit.updateReads = true;
								if(Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel < 0) {
									Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel = 0;
								}							
							}					
						}
						
						
						repaint();
					}
				}
			}
		}
	);
}

void drawCoverages(SplitClass split) {
	
	prepixel = -1;
	try {
		maxcovtemp = 0;
		
	for(int j = 0; j<readHash.getCoverages().length; j++) {
		
		if(j+readHash.getCoverageStart() < split.start-1) {			
			continue;
		}							
		if(j+readHash.getCoverageStart() > split.end) {
			break;
		}
		if(maxcovtemp < readHash.getCoverages()[j][0]) {
			maxcovtemp = (int)readHash.getCoverages()[j][0];
		}
		
		if(prepixel < 0 || prepixel < (int)((j*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel))) {
			
			split.getReadBuffer().fillRect((int)((j*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)),(int)(coveragebottom-(readHash.getCoverages()[j][0]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider()))+1), (int)split.pixel+1,(int)(readHash.getCoverages()[j][0]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider())));
			prepixel = (int)((j*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel));
		}
			
		if(readHash.getCoverages()[j][1] >= Settings.coverageAlleleFreq || readHash.getCoverages()[j][2] >= Settings.coverageAlleleFreq || readHash.getCoverages()[j][3] >= Settings.coverageAlleleFreq || readHash.getCoverages()[j][4] >= Settings.coverageAlleleFreq || readHash.getCoverages()[j][5] >= Settings.coverageAlleleFreq || readHash.getCoverages()[j][6] >= Settings.coverageAlleleFreq || readHash.getCoverages()[j][7] >= Settings.coverageAlleleFreq) {
			coveragevariantheight = 0;
			coveragevariant = j;
			
		}
		if(moveX-split.offset >= (int)(((j)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)) && moveX-split.offset < (int)(((j)*split.pixel)+((readHash.getCoverageStart()-split.start+1)*split.pixel))) {
			overlapCoverage = (int)readHash.getCoverages()[j][0];
		}
		if(coveragevariant != -1) { //&& (int)(((j+1)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)) > (int)(((j)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel))) {
			if(readHash.getCoverages()[coveragevariant][1] >= Settings.coverageAlleleFreq) {
				
				split.getReadBuffer().setColor(Color.green);	
				coveragevariantheight = (int)(readHash.getCoverages()[coveragevariant][1]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider()));
				split.getReadBuffer().fillRect((int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)),coveragebottom-coveragevariantheight, (int)split.pixel+1,(int)(readHash.getCoverages()[coveragevariant][1]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider())));
				
				if(moveX-split.offset >= (int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)) && moveX-split.offset < (int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start+1)*split.pixel)) && moveY-Main.drawScroll.getVerticalScrollBar().getValue() < coveragebottom) {
					
					if(!overlapText.equals("A : " +(int)readHash.getCoverages()[coveragevariant][1] +" ("+(int)(MethodLibrary.round(readHash.getCoverages()[coveragevariant][1]/(double)readHash.getCoverages()[coveragevariant][0],2)*100) +"%)")) {
						overlapText = "A : " +(int)readHash.getCoverages()[coveragevariant][1] +" (" +(int)(MethodLibrary.round(readHash.getCoverages()[coveragevariant][1]/(double)readHash.getCoverages()[coveragevariant][0],2)*100) +"%)";
						
					}										
					overlap = true;	
					overlapColor = Color.green;
				}
				
				
			}
			if(readHash.getCoverages()[coveragevariant][2] >= Settings.coverageAlleleFreq) {
				if(moveX-split.offset >= (int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)) && moveX-split.offset < (int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start+1)*split.pixel)) && moveY-Main.drawScroll.getVerticalScrollBar().getValue() < coveragebottom-coveragevariantheight) {
					if(!overlapText.equals("C : " +(int)readHash.getCoverages()[coveragevariant][2] +" ("+(int)(MethodLibrary.round(readHash.getCoverages()[coveragevariant][2]/(double)readHash.getCoverages()[coveragevariant][0],2)*100) +"%)")) {
						overlapText = "C : "+(int)readHash.getCoverages()[coveragevariant][2] +" (" +(int)(MethodLibrary.round(readHash.getCoverages()[coveragevariant][2]/(double)readHash.getCoverages()[coveragevariant][0],2)*100) +"%)";
					}										
					overlap = true;		
					overlapColor = Color.cyan;
				}
				coveragevariantheight += (int)(readHash.getCoverages()[coveragevariant][2]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider()));
				split.getReadBuffer().setColor(Color.cyan);								
				split.getReadBuffer().fillRect((int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)),coveragebottom-coveragevariantheight, (int)split.pixel+1,(int)(readHash.getCoverages()[coveragevariant][2]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider())));
					
				
			}
			if(readHash.getCoverages()[coveragevariant][3] >= Settings.coverageAlleleFreq) {
				if(moveX-split.offset >= (int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)) && moveX-split.offset < (int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start+1)*split.pixel)) && moveY-Main.drawScroll.getVerticalScrollBar().getValue() < coveragebottom-coveragevariantheight) {
					
					if(!overlapText.equals("G : " +(int)readHash.getCoverages()[coveragevariant][3] +" ("+(int)(MethodLibrary.round(readHash.getCoverages()[coveragevariant][3]/(double)readHash.getCoverages()[coveragevariant][0],2)*100) +"%)")) {
						overlapText = "G : " +(int)readHash.getCoverages()[coveragevariant][3] +" ("+(int)(MethodLibrary.round(readHash.getCoverages()[coveragevariant][3]/(double)readHash.getCoverages()[coveragevariant][0],2)*100) +"%)";
					}										
					overlap = true;		
					overlapColor = Color.orange;
				}
				split.getReadBuffer().setColor(Color.orange);
				
				coveragevariantheight += (int)(readHash.getCoverages()[coveragevariant][3]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider()));
				
												
				split.getReadBuffer().fillRect((int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)),coveragebottom-coveragevariantheight, (int)split.pixel+1,(int)(readHash.getCoverages()[coveragevariant][3]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider())));
				
			}
			if(readHash.getCoverages()[coveragevariant][4] >= Settings.coverageAlleleFreq) {
				
				if(moveX-split.offset >= (int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)) && moveX-split.offset < (int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start+1)*split.pixel)) && moveY-Main.drawScroll.getVerticalScrollBar().getValue() < coveragebottom-coveragevariantheight) {
					
					if(!overlapText.equals("T : " +(int)readHash.getCoverages()[coveragevariant][4] +" ("+(int)(MethodLibrary.round(readHash.getCoverages()[coveragevariant][4]/(double)readHash.getCoverages()[coveragevariant][0],2)*100) +"%)")) {
						
						overlapText = "T : " +(int)readHash.getCoverages()[coveragevariant][4] +" ("+(int)(MethodLibrary.round(readHash.getCoverages()[coveragevariant][4]/(double)readHash.getCoverages()[coveragevariant][0],2)*100) +"%)";
					}										
					overlap = true;		
					overlapColor = Color.red;
				}	
				coveragevariantheight += (int)(readHash.getCoverages()[coveragevariant][4]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider()));
				
				split.getReadBuffer().setColor(Color.red);								
				split.getReadBuffer().fillRect((int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)),coveragebottom-coveragevariantheight, (int)split.pixel+1,(int)(readHash.getCoverages()[coveragevariant][4]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider())));
				
			}
			if(readHash.getCoverages()[coveragevariant][5] >= Settings.coverageAlleleFreq) {
				
				if(moveX-split.offset >= (int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)) && moveX-split.offset < (int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start+1)*split.pixel)) && moveY-Main.drawScroll.getVerticalScrollBar().getValue() < coveragebottom-coveragevariantheight) {
					
					if(!overlapText.equals("N : "+(int)readHash.getCoverages()[coveragevariant][5] +" (" +(int)(MethodLibrary.round(readHash.getCoverages()[coveragevariant][5]/(double)readHash.getCoverages()[coveragevariant][0],2)*100) +"%)")) {
						overlapText = "N : " +(int)readHash.getCoverages()[coveragevariant][5] +" ("+(int)(MethodLibrary.round(readHash.getCoverages()[coveragevariant][5]/(double)readHash.getCoverages()[coveragevariant][0],2)*100) +"%)";
					}										
					overlap = true;		
					overlapColor = Color.gray;
				}
				coveragevariantheight += (int)(readHash.getCoverages()[coveragevariant][5]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider()));
				
				split.getReadBuffer().setColor(Color.gray);								
				split.getReadBuffer().fillRect((int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)),coveragebottom-coveragevariantheight, (int)split.pixel+1,(int)(readHash.getCoverages()[coveragevariant][5]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider())));
				
			}
			if(readHash.getCoverages()[coveragevariant][6] >= Settings.coverageAlleleFreq) {
				
				if(moveX-split.offset >= (int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)) && moveX-split.offset < (int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start+1)*split.pixel)) && moveY-Main.drawScroll.getVerticalScrollBar().getValue() < coveragebottom-coveragevariantheight) {
					
					if(!overlapText.equals("INS : "+(int)readHash.getCoverages()[coveragevariant][6] +" (" +(int)(MethodLibrary.round(readHash.getCoverages()[coveragevariant][6]/(double)readHash.getCoverages()[coveragevariant][0],2)*100) +"%)")) {
						overlapText = "INS : " +(int)readHash.getCoverages()[coveragevariant][6] +" ("+(int)(MethodLibrary.round(readHash.getCoverages()[coveragevariant][6]/(double)readHash.getCoverages()[coveragevariant][0],2)*100) +"%)";
					}					
					
					overlap = true;		
					overlapColor = Color.white;
				}
				coveragevariantheight += (int)(readHash.getCoverages()[coveragevariant][6]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider()));
				
				split.getReadBuffer().setColor(Color.white);								
				split.getReadBuffer().fillRect((int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)),coveragebottom-coveragevariantheight, (int)split.pixel+1,(int)(readHash.getCoverages()[coveragevariant][6]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider())));
				
			}
			
			if(readHash.getCoverages()[coveragevariant][7] >= Settings.coverageAlleleFreq) {
				if(moveX-split.offset >= (int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)) && moveX-split.offset < (int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start+1)*split.pixel)) && moveY-Main.drawScroll.getVerticalScrollBar().getValue() < coveragebottom-((readHash.getCoverages()[coveragevariant][0])/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider()))) {
					if(!overlapText.equals("DEL : "+(int)readHash.getCoverages()[coveragevariant][7] +" (" +(int)(MethodLibrary.round(readHash.getCoverages()[coveragevariant][7]/(double)(readHash.getCoverages()[coveragevariant][0]+readHash.getCoverages()[coveragevariant][7]),2)*100) +"%)")) {
						overlapText = "DEL : " +(int)readHash.getCoverages()[coveragevariant][7] +" ("+(int)(MethodLibrary.round(readHash.getCoverages()[coveragevariant][7]/(double)(readHash.getCoverages()[coveragevariant][0]+readHash.getCoverages()[coveragevariant][7]),2)*100) +"%)";
					}					
					
					overlap = true;		
					overlapColor = Color.white;
				}
				//coveragevariantheight += (int)(readHash.getCoverages()[coveragevariant][6]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider()));
				
			//	split.getReadBuffer().setColor(Color.white);								
			//	split.getReadBuffer().fillRect((int)(((coveragevariant)*split.pixel)+((readHash.getCoverageStart()-split.start)*split.pixel)),coveragebottom-coveragevariantheight, (int)split.pixel+1,(int)(readHash.getCoverages()[coveragevariant][6]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider())));
				
			}
			
			
			split.getReadBuffer().setColor(Color.black);
			coveragevariant = -1;
		}
	//	split.getReadBuffer().fillRect((int)((j*split.pixel)+((readHash.getReadStart()-split.start)*split.pixel)),(int)((readsample.getIndex()*drawVariables.sampleHeight)+(drawVariables.sampleHeight/split.getDivider())-(int)(coverages[j][1]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider())))-Main.drawScroll.getVerticalScrollBar().getValue(), (int)split.pixel+1,(int)(coverages[j][1]/(double)readHash.getMaxcoverage()*(drawVariables.sampleHeight/split.getDivider())));
	
	}
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	if((int)readHash.getMaxcoverage() != maxcovtemp) {
		readHash.setMaxcoverage(maxcovtemp);
	}
}

void drawCoverage(SplitClass split) {
	
	if(split.viewLength <= Settings.coverageDrawDistance && drawVariables.sampleHeight > 100) {		
	
	if(!updateReads && !FileRead.cancelreadread) {
	
		/*if(splitIndex != splits.indexOf(split)) {
			splitIndex = splits.indexOf(split);
		}
		*/
		for(int i =drawVariables.visiblestart; i<drawVariables.visiblestart+drawVariables.visiblesamples;i++ ) {
			
			try {	
				if(split.getDivider() == null || split.getDivider() != 1 && split.viewLength > Settings.readDrawDistance) {
					split.setDivider(1.0);
				}
				else if(split.getDivider() != 5.0 && split.viewLength <= Settings.readDrawDistance) {
					split.setDivider(5.0);
				}
				if(i >= Main.samples) {
					break;
				}
				if(sampleList.get(i).samFile == null) {			
					continue;
				}
							
			
				if(i < (int)(Main.drawScroll.getVerticalScrollBar().getValue()/drawVariables.sampleHeight)-1) {				
					continue;
				}
				if(i > (int)(Main.drawScroll.getVerticalScrollBar().getValue()/drawVariables.sampleHeight)+Main.drawScroll.getViewport().getHeight()/drawVariables.sampleHeight) {
					break;
				}		
				readsample = sampleList.get(i);
				
				split.getReadBuffer().setComposite( composite);					
			//	split.getReadBuffer().fillRect(0,0, split.getReadImage().getWidth(),50);	
				split.getReadBuffer().fillRect(0, (int)((readsample.getIndex())*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue(), this.drawWidth, (int)(drawVariables.sampleHeight/split.getDivider()));			
				split.getReadBuffer().setComposite(split.getBackupr());	
				split.getReadBuffer().setColor(Color.black);
				coveragebottom = (int)((readsample.getIndex()*drawVariables.sampleHeight)+(drawVariables.sampleHeight/split.getDivider()))-Main.drawScroll.getVerticalScrollBar().getValue();
				split.getReadBuffer().fillRect(0, coveragebottom, this.drawWidth,2);
				if(split.getDivider() != 1) {
					split.getReadBuffer().setColor(backTransparent );		
					split.getReadBuffer().fillRect(0, (int)((readsample.getIndex())*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue(), this.drawWidth, (int)(drawVariables.sampleHeight/split.getDivider())-2);	
				}
				if(readsample.getreadHash() != null && readsample.getreadHash().size() < splits.size()) {
					readsample.resetreadHash();
				}
				readHash = readsample.getreadHash().get(split);
				if(readHash != null && readHash.getCoverages() != null) {					
					
					split.getReadBuffer().setColor(Color.black);
					if(split.viewLength > Settings.readDrawDistance) {
						if((lineZoomer && !split.equals(selectedSplit)) || !lineZoomer) {
							
							for(int j = (int)((readHash.getCoverageStart()-split.start)*split.pixel); j<readHash.getCoverages().length; j++) {
								if(j < 0) {
									continue;
								}
								if((int)((readHash.getCoverageStart()-split.start)*split.pixel)+j > this.getWidth()) {								
									break;
								}		
								split.getReadBuffer().fillRect(j+(int)((readHash.getCoverageStart()-split.start)*split.pixel),coveragebottom-(int)((readHash.getCoverages()[j][0]/(double)readHash.getMaxcoverage())*(drawVariables.sampleHeight/split.getDivider())), 1,(int)((readHash.getCoverages()[j][0]/(double)readHash.getMaxcoverage())*(drawVariables.sampleHeight/split.getDivider())));
																
							}
							
							if(moveY > readsample.getIndex()*drawVariables.sampleHeight && moveY < readsample.getIndex()*drawVariables.sampleHeight+(drawVariables.sampleHeight/split.getDivider()) ) {
								if(!split.getReadBuffer().getFont().equals(Main.menuFont)) {
									split.getReadBuffer().setFont(Main.menuFont);	
								}
								split.getReadBuffer().setColor(Color.white);
								split.getReadBuffer().drawString(""+(int)((readsample.getIndex()*drawVariables.sampleHeight+(drawVariables.sampleHeight/split.getDivider()-2)-moveY)*(readHash.getMaxcoverage()*split.pixel/(drawVariables.sampleHeight/split.getDivider()))+1), moveX-split.offset+20, moveY-Main.drawScroll.getVerticalScrollBar().getValue()+10);
							}					
						}
					}
					else {
					
						coveragevariant = -1;
						overlap = false;
						
						drawCoverages(split);
						
						if(moveX > split.offset && moveY > readsample.getIndex()*drawVariables.sampleHeight && moveY < readsample.getIndex()*drawVariables.sampleHeight+(drawVariables.sampleHeight/split.getDivider())) {
							
							split.getSelectbuf().setColor(Color.black);
							split.getSelectbuf().fillRoundRect(moveX-split.offset+20, moveY-Main.drawScroll.getVerticalScrollBar().getValue(), (int)(Main.chromDraw.bounds*10+(Math.log10(overlapCoverage)*Main.chromDraw.bounds)), Main.defaultFontSize*4,6,6);
							split.getSelectbuf().setColor(Color.white);
							
							split.getSelectbuf().drawString("Pointer: "+(int)((readsample.getIndex()*drawVariables.sampleHeight+(drawVariables.sampleHeight/split.getDivider()-2)-moveY)*(readHash.getMaxcoverage()/(drawVariables.sampleHeight/split.getDivider()))+1), moveX-split.offset+26, moveY-Main.drawScroll.getVerticalScrollBar().getValue()+Main.defaultFontSize*3+8);
							split.getSelectbuf().drawString("Coverage: " +overlapCoverage, moveX-split.offset+26, moveY-Main.drawScroll.getVerticalScrollBar().getValue()+Main.defaultFontSize);
							
							if(overlap) {
								
								split.getSelectbuf().setColor(overlapColor);
								
								split.getSelectbuf().drawString(this.overlapText, moveX-split.offset+26, moveY-Main.drawScroll.getVerticalScrollBar().getValue()+Main.defaultFontSize*2+4);
								
							}
						}
						
					}
				}
			}
			catch(Exception e) {
				e.printStackTrace();
			}			
		
		}
	}
	}
	else {
		if(!split.clearedReads) {
			if(sampleList.size() > 0) {
				/*if(splitIndex != splits.indexOf(split)) {
					splitIndex = splits.indexOf(split);
				}
				*/
				for(int i =drawVariables.visiblestart; i<drawVariables.visiblestart+drawVariables.visiblesamples;i++ ) {
					if(sampleList.get(i).removed) {
						continue;
					}
					readsample = sampleList.get(i);					
					readHash = readsample.getreadHash().get(split);
					
					if(readHash.getCoverageStart()!= 0) {
						readHash.setCoverages(null);
						readHash.setCoverageEnd(0);
						readHash.setCoverageStart(Integer.MAX_VALUE);
					}					
				}
				split.clearedCoverages = true;
				}
		}		
	}
}

/*
public void drawInserts() {	
	
	if(this.insertList.size() > 0) {
		
		buf.setStroke(basicStroke);
		for(int i =0;i<this.insertList.size(); i++) {
			buf.setColor(Color.red);
			if(this.insertList.get(i)[0] > -1) buf.drawLine((int)((this.insertList.get(i)[0]-start)*selectedSplit.pixel),0, (int)((this.insertList.get(i)[0]-start)*selectedSplit.pixel), getHeight());
			
			if(this.insertList.get(i)[1] > -1) buf.drawLine((int)((this.insertList.get(i)[1]-start)*selectedSplit.pixel),0, (int)((this.insertList.get(i)[1]-start)*selectedSplit.pixel), getHeight());
			
			if(this.insertList.get(i)[0] > -1 && this.insertList.get(i)[1] > -1){
				buf.setColor(zoomColor);
				if(this.insertList.get(i)[0] < this.insertList.get(i)[1]) {
					buf.fillRect((int)((this.insertList.get(i)[0]-start)*selectedSplit.pixel), 0,(int)((((this.insertList.get(i)[1]-start))-((this.insertList.get(i)[0]-start)))*selectedSplit.pixel), getHeight() );
			
				}
				else {
					buf.fillRect((int)((this.insertList.get(i)[1]-start)*selectedSplit.pixel), 0,(int)((((this.insertList.get(i)[0]-start))-((this.insertList.get(i)[1]-start)))*selectedSplit.pixel), getHeight() );
					
				}
			}
		}
		
	}
	buf.setColor(Color.black);
}
*/
public int getDrawWidth() {
	return this.drawWidth;
}
public int getWidth() {
	return this.width;
}

public void drawVars(int offset) {	
	
	if(updatevars || splits.get(0).pixel > 3) {
		
	//	Long start = System.currentTimeMillis();
		if(!lineZoomer && !mouseDrag && (!loading || Main.drawCanvas.loadingtext.contains("Calling")) && !scrolldrag) {	
			if(splits.get(0).pixel > 3) {
				
			}
			else {
				calculateVars = true;
				varCalc = 0;	
				indelCalc = 0;
			}
		}
		else {			
			
			calculateVars = false;
		}
		
		for(int s = drawVariables.visiblestart; s < drawVariables.visiblestart+drawVariables.visiblesamples+removesamples+1; s++) {
			if(s < 0) {
				continue;
			}
			if(s>this.sampleList.size()-1) {
				break;
			}
			if(sampleList.get(s).removed) {
				continue;
			}
			if(calculateVars ) {
				this.sampleList.get(s).varcount = 0;
				this.sampleList.get(s).homozygotes = 0;
				this.sampleList.get(s).heterozygotes = 0;
				this.sampleList.get(s).indels = 0;
				this.sampleList.get(s).sitioRate = 0;
				this.sampleList.get(s).versioRate = 0;
			}
			this.sampleList.get(s).prepixel = -1;			
		}
		
		g2v.setComposite( composite);		
		g2v.fillRect(0,0, this.getWidth(), Main.drawScroll.getViewport().getHeight());	
		g2v.setComposite(backupv);
		
		if(current == null && Main.chromDraw.vardraw != null) {
			Main.chromDraw.vardraw = null;
		}
		
		varOverLap = null;
		baseHover = null;
		
		if(FileRead.head.getNext() != null && current != null) {
			
			try {
				
				if(current.getPosition()+1 < this.splits.get(0).start && current.getNext() != null) {
					if(current.getNext().getPosition()+1 < this.splits.get(0).start) {
						while(current.getPosition()+1 < this.splits.get(0).start && current.getNext() != null) {
							if(FileRead.novars) {								
								if(current.getNext().getChrom().equals(splits.get(0).chrom)) {
									current = current.getNext();
								}
								else {
									break;
								}
							}
							else {
								current = current.getNext();
							}													
						}
					}
				}
				else {
					while(current.getPosition()-1 > this.splits.get(0).start) {		
						if(FileRead.novars) {
							if(current.getPrev().getChrom() != null && current.getPrev().getChrom().equals(splits.get(0).chrom)) {
								current = current.getPrev();
							}
							else {
								break;
							}
						}
						else {
							current = current.getPrev();
						}
						if(current == null) {
							break;
						}
					}
				}		
				
			}
			catch(Exception e) {
				ErrorLog.addError(e.getStackTrace());
				e.printStackTrace();
			}
			
			if(current == null) {
				updatevars = false;
				return;
			}
			
			currentDraw = current;
		
			if(Main.chromDraw.varnode == null || !Main.chromDraw.varnode.equals(current)) {			
				Main.chromDraw.varnode = current;			
			}
		try {
				heightFrac = 3*drawVariables.sampleHeight/(double)4;
				hei = drawVariables.sampleHeight/(double)4;
				verticalScrollValue = Main.drawScroll.getVerticalScrollBar().getValue();
				while(currentDraw != null) {
					
					
					if(currentDraw.getPosition() < this.splits.get(0).start) {
						
						currentDraw = currentDraw.getNext();	
						continue;
					}
					if(currentDraw.getPosition() > this.splits.get(0).end) {
						break;
					}		
					
					if(hideNode(currentDraw)) {						
						currentDraw = currentDraw.getNext();
						continue;
					}
					
					if(currentDraw.getExons() != null) {			
						
						if(!currentDraw.coding) {
							if(!g2v.getColor().equals(intronColor)) {
								varColor = intronColor;
							}									
						}
						else {							
							if(!g2v.getColor().equals(redColor)) {
								varColor = redColor;					
							}						
						}										
					}
					else {		
						if(currentDraw.isInGene()) {
							if(!g2v.getColor().equals(intronColor)) {
								varColor = intronColor;
							}							
						}
						else {
							if(!g2v.getColor().equals(interColor)) {
								varColor = interColor;
							}							
						}
					}
					if(!g2v.getColor().equals(varColor)) {
						g2v.setColor(varColor);
					}					
						
						varpos = (int)((currentDraw.getPosition()+1-this.splits.get(0).start)*splits.get(0).pixel);
						
						for(int v = 0; v<currentDraw.vars.size(); v++) {
							if(currentDraw == null) {
								break;
							}
							entry = currentDraw.vars.get(v);
							
							if(hideNodeVar(currentDraw,entry)) {
								continue;
							}
						
							for(int i = 0; i<entry.getValue().size(); i++) {
								
								try {
									if(entry.getValue().get(i).allelenumber != null) {
										break;
									}
									sample = entry.getValue().get(i).getSample();
									/*
									if(sample == null) {
										continue;
									}*/
									
									if(sample.getIndex() < drawVariables.visiblestart) {											
										continue;
									}
								
									if(sample.getIndex() > drawVariables.visiblestart+drawVariables.visiblesamples+1) {										
										break;
									}
																		
									
									if(hideVar(entry.getValue().get(i), entry.getKey().length() > 1)) {									
										continue;
									}
									samplenode = entry.getValue().get(i);	
									if(currentDraw.indel && currentDraw.getExons() != null) { //&& currentDraw.getPosition()+1 > currentDraw.getExons().get(0).getTranscript().getCodingStart() && currentDraw.getPosition()+1 < currentDraw.getExons().get(0).getTranscript().getCodingEnd()) {
										
										if(!currentDraw.coding) {
											if(!g2v.getColor().equals(intronColor)) {
												varColor = intronColor;
											}
											
										}
										else {
											if( entry.getKey().length() > 1) {
												if(!g2v.getColor().equals( greenColor)) {
													varColor = greenColor;
												}
											}
											else {
												if(!g2v.getColor().equals( Color.red)) {
													varColor = Color.red;
												}
											}
											
										}
										if(!g2v.getColor().equals(varColor)) {
											g2v.setColor(varColor);
										}
									}
								
									if(calculateVars) {
										
										if(VariantHandler.indelFilters.isSelected()) {
											
											if(entry.getKey().length() > 1) {											
												indelCalc++;
											}
											else {
												varCalc++;	
											}
										}
										else {
											varCalc++;	
											if(entry.getKey().length() > 1) {											
												indelCalc++;
											}
										}
										if(drawVariables.sampleHeight > 30) {
											sample.varcount++;
											if(entry.getKey().length() > 1) {
												sample.indels++;
												indelCalc++;
											}
											if(drawVariables.sampleHeight > 45) {
												if(samplenode.isHomozygous()) {
													sample.homozygotes++;
												}
												else {
													sample.heterozygotes++;
												}
												if(drawVariables.sampleHeight > 75) {
													if(!currentDraw.indel || entry.getKey().length() < 2) {
														if(((char)currentDraw.getRefBase() == ('A') && entry.getKey().equals("G")) || ((char)currentDraw.getRefBase() == 'G' && entry.getKey().equals("A")) || ((char)currentDraw.getRefBase() == 'C' && entry.getKey().equals("T")) || ((char)currentDraw.getRefBase() == 'T' && entry.getKey().equals("C"))) {
															sample.sitioRate++;
														}
														else {
															sample.versioRate++;
														}														
													}														
												}													
											}												
										}
									}
								
									if(splits.get(0).viewLength > Settings.readDrawDistance) {
										
											if(sample.getMaxCoverage() < VariantHandler.maxSlideValue) {
												
												if(varpos != sample.prepixel/* && currentDraw.getExons() == null*/) {													
													g2v.drawLine(varpos, (int)(((sample.getIndex()+1)*drawVariables.sampleHeight-hei)-samplenode.heightValue*heightFrac/sample.getMaxCoverage())-verticalScrollValue, varpos, (int)(((sample.getIndex())*drawVariables.sampleHeight)+(drawVariables.sampleHeight))-verticalScrollValue);
												}											
												/*
												if(currentDraw.getExons() != null) {
													g2v.drawLine(varpos, (int)(((sample.getIndex()+1)*drawVariables.sampleHeight-hei)-samplenode.heightValue*heightFrac/sample.getMaxCoverage())-verticalScrollValue, varpos,  (int)(((sample.getIndex())*drawVariables.sampleHeight)+(drawVariables.sampleHeight))-verticalScrollValue);
													
												}*/
											}
											else {
												
												if(varpos != sample.prepixel/* && currentDraw.getExons() == null*/) {													
													g2v.drawLine(varpos, (int)(((sample.getIndex()+1)*drawVariables.sampleHeight-drawVariables.sampleHeight/(double)4)-samplenode.heightValue*((3*drawVariables.sampleHeight/(double)4)/VariantHandler.maxSlideValue))-verticalScrollValue, varpos, (int)(((sample.getIndex())*drawVariables.sampleHeight)+(drawVariables.sampleHeight))-verticalScrollValue);
													
												}											
												
												/*if(currentDraw.getExons() != null) {
													
													g2v.drawLine(varpos, (int)(((sample.getIndex()+1)*drawVariables.sampleHeight-drawVariables.sampleHeight/(double)4)-samplenode.heightValue*((3*drawVariables.sampleHeight/(double)4)/VariantHandler.maxSlideValue))-verticalScrollValue, varpos,  (int)(((sample.getIndex())*drawVariables.sampleHeight)+(drawVariables.sampleHeight))-verticalScrollValue);
													
												}*/
											}
											if(varpos != sample.prepixel) {
												sample.prepixel = (short)varpos;
												
											}
									}
									else {
										 
										  if(entry.getKey().length() > 1) {
											  varpos = (int)((currentDraw.getPosition()+2-this.splits.get(0).start)*splits.get(0).pixel);
										  }
											
										  if(moveX <=varpos+(int)splits.get(0).pixel+1+splits.get(0).offset && moveX >=varpos+splits.get(0).offset) {
											  if(moveY >=(int)(((sample.getIndex())*drawVariables.sampleHeight)) && moveY <=(int)(((sample.getIndex())*drawVariables.sampleHeight))+(int)drawVariables.sampleHeight) {
												  sampleOverLap = samplenode;
												  baseHover = entry.getKey();
												  varOverLap = currentDraw;
												  g2v.setColor(Color.white);
											  }
										  }
										  if(varpos == sample.prepixel) {
											  samePosCount++;
											 // g2v.fillRect(varpos, (int)(((sample.getIndex())*drawVariables.sampleHeight)-drawVariables.sampleHeight/2)-verticalScrollValue, (int)splits.get(0).pixel+1, (int)drawVariables.sampleHeight/4);
											 
										  }
										  else {
											  samePosCount = 1;
											  g2v.fillRect(varpos, (int)(((sample.getIndex())*drawVariables.sampleHeight))-verticalScrollValue, (int)splits.get(0).pixel+1, (int)drawVariables.sampleHeight);
										  }
									      if(this.splits.get(0).viewLength < 200) {
									    	g2v.setColor(Color.white);
									    	g2v.setFont(ChromDraw.seqFont);		
									    	 if(varpos == sample.prepixel) {
									    		
									    		 g2v.drawLine(varpos,(int)(((sample.getIndex())*drawVariables.sampleHeight)+drawVariables.sampleHeight/samePosCount)-verticalScrollValue,(int)(varpos+splits.get(0).pixel+1), (int)(((sample.getIndex())*drawVariables.sampleHeight)+drawVariables.sampleHeight/samePosCount)-verticalScrollValue);
									    		 g2v.drawString(entry.getKey(), varpos, (int)(((sample.getIndex())*drawVariables.sampleHeight+Main.defaultFontSize+2)+drawVariables.sampleHeight/samePosCount)-verticalScrollValue);
											    	
									    	 }
									    	 else {
									    		 g2v.drawString(entry.getKey(), varpos, (int)(((sample.getIndex())*drawVariables.sampleHeight+Main.defaultFontSize+2))-verticalScrollValue);
											    	
									    	 }
									    	if(currentDraw.coding) {
									    		if(!g2v.getColor().equals(Color.red)) {
													varColor = Color.red;
												}
											}
											else {
												
												if(currentDraw.isInGene()) {	
													if(!g2v.getColor().equals(intronColor)) {
														varColor = intronColor;
													}
												
												}
												else {		
													if(!g2v.getColor().equals(interColor)) {
														varColor = interColor;
													}
													
												}
											}
									    	if(!g2v.getColor().equals(varColor)) {
												g2v.setColor(varColor);
											}
									    	g2v.setFont(defaultFont);	
									    	
									    }
									      
								        if(entry.getKey().length() > 1) {
										   varpos = (int)((currentDraw.getPosition()+1-this.splits.get(0).start)*splits.get(0).pixel);
									    }
								        if(g2v.getColor().equals(Color.white)) {
								        	g2v.setColor(varColor);
								        	
								        }
								        if(varpos != sample.prepixel) {
											sample.prepixel = (short)varpos;
											
										}
									}											
								}
								catch(Exception r) {
									ErrorLog.addError(r.getStackTrace());
									r.printStackTrace();
									break;
								}									
							}
						}
						
					currentDraw = currentDraw.getNext();					
					
				}
			}
			catch(Exception e) {
				updatevars = false;
				ErrorLog.addError(e.getStackTrace());
				e.printStackTrace();
			}
		
				if(calculateVars) {	
					if(VariantHandler.indelFilters.isSelected()) {
						VariantHandler.totalVars.setText("SNV count on screen: " +varCalc);
						VariantHandler.totalVars2.setText("SNV count on screen: " +varCalc);
						VariantHandler.totalVarsIndel.setText("Indel count on screen: " +indelCalc);
						
					}
					else {
						VariantHandler.totalVars.setText("Variant count on screen: " +varCalc);
						VariantHandler.totalVars2.setText("Variant count on screen: " +varCalc);
						VariantHandler.totalVarsIndel.setText("Indel count on screen: " +indelCalc);
					}
					/*if(varCalc == 0 && variantsEnd != 1 && !VariantHandler.allChroms.isSelected()) {
						VariantHandler.varcalc.setEnabled(false);
						VariantHandler.varcalc.setToolTipText("No variants on screen.");
						
					}
					else {
						VariantHandler.varcalc.setEnabled(true);
						VariantHandler.varcalc.setToolTipText("Annotate variants.");
						Main.manage.setToolTipText("Manage variants.");
						Main.manage.setEnabled(true);
					}*/
				}
				if(varOverLap != null) {
					if(getCursor().getType() != Cursor.HAND_CURSOR) {					
						setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
					}
				}
				else {
					if(getCursor().getType() != Cursor.DEFAULT_CURSOR) {					
						setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
					}
				}
				
				calculateVars = false;
				updatevars = false;
				clusterCalc = false;
		}
		else {
			updatevars = false;
		}
	//	Long end = System.currentTimeMillis();
	//	System.out.println(end - start);
	}
	
	
	
	if(!lineZoomer && !sampleZoomer) {
		
		buf.drawImage(varbuffer, offset, 0, this);
	}
	else {
		if(splits.get(0).viewLength < Settings.readDrawDistance) {
			buf.drawImage(varbuffer, offset, 0, this);
		}
		else {
			//TODO venyta ruutua
		}
	}
}

boolean hideNodeTotal(VarNode node) {
	if(hideNode(node)) {
		
		return true;
	}
	
	foundvisible = false;
	if(!Main.drawCanvas.annotationOn) {
		for(int i=0;i<node.vars.size(); i++) {
			if(!hideNodeVar(node, node.vars.get(i))) {
				for(int s = 0; s <node.vars.get(i).getValue().size(); s++ ) {
					
					if(!hideVar(node.vars.get(i).getValue().get(s),node.vars.get(i).getKey().length() > 1)) {
						foundvisible = true;
						break;
					}
				}			
			}
			if(foundvisible) {
				break;
			}
		}
	}
	else {
		for(int i=0;i<node.vars.size(); i++) {
			if(!hideNodeVar(node, node.vars.get(i))) {
				for(int s = 0; s <node.vars.get(i).getValue().size(); s++ ) {
					if(node.vars.get(i).getValue().get(s).getSample() != null && node.vars.get(i).getValue().get(s).getSample().annotation) {
						continue;
					}
					if(!hideVar(node.vars.get(i).getValue().get(s),node.vars.get(i).getKey().length() > 1)) {
						foundvisible = true;
						break;
					}
				}			
			}
			if(foundvisible) {
				break;
			}
		}
	}
	if(foundvisible) {
		
		return false;
	}
	return true;
	
}
boolean hideNode(VarNode node) {
	
	if(node == null) {
		return true;
	}
	if(VariantHandler.hidenoncoding.isSelected()) {
		
		if(!node.coding) {			
			return true;
		}	
	}
	
	if(VariantHandler.rscode.isSelected()) {
		if(node.isRscode() != null) {
			return true;
		}
	}
	
	if(!FileRead.readFiles && VariantHandler.commonSlider.getValue() >= 1) {
		
		if(VariantHandler.clusterSize > 0) {
			if(VariantHandler.commonSlider.getUpperValue() < Main.varsamples) {
				if(node.incluster) {
					return true;
				}
			}
			else if(!node.incluster) {
				return true;
			}
			
		}		
	}
	
	if(annotationOn) {
		if(node.annotationOnly) {
			return false;
		}		
	}
	
	if(Main.bedCanvas.bedOn) {	
		
		if(!node.bedhit) {
			
			if(!FileRead.readFiles/* && !Main.bedCanvas.annotator*/) {
				
				return true;
			}		
		}
	}
	
	return false;
}
boolean hideNodeCluster(VarNode node) {
	
	if(node == null) {
		return true;
	}
	
	
	if(VariantHandler.hidenoncoding.isSelected()) {
		
		if(!node.coding) {			
			return true;
		}	
	}
	
	if(Main.bedCanvas.bedOn && Main.bedCanvas.intersected && !node.bedhit) {	
		
		if(!FileRead.readFiles && !Main.bedCanvas.annotator) {
			
			return true;
		}		
	}
	
	if(VariantHandler.rscode.isSelected()) {
		if(node.isRscode() != null) {
			return true;
		}
	}
	foundvisible = false;
	
	for(int i = 0; i<node.vars.size(); i++) {
		
		if(!hideNodeVar(node, node.vars.get(i))) {
			for(int j = 0 ; j<node.vars.get(i).getValue().size(); j++) {
				if(!hideVar(node.vars.get(i).getValue().get(j), (node.vars.get(i).getKey().length() > 1))) {
					foundvisible = true;
					
					break;
				}			
			}
		}
		if(foundvisible) {
			break;
		}
	}
	if(!foundvisible) {
		return true;
	}
	
	return false;
}
void calcClusters(VarNode node, int i) {
	clusterId = i;
	calcClusters(node);
}
void calcClusters(VarNode node) {

		try {
			
			if(VariantHandler.commonSlider.getValue() < 1 || VariantHandler.clusterSize < 1 || node == null) {
				return;
			}
			
			boolean first = true;
			
			if(node.getNext()!=null) {
				tempnode = node.getNext();
				
				while(tempnode != null) {
					tempnode.inVarList = false;
					tempnode.incluster = false;
					tempnode.clusterId = -1;
					tempnode = tempnode.getNext();
				}				
				tempnode = null;				
			}
			node = node.getNext(node);		
			
			while(node != null) {
				tempnode = node;
				if(!node.incluster) {
					clustersum = 0;
					first = true;
					
					if(first && tempnode.getNext(tempnode) == null) {	
						
						for(int v = 0; v<tempnode.vars.size(); v++) {							
							if(hideNodeVar(tempnode, tempnode.vars.get(v))) {
								continue;
							}
							for(int i = 0; i<tempnode.vars.get(v).getValue().size(); i++) {
								if(!hideVar(tempnode.vars.get(v).getValue().get(i), (tempnode.vars.get(v).getKey().length() > 1))) {									
									clustersum++;														
								}							
							}
						}						
					}
					if(tempnode.getNext(tempnode) != null && tempnode.getNext(tempnode).getPosition()-tempnode.getPosition() > VariantHandler.clusterSize) {
						for(int v = 0; v<tempnode.vars.size(); v++) {	
							
							if(hideNodeVar(tempnode, tempnode.vars.get(v))) {
								continue;
							}
							for(int i = 0; i<tempnode.vars.get(v).getValue().size(); i++) {
								if(!hideVar(tempnode.vars.get(v).getValue().get(i), (tempnode.vars.get(v).getKey().length() > 1))) {									
									clustersum++;															
								}							
							}
						}
						if(clustersum > 1 && clustersum >= VariantHandler.commonSlider.getValue()) {
							
							tempnode.incluster = true;								
							tempnode.clusterId = clusterId;
							
							clusterId++;						
							
						}
						node = tempnode.getNext(tempnode);
						continue;
					}
					
					while(tempnode.getNext(tempnode) != null && tempnode.getNext(tempnode).getPosition()-tempnode.getPosition() <= VariantHandler.clusterSize) {
						if(first) {
							clustersum = 0;
							for(int v = 0; v<tempnode.vars.size(); v++) {	
								
								if(hideNodeVar(tempnode, tempnode.vars.get(v))) {
									continue;
								}
								
								for(int i = 0; i<tempnode.vars.get(v).getValue().size(); i++) {
									if(!hideVar(tempnode.vars.get(v).getValue().get(i), (tempnode.vars.get(v).getKey().length() > 1))) {	
										tempnode.clusterId = clusterId;							
										tempnode.incluster = true;		
										clustersum++;																
									}							
								}								
							}
								tempnode = tempnode.getNext(tempnode);	
								
							for(int v = 0; v<tempnode.vars.size(); v++) {	
								if(hideNodeVar(tempnode, tempnode.vars.get(v))) {
									continue;
								}
								for(int i = 0; i<tempnode.vars.get(v).getValue().size(); i++) {
									if(!hideVar(tempnode.vars.get(v).getValue().get(i), (tempnode.vars.get(v).getKey().length() > 1))) {								
										tempnode.incluster = true;		
										tempnode.clusterId = clusterId;		
										clustersum++;								
									}							
								}					
								
								first = false;
							}
							
						}
						else {
							tempnode = tempnode.getNext(tempnode);	
								for(int v = 0; v<tempnode.vars.size(); v++) {	
									if(hideNodeVar(tempnode, tempnode.vars.get(v))) {
										continue;
									}
									
									for(int i = 0; i<tempnode.vars.get(v).getValue().size(); i++) {
										if(!hideVar(tempnode.vars.get(v).getValue().get(i), (tempnode.vars.get(v).getKey().length() > 1))) {	
											tempnode.clusterId = clusterId;			
											tempnode.incluster = true;		
											clustersum++;
																	
										}							
									}
							}
						
						}
						
					}
				}	
				else {
					node = node.getNext(node);
				}				
				
				if(clustersum < VariantHandler.commonSlider.getValue()) {
					
					tempnode.incluster = false;
					while(tempnode != null && !tempnode.equals(node)) {
						tempnode.incluster = false;
						tempnode = tempnode.getPrev();
					}
					node.incluster = false;
					
					
					node = node.getNext(node);
				}	
			
				else {
					
					node.incluster = true;
					node.clusterId = clusterId;
				
					clusterId++;					
					node = tempnode.getNext(node);
				}
				
			}
			tempnode = null;
			node = null;
			updatevars = true;
			repaint();
		}
		catch(Exception e ) {
			ErrorLog.addError(e.getStackTrace());
			e.printStackTrace();
		}
		
}

boolean hideNodeVar(VarNode node, Entry<String, ArrayList<SampleNode>> entry ) {
	
	if(VariantHandler.hideIndels.isSelected()) {
		if(entry.getKey().length() > 1) {			
			return true;
		}
	}
	if(VariantHandler.hideSNVs.isSelected()) {
		if(entry.getKey().length() == 1) {
			return true;
		}
	}		
	
	if(annotationOn) {
		
		if(node.annotationOnly) {
			return false;
		}		
		
		if(intersect) {
			boolean found = false;
			
			for(int i = entry.getValue().size()-1; i>= 0;i-- ) {	
				if(entry.getValue().get(i).alleles != null) {
					continue;
				}
				if(entry.getValue().get(i).getSample().annotation && entry.getValue().get(i).getSample().intersect) {
					
					found = true;
					break;
				}			
			}
			if(!found) {
				
				return true;
			}
		}
		
	}
	
	if(!FileRead.readFiles && (VariantHandler.commonSlider.getValue() >= 1 || VariantHandler.commonSlider.getUpperValue() < Main.varsamples)) {
		
		if(VariantHandler.clusterSize == 0) {
			
			int hidecount = 0;
			
			for(int i = 0; i<entry.getValue().size(); i++) {
				if(entry.getValue().get(i).getSample() == null) {
					continue;
				}
				
				if(entry.getValue().get(i).getSample().annotation) {
					continue;
				}
				
				if(!hideVarCommon(entry.getValue().get(i), (entry.getKey().length() > 1))) {	
					
					hidecount++;
				}			
				
			}		
			if(annotationOn) {
				//if(intersect) {
				if (hidecount >= VariantHandler.commonSlider.getValue()) {
					for(int i = 0; i<entry.getValue().size(); i++) {
						if(entry.getValue().get(i).getSample()==null) {
							continue;
						}
						if(entry.getValue().get(i).getSample().annotation) {
							continue;
						}
						if(!hideVarCommon(entry.getValue().get(i), (entry.getKey().length() > 1))) {					
							entry.getValue().get(i).common = true;
						}				
						else {
							entry.getValue().get(i).common = false;
						}
					}		
				}
				else {
					for(int i = 0; i<entry.getValue().size(); i++) {
						if(entry.getValue().get(i).getSample()== null) {
							continue;
						}
						if(entry.getValue().get(i).getSample().annotation) {
							continue;
						}
						
						entry.getValue().get(i).common = false;
						
					}		
				}
				
				/*}
				else {*/
				/*	for(int i = entry.getValue().size()-1; i>= 0;i-- ) {	
						if(entry.getValue().get(i).alleles != null) {
							continue;
						}
						if(entry.getValue().get(i).getSample().annotation) {
							return false;					
						}			
					}*/
				//}
			}
			else {
				
				if (hidecount < VariantHandler.commonSlider.getValue()) {
					return true;
				}
				if(hidecount > VariantHandler.commonSlider.getUpperValue()) {
					return true;
				}			
				
			}
		}
	}
	
	if(node.controlled) {		
		for(int i = entry.getValue().size()-1; i> 0;i-- ) {			
			if(entry.getValue().get(i).alleles == null) {	
				
				break;		
			}
			else {		
				if(!entry.getValue().get(i).getControlSample().controlOn) {
					continue;
				}				
				
				/*if(entry.getValue().get(i).getControlSample().alleleFreq > 1) {
					if(entry.getValue().get(i).alleles > entry.getValue().get(i).getControlSample().alleleFreq) {						
						return true;
					}
				}
				else {*/
					
					if(entry.getValue().get(i).alleles/(double)entry.getValue().get(i).allelenumber > entry.getValue().get(i).getControlSample().alleleFreq) {						
						return true;
					}
				//}
			}			
		}					
	}	
	
	
	return false;
}

boolean hideVar(SampleNode sample, boolean indel) {
	
	if(sample.alleles != null) {		
		return true;
	}
	
	if(sample.getSample() != null) {
		if(sample.getSample().removed) {
			
			return true;
		}
		
	}
	
	if(annotationOn) {
		
		if(sample.getSample().annotation) {			
			return false;
		}
		//if(intersect) {
			if(!FileRead.readFiles && (VariantHandler.commonSlider.getValue() >= 1 || VariantHandler.commonSlider.getUpperValue() < Main.varsamples)) {
			
				if(!sample.common) {
					return true;
				}
			}
		//}
	}
	
	if(!sample.getSample().annoTrack) {
		
		if(!VariantHandler.indelFilters.isSelected()) {	
			if(sample.getQuality() != null && VariantHandler.qualitySlider.getValue() > sample.getQuality()) {		
				return true;
			}
			if(sample.getGQ() != null && VariantHandler.gqSlider.getValue() > sample.getGQ()) {		
				return true;
			}
			if(VariantHandler.coverageSlider.getValue() > sample.getCoverage()) {	
				return true;
			}
			if(VariantHandler.maxCoverageSlider.getValue() < sample.getCoverage()) {		
				return true;
			}
			if(VariantHandler.callSlider.getValue()/100.0 > sample.getAlleleFraction() || sample.getAlleleFraction() == 0.0) {		
				return true;
			}
			if(VariantHandler.callSlider.getUpperValue() != 100) {
				if(VariantHandler.callSlider.getUpperValue()/100.0 < sample.getAlleleFraction()) {		
					return true;
				}
			}
			if(VariantHandler.readSlider.getValue() > sample.getCalls()) {		
				return true;
			}
		}
		else {
			if(indel) {
				if(sample.getQuality() != null && VariantHandler.qualitySliderIndel.getValue() > sample.getQuality()) {		
					return true;
				}
				if(sample.getGQ() != null && VariantHandler.gqSliderIndel.getValue() > sample.getGQ()) {		
					return true;
				}
				if(VariantHandler.coverageSliderIndel.getValue() > sample.getCoverage()) {	
					return true;
				}
				if(VariantHandler.maxCoverageSliderIndel.getValue() < sample.getCoverage()) {		
					return true;
				}
				if(VariantHandler.callSliderIndel.getValue()/100.0 > sample.getAlleleFraction() || sample.getAlleleFraction() == 0.0) {			
					return true;
				}
				if(VariantHandler.readSliderIndel.getValue() > sample.getCalls()) {		
					return true;
				}
			}
			else {
				if(sample.getQuality() != null && VariantHandler.qualitySlider.getValue() > sample.getQuality()) {		
					return true;
				}
				if(sample.getGQ() != null && VariantHandler.gqSlider.getValue() > sample.getGQ()) {		
					return true;
				}
				if(VariantHandler.coverageSlider.getValue() > sample.getCoverage()) {	
					return true;
				}
				if(VariantHandler.maxCoverageSlider.getValue() < sample.getCoverage()) {		
					return true;
				}
				if(VariantHandler.callSlider.getValue()/100.0 > sample.getAlleleFraction() || sample.getAlleleFraction() == 0.0) {		
					return true;
				}
				if(VariantHandler.readSlider.getValue() > sample.getCalls()) {		
					return true;
				}
			}
			
		}
		
	
	if(VariantHandler.hideHomos.isSelected()) {
		if(sample.isHomozygous()) {
			return true;
		}
	}
	}
	
	/*if(drawVariables.advQDraw != null) {
		if(sample.advQualities != null) {
			
			for(int i = 0 ; i<drawVariables.advQDraw.size(); i++) {				
				if(sample.advQualities.containsKey(drawVariables.advQDraw.get(i).key)) {				
					if(sample.advQualities.get(drawVariables.advQDraw.get(i).key) < drawVariables.advQDraw.get(i).value) {						
						return true;
					}
				}			
			}
		}
	}*/
	return false;
}
boolean hideVarCommon(SampleNode sample, boolean indel) {
	
	if(sample.alleles != null) {		
		return true;
	}
	if(sample.getSample() != null) {
		if(sample.getSample().removed) {
			
			return true;
		}
	}
	
	if(!VariantHandler.indelFilters.isSelected()) {	
		if(sample.getQuality() != null && VariantHandler.qualitySlider.getValue() > sample.getQuality()) {		
			return true;
		}
		if(sample.getGQ() != null && VariantHandler.gqSlider.getValue() > sample.getGQ()) {		
			return true;
		}
		if(VariantHandler.coverageSlider.getValue() > sample.getCoverage()) {	
			return true;
		}
		if(VariantHandler.maxCoverageSlider.getValue() < sample.getCoverage()) {		
			return true;
		}
		if(VariantHandler.callSlider.getValue()/100.0 > sample.getAlleleFraction()) {		
			return true;
		}
		if(VariantHandler.callSlider.getUpperValue() != 100) {
			if(VariantHandler.callSlider.getUpperValue()/100.0 < sample.getAlleleFraction()) {		
				return true;
			}
		}
	}
	else {
		if(indel) {
			if(sample.getQuality() != null && VariantHandler.qualitySliderIndel.getValue() > sample.getQuality()) {		
				return true;
			}
			if(sample.getGQ() != null && VariantHandler.gqSliderIndel.getValue() > sample.getGQ()) {		
				return true;
			}
			if(VariantHandler.coverageSliderIndel.getValue() > sample.getCoverage()) {	
				return true;
			}
			if(VariantHandler.maxCoverageSliderIndel.getValue() < sample.getCoverage()) {		
				return true;
			}
			if(VariantHandler.callSliderIndel.getValue()/100.0 > sample.getAlleleFraction()) {		
				return true;
			}
		}
		else {
			if(sample.getQuality() != null && VariantHandler.qualitySlider.getValue() > sample.getQuality()) {		
				return true;
			}
			if(sample.getGQ() != null && VariantHandler.gqSlider.getValue() > sample.getGQ()) {		
				return true;
			}
			if(VariantHandler.coverageSlider.getValue() > sample.getCoverage()) {	
				return true;
			}
			if(VariantHandler.maxCoverageSlider.getValue() < sample.getCoverage()) {		
				return true;
			}
			if(VariantHandler.callSlider.getValue()/100.0 > sample.getAlleleFraction()) {		
				return true;
			}
		}
		
	}
	
	
	if(VariantHandler.hideHomos.isSelected()) {
		if(sample.isHomozygous()) {
			return true;
		}
	}
	
	/*if(drawVariables.advQDraw != null) {
		if(sample.advQualities != null) {
			
			for(int i = 0 ; i<drawVariables.advQDraw.size(); i++) {				
				if(sample.advQualities.containsKey(drawVariables.advQDraw.get(i).key)) {				
					if(sample.advQualities.get(drawVariables.advQDraw.get(i).key) < drawVariables.advQDraw.get(i).value) {						
						return true;
					}
				}			
			}
		}
	}*/
	return false;
}
void drawZoom() {
	
	if(Main.chromDraw.lineZoomer || Main.bedCanvas.lineZoomer || Main.controlDraw.lineZoomer) {
		return;
	}
	buf.setColor(Color.black);
	
	
	if(lineZoomer) {
		buf.setStroke(dashed);
		buf.drawLine(pressX, pressY-Main.drawScroll.getVerticalScrollBar().getValue(), mouseX, mouseY-Main.drawScroll.getVerticalScrollBar().getValue());			
	}	
	else if(!sampleZoomer && zoomDrag && !lineZoomer) { //|| mouseX-pressX >= 0 && mouseY-pressY > -10)) {
		
		buf.setColor(Color.white);
	//	buf.drawRect(pressX, pressY-Main.drawScroll.getVerticalScrollBar().getValue(), mouseX-pressX, mouseY-pressY);		
		buf.drawRect(pressX, 0, mouseX-pressX-1, Main.drawScroll.getViewport().getHeight());	
		buf.setColor(zoomColor);
	//	buf.fillRect(pressX, pressY-Main.drawScroll.getVerticalScrollBar().getValue(), mouseX-pressX, mouseY-pressY);
		buf.fillRect(pressX, 0, mouseX-pressX, Main.drawScroll.getViewport().getHeight());
		if(this.mouseX-this.pressX > 0) {
			Main.widthLabel.setText(MethodLibrary.formatNumber((int)((this.mouseX-this.pressX)/selectedSplit.pixel)) +"bp  (Right click to cancel zoom)");
		}
	}
	if (!sampleZoomer && mouseX-pressX < 0 && mouseY-pressY < 0) {
		
		lineZoomer = true;		
		
	//	buf.drawLine(pressX, pressY, mouseX, mouseY);
		
	}

	buf.setStroke(basicStroke);
	buf.setColor(Color.black);
}

void checkSampleZoom() {
//	if(drawVariables.sampleHeight < defaultSampleHeight) {
	
		if(drawVariables.visiblesamples > Main.drawScroll.getViewport().getHeight()) {
			drawVariables.visiblesamples = (short)(Main.drawScroll.getViewport().getHeight());
		}
		else if(drawVariables.visiblesamples > Main.samples) {
			drawVariables.visiblesamples = (short)Main.samples;
		}
		else if(drawVariables.visiblestart > 0 && drawVariables.visiblestart+drawVariables.visiblesamples > Main.samples) {
			drawVariables.visiblestart--;
		}
//	}
//	else if(drawVariables.visiblesamples >  Main.drawScroll.getViewport().getHeight()/(double)(defaultSampleHeight)) {
//		drawVariables.visiblesamples = (short)(Main.drawScroll.getViewport().getHeight()/(double)(defaultSampleHeight));
//	}
}
void zoomSamples() {
	
	if(prezoom != (pressX-mouseX)/20 ) {
		
		if(prezoom < (pressX-mouseX)/20 ) {
			//shrink
			if(drawVariables.visiblesamples >= Main.samples) {
				drawVariables.visiblesamples =  (short)Main.samples;				
				
			}
			else {
				
				drawVariables.visiblesamples+=1+(drawVariables.visiblesamples/10);
				
				if(drawVariables.visiblesamples > Main.samples) {				
					drawVariables.visiblesamples =  (short)Main.samples;
				}
				if(drawVariables.visiblestart+drawVariables.visiblesamples > Main.samples) {
					drawVariables.visiblestart=(short)(Main.samples-drawVariables.visiblesamples);
					if(drawVariables.visiblestart<0) {
						drawVariables.visiblestart=0;
						
					}
				}
			}
		}
		else {		
			
			drawVariables.visiblesamples-=1+(drawVariables.visiblesamples/10);
			
			if(drawVariables.visiblesamples < 1) {
				drawVariables.visiblesamples = 1;
			}			
		}
		
		drawVariables.sampleHeight = Main.drawScroll.getViewport().getHeight()/(double)(drawVariables.visiblesamples);
		if(drawVariables.sampleHeight < 1) {
			drawVariables.sampleHeight = 1;
		}
		
		checkSampleZoom();		
		
		//setScrollbar((int)(drawVariables.visiblestart*drawVariables.sampleHeight));	 
		this.height = (int)(drawVariables.sampleHeight*Main.samples);		
		this.setPreferredSize(new Dimension(this.getWidth(),this.height));		
		this.revalidate();
		setScrollbar((int)(drawVariables.visiblestart*drawVariables.sampleHeight));	 		
		pressY = (int)(this.height*pressYrelative);		
		
	}
	prezoom = (pressX-mouseX)/20;
	
}

void drawSidebar() {	
	
	sidebuf.setColor(sidecolor);
	sidebuf.fillRect(0, 0, Main.sidebarWidth, Main.drawScroll.getViewport().getHeight());
	
	if(!sidebuf.getFont().equals(defaultFont)) {
		sidebuf.setFont(defaultFont);
	}
	if(selectedSampleIndex != -1) {
		sidebuf.setColor(Color.white);
		sidebuf.fillRect(0, (int)(selectedSampleIndex*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue()+2, Main.defaultFontSize*4, (int)drawVariables.sampleHeight);
		sidebuf.setColor(Color.black);
		sidebuf.setColor(Color.gray);
		sidebuf.drawLine(0,(int)((sidebarSample.getIndex()+1)*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue(), Main.sidebarWidth-1, (int)((sidebarSample.getIndex()+1)*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue());
		sidebuf.setColor(Color.lightGray);
		sidebuf.drawLine(0,(int)((sidebarSample.getIndex()+1)*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue()+1, Main.sidebarWidth-1, (int)((sidebarSample.getIndex()+1)*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue()+1);
		sidebuf.setColor(Color.black);
	}
	if(sidebar) {		
		sidebuf.setColor(Color.white);
		
		sidebuf.fillRect(0, (int)(selectedIndex*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue(), Main.defaultFontSize*4, (int)drawVariables.sampleHeight);
		if(hoverIndex > -1) {
			sidebuf.fillRect(0, (int)(hoverIndex*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue(), Main.defaultFontSize*4, (int)drawVariables.sampleHeight);
			
		}
	}
	sidebuf.setColor(Color.black);
	if(drawVariables.sampleHeight < Main.defaultFontSize) {
		sidebuf.drawString("Tracks: " +Main.samples, 10, Main.defaultFontSize + 5);
		if(selectedSample != null && !selectedSample.removed) {
			sidebuf.drawString("Selected: " +(selectedSample.getIndex()+1) +". " +selectedSample.getName(), 10, Main.defaultFontSize*3);
		}
		if(sidebar && selectedIndex > -1  && removeSample == null) {
			/*if(selectedSample == null) {
				sidebuf.drawString((selectedIndex+1) +". " +sampleList.get(selectedIndex).getName(), 10, Main.defaultFontSize*5);
			}
			else {*/
				sidebuf.drawString((selectedIndex+1) +". " +sampleList.get(selectedIndex).getName(), 10, Main.defaultFontSize*5);
			//}
		}
	}
	sampleInfo = false;
	bamHover = false;
	removeSample = null;

	if( drawVariables.visiblestart < 0) {
		 drawVariables.visiblestart = 0;
	}

	for(int i = drawVariables.visiblestart; i<drawVariables.visiblestart+drawVariables.visiblesamples+2+removesamples; i++) {
		if(i > sampleList.size()-1) {
			break;
		}
		
		if(sampleList.get(i).removed) {
			removesamples++;
			continue;
		}
		
		sidebarSample = sampleList.get(i);
		sampleYpos = (int)(drawVariables.sampleHeight*sidebarSample.getIndex())-Main.drawScroll.getVerticalScrollBar().getValue();
		bottomYpos = (int)(drawVariables.sampleHeight*(sidebarSample.getIndex()+1))-Main.drawScroll.getVerticalScrollBar().getValue();
		if(sidebarSample.familyColor != null) {			
			sidebuf.setColor(sidebarSample.familyColor);
			sidebuf.fillRect(0, sampleYpos+2, Main.defaultFontSize*2, (int)drawVariables.sampleHeight);
			sidebuf.setColor(Color.black);
		}
		if(sampleList.get(i).affected != null && sampleList.get(i).affected) {
			sidebuf.setColor(Color.red);
			sidebuf.fillOval(Main.sidebarWidth-Main.defaultFontSize-8, sampleYpos+4, Main.defaultFontSize, Main.defaultFontSize);
			sidebuf.setColor(Color.black);
		}
		if(sampleList.get(i).annotation) {
			sidebuf.setColor(Color.orange);
			sidebuf.fillOval(Main.sidebarWidth-Main.defaultFontSize-8, sampleYpos+4, Main.defaultFontSize, Main.defaultFontSize);
			sidebuf.setColor(Color.black);
		}
		if(/*((i == selectedIndex || i== hoverIndex) && sidebar) ||*/ sidebar && (selectedSampleIndex == i ||selectedIndex == i)) {		
			
			
			if((this.remoBox.getBounds().x != Main.sidebarWidth-(Main.defaultFontSize+10) || this.remoBox.getBounds().y != sampleYpos +(int)drawVariables.sampleHeight -(Main.defaultFontSize+6))) {
				this.remoBox.setBounds(Main.sidebarWidth-(Main.defaultFontSize+10), sampleYpos +(int)drawVariables.sampleHeight -(Main.defaultFontSize+6), Main.defaultFontSize+3, Main.defaultFontSize+3);
			}
			if(drawVariables.sampleHeight > Main.defaultFontSize) {
				
				
				/*sidebuf.drawString(sidebarSample.getName(), 50, sampleYpos+Main.defaultFontSize);
				sidebuf.drawString(""+(sidebarSample.getIndex()+1), 10, sampleYpos+Main.defaultFontSize);
		//		sidebuf.drawLine(0,(int)((sidebarSample.getIndex())*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue(), Main.sidebarWidth-1, (int)((sidebarSample.getIndex())*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue());
				sidebuf.setColor(Color.gray);
				sidebuf.drawLine(0,bottomYpos, Main.sidebarWidth-1, bottomYpos);
				sidebuf.setColor(Color.lightGray);
				sidebuf.drawLine(0,bottomYpos+1, Main.sidebarWidth-1, bottomYpos+1);
				
				sidebuf.setColor(Color.black);
				if(moveY >= sampleYpos && moveY <= sampleYpos+ Main.defaultFontSize+2) {
					sidebuf.setFont(defaultFont);
				}*/
			}
			/*if(drawVariables.sampleHeight > 4) {
				sidebuf.drawLine(0,(int)((sidebarSample.getIndex())*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue(), Main.sidebarWidth-1, (int)((sidebarSample.getIndex())*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue());
			}*/
			if(sidebar && mouseRect.intersects(this.remoBox)) {				
				sidebuf.setStroke(doubleStroke);
				removeSample = sidebarSample;
			}
			
			
			
			sidebuf.setColor(Color.white);
			sidebuf.fillRect(this.remoBox.x-2, this.remoBox.y-1, this.remoBox.width+2, this.remoBox.height+2);
			sidebuf.setColor(Color.black);
			sidebuf.drawRect(this.remoBox.x, this.remoBox.y, this.remoBox.width, this.remoBox.height);
			sidebuf.drawLine(this.remoBox.x, this.remoBox.y, this.remoBox.x+this.remoBox.width, this.remoBox.y+(int)this.remoBox.getHeight());
			sidebuf.drawLine(this.remoBox.x, this.remoBox.y+(int)this.remoBox.getHeight(), this.remoBox.x+this.remoBox.width,this.remoBox.y);
			if(sidebuf.getStroke().equals(doubleStroke)) {
				sidebuf.setStroke(basicStroke);
			}
		}		
	//	else {
			
		/*	if(sidebarSample.getName().length() > 15) {
				sidebuf.drawString(sidebarSample.getName().substring(0,15) +"...", 50, (int)(drawVariables.sampleHeight*sidebarSample.getIndex())-Main.drawScroll.getVerticalScrollBar().getValue()+13);
			}
			else {*/
			if(drawVariables.sampleHeight > Main.defaultFontSize) {
				if(sidebar) {
				
					if(moveX > 50 && moveY >= sampleYpos && moveY <= sampleYpos+ Main.defaultFontSize+2) {
						sidebuf.setFont(Main.menuFontBold);
						setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
						sampleInfo = true;
					}
				}
				
				sidebuf.drawString(sidebarSample.getName(), 50, sampleYpos+Main.defaultFontSize);
				sidebuf.drawString(""+(sidebarSample.getIndex()+1), 10, sampleYpos+Main.defaultFontSize);
				sidebuf.setColor(Color.gray);
				sidebuf.drawLine(0,bottomYpos, Main.sidebarWidth-1, bottomYpos);
				sidebuf.setColor(Color.lightGray);
				sidebuf.drawLine(0,bottomYpos+1, Main.sidebarWidth-1, bottomYpos+1);
				
				sidebuf.setColor(Color.black);
				if(sidebar) {
					if(sampleInfo) {
						sidebuf.setFont(defaultFont);
						
					}
				}
			}
	//	}
		if(sidebarSample.getTabixFile() == null && !sidebarSample.calledvariants && !sidebarSample.multipart && drawVariables.sampleHeight > Main.defaultFontSize*2) {
			if(sidebarSample.samFile != null) {
				
				sidebuf.setColor(ChromDraw.exonBarColor);
				if(sidebarSample.readString != null) {
					sidebuf.drawString(sidebarSample.readString, Main.defaultFontSize*5, (int)(drawVariables.sampleHeight*sidebarSample.getIndex())-Main.drawScroll.getVerticalScrollBar().getValue()+Main.defaultFontSize*2);
				}
			}
			sidebuf.setColor(Color.red);
			sidebuf.drawString("No VCF", 10, sampleYpos+Main.defaultFontSize*2);
			
			sidebuf.setColor(Color.black);
			continue;
		}
		if((sidebarSample.getTabixFile() != null || sidebarSample.calledvariants || sidebarSample.multipart) && drawVariables.sampleHeight > Main.defaultFontSize*3) {
			if(sidebarSample.indels == 1) {
				sidebuf.drawString("variants: "+sidebarSample.varcount +" (" +sidebarSample.indels +" indel)", 10, (int)(drawVariables.sampleHeight*sidebarSample.getIndex())-Main.drawScroll.getVerticalScrollBar().getValue()+Main.defaultFontSize*3);
				
			}
			else {
				sidebuf.drawString("variants: "+sidebarSample.varcount +" (" +sidebarSample.indels +" indels)", 10, (int)(drawVariables.sampleHeight*sidebarSample.getIndex())-Main.drawScroll.getVerticalScrollBar().getValue()+Main.defaultFontSize*3);
				
			}
			if(drawVariables.sampleHeight > Main.defaultFontSize*4) {
				sidebuf.drawString("het: "+sidebarSample.heterozygotes, 10, sampleYpos+Main.defaultFontSize*4);
				if(drawVariables.sampleHeight > Main.defaultFontSize*5) {
					sidebuf.drawString("hom: "+sidebarSample.homozygotes, 10, sampleYpos+Main.defaultFontSize*5);
					if(drawVariables.sampleHeight > Main.defaultFontSize*6) {
						sidebuf.drawString("Ts/Tv: "+MethodLibrary.round(sidebarSample.sitioRate/(double)sidebarSample.versioRate,2), 10, sampleYpos+Main.defaultFontSize*6);
					
						if(drawVariables.sampleHeight > Main.defaultFontSize*7) {
							if(sidebarSample.getTabixFile() !=null || sidebarSample.multipart) {
								sidebuf.setColor(ChromDraw.exonBarColor);
								sidebuf.drawString("VCF", 10, sampleYpos+Main.defaultFontSize*7);
								
							}
							else {
								sidebuf.setColor(Color.red);
								sidebuf.drawString("No VCF", 10, sampleYpos+Main.defaultFontSize*7);
								
							}
							if(sidebarSample.samFile != null) {
								
								sidebuf.setColor(ChromDraw.exonBarColor);
								if(sidebarSample.readString != null) {
									sidebuf.drawString(sidebarSample.readString, Main.defaultFontSize*4, sampleYpos+Main.defaultFontSize*7);
								}
							/*	if(sidebarSample.CRAM) {
									sidebuf.drawString("CRAM", 40, (int)(drawVariables.sampleHeight*sidebarSample.getIndex())-Main.drawScroll.getVerticalScrollBar().getValue()+90);							
									
								}
								else {
									sidebuf.drawString("BAM", 40, (int)(drawVariables.sampleHeight*sidebarSample.getIndex())-Main.drawScroll.getVerticalScrollBar().getValue()+90);							
								}*/
							}
							else {
								sidebuf.setColor(Color.red);								
								if(!FileRead.searchingBams && sidebarSample.readString != null) {
									if(sidebar) {
										//System.out.println(moveY +" " +(sampleYpos+Main.defaultFontSize*6));
										if(moveY-Main.drawScroll.getVerticalScrollBar().getValue() >= sampleYpos+Main.defaultFontSize*6 && moveY-Main.drawScroll.getVerticalScrollBar().getValue() <= sampleYpos+Main.defaultFontSize*6+ Main.defaultFontSize+2) {
											if(moveX >= Main.defaultFontSize*4) {												
											
												sidebuf.setFont(Main.menuFontBold);
												setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
												bamHover = true;
												
											}
										}
									}
									sidebuf.drawString(sidebarSample.readString, Main.defaultFontSize*4, sampleYpos+Main.defaultFontSize*7);
									if(sidebar) {
										if(bamHover) {
											sidebuf.setFont(defaultFont);
											
										}
									}
								}
								else {
									sidebuf.drawString("Searching reads...", Main.defaultFontSize*4, sampleYpos+Main.defaultFontSize*7);
									
								}
							}							
							sidebuf.setColor(Color.black);		
							if(drawVariables.sampleHeight > Main.defaultFontSize*8) {															
								sidebuf.drawString("Var height by " +Settings.varDrawList.getSelectedItem().toString() +" (max "+(int)sidebarSample.getMaxCoverage() +")", 10, sampleYpos+Main.defaultFontSize*8);
														
							}							
						}						
					}
				}
			}			
		}
	}
	
	if(sidebar) {
		if(!sampleInfo && !bamHover) {			
			setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));			
		}
	}
	sidebuf.setComposite( composite);		
	sidebuf.fillRect(Main.sidebarWidth,0, this.getWidth(),Main.drawScroll.getViewport().getHeight());	
	sidebuf.setComposite(backupside);		
	sidebuf.setStroke(doubleStroke);
	sidebuf.setColor(Color.gray);
	sidebuf.drawLine(Main.sidebarWidth-1, 0, Main.sidebarWidth-1, Main.drawScroll.getViewport().getHeight());
	sidebuf.drawLine(1, 0, 1, Main.drawScroll.getViewport().getHeight());
	sidebuf.setColor(Color.lightGray);
	sidebuf.drawLine(3, 0,3, Main.drawScroll.getViewport().getHeight());
	sidebuf.drawLine(Main.sidebarWidth-3, 0, Main.sidebarWidth-3, Main.drawScroll.getViewport().getHeight());
	sidebuf.setStroke(basicStroke);
	buf.drawImage(sidebuffer, 0, 0, null);
}

void resizeAllCanvas(int width) {
	
	backbuffer = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
	buf = (Graphics2D)backbuffer.getGraphics();
	buf.setRenderingHints(rh);
	
	backupb = buf.getComposite();
		
	varbuffer = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
	g2v = (Graphics2D)varbuffer.getGraphics();
	
	backupv = g2v.getComposite();
	
	readbuffer = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
	rbuf = (Graphics2D)readbuffer.getGraphics();
	backupr = rbuf.getComposite();
	rbuf.setRenderingHints(rh);
	coveragebuffer = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
	cbuf = (Graphics2D)coveragebuffer.getGraphics();
	backupc = cbuf.getComposite();
	
	sidebuffer = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));
	sidebuf = (Graphics2D)sidebuffer.getGraphics();
	backupside = sidebuf.getComposite();
	
	Main.chromDraw.chromImage = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
	Main.chromDraw.chromImageBuffer = (Graphics2D)Main.chromDraw.chromImage.getGraphics();
	//Main.chromDraw.exonImage = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
	//Main.chromDraw.exonImageBuffer = (Graphics2D)Main.chromDraw.exonImage.getGraphics();
	Main.chromDraw.selectImage = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));
	Main.chromDraw.selectImageBuffer = (Graphics2D)Main.chromDraw.selectImage.getGraphics();
	
//	Main.chromDraw.backupe = Main.chromDraw.exonImageBuffer.getComposite();	
	
	Main.bedCanvas.bufImage = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
	Main.bedCanvas.buf = (Graphics2D)Main.bedCanvas.bufImage.getGraphics();
	
	Main.bedCanvas.nodeImage = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));
	Main.bedCanvas.nodebuf = (Graphics2D)Main.bedCanvas.nodeImage.getGraphics();
	
	Main.bedCanvas.backupComposite = Main.bedCanvas.nodebuf.getComposite();
	
	for(int i = 0; i<splits.size(); i++) {
		splits.get(i).getSplitDraw().resizeImages(width);
	
	}
}

public void drawScreen(Graphics g) {
		if(Settings.wallpaper == null) {
			buf.setColor(backColor);		
			buf.fillRect(Main.sidebarWidth, 0,this.getWidth()-Main.sidebarWidth, Main.drawScroll.getViewport().getHeight());	
		}
		else {
			buf.drawImage(Settings.wallpaper, 0, 0, this.getWidth(), Main.drawScroll.getViewport().getHeight(),this);
			buf.setColor(backColor);	
			buf.fillRect(Main.sidebarWidth, 0,this.getWidth()-Main.sidebarWidth, Main.drawScroll.getViewport().getHeight());
		}
		if(Main.drawScroll.getVerticalScrollBar().getUnitIncrement() != (int)drawVariables.sampleHeight) {
			Main.drawScroll.getVerticalScrollBar().setUnitIncrement((int)drawVariables.sampleHeight);		
			
		}
		
		//System.out.println(Main.drawScroll.getVerticalScrollBar().getValue() +" " +(int)(drawVariables.visiblestart*drawVariables.sampleHeight));
		if(Main.chromDraw.timer > 0) {			
			if(System.currentTimeMillis()-Main.chromDraw.timer > 2000) {
				Main.chromDraw.timer = 0;
				Main.chromDraw.repaint();
			} 				
		}
		
		if((Main.varsamples > 0 || annotationOn) && splits.get(0).viewLength <= 1000) {			
			drawVars(this.splits.get(0).offset);		
		}
		
		
		if(Main.samples > 0 && drawVariables.sampleHeight > 100) {
			
			if(mouseDrag) {			
				if(splits.size() > 1) {			
					for(int s = 0; s<Main.drawCanvas.splits.size(); s++) {						
						if(Main.drawCanvas.splits.get(s)!=selectedSplit){
							rbuf.setComposite( composite);			
							rbuf.fillRect(Main.drawCanvas.splits.get(s).offset,0, this.drawWidth,Main.drawScroll.getViewport().getHeight());	
							rbuf.setComposite(backupr);		
							rbuf.drawImage(Main.drawCanvas.splits.get(s).getReadImage(), Main.drawCanvas.splits.get(s).offset, 0,this);	
							rbuf.drawImage(Main.drawCanvas.splits.get(s).getSelectbuffer(), this.splits.get(s).offset, 0, this);
						}
						else {
							rbuf.setComposite( composite);			
							rbuf.fillRect(selectedSplit.offset,0, this.drawWidth,Main.drawScroll.getViewport().getHeight());	
							rbuf.setComposite(backupr);	
							rbuf.drawImage(selectedSplit.getReadImage(),selectedSplit.offset,0,selectedSplit.offset+Main.drawScroll.getViewport().getWidth(), this.height, -(moveX-pressX), 0,-(moveX-pressX)+Main.drawScroll.getViewport().getWidth(), this.getHeight(),this);
							rbuf.drawImage(selectedSplit.getSelectbuffer(), selectedSplit.offset+(moveX-pressX), 0, this);
						}				
					}
				}
				else {					
					rbuf.setComposite(composite);		
					rbuf.fillRect(0,0, this.getWidth(),Main.drawScroll.getViewport().getHeight());					
					rbuf.setComposite(backupr);		
					rbuf.drawImage(selectedSplit.getReadImage(),selectedSplit.offset,0,selectedSplit.offset+Main.drawScroll.getViewport().getWidth(), this.height, -(moveX-pressX), 0,Main.drawScroll.getViewport().getWidth()-(moveX-pressX), this.getHeight(),this);
				}	
			}			
			else {		
				if(scrolldrag) {									
					rbuf.setComposite( composite);				
					rbuf.fillRect(0,0, selectedSplit.getReadImage().getWidth(),Main.drawScroll.getViewport().getHeight());	
					rbuf.setComposite(backupr);		
				}
				else if((!sampleZoomer) && (!zoomDrag || lineZoomer)) {					
					
					rbuf.setComposite( composite);				
					rbuf.fillRect(0,0, selectedSplit.getReadImage().getWidth(),Main.drawScroll.getViewport().getHeight());	
					rbuf.setComposite(backupr);						
					drawReads(selectedSplit);					
					drawCoverage(selectedSplit);					
					
					for(int s = 0; s<Main.drawCanvas.splits.size(); s++) {				
						if(Main.drawCanvas.splits.get(s) != selectedSplit && splits.get(s).updateReads) {						
							drawReads(this.splits.get(s));							
						}
						if(Main.drawCanvas.splits.get(s) != selectedSplit) {						
							drawCoverage(this.splits.get(s));
						}
						rbuf.drawImage(Main.drawCanvas.splits.get(s).getReadImage(), this.splits.get(s).offset, 0,this);	
						rbuf.drawImage(Main.drawCanvas.splits.get(s).getSelectbuffer(), this.splits.get(s).offset, 0, this);		
						rbuf.setComposite( composite);					
						rbuf.fillRect(Main.drawCanvas.splits.get(s).offset+this.getDrawWidth(),0, this.getDrawWidth(),Main.drawScroll.getViewport().getHeight());	
						rbuf.setComposite(backupr);			
					}					
				}				
			}
		}
		else {
			
			zoomNote = true;
			if(lineZoomer || sampleZoomer) {				
				//TODO clearaa coveraget paremmin
				splitList.clear();
				selectedRead = null;
				
				rbuf.setComposite( composite);				
				rbuf.fillRect(0,0, selectedSplit.getReadImage().getWidth(),Main.drawScroll.getViewport().getHeight());	
				rbuf.setComposite(backupr);		
				
			}
			
			splitList.clear();
			selectedRead = null;
			clickedRead = null;
			rbuf.setComposite( composite);				
			rbuf.fillRect(0,0, selectedSplit.getReadImage().getWidth(),Main.drawScroll.getViewport().getHeight());	
			rbuf.setComposite(backupr);		
			
		}
		if(Main.drawCanvas.variantsEnd > 0 && Main.varsamples > 0) {
			drawVarCurtain();
		}
		if(Main.samples > 0 && !sampleZoomer) {			
			buf.drawImage(readbuffer, 0, 0, this);		
		}
	
		if((Main.varsamples > 0 || annotationOn) && splits.get(0).viewLength > 1000) {	
			
			drawVars(this.splits.get(0).offset);		
		}
		drawSidebar();		
		
		buf.setColor(Color.black);
		buf.setStroke(basicStroke);
		
		if(Main.samples > 0 && drawVariables.sampleHeight > 10) {
			for(int i = drawVariables.visiblestart ; i<drawVariables.visiblestart+drawVariables.visiblesamples+1; i++) {
				if(i<Main.samples) {
					sampleList.get(i).prepixel = -1;
				}
			if(i != 0) {	
				buf.setColor(Color.lightGray);
				buf.drawLine(Main.sidebarWidth,(int)((i)*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue(), this.getWidth(), (int)((i)*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue());
				
				buf.setColor(Color.gray);
				buf.drawLine(Main.sidebarWidth,(int)((i)*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue()+1, this.getWidth(), (int)((i)*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue()+1);
				buf.setColor(Color.black);
			}
				if(i == selectedSampleIndex) {
					buf.setColor(Color.white);
					buf.drawLine(Main.sidebarWidth,(int)((i)*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue(), this.getWidth(), (int)((i)*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue());
					buf.drawLine(Main.sidebarWidth,(int)((i+1)*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue()-1, this.getWidth(), (int)((i+1)*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue()-1);
					
					buf.setColor(Color.black);
				}
			}
		}
		if(selectedSplit.viewLength < 300 && Main.genomehash.size() > 0) {
			buf.setColor(Color.black);			
			buf.setStroke(dashed);		
			buf.drawLine((this.getDrawWidth())/2+selectedSplit.offset, 0, (this.getDrawWidth())/2+selectedSplit.offset, Main.drawScroll.getViewport().getHeight());
			buf.drawLine((int)((this.getDrawWidth())/2+selectedSplit.pixel)+selectedSplit.offset, 0, (int)((this.getDrawWidth())/2+selectedSplit.pixel)+selectedSplit.offset, Main.drawScroll.getViewport().getHeight());
			buf.setStroke(basicStroke);
					
		}
		
		if(splits.size() > 1 && (selectedIndex < sampleList.size() || sampleList.size() == 0)) {
				
				buf.setStroke(strongStroke);
				if(sampleList.size() > 0 && sampleList.get(selectedIndex).getreadHash() == null) {
					sampleList.get(selectedIndex).resetreadHash();
				}
				
				if(mouseRect.x > this.drawWidth-35) {				
					if(mouseRect.x < this.drawWidth-5) {
						if(sampleList.size() > 0) {
							
							if(mouseRect.y/*-Main.drawScroll.getVerticalScrollBar().getValue()+sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel*/ > 5 && mouseRect.y/*-Main.drawScroll.getVerticalScrollBar().getValue()+sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel*/ < 35) {
								buf.setColor(Color.white);
								
								splitremove = true;
							}
							else {
								splitremove = false;
							}
						}
						else {
							if(mouseRect.y/*-Main.drawScroll.getVerticalScrollBar().getValue()*/ > 5 && mouseRect.y/*-Main.drawScroll.getVerticalScrollBar().getValue()*/ < 35) {
								buf.setColor(Color.white);
								splitremove = true;
							}
							else {
								splitremove = false;
							}
						}							
					}	
					else {
						splitremove = false;
					}
				}
				else {
					splitremove = false;
				}
				if(buf.getColor().equals(Color.black)) {
					buf.setColor(Color.white);
					buf.fillRect(selectedSplit.offset+this.drawWidth-35, 5, 30, 30);
					buf.setColor(Color.black);
				}
				else if(buf.getColor().equals(Color.white)) {
					buf.setColor(Color.black);
					buf.fillRect(selectedSplit.offset+this.drawWidth-35, 5, 30, 30);
					buf.setColor(Color.white);
				}
				buf.drawLine(selectedSplit.offset+this.drawWidth-33, 33, selectedSplit.offset+this.drawWidth-7, 7);
				buf.drawLine(selectedSplit.offset+this.drawWidth-7, 33, selectedSplit.offset+this.drawWidth-33, 7);
				
				buf.setColor(Color.black);
				
				buf.drawRect(selectedSplit.offset+this.drawWidth-35, 5, 30, 30);
				for(int i = 1; i<splits.size(); i++) {
					buf.setColor(Color.gray);
					buf.fillRect(splits.get(i).offset-3, 0, 5, this.getHeight());
					buf.setColor(Color.lightGray);
					buf.fillRect(splits.get(i).offset-1, 0, 2, this.getHeight());
					
				}
				buf.setStroke(basicStroke);
				
		}
		drawClickedRead();	
		if(zoomDrag || lineZoomer) {
			
			drawZoom();
		}
		else if (sampleZoomer) {
			
			 zoomSamples();
		}
		
		if(!lineZoomer && moveX > Main.sidebarWidth-3 && moveX < Main.sidebarWidth+3) {
		
			if(getCursor().getType() != Cursor.W_RESIZE_CURSOR) {
				
				sampleZoomer = false;
				setCursor(Cursor.getPredefinedCursor(Cursor.W_RESIZE_CURSOR));
			}
		}
		else {
			
			if(getCursor().getType() == Cursor.W_RESIZE_CURSOR) {
				//resizeSidebar = false;
				
				setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
			}
		}		

		//buf.setStroke(doubleStroke);
		
	//	buf.setColor(Color.white);
	//	buf.drawLine(Main.sidebarWidth, 0, Main.sidebarWidth, Main.drawScroll.getViewport().getHeight());
		if(mousePress) {
			drawNavbar();			
		}
	
		if(mouseWheel) {
			mouseWheel = false;
		}
		
		if(bam && zoomNote && splits.get(0).viewLength > Settings.coverageDrawDistance  && drawVariables.sampleHeight >= 100) {
			if(!buf.getColor().equals(Color.lightGray)) {
				buf.setColor(Color.lightGray);				
			}			
			
			buf.drawString("Zoom in closer than " +MethodLibrary.formatNumber(Settings.coverageDrawDistance) +" bp to see reads and coverages", Main.sidebarWidth +20, Main.defaultFontSize*2 +6);
			
		}	
		else if (bam && zoomNote && drawVariables.sampleHeight < 100) {
			if(!buf.getColor().equals(Color.lightGray)) {
				buf.setColor(Color.lightGray);				
			}	
			buf.drawString("Expand sample tracks to see reads and coverages", Main.sidebarWidth +20, 20);
			
		}
		else {
			if(FileRead.cancelreadread) {
				if(!buf.getColor().equals(Color.lightGray)) {
					buf.setColor(Color.lightGray);				
				}
				buf.drawString("Double click sample tracks to load reads again.", Main.sidebarWidth +20, 20);
			//	Main.putMessage("Double click sample tracks to load reads.");
			}
		}
		
		if(sampleDrag) {
			buf.setColor(softLight);
			buf.fillRect(moveX, moveY, drawWidth+Main.sidebarWidth, (int)drawVariables.sampleHeight);
			
		}
	//	buf.setColor(Color.white);
	//	buf.setFont(defaultFont);
		//buf.drawString("Font size: "+Main.defaultFontSize, 300, 100);
		
		g.drawImage(backbuffer, 0, Main.drawScroll.getVerticalScrollBar().getValue(), null);		
		
}

void drawVarCurtain() {	
	
	buf.setColor(curtainColor);
	
	if(this.varOverLap == null && !sidebar) {
		if((mouseRect.intersects(varStartRect) || mouseRect.intersects(varEndRect))) {
			
			buf.setFont(Main.menuFontBold);
			if(getCursor().getType() != Cursor.HAND_CURSOR) {
				
				setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
				getMoreVariants = true;
			}
			
		}
		else {		
			
			if(getCursor().getType() ==	Cursor.HAND_CURSOR ) {				
				setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
			}
				getMoreVariants = false;	
		}	
		
	}
	
	if(this.variantsEnd > 1 && this.variantsEnd < splits.get(0).chromEnd) {
		tempwidth = (this.variantsStart-splits.get(0).start)*splits.get(0).pixel;
		
		if(tempwidth > 0) {		
			
			buf.fillRect(splits.get(0).offset, 0, (int)tempwidth, this.getHeight());				
			varStartRect.setBounds((int)(tempwidth-this.varStringLen-10), 0, this.varStringLen, 20);
	
		}
		tempwidth = (splits.get(0).end - this.variantsEnd)*splits.get(0).pixel;
		if(tempwidth > 0) {			
			buf.fillRect((int)((this.variantsEnd+1-splits.get(0).start)*splits.get(0).pixel)+splits.get(0).offset, 0, (int)tempwidth, this.getHeight());
			varEndRect.setBounds((int)((this.variantsEnd+1-splits.get(0).start)*splits.get(0).pixel)+10, 0, this.varStringLen, 20);
		}		
	}
	else {		
		if(varEndRect.x != -1) {
			varEndRect.setBounds(-1, -1, 0,0);		
		}
	}
	
	buf.setColor(Color.white);
	
	if(this.variantsStart > splits.get(0).start) {
		
		if(varStartRect.x+ varStartRect.width +10 < this.drawWidth) {
		
			buf.drawString(this.varloadString, varStartRect.x+splits.get(0).offset, Main.defaultFontSize+2);
		}
		else {
			buf.drawString(this.varloadString, this.drawWidth+splits.get(0).offset - (varStartRect.width +10), Main.defaultFontSize+2);
			varStartRect.setBounds(this.drawWidth - varStartRect.width+10 , 0, this.varStringLen, Main.defaultFontSize*2);
		}
	}
	else {
		
		if(varStartRect.x != -1) {
			varStartRect.setBounds(-1, -1, 0,0);			
		}
	}
	
	if(variantsEnd == 1) {
		
		varEndRect.setBounds((int)splits.get(0).offset+10-Main.sidebarWidth,0, this.varStringLen, 20);
		buf.drawString(this.varloadString, splits.get(0).offset+10, Main.defaultFontSize+2);		
		buf.drawString(this.varloadhint, (int)(splits.get(0).offset+this.varloadString.length()*Main.chromDraw.bounds), Main.defaultFontSize+2);
		
		if(!Main.manage.isEnabled()) {
			Main.manage.setEnabled(true);
		}
		
	}	
	
	if(variantsEnd < splits.get(0).end) {
		
		if(varEndRect.x > 10) {			
			buf.drawString(this.varloadString, varEndRect.x+splits.get(0).offset, Main.defaultFontSize+2);
			
		}
		else {			
			buf.drawString(this.varloadString, splits.get(0).offset+10, Main.defaultFontSize+2);
			varEndRect.setBounds(10, 0, this.varStringLen, Main.defaultFontSize*2);
		}
	}
	if(!buf.getFont().equals(defaultFont)) {
		buf.setFont(defaultFont);
	}
	
}
void drawNavbar() {
	//buf.setFont(ChromDraw.seqFont);
	buf.setColor(new Color(255,255,255,100));
	buf.fillRect(pressX-50, pressY-50-Main.drawScroll.getVerticalScrollBar().getValue(), 50, 50);
	buf.fillRect(pressX, pressY-Main.drawScroll.getVerticalScrollBar().getValue(), 50, 50);
	buf.setColor(loadColor);
	buf.fillOval(pressX-50, pressY-50-Main.drawScroll.getVerticalScrollBar().getValue(), 100, 100);	
	buf.setColor(Color.black);
	buf.drawString("Zoom",pressX-50 , pressY-35-Main.drawScroll.getVerticalScrollBar().getValue());
	buf.drawString("out",pressX-50 , pressY-20-Main.drawScroll.getVerticalScrollBar().getValue());
	buf.drawString("Zoom",pressX+20 , pressY+25-Main.drawScroll.getVerticalScrollBar().getValue());
	buf.drawString("in",pressX+20 , pressY+40-Main.drawScroll.getVerticalScrollBar().getValue());
	
	if(Main.samples > 1) {
		if(drawVariables.visiblesamples > 1) {
			buf.drawString("Expand",pressX+5 , pressY-10-Main.drawScroll.getVerticalScrollBar().getValue());
		}
		if(drawVariables.sampleHeight > defaultSampleHeight && drawVariables.visiblesamples < Main.samples) {
			
			buf.drawString("Shrink",pressX-40 , pressY+15-Main.drawScroll.getVerticalScrollBar().getValue());
		}
	}
	
}


void searchMate(ReadNode readnode) {
	if(readnode.getMatePos() == -1 || clickedReadSample.getreadHash() == null || clickedReadSample.getreadHash().get(readnode.split) == null) {
		return;
	}
	
	if(readnode.getMates() == null) {
		
		for(int j = 0; j<clickedReadSample.getreadHash().get(readnode.split).getReads().size();j++) {
			
			
			read = clickedReadSample.getreadHash().get(readnode.split).getReads().get(j);
			
			
			while(read != null) {							
				
				if(read.getName().equals(readnode.getName()) && !read.equals(readnode)) {
					
					if(readnode.getMates() == null) {
						readnode.mates = new ArrayList<ReadNode>();
						
					}
					if(!readnode.getMates().contains(read)) {
						readnode.getMates().add(read);
						readnode.getMates().add(readnode);
						read.mates = readnode.getMates();					
					}
					
					read = null;
					return;
				}
				read = read.getNext();
			}
			
		}
	}
}

void drawClickedRead() {
	
	if(clickedRead != null) {
		/*if(clickedRead.getMismatches() != null) {			
			System.out.println("--");
			for(int i = 0; i<clickedRead.getMismatches().size(); i++) {
				System.out.println(clickedRead.getMismatches().get(i).getKey() +" " +Main.getBase.get(clickedRead.getMismatches().get(i).getValue()));
			}
		}*/
		//System.out.println(clickedRead.getM);
		clickedRect = clickedRead.getRect();
		if(!buf.getStroke().equals(strongStroke)) {
			buf.setStroke(strongStroke);
		}
	//	offset = clickedReadSplit.offset;	
		if(clickedReadSample.getreadHash() == null || clickedReadSample.getreadHash().get(clickedRead.split) == null) {
			clickedRead = null;
			saReads = null;
			return;
		}
		if(readwheel != clickedReadSample.getreadHash().get(clickedRead.split).readwheel) {
			readwheel = clickedReadSample.getreadHash().get(clickedRead.split).readwheel;
		}
		if(!buf.getFont().equals(defaultFont)) {
			buf.setFont(defaultFont);
			
		}
		if(readtextWidth != buf.getFontMetrics().stringWidth(clickedReadInfo[10])) {
			readtextWidth = buf.getFontMetrics().stringWidth(clickedReadInfo[10]);			
		}
		/*
		if(textWidth != null) {
			fm = buf.getFontMetrics();	
			if(textWidth.getWidth() != fm.getStringBounds(clickedRead.getName(), buf).getWidth()) {
				
				textWidth = fm.getStringBounds(clickedRead.getName(), buf);
				this.readInfoWidth = (int)textWidth.getWidth() +10;
			}
		}
		else {
			fm = buf.getFontMetrics();	
			textWidth = fm.getStringBounds(clickedRead.getName(), buf);
			this.readInfoWidth = (int)textWidth.getWidth() +10;
		}
		*/
		
	/*	if(clickedRead.getPrimary()) {
			if(!buf.getColor().equals(Color.white)) {
				buf.setColor(Color.white);
			}
		}
		else {
			if(!buf.getColor().equals(Color.orange)) {
				
			}
			
		}*/
		if(hoverMate != null) {
			if(hoverMate.equals(clickedRead)) {
				buf.setColor(Color.white);
			
			}
			else {
				buf.setColor(Color.orange);
			}
		}
		else {
			buf.setColor(Color.orange);
		}
		if(mouseDrag) {
			clickedRead.setRectBounds((int)((clickedRead.getPosition()-clickedRead.split.start)*clickedRead.split.pixel), clickedRead.getRect().y, (int)(clickedRead.getWidth()*clickedRead.split.pixel), clickedReadSample.getreadHash().get(clickedRead.split).readHeight);
		}
		if(clickedRect.x < 5 && clickedRect.x+clickedRect.width <= drawWidth) {
			buf.drawRect(4+clickedRead.split.offset, clickedRect.y-3, clickedRect.width+clickedRead.getRect().x-4 , clickedRect.height+5);
		}
		else if(clickedRect.x < 5 && clickedRect.x+clickedRect.width > drawWidth) {
			buf.drawRect(4+clickedRead.split.offset, clickedRect.y-3, drawWidth-8 , clickedRect.height+5);
			
		}
		else if(clickedRect.x+clickedRect.width > drawWidth) {
			if(clickedRead.isForward()) {
				buf.drawRect(clickedRect.x+clickedRead.split.offset-3, clickedRect.y-3, clickedRect.width+(drawWidth-(clickedRect.x+clickedRect.width))-3 , clickedRect.height+5);
			}
			else {
				buf.drawRect(clickedRect.x+clickedRead.split.offset+2, clickedRect.y-3, clickedRect.width+(drawWidth-(clickedRect.x+clickedRect.width))-3 , clickedRect.height+5);
				
			}
		}
		else {
			if(clickedRead.isForward()) {
				buf.drawRect(clickedRect.x+clickedRead.split.offset-3, clickedRect.y-3, clickedRect.width+1, clickedRect.height+5);
				
			}
			else {
				buf.drawRect(clickedRect.x+clickedRead.split.offset+2, clickedRect.y-3, clickedRect.width, clickedRect.height+5);
				
			}
		}
		
		drawTriangles(clickedRead, this.clickedReadSample,clickedRead.split,buf.getColor());
	/*	
		if( clickmates != null) {
		for(int i = 0 ; i<clickmates.size(); i++) {
			
			if(splits.indexOf(clickmates.get(i).split) < 0) {				
				
				clickmates.remove(i);	
				i--;
				continue;
			}				
			if(clickmates.get(i).equals(clickedRead)) {					
				continue;
			}
			
			clickedMateRect = clickedRead.getMates().get(i).getRect();
			clickedMateRect.setBounds((int)((clickedRead.getMates().get(i).getPosition()-clickedRead.getMates().get(i).split.start)*clickedRead.getMates().get(i).split.pixel),clickedRect.y+clickedReadSample.getreadHash().get(clickmates.get(i).split).readHeight+2,clickedMateRect.width, clickedMateRect.height);
			
			if(clickedMateRect.x < 5 && clickedMateRect.x+clickedMateRect.width <= drawWidth) {
				buf.drawRect(4+clickmates.get(i).split.offset, clickedMateRect.y-Main.drawScroll.getVerticalScrollBar().getValue(), clickedMateRect.width+clickedMateRect.x-4 , clickedMateRect.height);
			}
			else if(clickedMateRect.x < 5 && clickedMateRect.x+clickedMateRect.width > drawWidth) {
				buf.drawRect(4+clickmates.get(i).split.offset, clickedMateRect.y-Main.drawScroll.getVerticalScrollBar().getValue(), drawWidth-8 , clickedMateRect.height);
				
			}
			else if(clickedRect.x+clickedRect.width > drawWidth) {
				buf.drawRect(clickedMateRect.x+clickmates.get(i).split.offset, clickedMateRect.y-Main.drawScroll.getVerticalScrollBar().getValue(), clickedMateRect.width+(drawWidth-(clickedMateRect.x+clickedMateRect.width))-4 , clickedMateRect.height);
				
			}
			else {
				buf.drawRect(clickedMateRect.x+clickmates.get(i).split.offset, clickedMateRect.y-Main.drawScroll.getVerticalScrollBar().getValue(), clickedMateRect.width, clickedMateRect.height);
			}
			}
		}
		*/
		drawInfoBox(clickedReadInfo, 0, clickedRead, readtextWidth);
		
		buf.setStroke(strongStroke);		
	/*	if(clickedRect.width < 2) {
			
				clickedRead.setRectBounds((int)((clickedRead.getPosition()-clickedRead.split.start)*clickedRead.split.pixel), clickedRead.getRect().y, (int)(clickedRead.getWidth()*clickedRead.split.pixel), clickedReadSample.getreadHash().get(clickedRead.split).readHeight);
				buf.drawRect(clickedRect.x+clickedRead.split.offset, clickedRect.y+readwheel-Main.drawScroll.getVerticalScrollBar().getValue(), 2, clickedRect.height);
		}
		else {			
				clickedRead.setRectBounds((int)((clickedRead.getPosition()-clickedRead.split.start)*clickedRead.split.pixel), clickedRead.getRect().y, (int)(clickedRead.getWidth()*clickedRead.split.pixel), clickedReadSample.getreadHash().get(clickedRead.split).readHeight);
				if(clickedRead.getRect().x < 5 && clickedRead.getRect().x+clickedRead.getRect().width <= drawWidth) {
					buf.drawRect(4+clickedRead.split.offset, clickedRect.y+readwheel-Main.drawScroll.getVerticalScrollBar().getValue(), clickedRect.width+clickedRead.getRect().x-4 , clickedRect.height);
				}
				else if(clickedRead.getRect().x < 5 && clickedRead.getRect().x+clickedRead.getRect().width > drawWidth) {
					buf.drawRect(4+clickedRead.split.offset, clickedRect.y+readwheel-Main.drawScroll.getVerticalScrollBar().getValue(), drawWidth-8 , clickedRect.height);
					
				}
				else if(clickedRead.getRect().x+clickedRead.getRect().width > drawWidth) {
					buf.drawRect(clickedRect.x+clickedRead.split.offset, clickedRect.y+readwheel-Main.drawScroll.getVerticalScrollBar().getValue(), clickedRect.width+(drawWidth-(clickedRead.getRect().x+clickedRead.getRect().width))-4 , clickedRect.height);
					
				}
				else {
					buf.drawRect(clickedRect.x+clickedRead.split.offset, clickedRect.y+readwheel-Main.drawScroll.getVerticalScrollBar().getValue(), clickedRect.width, clickedRect.height);
				}
		}		
		*/
		
		buf.setColor(Color.orange);
		
		if(clickedRead.getMates() != null) {
			
			for(int i = 0 ; i<clickedRead.getMates().size(); i++) {
			/*	if(clickedRead.equals(clickedRead.getMates().get(i))) {
					continue;
				}
				*/
				
				if(splits.indexOf(clickedRead.getMates().get(i).split) < 0) {					
					clickedRead.getMates().remove(i);	
					i--;
					continue;
				}		
				
				if(hoverMate != null) {
					if(hoverMate.equals(clickedRead.getMates().get(i))) {
						buf.setColor(Color.white);
					
					}
					else {
						buf.setColor(Color.orange);
					}
				}
				else {
					buf.setColor(Color.orange);
				}
				clickedMateRect = clickedRead.getMates().get(i).getRect();
				if(mouseDrag) {
					clickedRead.getMates().get(i).setRectBounds((int)((clickedRead.getMates().get(i).getPosition()-clickedRead.getMates().get(i).split.start)*clickedRead.getMates().get(i).split.pixel), clickedMateRect.y, (int)(clickedRead.getMates().get(i).getWidth()*clickedRead.getMates().get(i).split.pixel), clickedReadSample.getreadHash().get(clickedRead.getMates().get(i).split).readHeight);
				}
				drawTriangles(clickedRead.getMates().get(i), this.clickedReadSample,clickedRead.getMates().get(i).split,buf.getColor());
				
				if(clickedRead.getMates().get(i).split != null) {
					
				//	if(!clickedRead.equals(clickedRead.getMates().get(i))) {
						if(clickedMateRect.width < 2) {
							if(clickedRead.getMates().get(i).split.equals(selectedSplit)) {
					//			clickedRead.getMates().get(i).setRectBounds((int)((clickedRead.getMates().get(i).getPosition()-clickedRead.getMates().get(i).split.start)*clickedRead.getMates().get(i).split.pixel), clickedRead.getMates().get(i).getRect().y, (int)(clickedRead.getMates().get(i).getWidth()*clickedRead.getMates().get(i).split.pixel), clickedReadSample.getreadHash().get( clickedRead.getMates().get(i).split).readHeight);
							
								buf.drawRect(clickedMateRect.x+clickedRead.getMates().get(i).split.offset,clickedRead.getMates().get(i).yposition-3, 2, clickedReadSample.getreadHash().get( clickedRead.getMates().get(i).split).readHeight);
							}						
						}
						else {						
							clickedRead.getMates().get(i).setRectBounds((int)((clickedRead.getMates().get(i).getPosition()-clickedRead.getMates().get(i).split.start)*clickedRead.getMates().get(i).split.pixel), clickedRead.getMates().get(i).rect.y, (int)(clickedRead.getMates().get(i).getWidth()*clickedRead.getMates().get(i).split.pixel), clickedReadSample.getreadHash().get( clickedRead.getMates().get(i).split).readHeight);
							
							if(clickedMateRect.x < 5 && clickedMateRect.x+clickedMateRect.width <= drawWidth) {
								buf.drawRect(4+ clickedRead.getMates().get(i).split.offset, clickedRead.getMates().get(i).yposition-3, clickedMateRect.width+clickedMateRect.x-4 , clickedMateRect.height+5);
							}
							else if(clickedMateRect.x < 5 && clickedMateRect.x+clickedMateRect.width > drawWidth) {
								buf.drawRect(4+ clickedRead.getMates().get(i).split.offset, clickedRead.getMates().get(i).yposition-3, drawWidth-8 , clickedMateRect.height+5);
								
							}
							else if(clickedMateRect.x+clickedMateRect.width > drawWidth) {
								buf.drawRect(clickedMateRect.x+ clickedRead.getMates().get(i).split.offset, clickedRead.getMates().get(i).yposition-3, clickedMateRect.width+(drawWidth-(clickedMateRect.x+clickedMateRect.width))-4 , clickedMateRect.height+5);
								
							}
							else {
								if(clickedRead.isForward()) {
									buf.drawRect(clickedMateRect.x+ clickedRead.getMates().get(i).split.offset-3, clickedRead.getMates().get(i).yposition-3, clickedMateRect.width+1, clickedMateRect.height+5);
								}
								else {
									buf.drawRect(clickedMateRect.x+ clickedRead.getMates().get(i).split.offset+2, clickedRead.getMates().get(i).yposition-3, clickedMateRect.width, clickedMateRect.height+5);
									
								}
							}
						}		
					
					testRect = clickedMateRect;
					testRect.setBounds(testRect.x, clickedRead.getMates().get(i).yposition, testRect.width, testRect.height);
					
					if(mouseRect.intersects(testRect) && clickedRead.getMates().get(i).split.equals(selectedSplit)) {
					
						if(hoverMate == null || !hoverMate.equals(clickedRead.getMates().get(i))) {
							hoverMate = clickedRead.getMates().get(i);
							if(!clickedRead.equals(clickedRead.getMates().get(i))) {
								hoverMateInfo = createReadInfo(hoverMate);	
								if(saReads != null && saReads.reads != null) {
									for(int r = 0 ; r<saReads.reads.size(); r++) {
										if(saReads.reads.get(r)[SAread.SARead] == null) {
											if((int)saReads.reads.get(r)[SAread.pos] == hoverMate.getPosition()) {
												
												saReads.reads.get(r)[SAread.SARead] = hoverMate;
												break;
											}
										}
									}
								}
							}
						}
						if(!clickedRead.equals(clickedRead.getMates().get(i))) {
							drawInfoBox(hoverMateInfo, readtextWidth+30,hoverMate, buf.getFontMetrics().stringWidth(hoverMateInfo[10]));
							buf.setColor(Color.orange);
							buf.setStroke(strongStroke);	
						}
					}
		//			}
					
					if(i > 0) { // && !clickedRead.getMates().get(i-1).equals(clickedRead.getMates().get(i))) {
						
						if(clickedRead.getMates().get(i-1).split.equals(clickedRead.getMates().get(i).split)) {
							
								if(clickedRead.getMates().get(i-1).getRect().x+clickedRead.getMates().get(i-1).getRect().width > clickedRead.getMates().get(i).getRect().x) {
									
									preread = clickedRead.getMates().get(i-1).getRect().x + clickedRead.getMates().get(i-1).split.offset;	
									thisread = clickedRead.getMates().get(i).getRect().x + clickedRead.getMates().get(i).split.offset;
									
								}
								else {
									preread = clickedRead.getMates().get(i-1).getRect().x +clickedRead.getMates().get(i-1).getRect().width + clickedRead.getMates().get(i-1).split.offset;
									thisread = clickedRead.getMates().get(i).getRect().x + clickedRead.getMates().get(i).split.offset;
									
								}						
							
						}
						else {
							preread = clickedRead.getMates().get(i-1).getRect().x +clickedRead.getMates().get(i-1).getRect().width + clickedRead.getMates().get(i-1).split.offset;
							thisread = clickedRead.getMates().get(i).getRect().x + clickedRead.getMates().get(i).split.offset;
							
						}	
						
						if(!clickedRead.getMates().get(i-1).split.equals(clickedRead.getMates().get(i).split)) {
							if(sampleList.get(selectedIndex).getreadHash().get(clickedRead.getMates().get(i-1).split) == null) {
								clickedRead.getMates().remove(i-1);
								continue;
							}		
							
							curve.setCurve(preread, clickedRead.getMates().get(i-1).yposition, Math.abs(thisread-preread)/2+preread,   clickedRead.getMates().get(i-1).yposition-100, thisread, clickedRead.getMates().get(i).yposition);					
						}
						else {
						
							curve.setCurve(preread, clickedRead.getMates().get(i-1).yposition, Math.abs(thisread-preread)/2+preread,  clickedRead.getMates().get(i).yposition-100, thisread, clickedRead.getMates().get(i).yposition);					
							
						}
						buf.setColor(Color.white);
						buf.setStroke(doubledashed);
						buf.draw(curve);
						buf.setColor(Color.orange);
						buf.setStroke(strongStroke);	
						
					}	
				}			
			}
				
		}	
	}
	buf.setStroke(basicStroke);	
}

void drawInfoBox(String[] readinfo, int boxoffset, ReadNode read, int boxwidth) {
	if(boxwidth < 320) {
		boxwidth = 320;
	}
	//TODO SCALE
	buf.setColor(infoBoxColor);	
	if(saReads == null || boxoffset > 100) {
		buf.fillRect(10+Main.sidebarWidth+boxoffset,10,boxwidth+10,(readinfo.length)*21);
		
		buf.drawRect(10+Main.sidebarWidth+boxoffset,10,boxwidth+10,(readinfo.length)*21);
		boxheight = (readinfo.length)*21 +10;
	}
	else {
		buf.fillRect(10+Main.sidebarWidth+boxoffset,10,boxwidth+10,(readinfo.length)*21+saReads.reads.size()*21);
		
		buf.drawRect(10+Main.sidebarWidth+boxoffset,10,boxwidth+10,(readinfo.length)*21+saReads.reads.size()*21);
		boxheight = (readinfo.length)*21 +10+saReads.reads.size()*21;
	}
	buf.setStroke(dashed);
	
	if(read.getRect().x > drawWidth || read.getRect().x < 0) {
		
	}
	else if(read.getRect().x > 0) {
		if(read.getRect().x+read.split.offset-Main.sidebarWidth > boxoffset + boxwidth +20) {
			buf.drawLine((int)((read.getPosition()-read.split.start)*read.split.pixel)+read.split.offset, read.getRect().y, boxwidth+boxoffset+20+Main.sidebarWidth, 10);
		}
		else {
			buf.drawLine((int)((read.getPosition()-read.split.start)*read.split.pixel)+read.split.offset, read.getRect().y, boxwidth+boxoffset+20+Main.sidebarWidth, boxheight);
		}
		if(read.getRect().y > boxheight) {
			buf.drawLine((int)((read.getPosition()-read.split.start)*read.split.pixel)+read.split.offset, read.getRect().y, 10+boxoffset+Main.sidebarWidth, boxheight);
		}
		else {
			buf.drawLine((int)((read.getPosition()-read.split.start)*read.split.pixel)+read.split.offset, read.getRect().y, boxwidth+20+boxoffset+Main.sidebarWidth, boxheight);
		}
	}
	else {
		if(read.getRect().x+read.getRect().width+read.split.offset-Main.sidebarWidth > boxoffset + boxwidth +20) {
			buf.drawLine((int)((read.getPosition()-read.split.start)*read.split.pixel)+read.split.offset+read.getRect().width, read.getRect().y, boxwidth+boxoffset+20+Main.sidebarWidth, 10);
		}
		else {
			buf.drawLine((int)((read.getPosition()-read.split.start)*read.split.pixel)+read.split.offset+read.getRect().width, read.getRect().y, boxwidth+boxoffset+20+Main.sidebarWidth, boxheight);
		}
		if(clickedRect.y > boxheight) {
			buf.drawLine((int)((read.getPosition()-read.split.start)*read.split.pixel)+read.split.offset+read.getRect().width, read.getRect().y, 20+boxoffset+Main.sidebarWidth, boxheight);
		}
		else {
			buf.drawLine((int)((read.getPosition()-read.split.start)*read.split.pixel)+read.split.offset+read.getRect().width, read.getRect().y, boxwidth+boxoffset+20+Main.sidebarWidth, boxheight);
		}
	}
	buf.setStroke(basicStroke);
	buf.setColor(Color.black);	
	for(int i = 0; i<readinfo.length-1; i++) {
		buf.drawString(readinfo[i], 15+Main.sidebarWidth+boxoffset, 25+i*20);
	}
	if(saReads != null && boxoffset < 100) {
		
		//buf.fillRect(15+Main.sidebarWidth+boxoffset,readinfo.length*20-3, boxwidth, saReads.reads.size()*20);
		buf.setColor(Color.lightGray);
		
		buf.fillRect(15+Main.sidebarWidth+boxoffset,readinfo.length*20-20, boxwidth-10, 15);		
		buf.fillPolygon(makeTriangle(15+Main.sidebarWidth+boxoffset+(boxwidth-10), readinfo.length*20-23, 10,21,true));
		buf.setColor(Color.black);
		buf.drawRect(15+Main.sidebarWidth+boxoffset,readinfo.length*20-20, boxwidth-10, 15);		
		buf.drawPolygon(makeTriangle(15+Main.sidebarWidth+boxoffset+(boxwidth-10), readinfo.length*20-23, 10,21,true));
		buf.drawString("Fragment length: " +saReads.readlength +"bp", 15+Main.sidebarWidth+boxoffset+5,readinfo.length*20-8);
		for(int i = 0 ; i<saReads.reads.size(); i++) {
			buf.setColor(getChromColor((String)(saReads.reads.get(i)[SAread.chrom])));
			if((boolean)saReads.reads.get(i)[SAread.forward]) {
				saReadRect.setBounds(15+Main.sidebarWidth+boxoffset +(int)((double)saReads.reads.get(i)[SAread.relativepos]*boxwidth),readinfo.length*20+(6+i*21), (int)((double)saReads.reads.get(i)[SAread.relativelength]*boxwidth), 15);
				buf.fillRect(saReadRect.x,saReadRect.y, saReadRect.width-5, 15);
				buf.fillPolygon(makeTriangle(saReadRect.x+saReadRect.width-5, saReadRect.y-3, 5, 21, true));
					if(saReads.reads.get(i)[SAread.SARead] != null) {
						if(hoverMate != null && saReads.reads.get(i)[SAread.SARead].equals(hoverMate)) {
							
							buf.setColor(Color.white);
							buf.setStroke(doubleStroke);
						}
						else if(saReads.reads.get(i)[SAread.SARead].equals(saReads.read)) {
							buf.setColor(Color.orange);
							buf.setStroke(doubleStroke);
							
						}
						
						else {
							buf.setColor(Color.black);
							buf.setStroke(basicStroke);
						}
					buf.drawRect(saReadRect.x,saReadRect.y, saReadRect.width-5, 15);
					buf.drawPolygon(makeTriangle(saReadRect.x+saReadRect.width-5, saReadRect.y-3, 5, 21, true));
					buf.setStroke(basicStroke);
				}
				else {
					buf.setColor(Color.black);
					buf.drawRect(saReadRect.x,saReadRect.y, saReadRect.width-5, 15);
					buf.drawPolygon(makeTriangle(saReadRect.x+saReadRect.width-5, saReadRect.y-3, 5, 21, true));
				}
					saReadRect.setBounds(15+boxoffset +(int)((double)saReads.reads.get(i)[SAread.relativepos]*boxwidth),saReadRect.y, saReadRect.width, saReadRect.height);
					
					if(mouseRect.intersects(saReadRect)) {
						if(hoverMate == null || !hoverMate.equals((ReadNode)saReads.reads.get(i)[SAread.SARead])) {
							hoverMate = (ReadNode)saReads.reads.get(i)[SAread.SARead];
							hoverRect.setBounds(15+boxoffset +(int)((double)saReads.reads.get(i)[SAread.relativepos]*boxwidth),saReadRect.y, saReadRect.width, saReadRect.height);
							
						}
					
						if(hoverMate == null && clickedRead.getMates() != null) {
							for(int m = 0 ; m<clickedRead.getMates().size();m++) {						
								for(int r = 0 ; r<saReads.reads.size(); r++) {
									if((int)saReads.reads.get(r)[SAread.pos] == clickedRead.getMates().get(m).getPosition()) {										
										saReads.reads.get(r)[SAread.SARead] = clickedRead.getMates().get(m);
										break;
									}							
								}
							}
						}
					}
			}
			else {
				saReadRect.setBounds(20+Main.sidebarWidth+boxoffset +(int)((double)saReads.reads.get(i)[SAread.relativepos]*boxwidth),readinfo.length*20+(6+i*21), (int)((double)saReads.reads.get(i)[SAread.relativelength]*boxwidth), 15);
				
				buf.fillRect(saReadRect.x,saReadRect.y, saReadRect.width-5, 15);
				buf.fillPolygon(makeTriangle(saReadRect.x, saReadRect.y-3, 5, 21, false));
				
				if(saReads.reads.get(i)[SAread.SARead] != null) {
				
					if(hoverMate != null && saReads.reads.get(i)[SAread.SARead].equals(hoverMate)) {
						
						buf.setColor(Color.white);
						buf.setStroke(doubleStroke);
					}
					else if(saReads.reads.get(i)[SAread.SARead].equals(saReads.read)) {
						buf.setColor(Color.orange);
						buf.setStroke(doubleStroke);
						
					}
					
					else {
						buf.setColor(Color.black);
						buf.setStroke(basicStroke);
					}
					
					buf.drawRect(saReadRect.x,saReadRect.y, saReadRect.width-5, 15);
					buf.drawPolygon(makeTriangle(saReadRect.x, saReadRect.y-3, 5, 21, false));
					buf.setStroke(basicStroke);
				}
				else {
					buf.setColor(Color.black);
					buf.drawRect(saReadRect.x,saReadRect.y, saReadRect.width-5, 15);
					buf.drawPolygon(makeTriangle(saReadRect.x, saReadRect.y-3, 5, 21, false));
				}
				saReadRect.setBounds(15+boxoffset +(int)((double)saReads.reads.get(i)[SAread.relativepos]*boxwidth),saReadRect.y, saReadRect.width, saReadRect.height);
				
				if(mouseRect.intersects(saReadRect)) {
					if(hoverMate == null || !hoverMate.equals((ReadNode)saReads.reads.get(i)[SAread.SARead])) {
						hoverMate = (ReadNode)saReads.reads.get(i)[SAread.SARead];
						hoverRect.setBounds(15+boxoffset +(int)((double)saReads.reads.get(i)[SAread.relativepos]*boxwidth),saReadRect.y, saReadRect.width, saReadRect.height);
						
						
					}
				
					if(hoverMate == null && clickedRead.getMates() != null) {
						for(int m = 0 ; m<clickedRead.getMates().size();m++) {						
							for(int r = 0 ; r<saReads.reads.size(); r++) {
								if((int)saReads.reads.get(r)[SAread.pos] == clickedRead.getMates().get(m).getPosition()) {										
									saReads.reads.get(r)[SAread.SARead] = clickedRead.getMates().get(m);
									break;
								}							
							}
						}
					}
				}
			}			
		}
	}
}	

void eraseReads() {
	
	for(int i = 0;i<splits.size(); i++) {
	
		splits.get(i).getSelectbuf().setComposite( composite);				
		splits.get(i).getSelectbuf().fillRect(0,0, (int)Main.screenSize.getWidth(), (int)Main.screenSize.getHeight());	
		splits.get(i).getSelectbuf().setComposite(splits.get(i).getBackups());
		splits.get(i).getReadBuffer().setComposite( composite);		
		splits.get(i).getReadBuffer().fillRect(0,0, (int)Main.screenSize.getWidth(),  (int)Main.screenSize.getHeight());	
		splits.get(i).getReadBuffer().setComposite(splits.get(i).getBackupr());	
		rbuf.drawImage(splits.get(i).getReadImage(), splits.get(i).offset, 0,this);	
		rbuf.drawImage(splits.get(i).getSelectbuffer(), splits.get(i).offset, 0, this);
		buf.drawImage(readbuffer, 0, 0, this);	
		
	}
	
}
/*
void eraseReads(SplitClass split) {

	split.getSelectbuf().setComposite( composite);				
	split.getSelectbuf().fillRect(0,0, this.getDrawWidth(), Main.drawScroll.getViewport().getHeight());	
	split.getSelectbuf().setComposite(split.getBackups());
	split.getReadBuffer().setComposite( composite);		
	split.getReadBuffer().fillRect(0,0, (int)Main.screenSize.getWidth(),  (int)Main.screenSize.getHeight());	
	split.getReadBuffer().setComposite(split.getBackupr());	
	
}
*/
void clearReads(SplitClass split) {
	
	try {
		
		if(split.clearedReads) {
			return;
		}
		
		split.getSelectbuf().setComposite( composite);		
		split.getReadBuffer().setComposite( composite);	
		split.getReadBuffer().fillRect(0,0,(int)Main.screenSize.getWidth(), Main.drawScroll.getViewport().getHeight());	
		split.getSelectbuf().fillRect(0,0,(int)Main.screenSize.getWidth(), Main.drawScroll.getViewport().getHeight());	
		split.getSelectbuf().setComposite(split.getBackups());		
		split.getSelectbuf().setComposite(split.getBackupr());		
		Reads reads;
	
	for(int i = 0 ; i< this.sampleList.size(); i++) {
		if(sampleList.get(i).removed) {
			continue;
		}
		sample = this.sampleList.get(i);
		if(sample.getMates() != null) {
			sample.getMates().clear();
		}
		if(sample.getreadHash() == null) {
			break;
		}
		reads = sample.getreadHash().get(split);
		
		if(reads == null) {
			continue;
		}	       
	       
        if(reads.getReads() != null) {
			reads.getReads().clear();
			reads.getHeadAndTail().clear();
		}
		reads.setReadEnd(0);
		reads.setReadStart(Integer.MAX_VALUE);
		reads.setFirstRead(null);
		reads.setLastRead(null);
		reads.setCoverageStart(Integer.MAX_VALUE);
		reads.setCoverageEnd(0);
		reads.searchstart = Integer.MAX_VALUE;
		reads.searchend = 0;

	}	
	
	reads = null;	
	row = null;
	read = null;
	selectedRead = null;
	clickedRead = null;
	saReads = null;
	splitList.clear();
	selectedMate = null;
	split.clearedReads = true;

	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
	}
}
void clearReads() {
	
	try {
	
	for(int i = 0; i<splits.size();i++) {
		splits.get(i).getSelectbuf().setComposite( composite);		
		splits.get(i).getSelectbuf().fillRect(0,0,(int)Main.screenSize.getWidth(), Main.drawScroll.getViewport().getHeight());	
		splits.get(i).getSelectbuf().setComposite(splits.get(i).getBackups());
		
	}
	
	 Iterator<Map.Entry<SplitClass, Reads>> it;
	 Map.Entry<SplitClass, Reads> pair;
	 Reads reads;
	
	for(int i = 0 ; i< this.sampleList.size(); i++) {
		if(sampleList.get(i).removed) {
			continue;
		}
		sample = this.sampleList.get(i);
		if(sample.getreadHash()== null) {
			continue;
		}
		if(sample.getMates() != null) {
			sample.getMates().clear();
		}
		it = sample.getreadHash().entrySet().iterator();
	    while (it.hasNext()) {
	        pair = (Map.Entry<SplitClass, Reads>)it.next();
	        reads = pair.getValue();
	        if(reads.getReads() != null) {
				reads.getReads().clear();
				reads.getHeadAndTail().clear();
			}
	        reads.searchend = 0;
	        reads.searchstart = Integer.MAX_VALUE;
			reads.setReadEnd(0);
			reads.setReadStart(Integer.MAX_VALUE);
			reads.setFirstRead(null);
			reads.setLastRead(null);
			reads.setCoverageStart(Integer.MAX_VALUE);
			reads.setCoverageEnd(0);
//	        it.remove(); // avoids a ConcurrentModificationException
	    }
	/*	for(Reads reads : sample.getreadHash()) {
			if(reads.getReads() != null) {
				reads.getReads().clear();
				reads.getHeadAndTail().clear();
			}
			reads.setReadEnd(0);
			reads.setReadStart(Integer.MAX_VALUE);
			reads.setFirstRead(null);
			reads.setLastRead(null);
		
			
	//		reads.getHelpArray().clear();
	//		reads.getHelpHash().clear();	
	}*/
		
	//	sample.getreadHash().readNames.clear();
		/*
		if(sample.getreadHash().reads != null) {
			
			sample.getreadHash().reads.clear();
			sample.getreadHash().headAndTail.clear();
			sample.getreadHash().setFirstRead(null);
			sample.getreadHash().setLastRead(null);			
	//		this.sampleList.get(i).getreadHash().readNames.clear();
		}
		*/
	}
	
	pair = null;
	reads = null;
	it = null;
	row = null;
	read = null;
	selectedRead = null;
	clickedRead = null;
	saReads = null;
//	clickmates = null;
	splitList.clear();
	selectedMate = null;	
	selectedSplit.clearedReads = true;
	
	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
	}
}
void clearReads(Sample sample) {
	if(sample.getreadHash() == null ) {
		return;
	}
	Iterator<Map.Entry<SplitClass, Reads>> it;
	 Map.Entry<SplitClass, Reads> pair;
	 Reads reads;
	it = sample.getreadHash().entrySet().iterator();
    while (it.hasNext()) {
        pair = (Map.Entry<SplitClass, Reads>)it.next();
        reads = pair.getValue();
        if(reads.getReads() != null) {
			reads.getReads().clear();
			reads.getHeadAndTail().clear();
		}
		reads.setReadEnd(0);
		reads.setReadStart(Integer.MAX_VALUE);
		reads.setFirstRead(null);
		reads.setLastRead(null);
    //    reads.reference = null;
        it.remove(); // avoids a ConcurrentModificationException
    }
    selectedRead = null;
//	clickedRead = null;
//	clickmates = null;
	splitList.clear();
//	clickedReadSplit = null;
/*	for(Reads reads : sample.getreadHash()) {
		reads.getReads().clear();
		reads.getHeadAndTail().clear();
		reads.setFirstRead(null);
		reads.setLastRead(null);	
	//	reads.setCoverages(null);
	}
	selectedRead = null;
	clickedRead = null;
	clickedReadSplit = 0;
	/*
	if(sample.getreadHash().reads != null) {
		
		sample.getreadHash().reads.clear();
		sample.getreadHash().headAndTail.clear();
		sample.getreadHash().setFirstRead(null);
		sample.getreadHash().setLastRead(null);			
	}
	*/
}

Polygon makeTriangle(int x, int y, int width, int height, boolean right) {
	if(right) {
		p1 = new Point(x, y);
		p2 = new Point(x, y+height);
		p3 = new Point(x+width, (int)(y+(height/2.0)));						
		xs[0] = p1.x;
		xs[1] = p2.x;
		xs[2] = p3.x;						
		ys[0] = p1.y;
		ys[1] = p2.y;
		ys[2] = p3.y;
	}
	else {
		p1 = new Point(x, y);
		p2 = new Point(x, y+height);
		p3 = new Point(x-width, (int)(y+(height/2.0)));						
		xs[0] = p1.x;
		xs[1] = p2.x;
		xs[2] = p3.x;						
		ys[0] = p1.y;
		ys[1] = p2.y;
		ys[2] = p3.y;
	}
	return new Polygon(xs, ys, xs.length);
						
}

void drawTriangles(ReadNode read, Sample sample, SplitClass split, Color color) {
	
	if(read.equals(clickedRead) || read.equals(selectedRead) || (clickedRead != null && read.getMates() != null && read.getMates().equals(clickedRead.getMates()))) {
		
		
		buf.setColor(color);
		
	//	if(read.getMates() == null) {
			if(read.isForward()) {
				p1 = new Point((int)read.getRect().x+read.getRect().width+4, read.getRect().y-2);
				p2 = new Point((int)read.getRect().x+read.getRect().width+4, read.getRect().y+sample.getreadHash().get(split).readHeight+2);
				p3 = new Point((int)read.getRect().x+read.getRect().width+12, (int)(read.getRect().y+((sample.getreadHash().get(split).readHeight+4)/2.0)));						
				xs[0] = p1.x+split.offset;
				xs[1] = p2.x+split.offset;
				xs[2] = p3.x+split.offset;						
				ys[0] = p1.y;
				ys[1] = p2.y;
				ys[2] = p3.y;																
				triangle = new Polygon(xs, ys, xs.length);
			
				buf.fillPolygon(triangle);
										
			}
			else {
				p1 = new Point((int)read.getRect().x-2, read.getRect().y-2);
				p2 = new Point((int)read.getRect().x-2, read.getRect().y+sample.getreadHash().get(split).readHeight+2);
				p3 = new Point((int)read.getRect().x-10, (int)(read.getRect().y+((sample.getreadHash().get(split).readHeight+4)/2.0)));						
				xs[0] = p1.x+split.offset;
				xs[1] = p2.x+split.offset;
				xs[2] = p3.x+split.offset;						
				ys[0] = p1.y;
				ys[1] = p2.y;
				ys[2] = p3.y;																
				triangle = new Polygon(xs, ys, xs.length);							
				
				buf.fillPolygon(triangle);
				
			}
	//	}
	/*	else {
			for (int i = 0; i<read.getMates().size(); i++) {
				read = read.getMates().get(i);
				if(read.isForward()) {
					p1 = new Point((int)read.getRect().x+read.getRect().width+4, read.getRect().y-2);
					p2 = new Point((int)read.getRect().x+read.getRect().width+4, read.getRect().y+sample.getreadHash().get(split).readHeight+2);
					p3 = new Point((int)read.getRect().x+read.getRect().width+12, (int)(read.getRect().y+((sample.getreadHash().get(split).readHeight+4)/2.0)));						
					xs[0] = p1.x+read.split.offset;
					xs[1] = p2.x+read.split.offset;
					xs[2] = p3.x+read.split.offset;						
					ys[0] = p1.y;
					ys[1] = p2.y;
					ys[2] = p3.y;																
					triangle = new Polygon(xs, ys, xs.length);
					
					buf.fillPolygon(triangle);
												
				}
				else {
					p1 = new Point((int)read.getRect().x-2, read.getRect().y-2);
					p2 = new Point((int)read.getRect().x-2, read.getRect().y+sample.getreadHash().get(split).readHeight+2);
					p3 = new Point((int)read.getRect().x-10, (int)(read.getRect().y+((sample.getreadHash().get(split).readHeight+4)/2.0)));						
					xs[0] = p1.x+read.split.offset;
					xs[1] = p2.x+read.split.offset;
					xs[2] = p3.x+read.split.offset;						
					ys[0] = p1.y;
					ys[1] = p2.y;
					ys[2] = p3.y;																
					triangle = new Polygon(xs, ys, xs.length);							
					
					buf.fillPolygon(triangle);
					
				}
			}
		}*/
	}
	
	else if(read.isForward()) {
		p1 = new Point((int)read.getRect().x+drawRect.width-5, read.getRect().y);
		p2 = new Point((int)read.getRect().x+drawRect.width-5, read.getRect().y+sample.getreadHash().get(split).readHeight);
		p3 = new Point((int)read.getRect().x+read.getRect().width, (int)(read.getRect().y+(sample.getreadHash().get(split).readHeight/2.0)));						
		xs[0] = p1.x;
		xs[1] = p2.x;
		xs[2] = p3.x;						
		ys[0] = p1.y;
		ys[1] = p2.y;
		ys[2] = p3.y;																
		triangle = new Polygon(xs, ys, xs.length);							
		split.getReadBuffer().fillPolygon(triangle);
							
	}
	else {
		p1 = new Point((int)read.getRect().x+5, read.getRect().y);
		p2 = new Point((int)read.getRect().x+5, read.getRect().y+sample.getreadHash().get(split).readHeight);
		p3 = new Point((int)read.getRect().x, (int)(read.getRect().y+(sample.getreadHash().get(split).readHeight/2.0)));						
		xs[0] = p1.x;
		xs[1] = p2.x;
		xs[2] = p3.x;						
		ys[0] = p1.y;
		ys[1] = p2.y;
		ys[2] = p3.y;																
		triangle = new Polygon(xs, ys, xs.length);							
		split.getReadBuffer().fillPolygon(triangle);
	}
}

void drawReads(SplitClass split) {	
	
	
	try {	
		
		if(variantcalculator) {
			return;
		}
		
		if(bam && split.viewLength > Settings.coverageDrawDistance) {
			if(!zoomNote) {
				zoomNote = true;
			}
		}
		else {
			if(zoomNote) {
				zoomNote = false;
			}
		}
		
		if(scrollbar) {
			
			split.getReadBuffer().setComposite( composite);					
			split.getReadBuffer().fillRect(0,0, split.getReadImage().getWidth(),Main.drawScroll.getViewport().getHeight());	
			split.getReadBuffer().setComposite(split.getBackupr());		
			return;
		}
		
		if(lineZoomer || sampleZoomer) {
		
			split.getReadBuffer().setComposite( composite);					
			split.getReadBuffer().fillRect(0,0, split.getReadImage().getWidth(),Main.drawScroll.getViewport().getHeight());	
			split.getReadBuffer().setComposite(split.getBackupr());		
			rbuf.setComposite( composite);				
			rbuf.fillRect(0,0, this.getWidth(),Main.drawScroll.getViewport().getHeight());	
			rbuf.setComposite(backupr);		
		}
		if(!lineZoomer && (split.viewLength > Settings.readDrawDistance || drawVariables.sampleHeight <= 100) ) {		
			if(!split.clearedReads) {
				updateReads = false;			
				clearReads(split);
				
			}			
		}	
		
		if(split.viewLength >= Settings.readDrawDistance || drawVariables.sampleHeight <= 100) {			
			updateReads = false;			
		}
		
		if(split.viewLength < Settings.coverageDrawDistance && drawVariables.sampleHeight > 100) {
			
	
		if(Main.readsamples > 0 && !FileRead.cancelreadread) {			
			
		if(updateReads || split.updateReads || updateCoverages) {
		
		if(mouseWheel) {
		
			split.getReadBuffer().setComposite( composite);					
			split.getReadBuffer().fillRect(0,(int)(selectedIndex*drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue(), drawWidth,(int)drawVariables.sampleHeight);	
			split.getReadBuffer().setComposite(split.getBackupr());					
			split.getReadBuffer().setStroke(basicStroke);
			
		}
		else {
			split.getReadBuffer().setComposite( composite);					
			split.getReadBuffer().fillRect(0,0, split.getReadImage().getWidth(),Main.drawScroll.getViewport().getHeight());	
			split.getReadBuffer().setComposite(split.getBackupr());					
			split.getReadBuffer().setStroke(basicStroke);
			
		}
		
	//		firstreadsample = true;
		if(!split.splitRead) {
			
			if(previsiblestart != drawVariables.visiblestart || previsiblesamples != drawVariables.visiblesamples) {
				
				for(int i = 0; i<drawVariables.visiblestart; i++) {
					clearReads(sampleList.get(i));
				}
				for(int i = drawVariables.visiblestart+drawVariables.visiblesamples; i<sampleList.size(); i++) {
					clearReads(sampleList.get(i));
				}
			}
		}
		
		
		//boolean first = true;
		
		if(split.viewLength <= Settings.readDrawDistance && split.viewLength > 10){
			try {
				
				Main.chromDraw.getReadSeq(split);
			}
			catch(Exception ex) {
				ex.printStackTrace();
			}
			
		}
		
		if((!lineZoomer && !mouseDrag) || split.clearedReads ) {
			
			if(!sampleZoomer && !loading && !scrollbar) {
				
				if(prestart != (int)split.start || previsiblestart != drawVariables.visiblestart || previsiblesamples != drawVariables.visiblesamples || split.clearedReads) {
					if(!split.splitRead) {
						getReads(split);	
					}
					else {
						split.splitRead = false;
					}
				}
				prestart = (int)split.start;
				
			}
		}
		previsiblestart = drawVariables.visiblestart;
		previsiblesamples = drawVariables.visiblesamples;
		for(int i =drawVariables.visiblestart; i<drawVariables.visiblestart+drawVariables.visiblesamples;i++ ) {
			
			try {	
			if(i >= Main.samples) {
				break;
			}			
			
			if(!sidebar && mouseWheel && i != selectedIndex) {
				continue;
			}
			
			if(sampleList.get(i).samFile == null) {			
				continue;
			}
			readsample = sampleList.get(i);
			if(readsample.getreadHash().size() < splits.size()) {
				readsample.resetreadHash();
			}
			
			if(!scrolldrag && !FileRead.cancelreadread) { 
				
			if(split.viewLength <= Settings.coverageDrawDistance && split.viewLength > Settings.readDrawDistance) {
				
				if(!lineZoomer && !mouseDrag && ((readsample.getreadHash().get(split).getCoverageStart() > (int)split.start || readsample.getreadHash().get(split).getCoverageEnd() < (int)split.end) || (readsample.getreadHash().get(split).getCoverageEnd() > (int)split.end && readsample.getreadHash().get(split).getCoverageStart() < (int)split.start))) {
					//clearReads();					
					FileRead reader = new FileRead();
					readsample.getreadHash().get(split).setReadStart(Integer.MAX_VALUE);
					readsample.getreadHash().get(split).setReadEnd(0);
					reader.readClass = readsample.getreadHash().get(split);					
					//readsample.getreadHash().get(split).loading = true;
					reader.chrom = split.chrom;
					reader.start = (int)(split.start);
					reader.end = (int)(split.end);
					reader.viewLength = (int)(split.viewLength);
					reader.pixel = split.pixel;				
					reader.splitIndex = split;						
					reader.getreads = true;							
					reader.execute();
					updateCoverages = false;
					split.clearedCoverages = false;
				
				}
			}
			else if(split.viewLength <= Settings.readDrawDistance) {
				
			/*	if(readsample.getreadHash().get(split).searchend-readsample.getreadHash().get(split).searchstart > Settings.readDrawDistance && (split.start - readsample.getreadHash().get(split).searchstart > split.viewLength*2 || readsample.getreadHash().get(split).searchend - split.end > split.viewLength*2)) {
					
					clearReads();	
					
					updateReads = true;
					
					drawReads(split);
					return;
				}
				*/
				/*
				if(!readsample.getreadHash().get(split).loading && !lineZoomer && !mouseDrag && (readsample.getreadHash().get(split).searchend < (int)split.end || readsample.getreadHash().get(split).searchstart > (int)split.start || split.clearedReads) ) {
					
					if(!sampleZoomer && !loading && !scrollbar) {
						
						getReads(split);				
						
							
					}
				}*/
			}
			if(i >= Main.samples) {
				break;				
			}			
			if(readsample.getreadHash().get(split) == null) {
				readsample.resetreadHash();
			}
			readsample.getreadHash().get(split).setReadSize(readsample.getreadHash().get(split).getReads().size());
		
			for(int j = 0; j<readsample.getreadHash().get(split).getReads().size();j++) {
					if((int)(((readsample.getIndex()+1)*this.drawVariables.sampleHeight-(readsample.getreadHash().get(split).readHeight+2)*j))+readsample.getreadHash().get(split).readwheel < (readsample.getIndex()*this.drawVariables.sampleHeight + readsample.getreadHash().get(split).readHeight)+this.drawVariables.sampleHeight/split.getDivider()) {
		
						break;
					}	
					else if((int)(((readsample.getIndex()+1)*this.drawVariables.sampleHeight-(readsample.getreadHash().get(split).readHeight+2)*j))+readsample.getreadHash().get(split).readwheel > (readsample.getIndex()+1)*this.drawVariables.sampleHeight) {
					
						continue;
					}
			
					
				try {
					
					read = readsample.getreadHash().get(split).getReads().get(j);
					
					while(read.getPrev() != null && read.getPosition() > split.start-split.viewLength) {							
						read = read.getPrev();
					}
					while(read.getNext() != null && read.getNext().getPosition() < split.start-split.viewLength) {
						
							read = read.getNext();						
					}
					
					if(read.getPosition() > split.end) {
				//		readsample.getreadHash().get(split).setReadSize(j);
					//	if(readsample.getreadHash().get(split).getReadSize()*readsample.getreadHash().get(split).readHeight < readsample.getreadHash().get(split).readwheel) {
							//System.out.println("jou");
						//	readsample.getreadHash().get(split).readwheel = 0;
					//	}
					//	break;
					}
					try {
					readsample.getreadHash().get(split).getReads().set(j, read);
					}
					catch(Exception e) {
						e.printStackTrace();
						break;
					}
				while(read != null) {							
					
					if(read.getPosition() > split.end) {						
						break;
					}
					if(read.getPosition()+read.getWidth() < split.start && !read.equals(clickedRead)) {
						read = read.getNext();
						continue;
					}
			
					readY = (int)(((readsample.getIndex()+1)*drawVariables.sampleHeight-this.bottom-(j+1)*(readsample.getreadHash().get(split).readHeight+2))+this.sampleList.get(i).getreadHash().get(split).readwheel-Main.drawScroll.getVerticalScrollBar().getValue());
									
					read.setRectBounds((int)((read.getPosition()-split.start)*split.pixel), readY, (int)(read.getWidth()*split.pixel), readsample.getreadHash().get(split).readHeight);
										
						
					drawRead(read,readsample.getreadHash().get(split), readY);								
					read = read.getNext();
				}		
				
				
				}
				catch(Exception e) {
					ErrorLog.addError(e.getStackTrace());
					e.printStackTrace();
				}
				
			}
			if(clickedRead != null && clickedReadSample.equals(sampleList.get(i))) {
				
				if(clickedRead.getMates() != null) {
					
					mateadd = 0;
					for(int m = 0; m<clickedRead.getMates().size();m++) {
						read = clickedRead.getMates().get(m);	
					
						if(read.getRect().y >((int)(readsample.getIndex()*this.drawVariables.sampleHeight-drawVariables.visiblestart*this.drawVariables.sampleHeight)+this.drawVariables.sampleHeight/read.split.getDivider()) && read.getRect().y < ((int)(sampleList.get(i).getIndex()*this.drawVariables.sampleHeight-drawVariables.visiblestart*this.drawVariables.sampleHeight)+this.drawVariables.sampleHeight)) {
							if(read.getPosition()+read.getWidth() > read.split.start && read.getPosition() < read.split.end) {
								read.yposition = read.getRect().y;
								read.setRectBounds((int)((read.getPosition()-split.start)*split.pixel), read.getRect().y, (int)(read.getWidth()*split.pixel), readsample.getreadHash().get(split).readHeight);
								
								drawRead(read,clickedReadSample.getreadHash().get(read.split),read.getRect().y);	
							}
							
							continue;					
						}
						
						if(read.split == null || clickedReadSample.getreadHash().get(read.split) == null) {
							clickedRead.getMates().remove(m);
							m--;
							continue;
						}
						if(!clickedRead.equals(clickedRead.getMates().get(m))) {
							if(!clickedRead.split.equals(read.split)) {
								read.split.updateReads = true;
							}
							if(m > 0) {
								if(clickedRead.getMates().get(m-1).split.equals(read.split)) {									
									mateadd++;									
								}
								else {
									mateadd = 0;
								}
							}
							else {
								if(read.getRect().intersects(clickedRead.getRect())) {								
									mateadd++;								
								}
								else {
									mateadd = 0;
								}
							}
							
						}
						else {
							mateadd = 0;
						}
						
						yposition = clickedRead.getRect().y-(clickedReadSample.getreadHash().get(read.split).readHeight+20)*mateadd;
						read.yposition = yposition;
						
						if(read.getPosition()+read.getWidth() > read.split.start && read.getPosition() < read.split.end) {
							drawRead(read,clickedReadSample.getreadHash().get(read.split), yposition);					
						}
					}
					
				}
				else {
					if(clickedRead.getRect().y >((int)(sampleList.get(i).getIndex()*this.drawVariables.sampleHeight-drawVariables.visiblestart*this.drawVariables.sampleHeight)+this.drawVariables.sampleHeight/clickedRead.split.getDivider()) && clickedRead.getRect().y < ((int)(sampleList.get(i).getIndex()*this.drawVariables.sampleHeight-drawVariables.visiblestart*this.drawVariables.sampleHeight)+this.drawVariables.sampleHeight)) {
						
						clickedRead.yposition = clickedRead.getRect().y;
						drawRead(clickedRead,clickedReadSample.getreadHash().get(clickedRead.split),clickedRead.getRect().y);		
						continue;					
					}
				}
			}
			
			
			if(readsample.getreadHash().get(split).getReadSize()*(sampleList.get(i).getreadHash().get(split).readHeight+2) > this.drawVariables.sampleHeight-(this.drawVariables.sampleHeight/split.getDivider())) {
				readsidebar = true;
				sampleheight = this.drawVariables.sampleHeight-(this.drawVariables.sampleHeight/split.getDivider());
				totalheight = readsample.getreadHash().get(split).getReadSize()*(readsample.getreadHash().get(split).readHeight+2);
				
				split.getReadBuffer().setColor(Color.lightGray);
				yValue = (int)((readsample.getIndex()*this.drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue()+this.drawVariables.sampleHeight/split.getDivider());
				
				int bottom = (int)((readsample.getIndex()*this.drawVariables.sampleHeight)-Main.drawScroll.getVerticalScrollBar().getValue()+this.drawVariables.sampleHeight);
				barHeight = sampleheight*(sampleheight/totalheight);
				if(barHeight < 4) {
					barHeight = 4;
				}
				readsample.getreadHash().get(split).getScrollBar().setBounds(drawWidth-(Main.defaultFontSize+8), yValue, Main.defaultFontSize+8, (int)(this.drawVariables.sampleHeight-(this.drawVariables.sampleHeight/split.getDivider())));
				readsample.getreadHash().get(split).getScroller().setBounds(drawWidth-(Main.defaultFontSize+6), (int)(bottom-(readsample.getreadHash().get(split).readwheel*(sampleheight/totalheight))-barHeight), Main.defaultFontSize+6,(int)barHeight );
				split.getReadBuffer().fillRect(readsample.getreadHash().get(split).getScrollBar().x, yValue, 20, (int)this.drawVariables.sampleHeight);
				split.getReadBuffer().setColor(Color.gray);
				split.getReadBuffer().fillRect(readsample.getreadHash().get(split).getScroller().x,readsample.getreadHash().get(split).getScroller().y, 16, (int)barHeight);
				readsample.getreadHash().get(split).setReadScroll(true);
			}
			else {
				readsample.getreadHash().get(split).setReadScroll(false);
				readsidebar = false;
			}
			}
			}
			catch(Exception e) {
				ErrorLog.addError(e.getStackTrace());
				e.printStackTrace();
			}
		}
		
			split.getReadBuffer().setComposite( composite);					
			split.getReadBuffer().fillRect(this.drawWidth,0, this.drawWidth/2,Main.drawScroll.getViewport().getHeight());	
			split.getReadBuffer().setComposite(split.getBackupr());			
		
			split.updateReads = false;
			updateReads = false;
			updateCoverages = false;
			// drawSelect(split);
		}
		
		if(split.viewLength <= Settings.readDrawDistance && selectedSplit.equals(split)) {
			drawSelect(split);
		}		
		}		
	}	

	}
	catch(Exception ex) {
		ErrorLog.addError(ex.getStackTrace());
		ex.printStackTrace();
	}
	
}

void drawRead(ReadNode read, Reads reads, int ypos) {
	if(reads == null) {
		return;
	}
	if(read.isDiscordant() && read.SA == null) {						
		tempColor = getColor(read, read.split.chrom);					
	}
	else {
		if(read.getMappingQuality() < Settings.mappingQ) {
			tempColor = dark;
		}
		else if(read.isForward()) {
			tempColor = forwardColor;
		}
		else {
			tempColor = reverseColor;
		}									
	}
	
	read.split.getReadBuffer().setColor(tempColor);	
	/*if(read.getName().equals("ST-E00214:53:H3JLWCCXX:1:2122:7354:64440")) {
		read.split.getReadBuffer().setColor(Color.magenta);	
	}*/
	drawRect = read.getRect();
	
		if(read.getCigar() == null) {
			if(read.split.pixel < 0.1) {
				if(read.isForward()) {
					read.setRectBounds((int)((read.getPosition()-read.split.start)*read.split.pixel), (int)drawRect.y, (int)(read.getWidth()*read.split.pixel)+5, reads.readHeight);
				}
				else {
					read.setRectBounds((int)((read.getPosition()-read.split.start)*read.split.pixel)-5,  (int)drawRect.y, (int)(read.getWidth()*read.split.pixel)+5, reads.readHeight);
				}							
			}
			
			if(read.isForward()) {
				read.split.getReadBuffer().fillRect((int)drawRect.getX(), ypos , (int)drawRect.getWidth()-5, (int)drawRect.getHeight());
			}
			else {
				read.split.getReadBuffer().fillRect((int)drawRect.getX()+5, ypos , (int)drawRect.getWidth()-5, (int)drawRect.getHeight());
			}
			if(read.getMatePos() == -1) {
				read.split.getReadBuffer().setColor(Color.white);
				if(read.isForward()) {
					read.split.getReadBuffer().drawRect((int)drawRect.getX()-1, ypos, (int)drawRect.getWidth()-5, (int)drawRect.getHeight());
				}
				else {
					read.split.getReadBuffer().drawRect((int)drawRect.getX()+5, ypos, (int)drawRect.getWidth()-5, (int)drawRect.getHeight());
				}
			}				
		}
	
		else {
		
		readpos = read.getPosition();
		insertion = 0;
		clipLength = 0;
		for(int c = 0; c<read.getCigar().numCigarElements(); c++) {
			
			element = read.getCigar().getCigarElement(c);
			if(element.getOperator().compareTo(CigarOperator.HARD_CLIP) == 0) {
				continue;
			}
			if(element.getOperator().compareTo(CigarOperator.MATCH_OR_MISMATCH)== 0 || element.getOperator().name().equals("=")) {		
				if((c==0 || (c==1 && read.getCigar().getCigarElement(0).getOperator().compareTo(CigarOperator.HARD_CLIP) == 0)) && !read.isForward()) {
					read.split.getReadBuffer().fillRect((int)((readpos+insertion-read.split.start)*read.split.pixel)+10, ypos, (int)((element.getLength()-insertion)*read.split.pixel )-10, (int)drawRect.getHeight());
					
				}
				else if(c == read.getCigar().numCigarElements()-1 && read.isForward()) {
					read.split.getReadBuffer().fillRect((int)((readpos+insertion-read.split.start)*read.split.pixel), ypos, (int)((element.getLength()-insertion)*read.split.pixel )-10, (int)drawRect.getHeight());

				}	
				else if((c == read.getCigar().numCigarElements()-2 && read.getCigar().getCigarElement(c+1).getOperator().compareTo(CigarOperator.HARD_CLIP) == 0) && read.isForward()) {
					read.split.getReadBuffer().fillRect((int)((readpos+insertion-read.split.start)*read.split.pixel),ypos, (int)((element.getLength()-insertion)*read.split.pixel )-10, (int)drawRect.getHeight());
					
				}
				else {
					
					read.split.getReadBuffer().fillRect((int)((readpos+insertion-read.split.start)*read.split.pixel), ypos, (int)((element.getLength()-insertion)*read.split.pixel ), (int)drawRect.getHeight());
											
				}
				if(insertion > 0) {
					insertion = 0;
				}
				readpos+=element.getLength();
			}
			else if(element.getOperator().compareTo(CigarOperator.DELETION)== 0) {	
				read.split.getReadBuffer().setColor(Color.black);
				if(element.getLength()*read.split.pixel > 5) {
					read.split.getReadBuffer().drawLine((int)((readpos-read.split.start)*read.split.pixel),  ypos+reads.readHeight/2, (int)((readpos+element.getLength()-read.split.start)*read.split.pixel),  (int)drawRect.getY()+reads.readHeight/2);
					
				}
				else {
					
					read.split.getReadBuffer().fillRect((int)((readpos-read.split.start)*read.split.pixel), ypos, (int)(element.getLength()*read.split.pixel+1), reads.readHeight);
					
				}
				readpos+=element.getLength();
				read.split.getReadBuffer().setColor(tempColor);
			}
			else if(element.getOperator().compareTo(CigarOperator.INSERTION)== 0) {	
				if(read.split.pixel > 5) {
					read.split.getReadBuffer().setColor(Color.white);
					read.split.getReadBuffer().fillRect((int)((readpos-read.split.start)*read.split.pixel), ypos, (int)(read.split.pixel+1), reads.readHeight);
					read.split.getReadBuffer().setColor(Color.black);
					read.split.getReadBuffer().drawString("I"+element.getLength(),(int)((readpos-read.split.start)*read.split.pixel)+2, ypos+reads.readHeight);
					read.split.getReadBuffer().setColor(tempColor);
				}
				else {
					read.split.getReadBuffer().setColor(Color.white);
					read.split.getReadBuffer().fillRect((int)((readpos-read.split.start)*read.split.pixel), ypos, (int)(element.getLength()*read.split.pixel+1), reads.readHeight);
					read.split.getReadBuffer().setColor(tempColor);
				}
				
				insertion = 1;
				//readpos++;
			}
			else if(element.getOperator().compareTo(CigarOperator.SOFT_CLIP)== 0) {
				
				if(read.SA == null /*|| Settings.softClips*/) {
					if(read.isForward() && c == read.getCigar().numCigarElements()-1) {
						read.split.getReadBuffer().fillRect((int)((readpos+insertion-read.split.start)*read.split.pixel)-1, ypos, (int)((element.getLength()-insertion)*read.split.pixel )-8, (int)drawRect.getHeight());
						read.split.getReadBuffer().setColor(softColor);							
						read.split.getReadBuffer().fillRect((int)((readpos+insertion-read.split.start)*read.split.pixel)-1, ypos, (int)((element.getLength()-insertion)*read.split.pixel )-8, (int)drawRect.getHeight());
					}
					else if(!read.isForward() && c == 0) {
						read.split.getReadBuffer().fillRect((int)((readpos+insertion-read.split.start)*read.split.pixel)+10, ypos, (int)((element.getLength()-insertion)*read.split.pixel )-8, (int)drawRect.getHeight());
						read.split.getReadBuffer().setColor(softColor);							
						read.split.getReadBuffer().fillRect((int)((readpos+insertion-read.split.start)*read.split.pixel)+10, ypos, (int)((element.getLength()-insertion)*read.split.pixel )-8, (int)drawRect.getHeight());
					}
					else {
						read.split.getReadBuffer().fillRect((int)((readpos+insertion-read.split.start)*read.split.pixel)-1, ypos, (int)((element.getLength()-insertion)*read.split.pixel ), (int)drawRect.getHeight());
						read.split.getReadBuffer().setColor(softColor);							
						read.split.getReadBuffer().fillRect((int)((readpos+insertion-read.split.start)*read.split.pixel)-1, ypos, (int)((element.getLength()-insertion)*read.split.pixel ), (int)drawRect.getHeight());
				
					}			
					
					read.split.getReadBuffer().setColor(tempColor);
					readpos+=element.getLength();
					
				}
				else {
					read.split.getReadBuffer().setColor(Color.white);
					if(c== 0) {
						clipLength = element.getLength();
						read.split.getReadBuffer().fillRect((int)((readpos+insertion-read.split.start)*read.split.pixel)-1, ypos-2,2,reads.readHeight+4); 
					}
					else {
						read.split.getReadBuffer().fillRect((int)((readpos+insertion-read.split.start)*read.split.pixel), ypos-2,2,reads.readHeight+4); 
						
					}
					read.split.getReadBuffer().setColor(tempColor);
					
				}
				
				if(insertion > 0) {
					insertion = 0;
				}								
			}
			else if(element.getOperator().compareTo(CigarOperator.SKIPPED_REGION)== 0) {
				read.split.getReadBuffer().setColor(Color.black);								
				read.split.getReadBuffer().drawLine((int)((readpos-read.split.start)*read.split.pixel),  ypos+reads.readHeight/2, (int)((readpos+element.getLength()-read.split.start)*read.split.pixel),  (int)drawRect.getY()+reads.readHeight/2);
				read.split.getReadBuffer().setColor(tempColor);
				readpos+=element.getLength();							
			}
		}	
		if(read.isForward()) {
			read.setRectBounds((int)((read.getPosition()-read.split.start)*read.split.pixel), read.getRect().y, (int)((readpos-read.getPosition())*read.split.pixel)-5, reads.readHeight);
		}
		else {
			
			read.setRectBounds((int)((read.getPosition()-read.split.start)*read.split.pixel)+5, read.getRect().y, (int)((readpos-read.getPosition())*read.split.pixel)-5, reads.readHeight);
		}
	}
		
	if(read.split.pixel > 0.1) {
		drawTriangles(read, reads.sample,read.split,read.split.getReadBuffer().getColor());
	}
	if(read.split.viewLength < 400 && read.getMismatches() != null) {
		if(read.split.getReadBuffer().getFont().getSize() != reads.readfont.getSize()) {
			read.split.getReadBuffer().setFont(reads.readfont);
		}
		
		mismatches = read.getMismatches();
		
		for(Entry<Integer, Byte> entry : mismatches) {
			
			if(!read.isForward()) {
				if(Character.isLowerCase((char)(byte)entry.getValue())) {								
					read.split.getReadBuffer().setColor(reverseTextLow);	
				}
				else {								
					read.split.getReadBuffer().setColor(reverseText);	
				}
			}
			else {
				
				if(Character.isLowerCase((char)(byte)entry.getValue())) {							
					read.split.getReadBuffer().setColor(forwardTextLow);	
				}
				else {
					read.split.getReadBuffer().setColor(forwardText);
				}
			}
			
			if(!Main.getBase.containsKey(entry.getValue())) {
				read.split.getReadBuffer().setColor(Color.black);	
				read.split.getReadBuffer().drawString(""+(char)(byte)entry.getValue(), (int)((((read.getPosition()+entry.getKey())-read.split.start)*read.split.pixel)), ypos+reads.readHeight);
				
			}
			else {
			
				read.split.getReadBuffer().drawString(Main.getBase.get(entry.getValue()), (int)((((read.getPosition()-clipLength+entry.getKey())-read.split.start)*read.split.pixel)+read.split.pixel/3),ypos+reads.readHeight);
					
										
			}	
			
		}
//		read.split.getReadBuffer().setColor(Color.white);
//		read.split.getReadBuffer().setFont(defaultFont);
	
	}
	else if (read.split.viewLength < Settings.readDrawDistance) {
	//	read.split.getReadBuffer().setFont(ChromDraw.seqFont);
		
		if(read.getMismatches() != null) {
			mismatches = read.getMismatches();
			for(int m= 0; m< mismatches.size(); m++) {
				mismatchentry = mismatches.get(m);
				if(!read.isForward()) {
					if(Character.isLowerCase((char)(byte)mismatchentry.getValue())) {		
						
						read.split.getReadBuffer().setColor(reverseTextLow);	
					}
					else {								
						read.split.getReadBuffer().setColor(reverseText);	
					}
				}
				else {
					
					if(Character.isLowerCase((char)(byte)mismatchentry.getValue())) {							
						read.split.getReadBuffer().setColor(forwardTextLow);	
					}
					else {
						read.split.getReadBuffer().setColor(forwardText);
					}
				}
				
				if((char)(byte)mismatchentry.getValue() == ('I')) {
					read.split.getReadBuffer().setColor(Color.black);	
				}					
										
				read.split.getReadBuffer().fillRect((int)(((read.getPosition()+mismatchentry.getKey())-read.split.start)*read.split.pixel), ypos, (int)(read.split.pixel+1), reads.readHeight);
			
			}
		}
	}
	if(reads.sample.hasMates &&  read.getMatePos() == -1) {
		read.split.getReadBuffer().setColor(Color.white);	
		read.split.getReadBuffer().drawRect(read.getRect().x, read.getRect().y, read.getRect().width, read.getRect().height);
	}
}

void drawSelect(SplitClass split) {
	try {
		
		split.getSelectbuf().setComposite( composite);		
		split.getSelectbuf().fillRect(0,0, split.getSelectbuffer().getWidth(),Main.drawScroll.getViewport().getHeight());	
		split.getSelectbuf().setComposite(split.getBackups());			
		
		if((readsidebar && (moveX > this.getDrawWidth()-(Main.defaultFontSize+8)+split.offset ) )|| moveY < selectedIndex*drawVariables.sampleHeight + drawVariables.sampleHeight/split.getDivider()) {		
			coverageregion = true;
			return;
		}
		coverageregion = false;
		if(mouseWheel || sidebar) {		
			return;			
		}
	
		if(sampleList.get(selectedIndex).getreadHash() != null && sampleList.get(selectedIndex).getreadHash().size() > 0 && readLevel < sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getReads().size() && readLevel > -1) {
				row = sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getReads().get(readLevel);
				foundread = false;			
			
				while(row != null) {
					try {
							if(row.getPosition() > (int)(split.start+mouseRect.x/split.pixel)) {
								if(selectedMate != null) {
									selectedMate.split.getSelectbuf().setComposite( composite);		
									selectedMate.split.getSelectbuf().fillRect(0,0, (int)Main.screenSize.getWidth(),  (int)Main.screenSize.getHeight());	
									selectedMate.split.getSelectbuf().setComposite(selectedMate.split.getBackups());
								}
								break;
							}
							
							if(row.getPosition()+row.getWidth() < split.start) {
								row = row.getNext();
								continue;
							}
							if(row.getRect().x+row.getRect().width < mouseRect.x) {
								row = row.getNext();
								continue;
							}											
							
							foundread = true;						
							split.getSelectbuf().setColor(Color.white);
							
							split.getSelectbuf().drawRect(row.getRect().x, row.getRect().y, row.getRect().width, row.getRect().height);	
							
							if(selectedRead == null) {
								selectedRead = row;
							}
							else if(!selectedRead.equals(row)) {
								selectedRead = row;								
							}
							drawTriangles(row, sampleList.get(selectedIndex), split,split.getSelectbuf().getColor());
							
							if(selectedRead.getMates() != null) {
															
								for(int i = 0; i<selectedRead.getMates().size(); i++) {									
									selectedMate = selectedRead.getMates().get(i);
									if(selectedMate.equals(selectedRead)) {
										continue;
									}
									if(splits.indexOf(selectedMate.split) > -1 && selectedMate.getPosition()+selectedMate.getWidth() > selectedMate.split.start && selectedMate.getPosition() < selectedMate.split.end/* && !selectedRead.split.equals(selectedMate.split)*/) {
										if(splits.indexOf(selectedMate.split) > -1 && !selectedRead.split.equals(selectedMate.split)) {
											selectedMate.split.getSelectbuf().setComposite( composite);		
											selectedMate.split.getSelectbuf().fillRect(0,0, (int)Main.screenSize.getWidth(),  (int)Main.screenSize.getHeight());	
											selectedMate.split.getSelectbuf().setComposite(selectedMate.split.getBackups());
										}
										selectedMate.split.getSelectbuf().drawRect(selectedMate.getRect().x, selectedMate.getRect().y, selectedMate.getRect().width, selectedMate.getRect().height);						
										
									}										
								}//TODO nanopore reads								
							}
							else if(selectedRead.getMates() != null && !selectedRead.isDiscordant()) {
								
						
								if(splits.size() > 1 && selectedMate != null) {
									
									selectedMate.split.getSelectbuf().setComposite( composite);		
									selectedMate.split.getSelectbuf().fillRect(0,0, (int)Main.screenSize.getWidth(),  (int)Main.screenSize.getHeight());	
									selectedMate.split.getSelectbuf().setComposite(selectedMate.split.getBackups());
								}								
							}
							else if (selectedRead.getMates() != null && selectedRead.isDiscordant() ) {
								
								 for(int i = 0; i<selectedRead.getMates().size(); i++) {
									
									 	if(selectedRead.getMates().get(i).equals(selectedRead)) {
									 		continue;
									 	}
									 	 
									 	selectedMate = selectedRead.getMates().get(i);
									 	if(splits.indexOf(selectedMate.split) <0) {
									 		 selectedRead.getMates().remove(i);
											 continue;
									 	}										
										selectedMate.split.getSelectbuf().drawRect(selectedMate.getRect().x, selectedMate.getRect().y, selectedMate.getRect().width, selectedMate.getRect().height);						
								 }								
							  }													
							break;
				//		}						
					//	row = row.getNext();					
					}
					catch(Exception ex) {
						ErrorLog.addError(ex.getStackTrace());
						ex.printStackTrace();
						break;
					}
			}
				
			if(!foundread) {				
				selectedRead = null;
					
				if(splits.size() > 1) {
					for (int i = 0 ; i<splits.size(); i++) {
						splits.get(i).getSelectbuf().setComposite( composite);		
						splits.get(i).getSelectbuf().fillRect(0,0, (int)Main.screenSize.getWidth(),  (int)Main.screenSize.getHeight());	
						splits.get(i).getSelectbuf().setComposite(splits.get(i).getBackups());
					}
				}					
			}			
		}
		else {
			if(splits.size() > 1) {
				for (int i = 0 ; i<splits.size(); i++) {
					splits.get(i).getSelectbuf().setComposite( composite);		
					splits.get(i).getSelectbuf().fillRect(0,0, (int)Main.screenSize.getWidth(),  (int)Main.screenSize.getHeight());	
					splits.get(i).getSelectbuf().setComposite(splits.get(i).getBackups());
				}
			}		
		}
		}
		catch(Exception e ) {
			ErrorLog.addError(e.getStackTrace());
			e.printStackTrace();
		}
}

public void loading(String text) {	
	
	
	
	 
	 if (!loading) {
		
		loading = true;
    	Main.loading = new Loader(text);
    	Main.loading.execute();	
	}
	 
	readyQueue.add(text);
	loadingtext = text;
	setGlasspane(true);
      
}
public void ready(String text) {	
	if(text.equals("all")) {
		
		readyQueue.clear();
		loading = false;
	}
	if(readyQueue.size() > 0) {
		Main.drawCanvas.readyQueue.remove(text);		
		
		if(readyQueue.size() > 0) {
			loadingtext = readyQueue.get(readyQueue.size()-1);
		}
		
	}
	
	if(readyQueue.size() == 0) {
		
		Main.frame.setCursor( Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
		
		loading = false;
		loadingtext="";
		//Main.loading.text = "";	
		/*if(!FileRead.novars && FileRead.changing && VariantHandler.commonSlider.getValue() > 1 && VariantHandler.clusterSize > 0) {
			
			calcClusters(FileRead.head, 1);
		}*/
		
		setGlasspane(false);
		
	}
}
public void paint(Graphics g) {
	try {				
		drawScreen(g);
	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
	}
}
public void resizeSidebar(int move) {
	Main.sidebarWidth = move;
	if(Main.sidebarWidth < 3) {
		Main.sidebarWidth = 3;
	}
//	resizeCanvas(this.getHeight());
	repaint();
	for(int i = 0 ; i<splits.size(); i++) {
		 splits.get(i).offset = Main.sidebarWidth +(int)(getDrawWidth()*i);
		 splits.get(i).chromOffset = (int)(getDrawWidth()*i);
    	 splits.get(i).pixel = (this.getDrawWidth())/(double)(this.splits.get(i).end-this.splits.get(i).start);	    	
    //	 Main.chromDraw.updateExons = true;	    	
    //	 Main.chromDraw.drawExons(splits.get(i));	    	
    //	 Main.chromDraw.repaint();
	}
//	Main.chromDraw.repaint();

}
@Override
public void mouseDragged(MouseEvent event) {
	moveY = event.getY();
	moveX = event.getX();
	
	switch(event.getModifiers()) {	
	
		case InputEvent.BUTTON1_MASK: {	
			
			
			mouseX = event.getX();
			mouseY = event.getY();
			
			if(readScrollDrag) {
				
				if(moveY-Main.drawScroll.getVerticalScrollBar().getValue() <= this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getScrollBar().y+this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getScrollBar().height && moveY-Main.drawScroll.getVerticalScrollBar().getValue() >= this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getScrollBar().y) {
					double sampleheight = this.drawVariables.sampleHeight-(this.drawVariables.sampleHeight/selectedSplit.getDivider());
					double totalheight = this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getReadSize()*(this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2);
					this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel -=(moveY-moveYtemp)/(sampleheight/totalheight);
					if(this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel < 0) {
						this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel = 0;
						
					}
					
					if(( (Main.drawCanvas.drawVariables.sampleHeight-(Main.drawCanvas.drawVariables.sampleHeight/selectedSplit.getDivider()))) > (Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getReadSize()*(this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2)) - (Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel)) {
						Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel = (int)((Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getReadSize()*(this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2)) - ((Main.drawCanvas.drawVariables.sampleHeight-(Main.drawCanvas.drawVariables.sampleHeight/selectedSplit.getDivider()))));
						
					}					
					
					moveYtemp = moveY;
					selectedSplit.updateReads = true;
					repaint();					
				}
				
				break;
			}
			if(resizeSidebar) {				
				resizeSidebar(mouseX);
				break;
			}
			if(!sidebar) {
				
				if (!sampleZoomer && (mouseX-pressX < -5 && mouseY-pressY < -5 || mouseX-pressX > 5 && mouseY-pressY > -5) || lineZoomer) {
					
					mousePress = false;
					task.cancel();
					zoomDrag = true;
					if(lineZoomer) {	
						if(splits.get(0).viewLength < Settings.readDrawDistance) {
							updatevars = true;
						}
						Main.chromDraw.updateExons = true;
						 if (selectedSplit.start > 0 || selectedSplit.end < selectedSplit.chromEnd ) {					
							 zoom();	
							 Main.chromDraw.repaint();						
						 }
					}
					
				}
				else if(!zoomDrag && !lineZoomer && (mouseX-pressX > 5 && mouseY-pressY < -10 || mouseX-pressX < -5 && mouseY-pressY > 5)) {
					mousePress = false;
					task.cancel();
					
					if(drawVariables.visiblesamples == 1 && mouseY-pressY < -10) {
						break;
					}
					
					if(drawVariables.sampleHeight < 1  && mouseY-pressY > 5) { //defaultSampleHeight /*|| (drawVariables.visiblesamples == Main.samples && drawVariables.sampleHeight > defaultSampleHeight ) )*/ && mouseY-pressY > 5) {
						
						break;
						
					}
				
					sampleZoomer = true;		
				}
				
				repaint();
				
			}
			else {
				if(!bamHover && selectedIndex > -1) {
					if(getCursor().getType() != Cursor.HAND_CURSOR) {					
						setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
					}
					sampleDrag = true;
					if(hoverIndex != (int)((moveY)/drawVariables.sampleHeight)) {
						hoverIndex = (int)((moveY)/drawVariables.sampleHeight);		
						
					}
					repaint();		
				}
			}
			break;
		}
		case InputEvent.BUTTON3_MASK: {	
			if((int)selectedSplit.start == 0 && (int)selectedSplit.end ==selectedSplit.chromEnd) {
				break;
			}
			mouseDrag = true;
	    	drag(event.getX());
			
			break;
		}
		case 17: {
			if((int)selectedSplit.start == 0 && (int)selectedSplit.end ==selectedSplit.chromEnd) {
				break;
			}
			mouseDrag = true;
	    	drag(event.getX());
			
			break;
		}
	}
	
	
	
}
void drag(int dragX) {	
		
		starttemp = selectedSplit.start;
		endtemp = selectedSplit.end;
		
		if(starttemp-(dragX-tempDrag)/selectedSplit.pixel < 0 ) {
			starttemp = 0;
			endtemp = starttemp + selectedSplit.viewLength;
		}
		else if(endtemp-(dragX-tempDrag)/selectedSplit.pixel > selectedSplit.chromEnd) {
			endtemp = selectedSplit.chromEnd;
			starttemp = endtemp-selectedSplit.viewLength;
		}
		else {		
			
			starttemp-= (dragX-tempDrag)/selectedSplit.pixel;			
			endtemp  -= (dragX-tempDrag)/selectedSplit.pixel;				
		}

		setStartEnd(starttemp, endtemp);
		updatevars = true;
//		updateReads = true;
//		updateCoverages = true;
		
		Main.bedCanvas.repaint();
			
		Main.chromDraw.updateExons = true;
		Main.chromDraw.repaint();
		tempDrag = dragX;
		repaint();
	
	
}

void addSplit(String chrom, int start, int end) {
    SplitClass split = new SplitClass();
    split.chrom = chrom.replace("chr","");
    if(selectedSplit.chrom.equals(split.chrom)) {
    	split.setGenes(selectedSplit.getGenes());
    }
    else {
    	
    	split.setGenes(Main.fileReader.getExons(split.chrom));
    }
   
    if(Main.samples > 0) {
    	for(int i = 0; i<this.sampleList.size(); i++) {
    		if(sampleList.get(i).removed) {
				continue;
			}
    		if(sampleList.get(i).getreadHash() == null) {
    			continue;
    		}
    		Reads newReads = new Reads();
    		newReads.sample = this.sampleList.get(i);
    		this.sampleList.get(i).getreadHash().put(split, newReads);
    	}
    }
    try {
    	if(Main.chromIndex.containsKey(Main.refchrom +split.chrom)) {
    		split.chromEnd = Main.chromIndex.get(Main.refchrom +split.chrom)[1].intValue();
    	}
    	else {
    		Main.showError("Target chromosome (" +(Main.refchrom +split.chrom) +") is not available.\nTry to download more comprehensive reference FASTA-file.", "Note");    			
    		return;
    	}
    }
    catch(Exception e) {
    	split.chromEnd = end+(end-start);
    	System.out.println(split.chrom);
    	ErrorLog.addError(e.getStackTrace());
    	e.printStackTrace();
    }
    split.start = start;
    split.end = end;
    split.viewLength =split.end-split.start;
   
    this.splits.add(split);	
//    split.updateReads = true;
    
    drawVariables.visiblestart = (short)selectedIndex;
    
    drawVariables.visiblesamples = (short)1;
 //   this.drawReads(split);
    Main.drawCanvas.resizeCanvas(Main.drawCanvas.getWidth(), Main.drawCanvas.getHeight());
    setScrollbar((int)(selectedIndex*drawVariables.sampleHeight));	
    
   
//    for(int i= 0; i<drawCanvas.splits.size(); i++) {
		//drawCanvas.splits.get(i).cytoImage = null;
		//chromDraw.drawCyto(drawCanvas.splits.get(i));
     
	//	Main.chromDraw.updateExons = true;
	//	Main.chromDraw.repaint();
//	}
	for(int i = 0; i<splits.size(); i++) {
		Main.chromDraw.drawCyto(Main.drawCanvas.splits.get(i));	
		Main.chromDraw.drawExons(Main.drawCanvas.splits.get(i));	
		Main.chromDraw.updateExons = true;
		 Main.chromDraw.repaint();
    	 splits.get(i).updateReads = true;
    	 
    	 
     }
	repaint();
}	

void removeSplits() {
	try {
			
	//	clickedReadSplit = null;
		selectedSplit = splits.get(0);
		for (int i = splits.size()-1; i>0; i--) {
			splits.get(i).setCytoImage(null);
			splits.get(i).setGenes(null);
			splits.get(i).removeSplitDraw();
			try {
				if(Main.samples > 0) {
			    	for(int s = 0; s<sampleList.size(); s++) {
			    		if(sampleList.get(s).removed) {
							continue;
						}
			    		sampleList.get(s).getreadHash().remove(splits.get(i));				    		
			    	}
			    }
			}
			catch(Exception e) {
				e.printStackTrace();
			}
			splits.remove(i);
		/*	if(Main.samples > 0) {
		    	for(int s = 0; s<sampleList.size(); s++) {
		    		sampleList.get(s).getreadHash().remove(i);				    		
		    	}
		    }*/
			
		}
			
		splits.get(0).setCytoImage(null);
		Main.chromDraw.drawCyto(splits.get(0));
		Main.chromDraw.updateExons = true;
		Main.chromDraw.repaint();
	
		Main.drawCanvas.resizeCanvas(Main.drawCanvas.getWidth(), Main.drawCanvas.getHeight());
	
//	Main.drawCanvas.repaint();
	
	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
	}
}

void removeSplit(SplitClass remsplit) {
	try {
		if(clickedRead != null && clickedRead.split.equals(remsplit)) {
	//		clickedRead = null;
			
			splitList.clear();
		}
	//	clickmates = null;
		selectedSplit = splits.get(0);
	    String prechrom = splits.get(0).chrom;
	    remsplit.setCytoImage(null);
	    remsplit.removed = true;
		splits.remove(remsplit);
		if(Main.samples > 0) {
	    	for(int s = 0; s<sampleList.size(); s++) {
	    		if(sampleList.get(s).removed) {
					continue;
				}
	    		if(sampleList.get(s).getreadHash() == null) {
	    			continue;
	    		}
	    		sampleList.get(s).getreadHash().remove(remsplit);				    		
	    	}
	    }
		remsplit = null;
		if(!splits.get(0).chrom.equals(prechrom)) {
			FileRead.search = true;
			gotoPos(splits.get(0).chrom, splits.get(0).start, splits.get(0).end);
		}
		
	for(int i = 0; i< splits.size(); i++) {
		splits.get(i).setCytoImage(null);
		Main.chromDraw.drawCyto(splits.get(i));
	}
	
	
	Main.chromDraw.updateExons = true;
	Main.chromDraw.repaint();

	Main.drawCanvas.resizeCanvas(Main.drawCanvas.getWidth(), Main.drawCanvas.getHeight());
	
	for(int i = 0; i<splits.size(); i++) {
		Main.chromDraw.drawCyto(Main.drawCanvas.splits.get(i));	
		Main.chromDraw.drawExons(Main.drawCanvas.splits.get(i));	
		Main.chromDraw.updateExons = true;
		 Main.chromDraw.repaint();
    	 splits.get(i).updateReads = true;
    	 
    	// repaint();
     }
	repaint();

}
catch(Exception e) {
	ErrorLog.addError(e.getStackTrace());
	e.printStackTrace();
	remsplit = null;
}
	remsplit = null;
}
void resizeCanvas(int height) {	
	this.height = height;
	this.setPreferredSize(new Dimension(this.getWidth(),height));
	drawVariables.sampleHeight = this.height/(double)drawVariables.visiblesamples;
   // drawVariables.visibleend = (short)(drawVariables.visiblestart + this.height/drawVariables.sampleHeight);
	
	updatevars = true;	
	updateReads = true;
	updateCoverages = true;
	
	this.revalidate();
}
void resizeCanvas(int width, int height) {	
	try {
		if(splits.get(0) == null || splits.get(0).getReadImage() == null) {
			return;
		}
		if(splits.get(0).getReadImage().getWidth() < this.getWidth()) {			
			resizeAllCanvas((int)(Main.screenSize.getWidth()*2));
		}
		if(Main.samples == 0) {
			this.setPreferredSize(new Dimension(this.getWidth(), height));
			this.width = width;
		    this.height = height;
		    this.drawWidth = (int)((this.width - Main.sidebarWidth)/(double)splits.size());
		}
		else {
			
		//	if(drawVariables.visiblesamples > Main.samples-1) {
		//		drawVariables.visiblesamples =(short)( Main.samples-1);
		//	}
			
		//	if(drawVariables.visibleend-drawVariables.visiblestart == 0) {
		//		drawVariables.sampleHeight = Main.drawScroll.getViewport().getHeight();
		//	}
		//	else {
		//		drawVariables.sampleHeight = Main.drawScroll.getViewport().getHeight()/(((double)drawVariables.visibleend-drawVariables.visiblestart)+1);
				
				drawVariables.sampleHeight = Main.drawScroll.getViewport().getHeight()/(double)drawVariables.visiblesamples;
				
		//	}
			if(drawVariables.sampleHeight < 1)  {
				drawVariables.sampleHeight = 1;
	//			drawVariables.visiblesamples = (short)(Main.drawScroll.getViewport().getHeight()/(double)drawVariables.sampleHeight);
			//	drawVariables.visibleend = (short)(drawVariables.visiblestart + Main.drawScroll.getViewport().getHeight()/drawVariables.sampleHeight);
			}			
			
			this.setPreferredSize(new Dimension(this.getWidth(),(int)(drawVariables.sampleHeight*Main.samples)));
			
			this.width = width;
		    this.height = (int)(drawVariables.sampleHeight*Main.samples);
		    this.drawWidth = (int)((this.width - Main.sidebarWidth)/(double)splits.size());
		   
			updatevars = true;	
			updateReads = true;
			updateCoverages = true;
			this.revalidate();
			
		}		
		
			for (int i = 0; i<splits.size(); i++) {	
				
				 splits.get(i).offset = Main.sidebarWidth +(int)(getDrawWidth()*i);
				 splits.get(i).chromOffset = (int)(getDrawWidth()*i);
		    	 splits.get(i).pixel = (this.getDrawWidth())/(double)(this.splits.get(i).end-this.splits.get(i).start);	    	
		    	 Main.chromDraw.updateExons = true;	    	
		    	 if(i > 0) {
		    		 Main.chromDraw.drawExons(splits.get(i));	    
		    	 }
		    	 Main.chromDraw.repaint();
		    	/* if(drawVariables.sampleHeight < 100) {		    	
		    		
					 sampleZoomer = true;				 			
		 		}*/		    	
		    }
			
	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
	}
	
}

@Override
public void mouseMoved(MouseEvent event) {
	try {
		
		moveX = event.getX();
		moveY = event.getY();		
		if(sampleZoomer) {
			sampleZoomer = false;
		}
		if(Main.samples > 0) {
		
			if((moveX-Main.sidebarWidth)/(this.getDrawWidth()) > -1 && (moveX-Main.sidebarWidth)/(this.getDrawWidth()) < splits.size()) {
				if(selectedSplit != splits.get((moveX-Main.sidebarWidth)/(this.getDrawWidth()))) {
					selectedSplit = splits.get((moveX-Main.sidebarWidth)/(this.getDrawWidth()));
					if(selectedSplit == null) {
						selectedSplit = splits.get(0);
					}				
					Main.chromDraw.repaint();					
				}
			}
			if(moveX < Main.sidebarWidth-3) {
				sidebar = true;
				mouseRect.setBounds(moveX, moveY-Main.drawScroll.getVerticalScrollBar().getValue(), 1,1);
				selectedIndex = (int)((moveY)/drawVariables.sampleHeight);
				
			}
			else {
				sidebar = false;
				removeSample = null;
			}
			
			if(!sidebar) {
				if(selectedIndex != (int)((moveY)/drawVariables.sampleHeight) ) {
					selectedIndex = (int)((moveY)/drawVariables.sampleHeight);
					
					if(selectedIndex > Main.samples-1) {
						selectedIndex =  Main.samples-1;
					}
				}
				
				if (this.sampleList.get(selectedIndex).getreadHash() != null && splits.size() > 0 && this.sampleList.get(selectedIndex).getreadHash().size() > 0) {
						mouseRect.setBounds(moveX-selectedSplit.offset, moveY-Main.drawScroll.getVerticalScrollBar().getValue(), 1,1);
						
					
					if(readLevel != (int)(((drawVariables.sampleHeight)*(selectedIndex+1)-(moveY))+this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel)/(this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2)) {
						readLevel = ((int)(((drawVariables.sampleHeight)*(selectedIndex+1))-(moveY))+this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel)/(this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2);		
					
					}					
				}
				else {
					mouseRect.setBounds(moveX-selectedSplit.offset, moveY, 1,1);
					
				}
			}			
			repaint();
		}
		else {			
		
			if(splits.size() > 1) {
				mouseRect.setBounds(moveX-selectedSplit.offset, moveY, 1,1);
				if(splits.size() > (moveX-Main.sidebarWidth)/(this.getDrawWidth())) {
					selectedSplit = splits.get((moveX-Main.sidebarWidth)/(this.getDrawWidth()));
				}
			}
			repaint();
		}
		
	
	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
	}	
}

void createSplit(Sample sample, String chrom, int pos) {
	try {
		
	/*	
		Main.chromDraw.splitClass = new SplitClass(sample, chrom);
		
		Main.chromDraw.cytoImage = Main.chromDraw.createBands(this.chrom);
	//	SplitClass split = Main.chromDraw.splitClass;	
		if (clickedRead.SA != null) {
			Main.chromDraw.splitClass.splitRead = true;
		}
		Main.chromDraw.splitClass.cytoImage = Main.chromDraw.createBands(chrom);	
		
		Main.chromDraw.splitClass.transcripts = Main.fileReader.getExons(chrom);
		
		Main.chromDraw.splitClass.getreadHash() = new Reads();
		Main.chromDraw.splitClass.getreadHash().chr = sample.getreadHash().chr;
		Main.chromDraw.splitClass.getreadHash().samFile = sample.getreadHash().samFile;
		Main.chromDraw.splitClass.gotoPos(pos-viewLength/2,pos+viewLength/2 );
	//	if(sample.getreadHash().samFiles.length == 1) {
//			split.getreadHash().samFileReader = sample.getreadHash().samFileReader;
	//	}
	/*	else {
			for(int i = 0; i<28; i++) {				
				if(split.getreadHash().samFiles[i].getName().contains("_" +chrom +".")) {					
					split.getreadHash().samFileReader = new SAMFileReader(split.getreadHash().samFiles[i]);
					break;
				}
			}
		}*/
//		split.getreadHash().bamIterator = sample.getreadHash().bamIterator;
	/*	if(sample.getreadHash().complete) {							
			Main.chromDraw.splitClass.getreadHash().complete = true;
		}
		
		Main.chromDraw.splitClass.getreadHash().sample = sample;
		Main.chromDraw.split = true;
	*/
	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
	}
}

void showsampleMenu(final Sample sample) {
	
	
	SampleDialog dialog = new SampleDialog(sample);
	dialog.createAndShowGUI();
	/*
	final JPopupMenu sampleMenu = new JPopupMenu();
	sampleMenu.addPopupMenuListener(new PopupMenuListener() {
		 
	    
	 
	    
		@Override
		public void popupMenuCanceled(PopupMenuEvent arg0) {
			// TODO Auto-generated method stub
			
		}

		@Override
		public void popupMenuWillBecomeInvisible(PopupMenuEvent arg0) {
			
			// TODO Auto-generated method stub
			
		}

		@Override
		public void popupMenuWillBecomeVisible(PopupMenuEvent arg0) {
			// TODO Auto-generated method stub
			
		}
	    });
	final JTextField nameField = new JTextField();
	
	String[] sexes = {"-", "Female", "Male"};
	
	DefaultComboBoxModel<String> refModel = new DefaultComboBoxModel<String>(sexes);
	
	SteppedComboBox sex = new SteppedComboBox(refModel);
	nameField.addKeyListener(new KeyListener() {
		@Override
		public void keyPressed(KeyEvent arg0) {}
		@Override
		public void keyReleased(KeyEvent arg0) {}
		@Override
		public void keyTyped(KeyEvent arg0) {
			sample.setName(nameField.getText() +arg0.getKeyChar());
			repaint();
		}		
	});
	nameField.setBorder(BorderFactory.createLineBorder(Color.darkGray));
	//JScrollPane menuscroll = new JScrollPane();
	nameField.setFont(Main.menuFont);
	nameField.setText(sample.getName());
	//sampleMenu.add(menuscroll);		
	//sampleMenu.setPreferredSize(new Dimension(300, 300));
	

	nameField.setCaretPosition(0);
	nameField.revalidate();
	sampleMenu.add(nameField);
	sampleMenu.add(sex);
	sampleMenu.pack();
	sampleMenu.show(this, moveX+(int)selectedSplit.pixel, moveY);
	*/
}

void addBam(int sampleindex) {
	try {			
		
		 if(VariantHandler.frame != null) {
		  VariantHandler.frame.setState(Frame.ICONIFIED);
		 }
		 boolean cram = false;
		 FileDialog fc = new FileDialog(Main.frame, "Choose BAM file", FileDialog.LOAD);
 		  fc.setDirectory(Main.path);
 		  fc.setFile("*.bam;*.cram;");
 		 fc.setFilenameFilter(new FilenameFilter() {
	  			public boolean accept(File dir, String name) {
			        return name.toLowerCase().endsWith(".bam") || name.toLowerCase().endsWith(".cram");
			     }
			 });
 		  fc.setMultipleMode(false);
 		  fc.setVisible(true);
 		  String openfile = fc.getFile();			
       
        if (openfile != null) {
         if(openfile.endsWith(".bai") || openfile.endsWith(".crai")) {
        	 Main.showError("Please, select BAM or CRAM, not the index file.", "Error");
        	 return;
         }
         File file = new File(fc.getDirectory() +openfile);
         if(!file.exists()) {
        	 Main.showError("File does not exist.", "Error");
	  		  return;
         }
       	  Main.path = fc.getDirectory();
       	  Main.writeToConfig("DefaultDir=" +Main.path);
       	 
       	Main.drawCanvas.bam = true;
  	 	Main.readsamples++;
  	 	Main.drawCanvas.sampleList.get(sampleindex).samFile = file;	 
  	 	Main.drawCanvas.sampleList.get(sampleindex).resetreadHash();
  	 	
  	 	
  	 	if(Main.readsamples==1) {
  	 		FileRead.checkSamples();
  	 	}
  	 
	  	  if(file.getName().endsWith(".cram")){			  			
	  			cram = true;		  		
	  	  }
	  	  else if(!file.getName().endsWith(".bam")) {
	  		  Main.showError("Open BAM or CRAM file.", "Error");
	  		  return;
	  	  }
	  	 
	  	Main.drawCanvas.sampleList.get(sampleindex).CRAM = cram;
	  	if(cram) {
	  		Main.drawCanvas.sampleList.get(sampleindex).readString = "CRAM";
	  	}
	  	else {
	  		Main.drawCanvas.sampleList.get(sampleindex).readString = "BAM";
	  	}
       	
        }
       /* else {
       	 Main.showError("File does not exist.", "Error");
        }*/
		 
	}
	catch(Exception ex) {
		Main.showError(ex.getMessage(), "Error");
	}
}

@Override
public void mouseClicked(MouseEvent event) {
	
	switch(event.getModifiers()) {	
	
		case 17: {
			if(sidebar) {
				
				if(secondSample == -1) {
					
					secondSample = (int)(pressY/drawVariables.sampleHeight);					
					
					if(secondSample-selectedSampleIndex >=0) {						
						
						drawVariables.visiblestart = (short)selectedSampleIndex;
						if(drawVariables.visiblestart < 0) {
							drawVariables.visiblestart = 0;
						}
						drawVariables.visiblesamples = (short)(secondSample-selectedSampleIndex+1);
						
						updateReads = true;
						updateCoverages = true;
					
						this.resizeCanvas(this.getWidth(), (int)(Main.samples*drawVariables.sampleHeight));	
						
						setScrollbar((int)(selectedSampleIndex*drawVariables.sampleHeight));	
						
					}
					else {
					
						drawVariables.visiblestart = (short)secondSample;
						if(drawVariables.visiblestart < 0) {
							drawVariables.visiblestart = 0;
						}
						drawVariables.visiblesamples = (short)(selectedSampleIndex-secondSample+1);
						updateReads = true;
						updateCoverages = true;
						this.resizeCanvas(this.getWidth(), (int)(Main.samples*drawVariables.sampleHeight));						
						setScrollbar((int)(secondSample*drawVariables.sampleHeight));						
					}					
				}
				secondSample = -1;
				break;
			}				
		}
		case InputEvent.BUTTON1_MASK: {
			
			if(event.getClickCount() == 2 && sidebar) {
				if(drawVariables.visiblesamples == 1) {				
					setScrollbar(0);
					drawVariables.visiblestart = 0;
					drawVariables.visiblesamples = (short)(Main.samples);
					checkSampleZoom();					
					this.resizeCanvas(this.getWidth(), (int)(Main.samples*drawVariables.sampleHeight));
				}
				else {
					drawVariables.visiblestart = (short)selectedIndex;
					if(drawVariables.visiblestart < 0) {
						drawVariables.visiblestart = 0;
						selectedIndex = (short)0;
					}
					drawVariables.visiblesamples = (short)1;
					this.resizeCanvas(this.getWidth(), (int)(Main.samples*drawVariables.sampleHeight));
					setScrollbar((int)(selectedIndex*drawVariables.sampleHeight));						
				}			
				break;
			}
			else if(event.getClickCount() == 2) {
				if(!FileRead.cancelreadread && coverageregion && selectedSplit.viewLength < Settings.readDrawDistance && Main.samples > 0) {
					try {
						Reads readclass = sampleList.get(selectedIndex).getreadHash().get(selectedSplit);		
						groupMismatchReads(Main.chromDraw.getPosition(event.getX()-Main.sidebarWidth, selectedSplit), selectedSplit,sampleList.get(selectedIndex), readclass);
					}
					catch(Exception e) {
						//TODO check this
					}
					//System.out.println(Main.chromDraw.getPosition(event.getX()-Main.sidebarWidth, selectedSplit));
					return;
				}
				if(drawVariables.visiblestart > this.selectedIndex) {
					drawVariables.visiblestart = (short)this.selectedIndex;
				}
				else if(drawVariables.visiblestart+drawVariables.visiblesamples-1 < this.selectedIndex) {
					drawVariables.visiblesamples++;
				}
				
				if(bam && zoomNote && splits.get(0).viewLength > Settings.coverageDrawDistance) {
							
				}	
				else if (bam && zoomNote && drawVariables.sampleHeight < 100) {
						
				}
				else {					
					if(FileRead.cancelreadread) {
						FileRead.cancelreadread = false;
						Main.putMessage(null);
						updateReads = true;
						repaint();
					}
				}
			}
			
			if(!sidebar && coverageregion && !splitremove) {
				
				return;
			}
			
			if(sidebar && bamHover) {
				
				addBam(selectedIndex);
				break;
			}
			if(varOverLap != null) {	
				MethodLibrary.showVariantMenu(this, varOverLap, sampleOverLap, moveX+(int)selectedSplit.pixel, moveY, "");
				varOverLap = null;
				break;
				
			}
			if(splits.size() > 1 && splitremove) {
				removeSplit(selectedSplit);
			}
			else if(removeSample != null) {
				
						
				removeSample(removeSample);
				
			}
			else if(Main.samples > 0) {				
			
				if(sidebar) {
					splitList.clear();
					
					selectedRead = null;
					selectedSampleIndex = (int)(pressY/drawVariables.sampleHeight);
					
					if(sampleList.get(selectedSampleIndex).removed) {
						
						for(int i =selectedSampleIndex+1; i<sampleList.size(); i++) {
							if(!sampleList.get(i).removed) {
								selectedSampleIndex = i;
								break;
							}
						}
					}
					if(selectedSample!=null && selectedSample.equals(sampleList.get(selectedSampleIndex))) {
						selectedSampleIndex = -1;
						selectedSample = null;
					}
					else {
						selectedSample = sampleList.get(selectedSampleIndex);
					}
					if(sampleInfo) {
						showsampleMenu(sampleList.get(selectedIndex));
					}
					else {
						if(SampleDialog.frame.isVisible()) {
							showsampleMenu(sampleList.get(selectedIndex));
						}
					}
					updatevars = false;
					repaint();
					break;					
				}
			
			if(event.getClickCount() == 2 && clickedRead != null && mouseRect.intersects(clickedRect)) {
				
				if(clickedRead.isDiscordant()) {
					selectedSplit.splitRead = true;
					if(!clickedRead.getMateChrom().equals(selectedSplit.chrom) || clickedRead.getMatePos() < readsample.getreadHash().get(selectedSplit).getReadStart() || clickedRead.getMatePos() > readsample.getreadHash().get(selectedSplit).getReadEnd()) {
						ArrayList<ReadNode[]> headAndTail = sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getHeadAndTail();
						ReadNode addread;						
						
						for(int j = 0; j<headAndTail.size();j++) {
							
							if(headAndTail.get(j).length == 0) {
								continue;
							}
							addread = headAndTail.get(j)[FileRead.headnode];
							
							while(addread != null) {
								if(!addread.isDiscordant()) {
									addread = addread.getNext();
								
									continue;
								}
								else if(!addread.getMateChrom().equals(clickedRead.getMateChrom())) {
									addread = addread.getNext();
									continue;
								}
								else {		
									if(sampleList.get(selectedIndex).getMates() == null) {
										sampleList.get(selectedIndex).setMates();
									}
								
									if(sampleList.get(selectedIndex).getMates().containsKey(addread.getName())) {
										
										if(!sampleList.get(selectedIndex).getMates().get(addread.getName()).contains(addread)) {
											
											sampleList.get(selectedIndex).getMates().get(addread.getName()).add(addread);	
											addread.mates = sampleList.get(selectedIndex).getMates().get(addread.getName());
											addread.split = selectedSplit;
										}
									}
									else {
										
										ArrayList<ReadNode> addlist = new ArrayList<ReadNode>();
										addlist.add(addread);
										addread.mates = addlist;
										sampleList.get(selectedIndex).getMates().put(addread.getName(), addlist);
										addread.split = selectedSplit;
										
									}
									
									addread = addread.getNext();
												
								}								
							}
						}
						
						headAndTail = null;
						addSplit(clickedRead.getMateChrom(), clickedRead.getMatePos()-(int)selectedSplit.viewLength/2, clickedRead.getMatePos()+(int)selectedSplit.viewLength/2);
					}
					else {
						
						ArrayList<ReadNode[]> headAndTail = sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getHeadAndTail();
						ReadNode addread;	
						if(sampleList.get(selectedIndex).getMates() == null) {
							sampleList.get(selectedIndex).setMates();
						}
						sampleList.get(selectedIndex).getMates().clear();
					
						for(int j = 0; j<headAndTail.size();j++) {
							
							addread = headAndTail.get(j)[FileRead.headnode];
							
							founddiscordant = false;
							while(addread != null) {
								
								if(!addread.isDiscordant()) {
									addread = addread.getNext();
									continue;
								}
								
								else if(!addread.getMateChrom().equals(clickedRead.getMateChrom())) {
									founddiscordant = true;
									
									addread = addread.getNext();
									continue;
								}
								else {		
									
									founddiscordant = true;
									if(sampleList.get(selectedIndex).getMates().containsKey(addread.getName())) {
										
										if(!sampleList.get(selectedIndex).getMates().get(addread.getName()).contains(addread)) {											
											
											addread.mates =sampleList.get(selectedIndex).getMates().get(addread.getName());
											addread.mates.add(addread);
											addread.split = selectedSplit;											
										}
									}
									else {
										
										ArrayList<ReadNode> addList = new ArrayList<ReadNode>();
										addList.add(addread);
										
										sampleList.get(selectedIndex).getMates().put(addread.getName(), addList);
										
										addread.mates = addList;
										addread.split = selectedSplit;											
									}										
									addread = addread.getNext();
								}	
								
							}
							if(!founddiscordant) {
								break;
							}
						}						
						headAndTail = null;
					}					
					
				/*	if(!Main.chromDraw.split || (Main.chromDraw.split && clickedRead.getRect().x < this.getDrawWidth()/2)) {
					//createSplit(sampleList.get(selectedSampleIndex), clickedRead.getMateChrom().replace("chr", ""), clickedRead.getMatePos());
					
			//		Main.drawCanvas.add(new Draw(this.getWidth()/2, this.getHeight()));
					
				//	this.drawWidth = (int)((this.width - Main.sidebarWidth)/(double)splits.size());
					
					//setStartEnd(clickedRead.getPosition()-this.viewLength/4, clickedRead.getPosition()-this.viewLength/4 + this.viewLength);
				}
				else {*/
					/*try {
						
						tempreadname = clickedRead.getName();
						
					//	setStartEnd(tempstart-this.viewLength/4, tempstart-this.viewLength/4 + this.viewLength);
					}
					catch(Exception e) {
						ErrorLog.addError(e.getStackTrace());
						e.printStackTrace();
					}
			//	}*/
					break;
				}
				else if (clickedRead.SA != null) {
					
					String[] SAs = clickedRead.SA.split(";");
					ArrayList<ReadNode[]> headAndTail = sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getHeadAndTail();
					ReadNode addread;		
					
					if(sampleList.get(selectedIndex).getMates() == null) {
						sampleList.get(selectedIndex).setMates();
					}
			//		sampleList.get(selectedIndex).getMates().clear();
					
					for(int j = 0; j<headAndTail.size();j++) {
						
						addread = headAndTail.get(j)[FileRead.headnode];
						
						while(addread != null) {
								if(addread.SA == null) {
									addread = addread.getNext();				
									continue;
								}
								if(addread.getName().contains("/")) {
									if(!sampleList.get(selectedIndex).getMates().containsKey(addread.getName().substring(0,addread.getName().lastIndexOf("/")))) {
										ArrayList<ReadNode> addList = new ArrayList<ReadNode>();
										addList.add(addread);
										sampleList.get(selectedIndex).getMates().put(addread.getName().substring(0,addread.getName().lastIndexOf("/")), addList);									
										addread.mates = addList;
										addread.split = selectedSplit;
									}
									else {					
										if(!sampleList.get(selectedIndex).getMates().get(addread.getName().substring(0,addread.getName().lastIndexOf("/"))).contains(addread)) {
											addread.mates =sampleList.get(selectedIndex).getMates().get(addread.getName().substring(0,addread.getName().lastIndexOf("/")));										
											addread.mates.add(addread);
											addread.split = selectedSplit;		
										}
										else {
											
											break;
										}
									}		
								}
								else {
									if(!sampleList.get(selectedIndex).getMates().containsKey(addread.getName())) {
										ArrayList<ReadNode> addList = new ArrayList<ReadNode>();
										addList.add(addread);
										sampleList.get(selectedIndex).getMates().put(addread.getName(), addList);									
										addread.mates = addList;
										addread.split = selectedSplit;
									}
									else {					
										if(!sampleList.get(selectedIndex).getMates().get(addread.getName()).contains(addread)) {
											addread.mates =sampleList.get(selectedIndex).getMates().get(addread.getName());										
											addread.mates.add(addread);
											addread.split = selectedSplit;		
										}
										else {
											
											break;
										}
									}		
								}
								addread = addread.getNext();														
						}
						
					}	
					
					addread = null;
					headAndTail = null;
					boolean found = false;
					int pos;
					ArrayList<splitTuple> splitList = new ArrayList<splitTuple>();					
					
					for(int i = 0; i<SAs.length; i++) {
						String[] sa = SAs[i].split(",");
						String chr = sa[0];
						found = false;
						pos = Integer.parseInt(sa[1]);		
						for(int s = 0 ; s<splits.size(); s++) {
							if(chr.equals(splits.get(s).chrom) && pos > readsample.getreadHash().get(splits.get(s)).getCoverageStart() && pos < readsample.getreadHash().get(splits.get(s)).getCoverageEnd()) {
								found = true;			
							}
						}
						if(!found) {
							splitList.add(new splitTuple(chr, pos));		
						}																		
					
					}
					if(splitList.size() > 1) {
						for(int i = 0; i<splitList.size()-1; i++) {
							for(int j = 1; j<splitList.size(); j++) {
								if(splitList.get(i).chr.equals(splitList.get(j).chr)) {
									if(Math.abs(splitList.get(i).pos - splitList.get(j).pos) < selectedSplit.viewLength+10000 ) {
										splitList.remove(j);
									}
								}
							}
						}
					}
					
					if(splitList.size() > 0) {
						for(int i = 0; i<splitList.size(); i++) {
							addSplit(splitList.get(i).chr, splitList.get(i).pos-(int)selectedSplit.viewLength/2, splitList.get(i).pos+(int)selectedSplit.viewLength/2);	
						}
					}
			//		saReads = new SAread(clickedRead);
					
				
				}
				if(clickedRead.getMates() != null) {
					Collections.sort(clickedRead.getMates(), mateSorter);
					/*if(saReads != null) {
						for(int i = 0 ; i<clickedRead.getMates().size();i++) {						
							for(int r = 0 ; r<saReads.reads.size(); r++) {
									if((int)saReads.reads.get(r)[SAread.pos] == clickedRead.getMates().get(i).getPosition()) {										
										saReads.reads.get(r)[SAread.SARead] = clickedRead.getMates().get(i);
								
									}							
							}
						}
					}*/					
				}
				
				repaint();
			}
			
			else if(selectedRead !=null) {
				saReads = null;
				splitList.clear();
				clickedRead = selectedRead;
				clickedReadSample = sampleList.get(selectedIndex);
				clickedReadInfo = createReadInfo(clickedRead);	
				
				searchMate(clickedRead);
				
				
				if(clickedRead.getMates() != null) {
					
					Collections.sort(clickedRead.getMates(), mateSorter);
					for(int i = 0 ; i<clickedRead.getMates().size()-1;i++) {
						
						if(clickedRead.getMates().get(i).getPosition() == clickedRead.getMates().get(i+1).getPosition()) {
							
							clickedRead.getMates().remove(i+1);
							i--;
						}
					}
					
		//			clickmates = (ArrayList<ReadNode>)clickedRead.getMates().clone();
				}		
				for(int i = 0; i <splits.size(); i++) {
					splits.get(i).updateReads = true;
					
				}
				if(clickedRead.SA != null ) {					
					saReads = new SAread(clickedRead);				
			}
				
				repaint();			
			}		
			}			
			
			break;
		}
		case InputEvent.BUTTON3_MASK: {
			if(coverageregion) {
				
				return;
			}
			if(clickedRead != null && mouseRect.intersects(clickedRead.getRect())) {
				//moveReadToBottom(clickedRead, selectedSplit, sampleList.get(selectedIndex), readLevel);
				return;
			}
			if(hoverMate != null && (mouseRect.intersects(hoverMate.getRect()) || (mouseRect.intersects(hoverRect)))) {				
				
				SAMRecord read = Main.fileReader.getRead(hoverMate.split.chrom, hoverMate.getPosition()-hoverMate.getStartOffset(), hoverMate.getPosition()+hoverMate.getWidth(),hoverMate.getName(), clickedReadSample.getreadHash().get(hoverMate.split));
				String myString = read.getReadString();
				StringSelection stringSelection = new StringSelection(myString);
				Clipboard clpbrd = Toolkit.getDefaultToolkit().getSystemClipboard();
				clpbrd.setContents(stringSelection, null);
				Main.putMessage("Read sequence copied to clipboard.");
				//Main.chromDraw.message = "Sequence copied to clipboard.";
				//Main.chromDraw.timer = System.currentTimeMillis();
				//Main.chromDraw.repaint();
			}
			else if (clickedRead != null && mouseRect.intersects(clickedRead.getRect())) {
				SAMRecord read = Main.fileReader.getRead(clickedRead.split.chrom, clickedRead.getPosition()-clickedRead.getStartOffset(),clickedRead.getPosition()+clickedRead.getWidth(), clickedRead.getName(), clickedReadSample.getreadHash().get(clickedRead.split));
				String myString = read.getReadString();
				StringSelection stringSelection = new StringSelection(myString);
				Clipboard clpbrd = Toolkit.getDefaultToolkit().getSystemClipboard();
				clpbrd.setContents(stringSelection, null);
				Main.putMessage("Read sequence copied to clipboard.");
			
			}
			else {
				clickedRead = null;
				saReads = null;
			//	clickmates = null;
				splitList.clear();
				for(int i = 0; i <splits.size(); i++) {
					splits.get(i).updateReads = true;
					
				}
			}
			repaint();
		}
	}
	
}
void groupMismatchReads(int position, SplitClass split, Sample readsample, Reads readClass) {
	if(readClass == null || readClass.getReads() == null) {
		return;
	}
	int addlevel = 0;
	
	ArrayList<ReadNode> removedReads = new ArrayList<ReadNode>();
	for(int j = 0; j<readClass.getReads().size();j++) {
			
			try {
				
			read = readClass.getReads().get(j);
			boolean moved = false;
			while(read != null) {			
				
				if(read.getPosition()+read.getWidth() < position) {
					read = read.getNext();
					continue;
				}
				if(read.getPosition() > position) {
					break;
				}
				if(read.getMismatches() == null) {
					read = read.getNext();
					continue;
				}
				
				for(int i = 0 ; i<read.getMismatches().size(); i++) {
					
					if(read.getPosition() +read.getMismatches().get(i).getKey() == position) {
						
						moveReadToBottom(read, split, readsample, j, addlevel, readClass, removedReads);
						
						addlevel++;
						moved = true;
						break;
					}
				}
				if (moved) {
					break;
				}
				read = read.getNext();
			}
			}
			catch(Exception e) {
				e.printStackTrace();
			}
	}
	

	for(int i = 0 ; i< removedReads.size(); i++) {
		//System.out.println(removedReads.get(i).getPosition());
		moveReadUp(removedReads.get(i), split, readsample, addlevel,readClass);
	}
	//System.out.println(removedReads.size());
	removedReads = null;
	Draw.updateReads = true;
	repaint();
}
void removeRead(ReadNode removeread) {
	if(removeread.getPrev() != null) {
		if(removeread.getNext() != null) {
			removeread.getPrev().setNext(removeread.getNext());
		}
		else {
			removeread.getPrev().setNext(null);
		}
	}
	if(removeread.getNext() != null) {
		if(removeread.getPrev() != null) {
			removeread.getNext().setPrev(removeread.getPrev());
		}
		else {
			removeread.getNext().setPrev(null);
		}
	}
	removeread.setNext(null);
	removeread.setPrev(null);
}
void moveReadUp(ReadNode moveread, SplitClass split, Sample readsample, int startrow, Reads readClass) {
	boolean found = false;
	
	
	for(int j = startrow; j<readClass.getReads().size();j++) {
		
		try {
			
			read = readClass.getReads().get(j);
		
		while(read != null) {
			
			if(read.getNext() != null) {
				if(read.getPosition()+read.getWidth() < moveread.getPosition() && read.getNext().getPosition() > moveread.getPosition()+moveread.getWidth()) {
					read.getNext().setPrev(moveread);
					moveread.setNext(read.getNext());
					read.setNext(moveread);
					
					moveread.setPrev(read);
					
					
					//System.out.println(j +" " +moveread.getName() +":" +moveread.getPosition() +" " +read.getName()+":" +read.getPosition());
					found = true;
					break;
				}
			}
			if(read.getPosition() > moveread.getPosition()+moveread.getWidth()) {
				break;
			}
			
			read = read.getNext();
					
		}
		if(found) {
			break;
		}
		
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	/*if(!found) {
	
		readClass.getReads().add(moveread);
		ReadNode[] addList = new ReadNode[2];
		addList[FileRead.headnode] = moveread;
		addList[FileRead.tailnode] = moveread;
		moveread.setPrev(null);
		moveread.setNext(null);
		readClass.getHeadAndTail().add(addList);
	}*/
}
void moveReadToBottom(ReadNode moveread, SplitClass split, Sample readsample, int readlevel, int addlevel, Reads readClass, ArrayList<ReadNode> removed) {
	clickedRead = null;
	selectedRead = null;
	
	if(readlevel == addlevel) {
		return;
	}
	removeRead(moveread);
	
	
	ReadNode read, nextRead;
	
	//boolean moved = false;
	
	try {	
	//for(int j = addlevel; j<readClass.getReads().size();j++) {		
		//if(j == readlevel) {
			//break;
	//	}
		
		read = readClass.getReads().get(addlevel);
		//moved = false;
		while(read != null) {			
			
			if(read.getPosition()+read.getWidth() < moveread.getPosition()) {
				read = read.getNext();
				continue;
			}
			if(read.getPosition() > moveread.getPosition()+moveread.getWidth()) {
				break;
			}
			
			
			//moveread.moved = true;
			if(read.equals(readClass.getReads().get(addlevel))) {
				readClass.getReads().set(addlevel,readClass.getHeadAndTail().get(addlevel)[0]);
				
			}
			nextRead = read.getNext();
			if(nextRead == null) {
				System.out.println(addlevel);
			}
			removed.add(read);
			
			removeRead(read);
			//moveReadUp(read, split, readsample, addlevel,readClass);
			read = nextRead;
			
			
		}		
		//if(moved) {
		if(read != null) {	
			if(read.getPrev() != null) {
				read.getPrev().setNext(moveread);
				moveread.setPrev(read.getPrev());
			}
			else {
				moveread.setPrev(null);
			}
			read.setPrev(moveread);
			moveread.setNext(read);
		}
			//	read = null;
		//	break;
	//	}
		
//	}
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	
}
private class splitTuple {
	String chr;
	int pos;
	
	public splitTuple(String chr, int pos) {
		this.chr = chr;
		this.pos = pos;
	}
}



void removeSample(Sample removeSample) {
	if(Main.samples < 1) {
		Main.samples = 0;
		return;
	}
	if(removeSample.multipart) {
		VarNode node = FileRead.head.getNext();
		
		while(node != null) {					
			node.removeSample(removeSample);					
			node = node.getNext();					
		}
		
		node = null;
		
		if(removeSample.samFile != null) {
			Main.readsamples--;
		}
		if(!removeSample.annotation) {
			Main.varsamples--;		
		}
		//sampleList.remove(removeSample);
		removeSample.removed = true;
		Main.samples--;
		
		if(Main.samples > 0) {
			if(removeSample.getIndex() < Main.samples) {	
				
				this.selectedSample = sampleList.get(removeSample.getIndex());
				this.selectedSampleIndex = removeSample.getIndex();
			}
			else {
				this.selectedSample = sampleList.get(removeSample.getIndex()-1);
				this.selectedSampleIndex = removeSample.getIndex()-1;
			}
		}	
		//TODO
	}
	else if(removeSample.multiVCF) {
		
		remindex = removeSample.getIndex()+1;
		
		Sample remsample;
		for(int i = remindex; i<Main.samples; i++) {
			remsample = sampleList.get(i);
			if(!remsample.multipart) {
				break;
			}
			VarNode node = FileRead.head.getNext();
			
			while(node != null) {					
				node.removeSample(remsample);					
				node = node.getNext();					
			}			
			if(remsample.samFile != null) {
				Main.readsamples--;
				remsample.samFile = null;
			}
			
			if(!removeSample.annotation) {
			Main.varsamples--;
			}
			sampleList.remove(remsample);
			Main.samples--;		
			
			
			node = null;
			i--;
		}
		remsample = null;
		Main.samples--;
		
		sampleList.remove(removeSample);
		
		if(Main.samples > 0) {
			if(removeSample.getIndex() < Main.samples) {
				this.selectedSample = sampleList.get(removeSample.getIndex());
				this.selectedSampleIndex = removeSample.getIndex();
			}
			else {
				this.selectedSample = sampleList.get(removeSample.getIndex()-1);
				this.selectedSampleIndex = removeSample.getIndex()-1;
			}
		}
		removeSample = null;
		this.removeSample = null;	
	}
	else {
		VarNode node = FileRead.head.getNext();
	
		while(node != null) {					
			node.removeSample(removeSample);					
			node = node.getNext();					
		}
		selectedIndex = 0;
		node = null;
		if(removeSample.calledvariants) {
			if(!removeSample.annotation) {
				Main.varsamples--;
			}
		}
		if(removeSample.samFile != null) {
			
			Main.readsamples--;
			
			removeSample.samFile = null;
		}
		if(removeSample.getTabixFile() != null) {
			if(!removeSample.annotation) {
				Main.varsamples--;
			}
			else {
				for(int i = 0 ; i <sampleList.size(); i++) {
					
				}
			}
		}
		sampleList.remove(removeSample);
		if(removeSample.annotation) {
			boolean found = false;
			for(int i = 0 ; i <sampleList.size(); i++) {
				if(sampleList.get(i).annotation) {
					found = true;
					break;
				}
			}
			annotationOn = found;
		}
		
		Main.samples--;
		if(Main.samples > 0) {
			if(removeSample.getIndex() < Main.samples) {		
				this.selectedSample = sampleList.get(removeSample.getIndex());
				this.selectedSampleIndex = removeSample.getIndex();
			}
			else {
				this.selectedSample = sampleList.get(removeSample.getIndex()-1);
				this.selectedSampleIndex = removeSample.getIndex()-1;
			}
		}	
	}
	VariantHandler.commonSlider.setMaximum(Math.max(Main.varsamples, 1));
	VariantHandler.geneSlider.setMaximum(Math.max(Main.varsamples, 1));
	if(VariantHandler.commonSlider.getUpperValue() > Main.varsamples) {
		VariantHandler.commonSlider.setUpperValue(Main.varsamples);
		
	}
	if(VariantHandler.commonSlider.getValue() > Main.varsamples) {
		VariantHandler.commonSlider.setValue(Math.max(Main.varsamples, 1));
	}
	if(VariantHandler.geneSlider.getValue() > Main.varsamples) {
		VariantHandler.geneSlider.setValue(Math.max(Main.varsamples, 1));
	}
	checkSampleIndices();
	FileRead.checkSamples();
	
	this.checkSampleZoom();
	updatevars = true;
	currentDraw = FileRead.head.getNext();
	current = FileRead.head.getNext();
	Main.chromDraw.vardraw = null;
	Main.chromDraw.varnode = null;
	removeSample = null;
	resizeCanvas(getWidth(), Main.drawScroll.getViewport().getHeight());
	sidebar = false;
	repaint();
}

void checkSampleIndices() {
	short index = 0;
	for(int j = 0; j<sampleList.size(); j++) {
		if(sampleList.get(j).removed) {
			continue;
		}
		sampleList.get(j).setIndex(index);
		index++;
	}
}

String[] createReadInfo(ReadNode read) {
	String[] temp = new String[11];
	temp[0] = "" +read.getName();
	if(read.isForward()) {
		temp[1] = "Position: chr" +read.split.chrom +":" +MethodLibrary.formatNumber(read.getPosition()) +" (+)";
	}
	else {
		temp[1] = "Position: chr" +read.split.chrom +":" +MethodLibrary.formatNumber(read.getPosition()) +" (-)";
		
	}
	temp[2] = "Mapping quality: " +read.getMappingQuality();
	temp[3] = "Insert size: " +MethodLibrary.formatNumber(read.getInsertSize()) +"bp";
	if(read.getCigar() == null) {
		temp[4] = "Cigar: " +read.getWidth() +"M";
		temp[5] = "Length: " +read.getWidth() +"bp";
	}
	else {
		if(read.getCigar().toString().length() > 50) {
			temp[4] = "Cigar: " +read.getCigar().toString().substring(0,20)  +" ... " +read.getCigar().toString().substring(read.getCigar().toString().length()-20,read.getCigar().toString().length());
		}
		else {
			temp[4] = "Cigar: " +read.getCigar().toString();
			
		}
		temp[5] = "Length: " +read.getWidth() +"bp";
	}				
	if(read.getMatePos() == -1) {
		
		temp[6] = "No mate.";
		temp[7] = "Primary alignment: " +read.getPrimary();
		temp[8] = "";								
		temp[9] = "";
	}
	else {
		if(clickedRead.isMateForward()) {
			temp[6] = "Mate position: chr" +read.getMateChrom()+":" +MethodLibrary.formatNumber(read.getMatePos()) +" (+)";
		}
		else {
			temp[6] = "Mate position: chr" +read.getMateChrom()+":" +MethodLibrary.formatNumber(read.getMatePos()) +" (-)";
			
		}
		if(read.isDiscordant()) {
			temp[7] = "Primary alignment: " +read.getPrimary();
			temp[8] = "Double click to see the mate.";								
			temp[9] = "";
		}
		else {
			temp[7] = "Primary alignment: " +read.getPrimary();
			temp[8] = "";								
			temp[9] = "";
		}
	}
	if (read.SA != null) {
		String[] SAs = read.SA.split(";");
		temp[6] = "Split position(s): ";
		for(int i = 0; i<SAs.length; i++) {
			String[] sa = SAs[i].split(",");
			sa[0] = sa[0].replace("chr","");
			temp[6] += "chr" +sa[0] +":" +MethodLibrary.formatNumber(Integer.parseInt(sa[1]));
			if(i < SAs.length-1) {
				temp[6] += ", ";
			}
		}
		temp[7] = "Primary alignment: " +read.getPrimary();
		temp[8] = "Double click to see the splitted read.";								
		temp[9] = "";
	}
	for(int i = 0; i<temp.length-1; i++) {
		if(temp[i].length() > temp[maxwidth].length()) {
			maxwidth = i;
		}
	}
	temp[10] = temp[maxwidth];
	return temp;
	
}
@Override
public void mouseEntered(MouseEvent event) {
	
}

@Override
public void mouseExited(MouseEvent event) {
	
}

@Override
public void mousePressed(MouseEvent event) {
	
	if(!this.hasFocus()) {
		this.requestFocus();
	}
	pressX = event.getX();
	pressY = event.getY();
	moveYtemp = pressY;
	
	tempDrag = pressX;	
	tempViewLength = selectedSplit.viewLength;
	switch(event.getModifiers()) {	
		case InputEvent.BUTTON1_MASK: {
			if(getMoreVariants) {
				FileRead.search = true;
				
				gotoPos(splits.get(0).chrom, splits.get(0).start, splits.get(0).end);
				getMoreVariants = false;
				
				break;
			}
			
			
			if(this.sampleList.size() > 0 && selectedIndex > -1 && this.sampleList.get(selectedIndex).getreadHash() != null && this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit) != null && this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).isReadScroll()) {
				
				if(mouseRect.intersects(this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getScrollBar())) {
					
					if(!mouseRect.intersects(this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getScroller())) {
						
						if(mouseRect.y < this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getScroller().y) {
							double sampleheight = this.drawVariables.sampleHeight-(this.drawVariables.sampleHeight/selectedSplit.getDivider());
							double totalheight = this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getReadSize()*(this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2);
							this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel += this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getScroller().height/(sampleheight/totalheight);
							
							if(( (Main.drawCanvas.drawVariables.sampleHeight-(Main.drawCanvas.drawVariables.sampleHeight/selectedSplit.getDivider()))) > (Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getReadSize()*(this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2)) - (Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel)) {
								Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel = (int)((Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getReadSize()*(this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2)) - ((Main.drawCanvas.drawVariables.sampleHeight-(Main.drawCanvas.drawVariables.sampleHeight/selectedSplit.getDivider()))));
								
							}
						}
						else {
							double sampleheight = this.drawVariables.sampleHeight-(this.drawVariables.sampleHeight/selectedSplit.getDivider());
							double totalheight = this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getReadSize()*(this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readHeight+2);
							
							this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel -= this.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).getScroller().height/(sampleheight/totalheight);
							
							if(Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel < 0) {
								Main.drawCanvas.sampleList.get(selectedIndex).getreadHash().get(selectedSplit).readwheel = 0;
								
							}
						}
						updateReads = true;
						repaint();
					}
					else {
						readScrollDrag = true;
					}
				}
				else {
					if(moveX > Main.sidebarWidth+3) {
						tempSample = (short)(pressY/drawVariables.sampleHeight);		
						task = new MyTimerTask();
						navTimer.scheduleAtFixedRate(task, 300, 1000);
					}
					if(!lineZoomer && moveX > Main.sidebarWidth-3 && moveX < Main.sidebarWidth+3) {
						resizeSidebar = true;
					}
					else {
						resizeSidebar = false;
					}
				}
				
			}
			else {
				if(moveX > Main.sidebarWidth+3) {
					tempSample = (short)(pressY/drawVariables.sampleHeight);		
					task = new MyTimerTask();
					navTimer.scheduleAtFixedRate(task, 300, 1000);
				}
				if(!lineZoomer && moveX > Main.sidebarWidth-3 && moveX < Main.sidebarWidth+3) {
					resizeSidebar = true;
				}
				else {
					resizeSidebar = false;
					
				}
			}
			break;
		}
		case InputEvent.BUTTON3_MASK: {
			if(sidebar) {
				if(selectedIndex < 0) {
					return;
				}
				showsampleMenu(sampleList.get(selectedIndex));
			}
		}
	}	
}


private class MyTimerTask extends TimerTask {
    public void run() {
        mousePress = true;
        Main.drawCanvas.repaint();
        
    }
}
public void gotoPos(double start, double end) {
//	clearReads();
/*	if(start > selectedSplit.chromEnd || end > selectedSplit.chromEnd) {
		return;
	}*/
	
	setStartEnd(start, end);
/*	if(!lineZoomer && !zoomDrag) {
		
		clearReads();
		Main.chromDraw.sequence = Main.chromDraw.getSeq(Main.chromosomeDropdown.getItemAt(Main.selectedChrom),(int)start-5000, (int)end+5000, Main.referenceFile);	
		Main.chromDraw.readSeqStart = (int)start-5000;
	}*/
	
//	Draw.updatevars = true;
//	Draw.updateReads = true;
//	updateCoverages = true;
//	splits.get(0).sequence = null;
	
//	Main.chromDraw.updateExons = true;
//	Main.chromDraw.repaint();
//	Main.bedCanvas.repaint();
//	FileRead.search = false;
//	repaint();
}

public void gotoPos(String chrom, double start, double end) {
	if(Draw.variantcalculator) {
		return;
	}
	Draw.updateReads = false;
	Draw.updatevars = false;	
	splits.get(0).transStart = 0;
	clearReads();
	if(FileRead.search) {
		splits.get(0).nullRef();
		
		removeSplits();
		if(Main.undoPointer < Main.undoList.size()-1) {
			Main.undoList.add(Main.undoPointer+1,chrom+":"+(int)start+"-"+(int)end);
			if(Main.undoPointer < Main.undoList.size()-1) {
				for(int i = Main.undoPointer+2 ; i< Main.undoList.size(); i++) {
					Main.undoList.remove(i);
					i--;
				}
			}
		}
		else {
			Main.undoList.add(chrom+":"+(int)start+"-"+(int)end);
		}
		Main.undoPointer = Main.undoList.size()-1;
		
		if(!FileRead.novars) {
			if(current != null && current.getPosition()+2 < start && current.getNext() != null) {
				if(current.getNext().getPosition()+2 < start) {
					while(current.getPosition()+2 < start) {				
						current = current.getNext();
						if(current == null) {
							break;
						}
					}
				}
			}
			else {
				while(current != null && !current.equals(FileRead.head) && current.getPosition()-1 > start) {				
					current = current.getPrev();
				}
			}
		}
		if(forcereload || (Main.chromosomeDropdown.getSelectedItem() != null &&!Main.chromosomeDropdown.getSelectedItem().toString().equals(chrom)) || current == null || current.equals(FileRead.head) || current.getNext() == null || variantsEnd < splits.get(0).end || variantsStart > splits.get(0).start) {
			if(Main.chromosomeDropdown.getSelectedItem() != null) {
			FileRead.searchStart = (int)start;
			FileRead.searchEnd = (int)end;
			if(FileRead.searchStart < 0) {
				FileRead.searchStart = 0;
			}
			if(Main.chromIndex.get(Main.refchrom + chrom) == null) {
				Main.showError("Select correct genome.", "Error.");
				return;
			}
			if(FileRead.searchEnd > Main.chromIndex.get(Main.refchrom + chrom)[1].intValue()) {
				FileRead.searchEnd = Main.chromIndex.get(Main.refchrom + chrom)[1].intValue();
			}
			if(!Main.chromosomeDropdown.getSelectedItem().toString().equals(chrom) ) {			
				for(int i = 0; i<Main.bedCanvas.bedTrack.size(); i++) {
					Main.bedCanvas.bedTrack.get(i).getHead().putNext(null);
				}				
			}
			
			
			Main.chromosomeDropdown.setSelectedItem(chrom);
			}
			
		}
		else {
			Main.drawCanvas.loading("Loading Tracks...");
			for(int i = 0; i<Main.bedCanvas.bedTrack.size(); i++) {
				if(!Main.bedCanvas.bedTrack.get(i).small) {
					Main.bedCanvas.bedTrack.get(i).getHead().putNext(null);
					
					Main.bedCanvas.getBEDfeatures(Main.bedCanvas.bedTrack.get(i),(int)start, (int)end);
					
				}
			}	
			Main.drawCanvas.ready("Loading Tracks...");
		}
		//Main.bedCanvas.repaint();
		if(FileRead.novars) {
			this.variantsStart = 0;
			this.variantsEnd = splits.get(0).chromEnd;
		}
	}
	
	setStartEnd(start, end);
	if(Main.samples > 0) {
		for(int i = drawVariables.visiblestart; i< drawVariables.visiblestart+drawVariables.visiblesamples; i++) {
			if(i> Main.samples-1) {
				break;
			}
			sampleList.get(i).prepixel = -1;
		}
	}
	/*if(Main.bedCanvas.bedTrack.size() > 0) {
		
		bedFeatureFetcher fetch = Main.bedCanvas.new bedFeatureFetcher();
		
		fetch.execute();
	}*/
	
	Main.chromDraw.updateExons = true;	
	Draw.updateReads = true;
	updateCoverages = true;
	Draw.updatevars = true;
	repaint();
	
	Main.chromDraw.repaint();	
	
	
}

void setStartEnd(double start, double end) {
	
	starttemp = start;
	this.endtemp = end;
	
	if(end -start < minzoom ){
		this.endtemp = start+((end-start)/2)+minzoom/2;
		starttemp = start+((end-start)/2)-minzoom/2;
	}
	if(end > selectedSplit.chromEnd) {
		this.endtemp = selectedSplit.chromEnd;
	}
	if(start < 0) {
		starttemp = 0;
	}	
	updateReads = true;
	updateCoverages = true;
//	SplitClass.updateReads = true;
	selectedSplit.start = starttemp;
//	selectedSplit.end = endtemp;
/*	if(Main.chromDraw.split) {
		
		selectedSplit.pixel = (this.getDrawWidth())/2/(double)(selectedSplit.end-selectedSplit.start);
		
		selectedSplit.viewLength = selectedSplit.end-selectedSplit.start;
	}
	else {
	*/	
	selectedSplit.end = endtemp;
	selectedSplit.viewLength = selectedSplit.end-selectedSplit.start;
	selectedSplit.pixel = (this.getDrawWidth())/(double)(selectedSplit.viewLength);
	
	if((int)splits.get(0).start == 0 && (int)splits.get(0).end == splits.get(0).chromEnd) {
		Main.zoomout.setEnabled(false);
	}
	else {
		if(!Main.zoomout.isEnabled()) {
			Main.zoomout.setEnabled(true);
		}
		
	}
		
//	}
	Main.updatePositions(selectedSplit.start, selectedSplit.end);
	
}


Color getColor(ReadNode read, String chrom) {
	
	if(!read.getMateChrom().equals("") && !read.getMateChrom().equals(chrom)) {	
		return getChromColor(read.getMateChrom());
		
	}
	
	if(read.getMatePos() != -1) {
		if(read.isDiscordant() && read.isForward() == read.isMateForward()) {
			return Color.blue;
		}
		
		if(Math.abs(read.getInsertSize()) > Settings.insertSize) {			
			
			return Color.green;
		}
		if(read.getInsertSize() < 0) {
			if(read.getMatePos() > read.getPosition()) {					
				return Color.red;
			}
			
		}
		else {
			if(read.getMatePos() < read.getPosition()) {					
				return Color.red;
			}
		}	
	}
	return Color.lightGray;
}
Color getChromColor(String chrom) {
	if(chrom.matches("\\d+")) {
		switch (Integer.parseInt(chrom)) {	
		
		//chr1
			case 1:
				readColor = chr1color;
				break;
				//chr2
			case 2:
				readColor = chr2color;
				break;
				//chr3
			case 3:														
				readColor = chr3color;
				break;
				//chr4
			case 4:
				readColor = chr4color;
				break;
				//chr5
			case 5:
				readColor = chr5color;
				break;
				//chr6
			case 6:
				readColor = chr6color;
				break;
				//chr7
			case 7:
				readColor = chr7color;
				break;
				//chr8
			case 8:
				readColor = chr8color;
				break;
				//chr9
			case 9:
				readColor = chr9color;
				break;
				//chr10
			case 10:
				readColor = chr10color;
				break;
				//chr11
			case 11:														
				readColor = chr11color;
				break;
				//chr12
			case 12:														
				readColor = chr12color;
				break;
				//chr13
			case 13:														
				readColor = chr13color;
				break;
				//chr14
			case 14:
				readColor = chr14color;
				break;
				//chr15
			case 15:
				readColor = chr15color;
				break;
				//chr16
			case 16:
				readColor = chr16color;
				break;
				//chr17
			case 17:
				readColor = chr17color;
				break;
				//chr18
			case 18:
				readColor = chr18color;
				break;
				//chr19
			case 19:
				readColor = chr19color;
				break;
				//chr20
			case 20:
				readColor = chr20color;
				break;
				//chr21
			case 21:
				readColor = chr21color;
				break;
				//chr22
			case 22:
				readColor = chr22color;
				break;
			}
		}												
	else {													
			switch (chrom.getBytes()[0]) {
					//chrX
				case 88:
					readColor = chrXcolor;
					break;
					//chrY
				case 89:
					readColor = chrYcolor;
					break;
					//chrM
				case 77:
					readColor = chrMTcolor;
					break;														
			}
			return Color.cyan;
											
	}		
	return readColor;
}
void zoom() {	
	
		if(lineZoomer) {			
			gotoPos(selectedSplit.start-(tempDrag-mouseX)/selectedSplit.pixel*2, selectedSplit.end+(tempDrag-mouseX)/selectedSplit.pixel*2);			
			tempDrag = mouseX;
		}
		else if(tempSample <= mouseY/drawVariables.sampleHeight && pressX < mouseX) {
				if(selectedSplit.viewLength < minzoom +2) {
					return;
				}
				if(Main.samples > 0) {								
					gotoPos((int)(selectedSplit.start+((pressX-selectedSplit.offset)/selectedSplit.pixel)), (int)(selectedSplit.start+((mouseX-selectedSplit.offset)/selectedSplit.pixel)));
					
				}
				else {					
					gotoPos((int)(selectedSplit.start+((pressX-selectedSplit.offset)/selectedSplit.pixel)), (int)(selectedSplit.start+((mouseX-selectedSplit.offset)/selectedSplit.pixel)));
				}
		}	
		else {
			
			if(drawVariables.sampleHeight*Main.samples > Main.drawScroll.getViewport().getHeight() ) {				
				//	TODO
				//	drawVariables.visiblestart = 0;
				//	drawVariables.visibleend = (short)(Main.samples-1);								
				//	this.resizeCanvas(this.getWidth(), Main.drawScroll.getViewport().getHeight());				
				//	setScrollbar(0);						
				//	repaint();
			}
		
	
		}
	if(selectedSplit.viewLength > Settings.readDrawDistance) {
		
	//	clickedRead = null;
	//	clickmates = null;
		splitList.clear();
	//	clickedReadSplit = null;
	}
	
}

void release() {
	Main.bedCanvas.repaint();
	task.cancel();
	timer = 0;
	zoomDrag = false;
	mouseDrag = false;	
	
	lineZoomer = false;
	mousePress = false;
	prezoom = 0;
	
	
	updateReads = true;
	updateCoverages = true;
	sampleZoomer = false;
	scrollbar = false;
}
void moveSample(int move, int where) {	
	
	Sample movesample = this.sampleList.get(move);
	this.sampleList.remove(move);
	this.sampleList.add(where, movesample);
	
	checkSampleIndices();
	if(selectedSampleIndex != -1) {
		selectedSampleIndex = where;
		selectedIndex = -1;
	}
	//repaint();
	
	VarNode node = FileRead.head.getNext();
	
	while(node != null) {					
		node.moveSample(movesample);					
		node = node.getNext();					
	}
	
	node = null;
	movesample = null;
}
public void getReads(SplitClass split) {
	if(split.viewLength <= Settings.readDrawDistance && drawVariables.sampleHeight > 100) {		
		
		split.clearedReads = false;		
		
		FileRead reader = new FileRead();
		reader.splitIndex = split;
		reader.getreads = true;						
		reader.execute();
	}
}
@Override
public void mouseReleased(MouseEvent event) {
	
	readScrollDrag = false;
	
	if(getCursor().getType() == Cursor.HAND_CURSOR) {					
		setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
	}
	if(sidebar) {
		if(sampleDrag) {
			//System.out.println(selectedIndex +" " +hoverIndex);
			moveSample(selectedIndex, hoverIndex);
			
			updatevars = true;
			repaint();
		}
		
	}
	hoverIndex = -1;
	sampleDrag = false;
	if(zoomDrag) {
		
		zoom();		
		updatevars = true;
		Main.chromDraw.updateExons = true;
		Main.chromDraw.repaint();
	}
	if(resizeSidebar) {
		
		resizeCanvas(this.getWidth(), this.getHeight());
		Main.bedCanvas.repaint();
		Main.controlDraw.repaint();
		Main.upPanel.setDividerLocation(Main.sidebarWidth-4);	
		resizeSidebar = false;	
	}
	if(tempViewLength > Settings.readDrawDistance && selectedSplit.viewLength < Settings.readDrawDistance && drawVariables.sampleHeight > 100) {
		
		
		for(int i = drawVariables.visiblestart; i<= drawVariables.visiblestart+drawVariables.visiblesamples; i++) {
			
			if(i > Main.samples-1 ) {
				break;
			}
			if(sampleList.get(i).samFile == null) {
				continue;
			}
			if(sampleList.get(i).getreadHash().get(selectedSplit) != null) {
				sampleList.get(i).getreadHash().get(selectedSplit).setCoverageStart(Integer.MAX_VALUE);
				sampleList.get(i).getreadHash().get(selectedSplit).setCoverages(null);
				sampleList.get(i).getreadHash().get(selectedSplit).getReads().clear();		
				sampleList.get(i).getreadHash().get(selectedSplit).getHeadAndTail().clear();	
			}
		}
		
	}
	else if(tempViewLength < Settings.readDrawDistance && selectedSplit.viewLength > Settings.readDrawDistance && drawVariables.sampleHeight > 100) {
		for(int i = drawVariables.visiblestart; i<= drawVariables.visiblestart+drawVariables.visiblesamples; i++) {
			if(i > Main.samples-1) {
				break;
			}
			if(sampleList.get(i).samFile == null) {
				continue;
			}
			try {
				if(sampleList.get(i).getreadHash().get(selectedSplit) != null) {
					sampleList.get(i).getreadHash().get(selectedSplit).setCoverageStart(Integer.MAX_VALUE);	
				}
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
	}
	if(sampleZoomer) {
		if(splits.size() > 1) {
			
			for(int i = 0; i<splits.size(); i++) {
				splits.get(i).updateReads = true; 
			}
		}		
	}
	
	release();
	Main.widthLabel.setText("" +MethodLibrary.formatNumber((int)(splits.get(0).end-splits.get(0).start)) +"bp");	
	
	switch(event.getModifiers()) {
		case InputEvent.BUTTON3_MASK: {	
			
			
		}
	}
	if(Main.bedCanvas.bedTrack.size() > 0) {
		bedFeatureFetcher fetch = Main.bedCanvas.new bedFeatureFetcher();
		fetch.execute();
	}
	
	if(!loading && drawVariables.sampleHeight > 100) {
		updatevars = true;
		calculateVars = true;
		
		
	}
	sampleDrag = false;
	repaint();
}	
	
static void setScrollbar(final int value) {
	
	SwingUtilities.invokeLater(new Runnable() {
	    public void run() {
	    		
		    	updatevars = true;
		    	for(int i = 0; i<Main.drawCanvas.splits.size(); i++) {
		    		Main.drawCanvas.splits.get(i).updateReads = true;
		    	}
		    	
		        updateCoverages = true;
		    	Main.drawScroll.getVerticalScrollBar().setValue(value);	      
		    	Main.drawScroll.getVerticalScrollBar().revalidate();
		    	Main.drawCanvas.scrollbar = false;
		    	if(Main.drawCanvas.sampleList.size() > 0) {
		    		
			    	for(int i = Main.drawCanvas.drawVariables.visiblestart; i<Main.drawCanvas.drawVariables.visiblestart+Main.drawCanvas.drawVariables.visiblesamples; i++) {
			    		
			    		if(i > Main.drawCanvas.sampleList.size()-1) {
			    			break;
			    		}
			    		if(Main.drawCanvas.sampleList.get(i).removed) {
							continue;
						}
			    		Main.drawCanvas.sampleList.get(i).prepixel = -1;
					}
			    	
					updatevars = true;
		    	}
		    	if(Main.drawCanvas.drawVariables.visiblestart != (short)(0.5+Main.drawScroll.getVerticalScrollBar().getValue()/Main.drawCanvas.drawVariables.sampleHeight)) {
		    		Main.drawCanvas.drawVariables.visiblestart = (short)(0.5+Main.drawScroll.getVerticalScrollBar().getValue()/Main.drawCanvas.drawVariables.sampleHeight);
				}
		    	Main.drawCanvas.repaint();	
		    	
	    }
	});	
}
static void setGlasspane(final boolean value) {
	SwingUtilities.invokeLater(new Runnable() {
	    public void run() {	    	
	    	
	    	Main.glassPane.setVisible(value);   	 
	    	Main.glassPane.repaint();
	    	
	    	
	    	
	    }
	});
}

}
