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
import java.awt.Composite;
import java.awt.Cursor;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.Graphics2D;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.tribble.annotation.Strand;
import base.BBfile.BedFeature;
import base.BBfile.BigBedIterator;
import base.BBfile.BigWigIterator;
import base.BBfile.WigItem;
import base.BBfile.ZoomDataRecord;
import base.BBfile.ZoomLevelIterator;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.index.Block;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;

import javax.swing.JPanel;
import javax.swing.SwingWorker;

import org.apache.commons.io.FilenameUtils;

public class BedCanvas extends JPanel implements MouseMotionListener, MouseListener, KeyListener {

private static final long serialVersionUID = 1L;
	static int testfield = 0;
	Color forwardColor = new Color(171,194,171, 255), reverseColor = new Color(194,171,171, 255);
	BufferedImage bufImage, nodeImage;
	Graphics2D buf, nodebuf;
	int width, height, bedwidth;
	Color zoomColor = new Color(145,145,130,50);
	ArrayList<String> infolist = new ArrayList<String>();
	Composite backupComposite;
	ArrayList<BedTrack> bedTrack = new ArrayList<BedTrack>();
	TabixReaderMod.Iterator iterator = null;
	TabixReaderMod.Iterator bedIterator = null;
	ArrayList<Double> trackDivider = new ArrayList<Double>();
	String scaletext = "";
	int scalewidth = 0;
	ArrayList<Double> tempDivider;	
	BedTrack track;
	static int counter = 0;
	Rectangle remoBox = new Rectangle();
	private int bedEndPos;
	private boolean foundlevel;
	private int[][] matrix;
	private int sum;
	private int prevY;
	private Color bedcolor;
	private int bedcounter;
	private int trackstart;
	FontMetrics fm;
	private Rectangle2D textWidth;
	private int lettercount;
	private int mouseY;
	private int hoverIndex;
	boolean bedOn = false, intersected = false;
	private Rectangle mouseRect = new Rectangle();
	private Rectangle sideMouseRect = new Rectangle();
	private boolean overlapping = false;
	private int selectedPlay = -1;
	private Rectangle testRect = new Rectangle();
	BedNode infoNode = new BedNode("",0,0, null), preInfoNode = new BedNode("",0,0, null);
	private Rectangle boxmetrics = new Rectangle();
	private boolean overlap = false;
	private boolean sidebar = false;
	private int removeTrack = -1;
	//public boolean graph = false;
	
	private int bedheight=0;
	private int mouseX =0;
	private boolean zoomDrag=false;
	int overlapindex = 0;
	boolean lineZoomer = false;
	private String zoomtext= "";
	private int zoompostemp=0;
	private int pressY = 0;
	//int nodeHeight = 10;
	BedTrack updateTrack = null;
	private String[] colorsplit;
	private Color addColor;
	private int selexheight = 50;
	public boolean resize = false;
	private boolean resizer = false;
	private int resizeDivider;
	private Integer trackheight;
	private boolean mouseDrag = false;
	private int preresize = 0;
	private double preresizer = 0.0;
	private boolean negative = false;
	private boolean negativelock = false,positivelock = false;
	boolean heightchanged = false;
	private int selexstart;
	public boolean annotator = false;
	private boolean found;
	
	static BedTrack annoTrack = null;
	
	BedCanvas(int width, int height) {	
		
		this.width = width;
		this.height = height;
		bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)Main.screenSize.getWidth(), (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
		buf = (Graphics2D)bufImage.getGraphics();
	//	buf.setRenderingHints(Draw.rh);
		nodeImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)Main.screenSize.getWidth(), (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));
		nodebuf = (Graphics2D)nodeImage.getGraphics();
	//	nodebuf.setRenderingHints(Draw.rh);
		backupComposite = nodebuf.getComposite();
		buf.setFont(Draw.defaultFont);
		//nodebuf.setFont(Draw.defaultFont);
		addMouseListener(this);
		addKeyListener(this);
		addMouseMotionListener(this);
		addMouseWheelListener(
				new MouseWheelListener() {
					public void mouseWheelMoved(MouseWheelEvent mwe) {							
						bedTrack.get(hoverIndex).mouseWheel -=mwe.getWheelRotation();
						if(bedTrack.get(hoverIndex).mouseWheel > 0) {
							bedTrack.get(hoverIndex).mouseWheel = 0;
						}							
						repaint();							
					}					
				}
			);
	}	

void drawScreen(Graphics g) {
	
	if(Main.readingbeds) {
		buf.setColor(Draw.backColor);
		buf.fillRect(Main.sidebarWidth-4, 0, this.getWidth(), nodeImage.getHeight());	
		buf.drawString("Loading tracks...", 10, Main.bedScroll.getViewport().getHeight());
		g.drawImage(bufImage, 0, 0, null);	
		return;
	}
	if(this.trackDivider.size() > 0 && this.trackDivider.get(this.trackDivider.size()-1) != 1.0) {	
		 for(int i = 0 ; i<Main.bedCanvas.trackDivider.size(); i++) {
			 Main.bedCanvas.trackDivider.set(i, ((i+1)*(this.getHeight()/(double)trackDivider.size())/this.getHeight()));
		 }		
	}
	
	drawSidebar();
		
	
	//if(Settings.wallpaper == null) {
		buf.setColor(Draw.backColor);
		buf.fillRect(Main.sidebarWidth-4, 0, this.getWidth(), this.getHeight());	
		/*}
	else {
		
		buf.drawImage(Settings.wallpaper, Main.sidebarWidth-4,0, this);
		buf.setColor(Draw.backColor);	
		buf.fillRect(Main.sidebarWidth-4, 0,this.getWidth(), this.getHeight());
	}*/
	//buf.setColor(Color.gray);
	
	if(!zoomDrag && !resize) {	
		
		try {
			
			drawNodes();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}	
	
	
	
	
	if(resizer && !mouseDrag) {
		resizer = false;
	}
	if(resize) {		
		buf.drawImage(nodeImage, Main.sidebarWidth-4, 0, nodeImage.getWidth(),(int)(Main.vardivider*Main.varPaneDivider.getY()), null);
	}
	else {		
		buf.drawImage(nodeImage, Main.sidebarWidth-4, 0, null);
	}
	for(int i = 0 ; i<bedTrack.size(); i++) {
		if(i <bedTrack.size()-1) {
			buf.setColor(Color.lightGray);
			buf.drawLine(0, (int)(trackDivider.get(i)*this.getHeight()), this.getWidth(), (int)(trackDivider.get(i)*this.getHeight()));
			
			buf.setColor(Color.gray);
			buf.drawLine(0, (int)(trackDivider.get(i)*this.getHeight())+1, this.getWidth(), (int)(trackDivider.get(i)*this.getHeight())+1);
			
				if(!lineZoomer && mouseY < (int)(trackDivider.get(i)*this.getHeight())+4 && mouseY > (int)(trackDivider.get(i)*this.getHeight())-4) {
				resizer = true;
				if(getCursor().getType() != Cursor.N_RESIZE_CURSOR) {
					resizeDivider = i;
					setCursor(Cursor.getPredefinedCursor(Cursor.N_RESIZE_CURSOR));					
				}	
			}			
		}
		if(bedTrack.get(i).graph && bedTrack.get(i).minvalue != Double.MAX_VALUE && bedTrack.get(i).getHead().getNext() != null) {
			if(!buf.getColor().equals(Color.white)) {
				buf.setColor(Color.white);				
			}
			if(bedTrack.get(i).getLogscale().isSelected()) {
				scaletext = "Log scale [" +MethodLibrary.round(bedTrack.get(i).minvalue, 2) +", " +MethodLibrary.round(bedTrack.get(i).maxvalue,2) +"]";
				
				scalewidth = buf.getFontMetrics().stringWidth(scaletext);
				buf.fillRoundRect(Main.sidebarWidth+5,  (int)(trackDivider.get(i)*this.getHeight())-5-(Main.defaultFontSize+4), scalewidth+4, Main.defaultFontSize+4, 4, 4);
				buf.setColor(Color.black);	
				buf.drawString(scaletext, Main.sidebarWidth+7, (int)(trackDivider.get(i)*this.getHeight())-9);
				
			}
			else {
				scaletext="Scale [" +MethodLibrary.round(bedTrack.get(i).minvalue, 2) +", " +MethodLibrary.round(bedTrack.get(i).maxvalue,2) +"]";
				scalewidth = buf.getFontMetrics().stringWidth(scaletext);
				buf.fillRoundRect(Main.sidebarWidth+5,  (int)(trackDivider.get(i)*this.getHeight())-5-(Main.defaultFontSize+4), scalewidth+4, Main.defaultFontSize+4, 4, 4);
				buf.setColor(Color.black);	
				buf.drawString(scaletext, Main.sidebarWidth+7, (int)(trackDivider.get(i)*this.getHeight())-9);
			}
			buf.setColor(Color.black);
		}
	}
	if(overlap ) {
		
		drawInfo();
	}
	if(!resizer && !overlapping) {		
		if(getCursor().getType() != Cursor.DEFAULT_CURSOR) {			
			setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
		}
	}
	if(Main.drawCanvas.splits.get(0).pixel > 1) {
		// Middle line		
		buf.setColor(Color.black);
		buf.setStroke(Draw.dashed);			
		buf.drawLine((Main.drawCanvas.getDrawWidth())/2+Main.sidebarWidth-2, 0, ((Main.drawCanvas.getDrawWidth()))/2+Main.sidebarWidth-2, Main.bedScroll.getViewport().getHeight());
		buf.drawLine((int)((Main.drawCanvas.getDrawWidth())/2+Main.drawCanvas.splits.get(0).pixel+Main.sidebarWidth-2), 0, (int)(((Main.drawCanvas.getDrawWidth()))/2+Main.drawCanvas.splits.get(0).pixel+Main.sidebarWidth-2), Main.bedScroll.getViewport().getHeight());
//		buf.setStroke(Draw.doubleStroke);
		buf.setStroke(Draw.basicStroke);	
	}	
	for(int i = 1; i<Main.drawCanvas.splits.size(); i++) {
		buf.setColor(Color.gray);
		buf.fillRect(Main.drawCanvas.splits.get(i).offset-3, 0, 5, this.getHeight());
		buf.setColor(Color.lightGray);
		buf.fillRect(Main.drawCanvas.splits.get(i).offset-1, 0, 2, this.getHeight());
		
	}
	if(getCursor().getType() != Cursor.N_RESIZE_CURSOR) {
		
		drawZoom();
	}
	g.drawImage(bufImage, 0, 0, null);	
}	

void drawSidebar() {
	
	buf.setColor(Draw.sidecolor.darker());
	buf.fillRect(0, 0, Main.sidebarWidth-4, this.getHeight());
//	buf.setColor(Draw.softColor);
//	buf.fillRect(0, 0, Main.sidebarWidth, this.getHeight());
	buf.setColor(Color.black);
	
	if(bedTrack.size() > 0) {
	//	buf.setStroke(Draw.doubleStroke);
		overlapping = false;
		for(int i = 0; i<bedTrack.size(); i++) {
			if(!buf.getColor().equals(Color.black)) {
				buf.setColor(Color.black);
			}
			if(i == 0) {
				trackstart = Main.defaultFontSize;
				trackheight = (int)(trackDivider.get(i)*this.getHeight());
			}
			else {
				
				trackstart = Main.defaultFontSize+(int)(trackDivider.get(i-1)*this.getHeight());
				trackheight = (int)(trackDivider.get(i)*this.getHeight())-(int)(trackDivider.get(i-1)*this.getHeight());
			}
			
			if(bedTrack.get(i).file != null) {
				
				buf.drawString(bedTrack.get(i).file.getName(), 10, trackstart);
			}
			else {
				
				buf.drawString(FilenameUtils.getName(bedTrack.get(i).url.getFile()), 10, trackstart);
			}
			
			if(trackheight > Main.defaultFontSize+Main.defaultFontSize+(int)(Main.defaultFontSize*1.4)) {
				
				if((int)bedTrack.get(i).playbox.y != trackstart+Main.defaultFontSize) {
					
					bedTrack.get(i).playbox.setBounds(10, trackstart+Main.defaultFontSize, (int)(Main.defaultFontSize*1.4), (int)(Main.defaultFontSize*1.4));
					bedTrack.get(i).playTriangle.reset();					
					bedTrack.get(i).playTriangle.addPoint(bedTrack.get(i).playbox.x+(Main.defaultFontSize/5), bedTrack.get(i).playbox.y+(Main.defaultFontSize/5));
					bedTrack.get(i).playTriangle.addPoint(bedTrack.get(i).playbox.x+(Main.defaultFontSize/5), bedTrack.get(i).playbox.y+(bedTrack.get(i).playbox.width-(Main.defaultFontSize/5)));
					bedTrack.get(i).playTriangle.addPoint(bedTrack.get(i).playbox.x+(bedTrack.get(i).playbox.width-(Main.defaultFontSize/5)),bedTrack.get(i).playbox.y +bedTrack.get(i).playbox.height/2);						
					bedTrack.get(i).graphBox.setBounds(bedTrack.get(i).playbox.x+(int)(Main.defaultFontSize*1.4)+Main.defaultFontSize, trackstart+Main.defaultFontSize, (int)(Main.defaultFontSize*1.4), (int)(Main.defaultFontSize*1.4));
					if(bedTrack.get(i).settingsButton == null) {
						bedTrack.get(i).settingsButton = new Rectangle();
					}
					bedTrack.get(i).settingsButton.setBounds(Main.sidebarWidth-(int)(this.remoBox.width*(1.8))-4,trackstart-Main.defaultFontSize , (int)(this.remoBox.width*(1.5))+2, (int)(this.remoBox.height*(1.5)));
										
				}
				if(bedTrack.get(i).settingsButton == null) {
					bedTrack.get(i).settingsButton = new Rectangle();
				}
				if(bedTrack.get(i).settingsButton.y == 0 || bedTrack.get(i).settingsButton.x != Main.sidebarWidth-(int)(this.remoBox.width*(1.8)) || bedTrack.get(i).settingsButton.width !=  (int)(this.remoBox.width*(1.5))+2) {
					bedTrack.get(i).settingsButton.setBounds(Main.sidebarWidth-(int)(this.remoBox.width*(1.8)),trackstart-Main.defaultFontSize-4 , (int)(this.remoBox.width*(1.5))+2, (int)(this.remoBox.height*(1.5)));
					
				}
				if(Main.varsamples > 0 || FileRead.caller) {
					if(bedTrack.get(i).getCurrent() != null || (!bedTrack.get(i).small || bedTrack.get(i).getZoomlevel() != null)) {
						buf.setColor(Color.white);
						buf.fillRoundRect(bedTrack.get(i).playbox.getBounds().x-1,bedTrack.get(i).playbox.getBounds().y-1,bedTrack.get(i).playbox.getBounds().width, bedTrack.get(i).playbox.getBounds().height,2,2 );
						buf.setColor(Color.gray);
						buf.fillRoundRect(bedTrack.get(i).playbox.getBounds().x+1,bedTrack.get(i).playbox.getBounds().y+1,bedTrack.get(i).playbox.getBounds().width, bedTrack.get(i).playbox.getBounds().height,2,2 );
						if(sideMouseRect.intersects(bedTrack.get(i).playbox)) {
							overlapping = true;
							if(getCursor().getType() != Cursor.HAND_CURSOR) {
								selectedPlay = i;
								setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
							}	
							buf.setColor(Color.white);
						//	buf.drawRect(bedTrack.get(i).playbox.getBounds().x-1,bedTrack.get(i).playbox.getBounds().y-1,bedTrack.get(i).playbox.getBounds().width+1, bedTrack.get(i).playbox.getBounds().height+1);
							buf.fillRoundRect(bedTrack.get(i).playbox.getBounds().x,bedTrack.get(i).playbox.getBounds().y,bedTrack.get(i).playbox.getBounds().width, bedTrack.get(i).playbox.getBounds().height,2,2 );
							
						}
						else {
							buf.setColor(Draw.sidecolor);
							buf.fillRoundRect(bedTrack.get(i).playbox.getBounds().x,bedTrack.get(i).playbox.getBounds().y,bedTrack.get(i).playbox.getBounds().width, bedTrack.get(i).playbox.getBounds().height,2,2 );
						}
						if(bedTrack.get(i).intersect) {
							buf.setColor(Draw.greenColor);
							buf.fillRect(bedTrack.get(i).playTriangle.getBounds().x, bedTrack.get(i).playTriangle.getBounds().y, bedTrack.get(i).playTriangle.getBounds().width,bedTrack.get(i).playTriangle.getBounds().height);
						}
						else {
							buf.setColor(Draw.redColor);
							buf.fillPolygon(bedTrack.get(i).playTriangle);
						}
					}						
				}
				
				if(bedTrack.get(i).hasvalues) {
					
					if(bedTrack.get(i).getCurrent() != null || !bedTrack.get(i).small) {
					buf.setColor(Color.lightGray);
					buf.fillRect(bedTrack.get(i).graphBox.getBounds().x,bedTrack.get(i).graphBox.getBounds().y,bedTrack.get(i).graphBox.getBounds().width,bedTrack.get(i).graphBox.getBounds().height );
					buf.setColor(Color.white);
					buf.fillRoundRect(bedTrack.get(i).graphBox.getBounds().x-1,bedTrack.get(i).graphBox.getBounds().y-1,bedTrack.get(i).graphBox.getBounds().width, bedTrack.get(i).graphBox.getBounds().height,2,2 );
					buf.setColor(Color.gray);
					buf.fillRoundRect(bedTrack.get(i).graphBox.getBounds().x+1,bedTrack.get(i).graphBox.getBounds().y+1,bedTrack.get(i).graphBox.getBounds().width, bedTrack.get(i).graphBox.getBounds().height,2,2 );
					if(sideMouseRect.intersects(bedTrack.get(i).graphBox)) {
						overlapping = true;
						if(getCursor().getType() != Cursor.HAND_CURSOR) {
							selectedPlay = i;
							setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
						}	
						buf.setColor(Color.white);
					//	buf.drawRect(bedTrack.get(i).graphBox.getBounds().x-1,bedTrack.get(i).graphBox.getBounds().y-1,bedTrack.get(i).graphBox.getBounds().width+1, bedTrack.get(i).graphBox.getBounds().height+1);
						buf.fillRoundRect(bedTrack.get(i).graphBox.getBounds().x,bedTrack.get(i).graphBox.getBounds().y,bedTrack.get(i).graphBox.getBounds().width, bedTrack.get(i).graphBox.getBounds().height,2,2 );
						
					}
					else {
						buf.setColor(Draw.sidecolor);
						buf.fillRoundRect(bedTrack.get(i).graphBox.getBounds().x,bedTrack.get(i).graphBox.getBounds().y,bedTrack.get(i).graphBox.getBounds().width, bedTrack.get(i).graphBox.getBounds().height,2,2 );
					}
					if(bedTrack.get(i).graph) {
						buf.setColor(Draw.greenColor);
						buf.drawLine(bedTrack.get(i).graphBox.getBounds().x, bedTrack.get(i).graphBox.getBounds().y+bedTrack.get(i).graphBox.getBounds().height, bedTrack.get(i).graphBox.getBounds().x +bedTrack.get(i).graphBox.getBounds().width/4, bedTrack.get(i).graphBox.getBounds().y+bedTrack.get(i).graphBox.getBounds().height/2);
						buf.drawLine(bedTrack.get(i).graphBox.getBounds().x+bedTrack.get(i).graphBox.getBounds().width/4, bedTrack.get(i).graphBox.getBounds().y+bedTrack.get(i).graphBox.getBounds().height/2, bedTrack.get(i).graphBox.getBounds().x+bedTrack.get(i).graphBox.getBounds().width/2, bedTrack.get(i).graphBox.getBounds().y+(int)(bedTrack.get(i).graphBox.getBounds().height*0.66));
						buf.drawLine(bedTrack.get(i).graphBox.getBounds().x+bedTrack.get(i).graphBox.getBounds().width/2, bedTrack.get(i).graphBox.getBounds().y+(int)(bedTrack.get(i).graphBox.getBounds().height*0.66), (int)(bedTrack.get(i).graphBox.getBounds().x+bedTrack.get(i).graphBox.getBounds().width*0.66), bedTrack.get(i).graphBox.getBounds().y+bedTrack.get(i).graphBox.getBounds().height/4);
						buf.drawLine((int)(bedTrack.get(i).graphBox.getBounds().x+bedTrack.get(i).graphBox.getBounds().width*0.66), bedTrack.get(i).graphBox.getBounds().y+bedTrack.get(i).graphBox.getBounds().height/4, bedTrack.get(i).graphBox.getBounds().x+bedTrack.get(i).graphBox.getBounds().width, bedTrack.get(i).graphBox.getBounds().y+bedTrack.get(i).graphBox.getBounds().height/2);
						
					}
					else {
						buf.setColor(Draw.redColor);
						buf.drawLine(bedTrack.get(i).graphBox.getBounds().x, bedTrack.get(i).graphBox.getBounds().y+bedTrack.get(i).graphBox.getBounds().height, bedTrack.get(i).graphBox.getBounds().x +bedTrack.get(i).graphBox.getBounds().width/4, bedTrack.get(i).graphBox.getBounds().y+bedTrack.get(i).graphBox.getBounds().height/2);
						buf.drawLine(bedTrack.get(i).graphBox.getBounds().x+bedTrack.get(i).graphBox.getBounds().width/4, bedTrack.get(i).graphBox.getBounds().y+bedTrack.get(i).graphBox.getBounds().height/2, bedTrack.get(i).graphBox.getBounds().x+bedTrack.get(i).graphBox.getBounds().width/2, bedTrack.get(i).graphBox.getBounds().y+(int)(bedTrack.get(i).graphBox.getBounds().height*0.66));
						buf.drawLine(bedTrack.get(i).graphBox.getBounds().x+bedTrack.get(i).graphBox.getBounds().width/2, bedTrack.get(i).graphBox.getBounds().y+(int)(bedTrack.get(i).graphBox.getBounds().height*0.66), (int)(bedTrack.get(i).graphBox.getBounds().x+bedTrack.get(i).graphBox.getBounds().width*0.66), bedTrack.get(i).graphBox.getBounds().y+bedTrack.get(i).graphBox.getBounds().height/4);
						buf.drawLine((int)(bedTrack.get(i).graphBox.getBounds().x+bedTrack.get(i).graphBox.getBounds().width*0.66), bedTrack.get(i).graphBox.getBounds().y+bedTrack.get(i).graphBox.getBounds().height/4, bedTrack.get(i).graphBox.getBounds().x+bedTrack.get(i).graphBox.getBounds().width, bedTrack.get(i).graphBox.getBounds().y+bedTrack.get(i).graphBox.getBounds().height/2);
						
					
					}
					/*if(sideMouseRect.intersects(bedTrack.get(i).graphBox)) {
						
						overlapping = true;
						selectedPlay = i;
						if(getCursor().getType() != Cursor.HAND_CURSOR) {
							
							setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
						}	
						buf.setColor(Color.white);
						buf.drawRect(bedTrack.get(i).graphBox.getBounds().x-1,bedTrack.get(i).graphBox.getBounds().y-1,bedTrack.get(i).graphBox.getBounds().width+1, bedTrack.get(i).graphBox.getBounds().height+1);
						
					}*/
					}
				}				
			}
			if(trackheight > bedTrack.get(i).settingsButton.height) {
				/*if(this.sideMouseRect.intersects(bedTrack.get(i).settingsButton)) {								
					buf.setColor(Color.white);
				}
				else {			
					
				}*/
				buf.setColor(Draw.sidecolor.darker());
				buf.fillRect(bedTrack.get(i).settingsButton.x, bedTrack.get(i).settingsButton.y,bedTrack.get(i).settingsButton.width,bedTrack.get(i).settingsButton.height);
				buf.drawImage(Main.settingsIcon.getImage(),Main.sidebarWidth-(int)(this.remoBox.width*(1.8))-4,trackstart-Main.defaultFontSize , (int)(this.remoBox.width*(1.5)), (int)(this.remoBox.height*(1.5)), this);
			}
		}
		
		if(!overlapping && !resizer) {			
			
			if(getCursor().getType() != Cursor.DEFAULT_CURSOR) {
				
				selectedPlay  = -1;
				setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
			}			
		}
		if(sidebar) {
			buf.setColor(Color.black);
			
		
			if(sidebar && sideMouseRect.intersects(this.remoBox)) {				
			//	buf.setStroke(Draw.doubleStroke);
				removeTrack  = hoverIndex;
				buf.setColor(Color.white);
				buf.fillRect(this.remoBox.x-2, this.remoBox.y-1, this.remoBox.width+2, this.remoBox.height+2);
			}
			else {				
				removeTrack = -1;
				if(this.remoBox.getBounds().x != Main.sidebarWidth-(Main.defaultFontSize+10)-4 || this.remoBox.getBounds().y != (int)(trackDivider.get(hoverIndex)*this.getHeight()) -(Main.defaultFontSize+6)) {
					this.remoBox.setBounds(Main.sidebarWidth-(Main.defaultFontSize+10)-4, (int)(trackDivider.get(hoverIndex)*this.getHeight()) -(Main.defaultFontSize+6), Main.defaultFontSize+3, Main.defaultFontSize+3);
				}
			/*	if(this.remoBox.getBounds().y != (int)(trackDivider.get(hoverIndex)*this.getHeight())-11|| this.remoBox.getBounds().x != Main.sidebarWidth-11) {
					this.remoBox.setBounds(Main.sidebarWidth-11, (int)(trackDivider.get(hoverIndex)*this.getHeight())-11, 8, 8);
				}*/
			}
		
			buf.setColor(Color.black);
			buf.drawRect(this.remoBox.x, this.remoBox.y, this.remoBox.width, this.remoBox.height);
			buf.drawLine(this.remoBox.x, this.remoBox.y, this.remoBox.x+this.remoBox.width, this.remoBox.y+(int)this.remoBox.getHeight());
			buf.drawLine(this.remoBox.x, this.remoBox.y+(int)this.remoBox.getHeight(), this.remoBox.x+this.remoBox.width,this.remoBox.y);
		
		}
		buf.setStroke(Draw.doubleStroke);
		buf.setColor(Color.gray);
		buf.drawLine(Main.sidebarWidth-5, 0, Main.sidebarWidth-5, this.getHeight());
		buf.drawLine(1, 0, 1, this.getHeight());
		buf.setColor(Color.lightGray);
		buf.drawLine(3, 0,3, this.getHeight());
		buf.drawLine(Main.sidebarWidth-7, 0, Main.sidebarWidth-7, this.getHeight());
		buf.setStroke(Draw.basicStroke);
	}
	
	
}

void drawZoom() {	
	
	if(lineZoomer ) {		
		buf.setColor(Color.black);
		buf.setStroke(Draw.dashed);		
		buf.drawLine(Main.drawCanvas.pressX, pressY, mouseX, mouseY);	
		
	}
	else if(zoomDrag) {

		buf.setColor(Color.white);

		if(this.mouseX-Main.drawCanvas.pressX >= 0) {
			buf.drawRect(Main.drawCanvas.pressX+1, 0, this.mouseX-Main.drawCanvas.pressX-2, this.getHeight());	
			if(Main.drawCanvas.getDrawWidth()-this.mouseX > 200) {
				buf.drawString("" +MethodLibrary.formatNumber((int)((this.mouseX-Main.drawCanvas.pressX)/Main.drawCanvas.selectedSplit.pixel)) +"bp", Main.drawCanvas.pressX+(this.mouseX-Main.drawCanvas.pressX)+4, this.mouseY-35);
				buf.drawString("Right click to cancel zoom" , Main.drawCanvas.pressX+(this.mouseX-Main.drawCanvas.pressX)+4, this.mouseY-6);
			}
			else {
				fm = buf.getFontMetrics();
				zoomtext = ""+MethodLibrary.formatNumber((int)((this.mouseX-Main.drawCanvas.pressX)/Main.drawCanvas.selectedSplit.pixel)) +"bp";
				textWidth = fm.getStringBounds(zoomtext, buf);
				this.zoompostemp = (int)((Main.drawCanvas.pressX+(this.mouseX-Main.drawCanvas.pressX))-textWidth.getWidth());
				buf.drawString(zoomtext, this.zoompostemp, this.mouseY-35);
				
				textWidth = fm.getStringBounds("Right click to cancel zoom", buf);
				this.zoompostemp = (int)((Main.drawCanvas.pressX+(this.mouseX-Main.drawCanvas.pressX))-textWidth.getWidth())-6;
				buf.drawString("Right click to cancel zoom" , this.zoompostemp, this.mouseY-6);
	
			}
			buf.setColor(Draw.zoomColor);
			buf.fillRect(Main.drawCanvas.pressX, 0, this.mouseX-Main.drawCanvas.pressX, this.getHeight());		
		}
		else {			
			lineZoomer = true;
			Main.drawCanvas.lineZoomer = true;
			zoomDrag = false;		
		}	
	}
}

void drawLogo(BedTrack track, BedNode node, int trackstart) {	
	try {
		
		
	if(node.id != null && Main.SELEXhash.containsKey(node.id)) {	
		matrix = Main.SELEXhash.get(node.id);
		
		if((trackstart +this.trackheight) - (trackstart+15+((node.level-1)*(selexheight+15))+(track.mouseWheel*((selexheight+15)*2))) < selexheight) {
			selexheight = (trackstart +this.trackheight) - (trackstart+15+((node.level-1)*(selexheight+15))+(track.mouseWheel*((selexheight+15)*2)));
		}		
		selexstart = (int)((node.getDrawPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel);
		
		if(matrix != null && (trackstart +this.trackheight) > trackstart+10+((node.level-1)*(40+15))+(track.mouseWheel*((40+15)*2))) {
			
			sum = 0;
			if(node.forward == null || node.forward) {
				for(int j = 0; j<matrix[0].length; j++) {
					prevY = trackstart+10+((node.level-1)*(40+15))+(track.mouseWheel*((40+15)*2));
					sum = matrix[0][j] + matrix[1][j] +matrix[2][j] +matrix[3][j];
					nodebuf.drawImage(Main.A, (int)(selexstart+(j*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[0][j]/(double)sum)), null);
					prevY += (int)(selexheight*(matrix[0][j]/(double)sum));
					nodebuf.drawImage(Main.C, (int)(selexstart+(j*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[1][j]/(double)sum)), null);
					prevY += (int)(selexheight*(matrix[1][j]/(double)sum));
					nodebuf.drawImage(Main.G, (int)(selexstart+(j*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[2][j]/(double)sum)), null);
					prevY += (int)(selexheight*(matrix[2][j]/(double)sum));
					nodebuf.drawImage(Main.T, (int)(selexstart+(j*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[3][j]/(double)sum)), null);
				}
			}
			else {
				bedcounter = 0;
				for(int j = matrix[0].length-1; j>=0; j--) {
					prevY = trackstart+10+((node.level-1)*(40+15))+(track.mouseWheel*((40+15)*2));
					sum = matrix[0][j] + matrix[1][j] +matrix[2][j] +matrix[3][j];
				//	System.out.print((100*(matrix[0][j]/(double)sum)) +"\t" +(100*(matrix[1][j]/(double)sum)) +"\t" +(100*(matrix[2][j]/(double)sum)) +"\t" +(100*(matrix[3][j]/(double)sum)) +"\n");
					nodebuf.drawImage(Main.A, (int)(selexstart+(bedcounter*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[3][j]/(double)sum)), null);
					prevY += (int)(selexheight*(matrix[3][j]/(double)sum));
					nodebuf.drawImage(Main.C, (int)(selexstart+(bedcounter*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[2][j]/(double)sum)), null);
					prevY += (int)(selexheight*(matrix[2][j]/(double)sum));
					nodebuf.drawImage(Main.G, (int)(selexstart+(bedcounter*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[1][j]/(double)sum)), null);
					prevY += (int)(selexheight*(matrix[1][j]/(double)sum));
					nodebuf.drawImage(Main.T, (int)(selexstart+(bedcounter*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[0][j]/(double)sum)), null);
					bedcounter++;
				}
			}
		}
		nodebuf.setColor(Color.black);
		nodebuf.fillRect(selexstart, trackstart+10+((node.level-1)*(40+15))+(track.mouseWheel*((40+15)*2)),(int)(matrix[0].length*Main.drawCanvas.splits.get(0).pixel), 1);
		
	}	
	if(selexheight != 40) {
		selexheight = 40;
	}
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}

void drawNodes() {
	try {
		
	nodebuf.setComposite( Main.drawCanvas.composite);	
	nodebuf.fillRect(0,0, bufImage.getWidth(),this.getHeight());	
	nodebuf.setComposite(this.backupComposite);	
	overlap = false;
	for(int i = 0; i<bedTrack.size(); i++) {
		
		track = bedTrack.get(i);
		
		if(updateTrack != null) {
			if(!updateTrack.equals(track)) {
				continue;
			}
			else {
				updateTrack = null;
			}
		}
		
		track.prepixel = -1;
		if(i ==0) {
			trackstart = 5;
			trackheight = (int)(trackDivider.get(i)*this.getHeight());
		}
		else {
			trackstart = 5+(int)(trackDivider.get(i-1)*this.getHeight());
			trackheight = (int)(trackDivider.get(i)*this.getHeight())-(int)(trackDivider.get(i-1)*this.getHeight());
		}
		if(trackheight < 22) {			
			continue;
		}
		
		if(!annotator) {
			if(Main.drawCanvas.splits.get(0).viewLength < Settings.windowSize && !track.small) {
				
			}
			else if(!track.small) {
				
				if(!track.cleared) {
				
					removeBeds(track);
				}
				nodebuf.setColor(zoomColor);
				if(trackheight < 70 ) {
					nodebuf.setFont(new Font("SansSerif", Font.BOLD, trackheight-20));
				}
				else {
					nodebuf.setFont(new Font("SansSerif", Font.BOLD, 50));
				}
				nodebuf.drawString("ZOOM IN < " +MethodLibrary.formatNumber(Settings.windowSize) +" bp", 100, trackstart+(trackheight-22));
				nodebuf.setFont(Draw.defaultFont);
				continue;
			}
		}
		
		
		
		if(track.getHead().getNext() == null) {
			
			continue;
		}
	
		try {
			if(track.getCurrent().getPrev()!=null &&!track.getCurrent().getPrev().equals(track.getHead()) && track.getCurrent().getPosition() >= Main.drawCanvas.splits.get(0).start-1) {
				
				while(track.getCurrent().getPosition() >= Main.drawCanvas.splits.get(0).start-1) {
					if(!track.getCurrent().equals(track.getHead())) {
						track.setCurrent(track.getCurrent().getPrev());
					}
					else {
						break;
					}
				}
			}
			else {
				try {
				
				while(track.getCurrent() != null &&(track.getCurrent().getPosition()+track.getCurrent().getLength() < Main.drawCanvas.splits.get(0).start-1 && (track.getCurrent().getNext() != null && track.getCurrent().getNext().getPosition()+track.getCurrent().getLength() < Main.drawCanvas.splits.get(0).start ))) {
					
					if(track.getCurrent().getNext() != null) {
						
						track.setCurrent(track.getCurrent().getNext());
					}
					else {
						break;
					}
				}
				}
				catch(Exception e) {
					System.out.println(track.getCurrent().getPosition());
				}
			}
			
		}
		catch(Exception ex) {
			ex.printStackTrace();
			track.setCurrent(track.getHead().getNext());
			
		}

	
	//	if(!track.getDrawNode().equals(track.getCurrent())) {
			track.setDrawNode(track.getCurrent());
	//	}
	
			
		while(track.getDrawNode() != null && track.getDrawNode().getPosition() < Main.drawCanvas.splits.get(0).end) {
			
			if(Main.drawCanvas.splits.get(0).viewLength > 100000 ) {
				if(!this.heightchanged) {
					if(track.nodeHeight != (short)Main.defaultFontSize/3) {
						track.nodeHeight = (short)(Main.defaultFontSize/3);
					}
				}
			}
			else {
				if(!this.heightchanged) {
					if(track.nodeHeight != (short)(Main.defaultFontSize)) {
						track.nodeHeight = (short)(Main.defaultFontSize);
					}
				}
			}
			
			bedwidth = (int)((track.getDrawNode().getLength())*Main.drawCanvas.splits.get(0).pixel);
			if(bedwidth < 1) {
				bedwidth = 1;
			}
			
			if(track.getCollapseBox().isSelected() && (bedwidth == 1 && (int)((track.getDrawNode().getDrawPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel) == track.prepixel)) {
				
				track.setDrawNode(track.getDrawNode().getNext());
				continue;
			}
			
			if(!track.getColors().isEmpty()) {
				
				bedcolor = track.getColors().get(track.getDrawNode().color);
			}
			else if(track.getDrawNode().forward == null) {
				bedcolor = forwardColor;
			}
			else if(track.getDrawNode().forward) {
				bedcolor = forwardColor;				
			}
			else {
				bedcolor= reverseColor;
			}
			nodebuf.setColor(bedcolor);
			
			if(track.graph) {
				if(!track.negatives) {
					if(track.getDrawNode().value == null || Double.isNaN(track.getDrawNode().value)) {
						track.prepixel = (int)((track.getDrawNode().getDrawPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel);
						testRect.setBounds(track.prepixel,trackstart,bedwidth,track.nodeHeight);						
						nodebuf.fillRect(track.prepixel, trackstart,bedwidth,track.nodeHeight);	
					}
					else {
						if(track.getDrawNode().value == null) {
							track.getDrawNode().value = 0.0;
						}
						track.prepixel = (int)((track.getDrawNode().getDrawPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel);
						if(track.getLogscale().isSelected()) {
							bedheight = (int)((Math.log(track.getDrawNode().value*this.trackheight)/(double)Math.log((track.scale*this.trackheight)))*this.trackheight-10);
							
						}
						else {
							try {
								if(track.getDrawNode().value != null) {
									bedheight = (int)((track.getDrawNode().value/(track.scale))*this.trackheight-10);
								}
							}
							catch(Exception e) {
								e.printStackTrace();
								
							}
						}
					
						testRect.setBounds(track.prepixel,(int)(trackDivider.get(i)*this.getHeight())-bedheight,bedwidth,bedheight);
						nodebuf.fillRect(track.prepixel,testRect.y,bedwidth,bedheight);	
						
					}
				}
				else if(track.getDrawNode().value!=null) {
					track.prepixel = (int)((track.getDrawNode().getDrawPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel);
					if(track.getDrawNode().value != null && !Double.isNaN(track.getDrawNode().value) && track.getDrawNode().value >= 0) {
						if(track.getLogscale().isSelected()) {
							bedheight = (int)((Math.log(track.getDrawNode().value*this.trackheight)/(Math.log((track.scale*2*this.trackheight))))*this.trackheight-10);
							
						}
						else {
							bedheight = (int)((track.getDrawNode().value/(track.scale*2))*this.trackheight-10);
						}
						if( track.getDrawNode().value == 0) {
							bedheight = 1;
							if(bedcolor != Color.gray) {
								nodebuf.setColor(Color.gray);
							}
						}
						else if(bedheight < 1) {
							bedheight = 1;
						}
							testRect.setBounds(track.prepixel,(int)((int)(trackDivider.get(i)*this.getHeight())-(bedheight+this.trackheight/2.0)),bedwidth,bedheight);
							nodebuf.fillRect(track.prepixel, testRect.y,bedwidth,bedheight);
						}
					else {
						if(bedcolor != reverseColor) {
							nodebuf.setColor(reverseColor);
						}
						
						if(track.getLogscale().isSelected()) {
							bedheight = (int)((Math.log(Math.abs(track.getDrawNode().value)*this.trackheight)/Math.log((track.scale*2)*this.trackheight))*this.trackheight-10);
						}
						else {
							bedheight = (int)((Math.abs(track.getDrawNode().value)/(track.scale*2))*this.trackheight-10);
						}
						if(track.getDrawNode().value == 0) {
							bedheight = 1;
							if(bedcolor != Color.gray) {
								nodebuf.setColor(Color.gray);
							}
						}
						else if(bedheight < 1) {
							bedheight = 1;
						}
							testRect.setBounds(track.prepixel, (int)((int)(trackDivider.get(i)*this.getHeight())-(this.trackheight/2.0)),bedwidth,bedheight);							
							nodebuf.fillRect(track.prepixel, testRect.y,bedwidth,bedheight);
						}
				}
				if(!sidebar && mouseRect.intersects(testRect)) {
					
					nodebuf.setColor(Color.white);
					nodebuf.fillRect(testRect.getBounds().x, testRect.getBounds().y,(int)testRect.getWidth(), (int)testRect.getHeight());
					if(infoNode != track.getDrawNode()) {
						infoNode = track.getDrawNode();
					}
					overlapindex = trackstart;
					overlap = true;					
				}				
			}
			else {
		/*		if(Main.drawCanvas.splits.get(0).viewLength > 100000) {
					if(track.prepixel == (int)((drawNode.getPosition()-Main.drawCanvas.selectedSplit.start)*Main.drawCanvas.selectedSplit.pixel)) {
						drawNode = drawNode.getNext();
						continue;
					}				
					
					track.prepixel = (int)((drawNode.getPosition()-Main.drawCanvas.selectedSplit.start)*Main.drawCanvas.selectedSplit.pixel);
					nodebuf.fillRect(track.prepixel, trackstart,bedwidth,10);				
					track.mouseWheel =0;
				
					if(hoverIndex < bedTrack.size() && track.equals(bedTrack.get(hoverIndex)) && mouseY < trackstart+10) {
			
						testRect.setBounds((int)((drawNode.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel), trackstart,bedwidth, 10);
						if(testRect.intersects(mouseRect)) {
							nodebuf.setColor(Color.white);
							nodebuf.fillRect(testRect.getBounds().x, testRect.getBounds().y,bedwidth, 10);	
							if(infoNode != drawNode) {
								infoNode = drawNode;
							}
							overlap = true;					
							
						}				
					}
					
				}
				else {
						*/
				
						track.prepixel = (int)((track.getDrawNode().getDrawPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel);
						bedEndPos = track.getDrawNode().getPosition()+track.getDrawNode().getLength();
						
						if(track.getBedLevels().isEmpty()) {
						
							track.getBedLevels().add(bedEndPos);
							track.getDrawNode().level = 1;
					//		level = 1;
						}
						else {
							
							foundlevel = false;
							for(int c=0;c<track.getBedLevels().size(); c++) {
								
								if(track.getBedLevels().get(c) <= track.getDrawNode().getPosition()) {
									
								//	level = c+1;
									track.getDrawNode().level = c+1;
									foundlevel = true;				
									track.getBedLevels().set(c, bedEndPos);
									break;
								}
							
							}
								
						if(!foundlevel) {			
							track.getBedLevels().add(bedEndPos);
							//	level = track.bedLevelMatrix.size();
							track.getDrawNode().level = track.getBedLevels().size();
							}
						}					
						
					if(Main.drawCanvas.splits.get(0).viewLength <= 300 && (track.selex || track.getDrawNode().id != null)) {
						
						if((track.getDrawNode().level-1)*(selexheight+15)+(track.mouseWheel*((selexheight+15)*2))+5 >= 5 && (track.getDrawNode().level-1)*(selexheight+15)+(track.mouseWheel*((selexheight+15)*2))+15 < this.trackheight) {
							testRect.setBounds((int)((track.getDrawNode().getDrawPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel),trackstart+((track.getDrawNode().level-1)*(selexheight+15))+(track.mouseWheel*((selexheight+15)*2)),bedwidth,10);
							if(!sidebar && testRect.intersects(mouseRect)) {
								nodebuf.setColor(Color.white);
								nodebuf.fillRect(testRect.x,testRect.y,bedwidth,10);
								if(infoNode != track.getDrawNode()) {
									infoNode = track.getDrawNode();
								}
								overlapindex = trackstart;
								overlap = true;
							}
							else {
								nodebuf.fillRect(testRect.x,testRect.y,bedwidth,10);
							}							
							
							if(track.getDrawNode().name != null) {
								fm = nodebuf.getFontMetrics();
								textWidth = fm.getStringBounds(track.getDrawNode().name, nodebuf);							
								nodebuf.setColor(Color.black);
								
								if(bedwidth > textWidth.getWidth()) {
									nodebuf.drawString(track.getDrawNode().name, (int)((track.getDrawNode().getDrawPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel),trackstart+(track.getDrawNode().level*(selexheight+15))-(selexheight+7)+(track.mouseWheel*((selexheight+15)*2))+1);
								}
								else {
									lettercount = (int)(bedwidth/fm.getStringBounds("R", nodebuf).getWidth());
									nodebuf.drawString(track.getDrawNode().name.substring(0, lettercount), (int)((track.getDrawNode().getDrawPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel),trackstart+(track.getDrawNode().level*(selexheight+15))-(selexheight+7)+(track.mouseWheel*((selexheight+15)*2))+1);
								}							
							}							
							drawLogo(track,track.getDrawNode(), trackstart);							
						}
					}
					else {
						if(((track.getDrawNode().level-1)*(track.nodeHeight+(track.nodeHeight/2))+(track.mouseWheel*75)+(track.nodeHeight*2) > this.trackheight)) {
							track.setDrawNode(track.getDrawNode().getNext());
							
							continue;
						}
						
						if((track.getDrawNode().level-1)*(track.nodeHeight+(track.nodeHeight/2))+(track.mouseWheel*75)+5 >= 5) {
							testRect.setBounds((int)((track.getDrawNode().getDrawPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel), trackstart+((track.getDrawNode().level-1)*(track.nodeHeight+(track.nodeHeight/2)))+(track.mouseWheel*75),bedwidth, track.nodeHeight);
							if(!sidebar && testRect.intersects(mouseRect)) {
								nodebuf.setColor(Color.white);
								nodebuf.fillRect(testRect.x, testRect.y,bedwidth, track.nodeHeight);	
								if(infoNode != track.getDrawNode()) {
									infoNode = track.getDrawNode();
								}
								overlapindex = trackstart;
								overlap = true;								
							}
							else {
								nodebuf.fillRect(testRect.x, testRect.y,bedwidth, track.nodeHeight);
							}
						
			//				nodebuf.fillRect((int)((drawNode.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel), trackstart+((drawNode.level-1)*15)+(track.mouseWheel*75),bedwidth, 10);
							
						if(track.getDrawNode().name != null && testRect.width > 10 && track.nodeHeight > 8 ) {
							
							fm = nodebuf.getFontMetrics();
							textWidth = fm.getStringBounds(track.getDrawNode().name, nodebuf);							
							nodebuf.setColor(Color.black);
							if( testRect.x < 0 && testRect.x+testRect.width > 0) {
								if(bedwidth > textWidth.getWidth()) {
									nodebuf.drawString(track.getDrawNode().name, 0,testRect.getBounds().y+(track.nodeHeight)-1);
								}
								else {
									lettercount = (int)(bedwidth/fm.getStringBounds("R", nodebuf).getWidth());
									nodebuf.drawString(track.getDrawNode().name.substring(0, lettercount),0,testRect.getBounds().y+(track.nodeHeight)-1);
									
								}
							}
							else {
								try {
									if(bedwidth > textWidth.getWidth()) {
										nodebuf.drawString(track.getDrawNode().name, (int)((track.getDrawNode().getDrawPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel),testRect.getBounds().y+(track.nodeHeight)-1);
									}
									else {
										lettercount = (int)(bedwidth/fm.getStringBounds("R", nodebuf).getWidth());
										nodebuf.drawString(track.getDrawNode().name.substring(0, lettercount), (int)((track.getDrawNode().getDrawPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel),testRect.getBounds().y+(track.nodeHeight)-1);
										
									}
								}
								catch(Exception e) {
									e.printStackTrace();
								}
							}
						}
					}
				}
	//		}
			}
			
			//System.out.println(drawNode.getNext().getPosition());	
			track.setDrawNode(track.getDrawNode().getNext());
		}
		
		
		if(Main.drawCanvas.splits.size() > 1) {
			nodebuf.setComposite( Main.drawCanvas.composite);		
			nodebuf.fillRect(Main.drawCanvas.splits.get(1).chromOffset,0, bufImage.getWidth(),this.getHeight());	
			nodebuf.setComposite(this.backupComposite);	
		}
		//drawNode = null;
		/*if(Main.drawCanvas.splits.get(0).viewLength <= 300 && track.selex) {
			this.setPreferredSize(new Dimension(this.getWidth(), ((track.bedLevelMatrix.size()+1) * 65)));
			if(nodeImage.getHeight() < ((track.bedLevelMatrix.size()+1) * 65)) {
				bufImage = new BufferedImage((int)Main.screenSize.getWidth(),((track.bedLevelMatrix.size()+1) * 65), BufferedImage.TYPE_INT_ARGB);
				buf = (Graphics2D)bufImage.getGraphics();
				nodeImage = new BufferedImage((int)Main.screenSize.getWidth(),((track.bedLevelMatrix.size()+1) * 65), BufferedImage.TYPE_INT_ARGB);	
				nodebuf = (Graphics2D)nodeImage.getGraphics();
				buf.setFont(Draw.defaultFont);
				nodebuf.setFont(Draw.defaultFont);
			}
		}
		else {
			this.setPreferredSize(new Dimension(this.getWidth(), ((track.bedLevelMatrix.size()+1) * 15)));	
			if(nodeImage.getHeight() < ((track.bedLevelMatrix.size()+1) * 15)) {
				bufImage = new BufferedImage((int)Main.screenSize.getWidth(),((track.bedLevelMatrix.size()+1) * 15), BufferedImage.TYPE_INT_ARGB);
				buf = (Graphics2D)bufImage.getGraphics();
				nodeImage = new BufferedImage((int)Main.screenSize.getWidth(),((track.bedLevelMatrix.size()+1) * 15), BufferedImage.TYPE_INT_ARGB);	
				nodebuf = (Graphics2D)nodeImage.getGraphics();
				buf.setFont(Draw.defaultFont);
				nodebuf.setFont(Draw.defaultFont);
			}
		}		
		
		this.revalidate();	*/	
		track.getBedLevels().clear();
		
		
	}	
	
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}

void drawInfo() {
	
//	if(boxmetrics.width != 0) {
	
		if(!infoNode.equals(preInfoNode)) {	
			
			preInfoNode = infoNode;
			infolist.clear();
			int maxlength = 0;
			
			if(infoNode.name != null && infoNode.name.length() > 0) {
				
				infolist.add(infoNode.name);
				if(infoNode.name.length() > infolist.get(maxlength).length()) {
					maxlength = infolist.size()-1;
				}
				if(infoNode.secondaryName != null && infoNode.secondaryName.length() > 0) {
					infolist.add(infoNode.secondaryName);
					if(infoNode.secondaryName.length() > infolist.get(maxlength).length()) {
						maxlength = infolist.size()-1;
					}
				}
			}
			
			else {
				if(infoNode.getTrack().file != null) {
					infolist.add(infoNode.getTrack().file.getName().substring(0, 10) +"...");
					if(13 > infolist.get(maxlength).length()) {
						maxlength = infolist.size()-1;
					}
				}
			}
			infolist.add("Position: " +Main.drawCanvas.selectedSplit.chrom +":" +(infoNode.getDrawPosition()));
			if(infolist.get(infolist.size()-1).length() > infolist.get(maxlength).length()) {
				maxlength = infolist.size()-1;
			}
			infolist.add("Length: " +infoNode.getLength());
			if(infolist.get(infolist.size()-1).length() > infolist.get(maxlength).length()) {
				maxlength = infolist.size()-1;
			}
			if(infoNode.getTrack().hasvalues) {
				infolist.add("Score: " +MethodLibrary.round(infoNode.value,2));
				if(infolist.get(infolist.size()-1).length() > infolist.get(maxlength).length()) {
					maxlength = infolist.size()-1;
				}
			}
			if(infoNode.forward != null) {
				if(infoNode.forward) {
					infolist.add("Strand: +");
				}
				else {
					infolist.add("Strand: -");
				}			
			}
			if(infoNode.info != null) {
				String[] split = infoNode.info.split(";");
				for(int i = 0; i<split.length;i++) {
					infolist.add(split[i]);
					if(infolist.get(infolist.size()-1).length() > infolist.get(maxlength).length()) {
						maxlength = infolist.size()-1;
					}
				}
			}
			if(fm== null) {
				fm = nodebuf.getFontMetrics();
			}
			boxmetrics = fm.getStringBounds(infolist.get(maxlength), nodebuf).getBounds();				
			
		}
		
		buf.setColor(Color.white);
		buf.fillRoundRect(Main.sidebarWidth +10, overlapindex+10, boxmetrics.width+8, boxmetrics.height*infolist.size(),10,10);
		buf.setColor(Color.black);
		for(int i = 0; i<infolist.size(); i++) {
			buf.drawString(infolist.get(i), Main.sidebarWidth +14, (overlapindex+6)+(boxmetrics.height*(i+1)));
		}	
	//}
}

public class Annotator extends SwingWorker<String, Object> {
 
	BedTrack track;
	
	public Annotator(BedTrack track) {
		this.track = track;
	}
	
	public void annotateVars() {
		try {
			
		if(FileRead.head.getNext() == null) {
			return;
		}
		track.annotator = true;
		track.used = false;
		Main.bedCanvas.annotator = true;
		counter = 0;
		VarNode node = FileRead.head.getNextVisible(FileRead.head), first = null, last = null;
		//VarNode node = FileRead.head.getNext(FileRead.head), first = null, last = null;
		boolean isfirst = true, cancelled = false; 		
		
		if(node != null) {
			
		int start = node.getPosition()-3;
		if(track.file != null && track.file.length() / 1048576 < Settings.settings.get("bigFile")) {
			start = 1;
		}
		int end = start +Settings.windowSize;
		BedReader reader = new BedReader(track, start, end);
		if(track.file != null && track.file.length() / 1048576 < Settings.settings.get("bigFile")) {
			first = node;
		}
		
			while(node != null && start < Main.drawCanvas.splits.get(0).chromEnd) {
				
				Main.bedCanvas.annotator = true;
				if(track.file != null) {
					Main.drawCanvas.loadingtext = "Annotating variants with " +track.file.getName();
				}
				else {
					Main.drawCanvas.loadingtext = "Annotating variants with " +FilenameUtils.getName(track.url.getFile());
				}
				
				if(!Main.drawCanvas.loading) {
					BedCanvas.getZoomlevel(track, start, end);
					annotator = false;
					track.used = false;
					track.intersect = false;
					removeBedhits(track);
					cancelled = true;
					repaint();
					track.annotator = false;					
					break;
				}
				if((track.file != null && track.file.length() / 1048576 >= Settings.settings.get("bigFile")) || !track.small) {
					if(last != null) {
						node = last.getNextVisible(last);
						if(node == null) {
							break;
						}
						start = node.getPosition()-3;				
						//end = start +Settings.windowSize;
						if(track.url != null) {
							end = node.getPosition()+3;
						}
						else {
							end = start +Settings.windowSize;
						}
					}
				}
				
				if((track.file != null && track.file.length() / 1048576 >= Settings.settings.get("bigFile"))  || !track.small) {
					isfirst = true;
					while(node != null && node.getPosition() < end) {
						
						if(!Main.drawCanvas.loading) {
							BedCanvas.getZoomlevel(track, start, end);
							annotator = false;
							track.used = false;
							annotator = false;
							track.annotator = false;
							cancelled = true;
							break;
						}
				
							if(isfirst ) {
								first = node;
								isfirst = false;
								
							}
						
							last = node;							
							node = node.getNextVisible(node);
					  }
				}
			
				Main.drawCanvas.loadBarSample = (int)((start/(double)(Main.drawCanvas.splits.get(0).chromEnd)*100));    
				Main.drawCanvas.loadbarAll = (int)((start/(double)(Main.drawCanvas.splits.get(0).chromEnd)*100));    
				if((track.file != null && track.file.length() / 1048576 >= Settings.settings.get("bigFile"))  || !track.small) {
					reader.start = first.getPosition()-3;				
				//	reader.end = last.getPosition()+3;
					if(track.url != null) {
						reader.end = first.getPosition()+3;
					}
					else {
						reader.end = last.getPosition()+3;
					}
				}
				
				
				removeBeds(track);				
				reader.getNodes();
				//System.out.println(start +" " +end);
				if(track.getHead().getNext() != null) {					
					annotate(track.getHead(), first);
				}			
				if((track.file != null && track.file.length() / 1048576 < Settings.settings.get("bigFile")) || !track.small) {
					start = end - 2000;
					end = start + Settings.windowSize;
					reader.start = start;
					reader.end = end;
					while(node != null && node.getPosition() < start)  {
						first = node;
						node = node.getNext();
					}
				}
			/*	Draw.updatevars = true;
				Main.drawCanvas.repaint();			
		*/
		}
		if(reader.tabixReader != null) {
			reader.tabixReader.close();
		}
		}
		/*
		while(start < Main.drawCanvas.splits.get(0).chromEnd) {						
			if(!Main.drawCanvas.loading ) {
				break;
			}
			
			reader.getNodes();
		
			Main.drawCanvas.loadBarSample = (int)((start/(double)(Main.drawCanvas.splits.get(0).chromEnd)*100));    
			Main.drawCanvas.loadbarAll = (int)((start/(double)(Main.drawCanvas.splits.get(0).chromEnd)*100));    
			if(track.getHead().getNext() != null) {
				
				annotate(track.getHead());
			}		
			
			Draw.updatevars = true;
			Main.drawCanvas.repaint();
			
				start = end;
				
				end=start +10000000;	
				reader.start = start;
				reader.end = end;
			}	
		*/
		
		if(track.annotator) {
			VarNode current = FileRead.head.getNext();
			Object[] t2 = new Object[Main.bedCanvas.bedTrack.size()];
			Object[] t1 = new Object[Main.bedCanvas.bedTrack.size()];
			track.annotator = false;
			for(int i = 0; i<t1.length; i++) {
				if(!Main.bedCanvas.bedTrack.get(i).intersect) {
					continue;
				}
				if(Main.bedCanvas.bedTrack.get(i).getIntersectBox().isSelected()) {			
					t1[i] = 1;		
				}				
			}	
			while(current != null) {	
				
				current.bedhit = checkIntersect(current, t1, t2);	
				current = current.getNext();
			}
			
		}
		
		if(VariantHandler.commonSlider.getValue() > 1 && VariantHandler.clusterSize > 0) {
	    	Main.drawCanvas.clusterId = 1;
			Main.drawCanvas.calcClusters(FileRead.head);
		}
		VariantHandler.clusterTable.repaint();
		Draw.updatevars = true;
		Draw.calculateVars = true;
		Main.drawCanvas.repaint();
		
		if(!cancelled) {
			track.used = true;
		}
		else {
			track.used = false;
			if(track.small) {
				Main.bedCanvas.getMoreBeds(track);
			}			
		}
		removeBeds(track);
		if(FileRead.bigcalc) {
			Main.drawCanvas.calcClusters(FileRead.head);
		}
		else {
			
			Main.drawCanvas.calcClusters(FileRead.head,1);
			
		}
		
		Main.bedCanvas.annotator = false;
		
		
		if(!cancelled) {
			bedOn = true;
		}
		
		Main.drawCanvas.calcClusters(FileRead.head.getNext(),1);
		BedCanvas.annoTrack = null;
		Draw.updatevars = true;
		Main.drawCanvas.repaint();
		//getZoomlevel(track, track.bedstart, track.bedend);
	
			getBEDfeatures(track, track.bedstart, track.bedend);
		
		//repaint();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	
	}
	protected String doInBackground() throws Exception {
		Main.drawCanvas.ready("Loading tracks");
		if(!Main.drawCanvas.loading) {
			
			Main.drawCanvas.loading("Annotating variants");
		}
		if(!track.small || track.file.getName().toLowerCase().endsWith(".bw") || track.file.getName().toLowerCase().endsWith(".bigwig")) {
			
			annotateVars();
			
		}
		else {
			
			try {
				
				annotate(track.getHead(), FileRead.head.getNext());
				BedCanvas.annoTrack = null;
				annotator = false;
							
				if(VariantHandler.commonSlider.getValue() > 1 && !VariantHandler.clusterBox.getText().equals("0") && !VariantHandler.clusterBox.getText().equals("")) {
			    	Main.drawCanvas.clusterId = 1;
					Main.drawCanvas.calcClusters(FileRead.head);
				}
				
				VariantHandler.clusterTable.repaint();
				Draw.updatevars = true;
				Draw.calculateVars = true;
			
				Main.drawCanvas.repaint();
			}
			catch(Exception e) {
				Main.drawCanvas.ready("Annotating variants");
				Main.drawCanvas.ready("Loading tracks");
				e.printStackTrace();
			}
		}
		Main.drawCanvas.ready("Annotating variants");
		Main.drawCanvas.ready("Loading tracks");
		Draw.updatevars = true;
		Draw.calculateVars = true;
		
		Main.drawCanvas.repaint();
		return null;
	}	
	
}

public class bedFeatureFetcher extends SwingWorker<String, Object> {	
	
	@Override
	protected String doInBackground() throws Exception {
		
		try {
			if(Main.drawCanvas.loading) {
				return null;
			}
			Main.drawCanvas.loading("Loading tracks");
		for(int i = 0 ; i<bedTrack.size(); i++) {
			if(bedTrack.get(i).graph) {
				calcScale(bedTrack.get(i));
				
			}
			
			getMoreBeds(bedTrack.get(i));
			
		}
		
		Main.drawCanvas.ready("Loading tracks");
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
		return null;
	}	
}
public static void getZoomlevel(BedTrack track, int start, int end) {
	if(track.getBBfileReader() == null) {
		return;
	}
	track.setZoomlevel(1);
	for(int i =2;i<track.getBBfileReader().getZoomLevels().getZoomHeaderCount();i++) {
		if(track.getBBfileReader().getZoomLevels().getZoomLevelHeader(i).getReductionLevel() < ((end-start)/(Main.drawCanvas.splits.get(0).pixel*(end-start)))) {
			track.setZoomlevel(i);
		}
		else {
			break;
		}
	}
}
public class BedReader extends SwingWorker<String, Object> {
	int start, end;
	BedTrack track;
	SeekableStream stream = null;
	TabixReaderMod tabixReader = null;
	BufferedReader bedreader = null;
	
	public BedReader(BedTrack track, int start, int end) {
		this.start = start;
		
		if(this.start < 0) {
			this.start = 0;
		}
		this.end = end;
		this.track = track;
		track.loading = true;
		if(track.file != null && track.file.getName().endsWith(".txt")) {
			
			getGeneTxt(track);
			track.loading = false;
			if(track.file == null) {
				Main.drawCanvas.ready("Loading track: " +track.url.toString().substring(track.url.toString().lastIndexOf("/")+1)); 
			}
			else {
				Main.drawCanvas.ready("Loading track " +track.file.getName());
			}
		}
	}
	
	public void getNodes() {
		
		try {			
			
			if(!Main.drawCanvas.loading) {
				return;
			}
			
			if(tabixReader != null) {
				tabixReader.close();
			}			
			
			if(track.file != null) {	
				if(track.file.getName().endsWith(".txt")) {					
					getGeneTxt(track);
					return;
				}
				
				else if(track.file.getName().toLowerCase().endsWith(".bb") || track.file.getName().toLowerCase().endsWith(".bigwig") || track.file.getName().toLowerCase().endsWith(".bigbed")|| track.file.getName().toLowerCase().endsWith(".bw")) {
					if(track.getBBfileReader() == null) {
						track.setBBfileReader(new BBFileReader(track.file.getCanonicalPath(), track));
					}
					 if(track.first) {
						 if(track.getBBfileReader().getChromosomeNames().get(0).contains("chr")) {
							 track.chr = "chr";
						 }					
						 track.first = false;
					 }						
				}
				else {
					if(track.file.getName().toLowerCase().endsWith(".bed")) {						
						bedreader = new BufferedReader(new FileReader(track.file.getCanonicalPath()));						
					}
					else {
						try {					
							
							tabixReader = new TabixReaderMod(track.file.getCanonicalPath());  
							
						}
						catch(Exception e) {
							bedreader = new BufferedReader(new FileReader(track.file.getCanonicalPath()));  
						}
 					}
				}
			}
			else {			
				if(stream != null) {
					stream.close();
				}
				
				stream = SeekableStreamFactory.getInstance().getStreamFor(track.url);		
				
				if(track.index != null) {
					if(tabixReader != null) {
						tabixReader.close();
					}
					
					tabixReader = new TabixReaderMod(track.url.toString(), track.index.toString(), stream);
					
				}
				else {
					try {						
						if(track.getBBfileReader() == null) {							
						 track.setBBfileReader(new BBFileReader(track.url.toString(), stream, track));
						}
						 if(track.first) {
							 
							 if(track.getBBfileReader().getChromosomeNames() != null && track.getBBfileReader().getChromosomeNames().get(0).contains("chr")) {
								 track.chr = "chr";
							 }
							
							 track.first = false;
						 }
					}
					catch(Exception ex) {
						ex.printStackTrace();
					}					
				}				
			}
			
			TabixReaderMod.Iterator bedIterator = null;
			BigWigIterator wigiter = null;
			BigBedIterator bediter = null;
			ZoomLevelIterator zoomiterator = null;
			
			String chrom ="";
					if(bedreader != null) {
						Index index = IndexFactory.loadIndex(track.file+".tbi");
						java.util.List<Block> blocks = index.getBlocks(Main.chromosomeDropdown.getSelectedItem().toString(), start, end);	
						chrom = Main.chromosomeDropdown.getSelectedItem().toString();						
					
						if(blocks.size() > 0) {							
							bedreader = new BufferedReader(new FileReader(track.file));
							bedreader.skip(blocks.get(0).getStartPosition());								
						}								
					}
					else if(tabixReader == null) {
						
						if(track.getBBfileReader().getChromosomeNames() != null && !track.getBBfileReader().getChromosomeNames().contains(track.chr+Main.chromosomeDropdown.getSelectedItem().toString())) {
							return;
						}
					 if(track.getBBfileReader().isBigWigFile()) {
						
						   track.setZoomlevel(1);
						
						if(!annotator) {
						   getZoomlevel(track, start, end);
						   
						}
						
						if( track.getZoomlevel() == 1) {
							
							if(track.file != null && track.file.length() / 1048576 >= Settings.settings.get("bigFile") && end-start > Settings.windowSize) {
								return;
							}
							/*
							else {*/
							
							/*if(track.file != null) {
								if(end-start > Settings.windowSize) {
									return;
								}
							 }
							*/
							
							 wigiter = track.getBBfileReader().getBigWigIterator(track.chr+Main.chromosomeDropdown.getSelectedItem().toString() , start, track.chr+Main.chromosomeDropdown.getSelectedItem().toString(), end, true);
							//}
						}
						else {
							
							 zoomiterator = track.getBBfileReader().getZoomLevelIterator(track.getZoomlevel(),track.chr+Main.chromosomeDropdown.getSelectedItem().toString() , start, track.chr+Main.chromosomeDropdown.getSelectedItem().toString(), end, true);
							
						}
							
						 
							
						//	 ZoomDataRecord features;
								/*
							   while(iterator.hasNext()) {
							    	
										//split = line.split("\t");
							    			    		
							    		features = iterator.next();
							    		System.out.println((features.getChromEnd()-features.getChromStart()) +" " +features.getMaxVal());
							   }*/
							 }
					 else {
						   track.setZoomlevel(1);
						   
							if(!annotator) {
							   for(int i =2;i<track.getBBfileReader().getZoomLevels().getZoomHeaderCount();i++) {
									if(track.getBBfileReader().getZoomLevels().getZoomLevelHeader(i).getReductionLevel() < (Main.drawCanvas.splits.get(0).viewLength/(Main.drawCanvas.splits.get(0).pixel*Main.drawCanvas.splits.get(0).viewLength))) {
										track.setZoomlevel(i);
									}
									else {
										break;
									}
								}
							}
							
							if(track.getZoomlevel() == 1) {
								 bediter  = track.getBBfileReader().getBigBedIterator(track.chr+Main.chromosomeDropdown.getSelectedItem().toString() , start, track.chr+Main.chromosomeDropdown.getSelectedItem().toString(), end, false);
									
							}
							else {
								 zoomiterator = track.getBBfileReader().getZoomLevelIterator(track.getZoomlevel(),track.chr+Main.chromosomeDropdown.getSelectedItem().toString() , start, track.chr+Main.chromosomeDropdown.getSelectedItem().toString(), end, false);
								 
							}
							//if(track.file.getName().toLowerCase().endsWith(".bb") || track.file.getName().toLowerCase().endsWith(".bigbed")) {
										
							/*}
							else if( track.getZoomlevel() == 1) {
								
								 bediter  = track.getBBfileReader().getBigBedIterator(track.chr+Main.chromosomeDropdown.getSelectedItem().toString() , start, track.chr+Main.chromosomeDropdown.getSelectedItem().toString(), end, false);
										
							}
							else {
								
								 zoomiterator = track.getBBfileReader().getZoomLevelIterator(track.getZoomlevel(),track.chr+Main.chromosomeDropdown.getSelectedItem().toString() , start, track.chr+Main.chromosomeDropdown.getSelectedItem().toString(), end, false);
								 
							}
							*/
						}				
					 
				 }
				
				 else {
					 
					 if(track.first) {
						
						 if(track.chr == null) {
							 track.chr = "";
						 }
						 String chr;
						 Iterator<String> chriterator = tabixReader.getChromosomes().iterator();
						 while(chriterator.hasNext()) {
							 chr = chriterator.next();
							 if(chr.matches("\\d+")) {
								 break;
							 }
							 else if(chr.startsWith("chr")) {
								 track.chr = "chr";
								 break;
							 }
						 }						
						
						 track.first = false;						
					 }					
					try {
						
						if(start == 0 && end == 0) {
							return;
						}
						
						bedIterator = tabixReader.query(track.chr +Main.chromosomeDropdown.getSelectedItem().toString() +":" +start +"-" +end);
						
					}
					catch(Exception e) {
						e.printStackTrace();
						System.out.println(track.chr +Main.chromosomeDropdown.getSelectedItem().toString() +" " +start +" " +end);
						track.getHead().putNext(null);						 
					}
				 }				
				/*
				  if((wigiter != null && !wigiter.hasNext() && zoomiterator == null) || bedIterator == null || bedIterator.next() == null) {
					  try {
						  	if(tabixReader == null && bedreader == null) {
						  		
						  		if(track.getBBfileReader().isBigWigFile()) {
						  			
						  			wigiter = track.getBBfileReader().getBigWigIterator(track.chr +Main.chromosomeDropdown.getSelectedItem().toString() , start, "chr" +Main.chromosomeDropdown.getSelectedItem().toString(), end, true);
						  			
						  		}
						  		else {
						  			
						  			 bediter  = track.getBBfileReader().getBigBedIterator(track.chr +Main.chromosomeDropdown.getSelectedItem().toString() , start,"chr" + Main.chromosomeDropdown.getSelectedItem().toString(), end, true);
										
						  		}
							 }
							 else {
								 bedIterator = tabixReader.query(track.chr +Main.chromosomeDropdown.getSelectedItem().toString()+":" +start +"-" +end );
							 }
					  }
					  catch(Exception ex) {
						  ex.printStackTrace();
						  ErrorLog.addError(ex.getStackTrace());
						
					  }
				  }
			 */
			
			if(track.file != null) {  
				if(track.file.getName().toLowerCase().endsWith("bed")) {
					
					iterateBEDfile(bedreader, track, chrom);
				}
				else if(track.file.getName().toLowerCase().endsWith("bed.gz") || track.file.getName().toLowerCase().endsWith("bedgraph.gz") || track.file.getName().toLowerCase().endsWith("bedgraph")) {
					
					iterateBED(bedIterator, track);
				}
				else if (track.file.getName().toLowerCase().endsWith("gff.gz") || track.file.getName().toLowerCase().endsWith("gff3.gz")){
					
					iterateGFF(bedIterator, track);
					
				}
				else if(track.file.getName().toLowerCase().endsWith(".bigwig") || track.file.getName().toLowerCase().endsWith(".bw")) {
					if(zoomiterator != null && !annotator) {

						
						iterateZoom(zoomiterator, track); 
					}
					else {
						
						iterateWig(wigiter, track);
					}
				}
				else if(track.file.getName().toLowerCase().endsWith(".bb") || track.file.getName().toLowerCase().endsWith(".bigbed")) {					
					iterateBigBed(bediter, track);					
					
				}
				else if(track.file.getName().toLowerCase().endsWith("tsv.gz") || track.file.getName().toLowerCase().endsWith("tsv.bgz")) {
					iterateTSV(bedIterator, track);
				}
			}
			else {
					if(track.url.toString().toLowerCase().endsWith("bed.gz") || track.url.toString().toLowerCase().endsWith("bedgraph.gz")) {
					  iterateBED(bedIterator, track);
					}
					else if(track.url.toString().toLowerCase().endsWith("tsv.gz") || track.url.toString().toLowerCase().endsWith("tsv.bgz")) {
						iterateTSV(bedIterator, track);
					}
					else if (track.url.toString().toLowerCase().endsWith("gff.gz")){
						iterateGFF(bedIterator, track);
					}
					else if(track.getBBfileReader().isBigBedFile()) 
					{
						if(zoomiterator != null) {
							
							iterateZoom(zoomiterator, track); 
						}
						else {
							iterateBigBed(bediter, track);
						}
						
					}
					else if(track.getBBfileReader().isBigWigFile()) {
						if(zoomiterator != null) {
							
							iterateZoom(zoomiterator, track); 
						}
						else {
							iterateWig(wigiter, track);
						}
					}
			}
			
			if( track.getHead().getNext() != null && (track.minvalue != Double.MAX_VALUE && (track.minvalue != 0 || track.maxvalue != 0))) {
				
				track.getLimitField().setEditable(true);
				track.hasvalues = true;
				 
				track.getLimitField().setEditable(true);
				//if(!track.selex) {
					/*if(!track.graph) {
						pressGraph(track);
					}
					*/
			//	}
	   		}
			/*if(!track.hasvalues) {				
		    	track.getLimitField().setVisible(false);		    	  
			}*/
			
		  if(!annotator) {
			  if(tabixReader != null) {
					tabixReader.close();
				}
		  }
	      if(stream != null) {
	    	  stream.close();
	      }
	      if(bedreader != null) {
	    	  bedreader.close();
	      }
	    
		
		if(Draw.variantcalculator && track.intersect && track.getIntersectBox().isSelected()) {
			if(track.getHead().getNext() == null && track.small && track.getBBfileReader() == null) {
				FileRead.nobeds = true;
			}
		}	
		
		track.setDrawNode(track.getHead());
		if(track.waiting) {		
		
			track.waiting = false;
			Main.bedCanvas.annotate(track.getHead(), FileRead.head.getNext());
			
			if(FileRead.bigcalc) {
				Main.drawCanvas.calcClusters(FileRead.head);
			}
			else {
				Main.drawCanvas.calcClusters(FileRead.head,1);
			}			
			Draw.updatevars = true;
		}
		}
		catch(Exception e) {
			ErrorLog.addError(e.getStackTrace());
			e.printStackTrace();
			
		//	Main.drawCanvas.ready("Reading BED-file");
		}		
		track.loading = false;
		updateTrack = track;
		
		
		
		//repaint();
		//Main.drawCanvas.ready("Reading BED-file");
		 
	}
	protected String doInBackground() {		
		
		if(track.file != null && track.file.getName().endsWith(".txt")) {			
			return "";
		}
		if(!Main.bedCanvas.annotator) {
			if(track.file == null) {
				Main.drawCanvas.loading("Loading track: " +track.url.toString().substring(track.url.toString().lastIndexOf("/")+1)); 
			}
			else {
				Main.drawCanvas.loading("Loading track " +track.file.getName());	
			}
			getNodes();
			
			if(track.file == null) {
				Main.drawCanvas.ready("Loading track: " +track.url.toString().substring(track.url.toString().lastIndexOf("/")+1)); 
			}
			else {
				Main.drawCanvas.ready("Loading track " +track.file.getName());
			}
		}
		
		return "";
	}
	
}
void getGeneTxt(BedTrack track) {
	try {
	 String line;
	 
	 track.getHead().putNext(null);
	 BedNode addNode = track.getHead();
	 addNode.putNext(null);
	 track.setCurrent(track.getHead());
	 int position = 0;
	 try {
		
	 BufferedReader reader = new BufferedReader(new FileReader(track.file.getCanonicalFile()));
	 String gene;
	 String[] result = {};
	 
	 while((line = reader.readLine()) != null) {
		
		 gene = line.replace(" ", "");	
		 if(Main.searchTable.containsKey(gene)) {
			 result = Main.searchTable.get(gene);
		 }
		 else if(Main.geneIDMap.containsKey(gene)) {
			 result = Main.searchTable.get(Main.geneIDMap.get(gene));
		 }
		 else {
			 ErrorLog.addError(gene +" not found.");
			 continue;
		 }		 
		
		 if(result != null && result.length == 3) {
			
			 //System.out.println(result[0] +"\t" +result[1] +"\t" +result[2] +"\t" +gene);
			 if(result[0].equals(Main.chromosomeDropdown.getSelectedItem().toString())){
				 position = Integer.parseInt(result[1]);
				 addNode = track.getHead();
				 while(addNode.getNext() != null && addNode.getNext().getPosition() < position) {
					 addNode = addNode.getNext();
				 }
				 
				 if(addNode.getNext() == null) {
					
					 addNode.putNext(new BedNode(result[0],position, Integer.parseInt(result[2])-Integer.parseInt(result[1]), track));				
					 addNode.getNext().putPrev(addNode);
					 addNode = addNode.getNext();
					 addNode.name = gene;
				 }
				 else {
				
					 BedNode adder = new BedNode(result[0],position, Integer.parseInt(result[2])-Integer.parseInt(result[1]), track);
					 addNode.getNext().putPrev(adder);
					 adder.putNext(addNode.getNext());
					 addNode.putNext(adder);
					 adder.putPrev(addNode);
					 adder.name = gene;
				 }
			 }
			
		 }	
		 else {
			 ErrorLog.addError(gene +" not found.");
			 continue;
		 }	
		
	 }
	 reader.close();
	 track.setCurrent(track.getHead());
	
	 }
	 catch(Exception e) {
		 
		 e.printStackTrace();
	 }
	}
	catch(Exception ex) {
		ex.printStackTrace();
	}
	}

void iterateBEDfile(BufferedReader reader, BedTrack track, String chrom) {
	
	try {
		String line;
		BedNode addNode = track.getHead();
		track.maxvalue = 0;
		track.minvalue = Double.MAX_VALUE;
		chrom = Main.chromosomeDropdown.getSelectedItem().toString();
		addNode = track.getHead();
		track.getHead().putNext(null);
		track.setCurrent(track.getHead());
		boolean first = true, bedgraph = false, firstrow = true;
		
		 BEDFeature features;
		 BEDCodecMod bedcodec = new BEDCodecMod();
		if(track.limitValue == null) {
			track.limitValue = (double)Integer.MIN_VALUE;
		}
		
	   while((line = reader.readLine()) != null) {			
		   if(!Main.drawCanvas.loading) {
			   
			   track.used = false;
				track.intersect = false;
				removeBedhits(track);
				
				repaint();
	    		break;
	    	}
	    	try {
	    		if(line.startsWith("#") || line.length() < 10 || line.startsWith("track")) {
	    			continue;
	    		}
	    		if(firstrow) {
	    			if(line.startsWith("chr")) {
	    				track.chr = "chr";
	    			}
	    			firstrow = false;
	    		}
	    		features = bedcodec.decode(line);
	    	
	    		if(!features.getContig().equals(track.chr +chrom)) {
	    			break;
	    		}
	    		if(track.selex && (features.getName() == null || !Main.SELEXhash.containsKey(features.getName()))) {	    			
					continue;
				}
	    		
	    		if(track.limitValue != (double)Integer.MIN_VALUE) {
	    			if(bedgraph) {
	    				if(!Double.isNaN(Double.parseDouble(features.getName()))) {
		    				if(Double.parseDouble(features.getName()) < track.limitValue) {
		    					continue;
		    				}
	    				}	
	    			}
	    			else { 
	    				if(!Double.isNaN((double)features.getScore())) {
		    				if(features.getScore() < track.limitValue) {
		    					continue;
		    				}
	    				}	    			
	    			}
	    		}
	    		
	    		if( (features.getEnd()-(features.getStart()-1)) == 1) {
	    			addNode.putNext(new BedNode(features.getContig(),features.getStart()-1-track.iszerobased, (features.getEnd()-(features.getStart()-1)), track));
	    		}
	    		else {
	    			addNode.putNext(new BedNode(features.getContig(),features.getStart()-1-track.iszerobased, (features.getEnd()-(features.getStart()-1))+track.iszerobased, track));	
	    		}
				addNode.getNext().putPrev(addNode);
				addNode = addNode.getNext();
				
				if(features.getName() != null) {					
					
					if(first) {						
						if(Main.SELEXhash.containsKey(features.getName().replace(".pfm", ""))) {						
							track.selex = true;
							track.getAffinityBox().setVisible(true);
							//track.iszerobased = 1;
							//track.getZerobased().setSelected(false);
						}						
					}
					if(track.selex) {
						addNode.id = features.getName().replace(".pfm", "");
						addNode.name = Main.factorNames.get(addNode.id);
						
					}
					else {
						addNode.name = features.getName();
					}					
				}
				
				if(bedgraph) {
					addNode.value = Double.parseDouble(features.getName());
					
				}
				else {
					addNode.value = (double)features.getScore();					
				}				
				
				if(!Double.isNaN(addNode.value)) {
					if(track.maxvalue < addNode.value) {
						track.maxvalue = addNode.value;
						
					}
					if(track.minvalue > addNode.value) {
						track.minvalue = addNode.value;
					}
				}
				if(features.getStrand() != null) {					
					addNode.forward = features.getStrand().equals(Strand.NEGATIVE) ? false : true;
				}
				
				if(features.getColor() != null) {	
					if(features.getColor().getRGB() != -16777216) {
						if(track.getColors().containsKey(features.getColor().getRGB())) {
							addNode.color = features.getColor().getRGB();
						}
						else {						
							track.getColors().put(features.getColor().getRGB(),features.getColor());
							addNode.color = features.getColor().getRGB();
						}
					}
				}			
	    	}
	    	catch(Exception ex) {
	    		ex.printStackTrace();
	    	}
	 
	    	first = false;
	     }     
	 		
   		
	   if(track.minvalue < 0) {
		   track.negatives = true;
		   track.scale = Math.max(Math.abs(track.maxvalue), Math.abs(track.minvalue));
	   }
	   else {
		   track.scale = track.maxvalue;
	   }
		   if(track.getBBfileReader() == null && ((track.minvalue == 0 && track.maxvalue == 0) || track.minvalue == Double.MAX_VALUE)) {
			   track.hasvalues = false;
		   }
		   else {
			   track.hasvalues = true;
		   }
	   iterator = null;
	     addNode = null;
	     track.cleared = false;
	     reader.close();
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	
}
void iterateBED(TabixReaderMod.Iterator iterator, BedTrack track) {
	 
	if(iterator == null) {
		return;
	}
	try {
		
	String line;
	 BedNode addNode = track.getHead();
	 track.maxvalue = 0;
		track.minvalue = Double.MAX_VALUE;

		track.getHead().putNext(null);
		
		 track.setCurrent(track.getHead());
		 boolean first = true, bedgraph = false;
		 
		 if(track.file != null && (track.file.getCanonicalPath().endsWith("bedgraph.gz") || track.file.getCanonicalPath().endsWith("bedgraph"))) {
			 bedgraph = true;			
		 }
		
		 BEDFeature features;
		 BEDCodecMod bedcodec = new BEDCodecMod();
		if(track.limitValue == null) {
			track.limitValue = (double)Integer.MIN_VALUE;
		}
		
	   while((line = iterator.next()) != null) {
		  
	    	if(!Main.drawCanvas.loading) {
	    		track.used = false;
	    		
				track.intersect = false;
				removeBedhits(track);
				
				repaint();
	    		break;
	    	}
		   try {
	    		features = bedcodec.decode(line);
	    		if(track.selex && (features.getName() == null || !Main.SELEXhash.containsKey(features.getName()))) {
	    			
					continue;
				}
	    		
	    		if(track.limitValue != (double)Integer.MIN_VALUE) {
	    			if(bedgraph) {
	    				if(!Double.isNaN(Double.parseDouble(features.getName()))) {
		    				if(Double.parseDouble(features.getName()) < track.limitValue) {
		    					continue;
		    				}
	    				}	
	    			}
	    			else { 
	    				if(!Double.isNaN((double)features.getScore())) {
		    				if(features.getScore() < track.limitValue) {
		    					continue;
		    				}
	    				}	    			
	    			}
	    		}
	    		if( (features.getEnd()-(features.getStart()-1)) == 1) {
	    			addNode.putNext(new BedNode(features.getContig(),features.getStart()-1-track.iszerobased, (features.getEnd()-(features.getStart()-1)), track));
	    		}
	    		else {
	    			addNode.putNext(new BedNode(features.getContig(),features.getStart()-1-track.iszerobased, (features.getEnd()-(features.getStart()-1))+track.iszerobased, track));	
	    		}
				addNode.getNext().putPrev(addNode);
				addNode = addNode.getNext();
				
				if(features.getName() != null) {					
					
					if(first) {						
						if(Main.SELEXhash.containsKey(features.getName().replace(".pfm", ""))) {						
							track.selex = true;
							track.getAffinityBox().setVisible(true);
							//track.iszerobased = 1;
							//track.getZerobased().setSelected(false);
						}						
					}
					if(track.selex) {
						addNode.id = features.getName().replace(".pfm", "");
						addNode.name = Main.factorNames.get(addNode.id);
						
					}
					else {
						addNode.name = features.getName();
					}					
				}
				
				if(bedgraph) {
					addNode.value = Double.parseDouble(features.getName());
					
				}
				else {
					addNode.value = (double)features.getScore();
					
				}				
				
				if(!Double.isNaN(addNode.value)) {
					
					if(track.maxvalue < addNode.value) {
						track.maxvalue = addNode.value;						
					}
					if(track.minvalue > addNode.value) {
						track.minvalue = addNode.value;
					}
					
				}
				
				if(features.getStrand() != null) {
					
					addNode.forward = features.getStrand().equals(Strand.NEGATIVE) ? false : true;
				}
				
				if(features.getColor() != null) {	
					if(features.getColor().getRGB() != -16777216) {
						
						if(track.getColors().containsKey(features.getColor().getRGB())) {
							addNode.color = features.getColor().getRGB();
						}
						else {
							
							track.getColors().put(features.getColor().getRGB(),features.getColor());
							addNode.color = features.getColor().getRGB();
						}
					}
				}
			
	    	}
	    	catch(Exception ex) {
	    		ex.printStackTrace();
	    	}
	    	
	    	first = false;
	     }     
	 		
	  
	   if(track.minvalue < 0) {
		   track.negatives = true;
		   track.scale = Math.max(Math.abs(track.maxvalue), Math.abs(track.minvalue));
	   }
	   else {
		   track.scale = track.maxvalue;
	   }
	   
	   if(track.getBBfileReader() == null && ((track.minvalue == 0 && track.maxvalue == 0) || track.minvalue == Double.MAX_VALUE)) {
		   track.hasvalues = false;
	   }
	   else {
		   track.hasvalues = true;
	   }
	  
	   iterator = null;
	     addNode = null;
	     track.cleared = false;
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}
void iterateBigBed(BigBedIterator iterator, BedTrack track) {
	try {
	
	 BedNode addNode = track.getHead();
	 addNode.putNext(null);
	 track.maxvalue = 0;
		track.minvalue = Double.MAX_VALUE;
		 track.setCurrent(track.getHead());
		
		 BedFeature features;
		String[] allfields = null;
	
	   while(iterator.hasNext()) {
		   if(!Main.drawCanvas.loading) {
			   track.used = false;
				track.intersect = false;
				removeBedhits(track);
				
				repaint();
	    		break;
	    	}
	    	try {
				//split = line.split("\t");
	    			    		
	    		features = iterator.next();
	    		allfields = features.getRestOfFields();
	    		
	    		if(track.limitValue != (double)Integer.MIN_VALUE) {	    			
    				if(!Double.isNaN(Double.parseDouble(allfields[1]))) {
	    				if(Double.parseDouble(allfields[1]) < track.limitValue) {
	    					continue;
	    				}
    				}    			
	    		}
	    //		position = features.getStart();
	    		if((features.getEndBase()-(features.getStartBase())) == 1) {
	    			addNode.putNext(new BedNode(features.getChromosome(),features.getStartBase()-track.iszerobased, (features.getEndBase()-(features.getStartBase())), track));
	    			
	    		}
	    		else {
	    			addNode.putNext(new BedNode(features.getChromosome(),features.getStartBase()-track.iszerobased, (features.getEndBase()-(features.getStartBase()))+track.iszerobased, track));
	    	    	
	    		}
				addNode.getNext().putPrev(addNode);
				addNode = addNode.getNext();
				
				if(allfields.length > 0) {
				//	if(!track.getNames().contains(allfields[0])) {
				//		track.getNames().put(allfields[0], allfields[0]);
						addNode.name = allfields[0];// track.getNames().get(allfields[0]);
				//	}
				//	else {
				//		addNode.name = track.getNames().get(allfields[0]);		
				//	}
					if(allfields.length > 1) {
						addNode.value = Double.parseDouble(allfields[1]);	
						if(!Double.isNaN(addNode.value)) {
							if(track.maxvalue < addNode.value) {
								track.maxvalue = addNode.value;
								
							}
							if(track.minvalue > addNode.value) {
								track.minvalue = addNode.value;
							}
						}
						if(allfields.length > 2) {
							addNode.forward = allfields[2].contains("-") ? false : true;
							if(allfields.length > 5) {
								if(allfields[4].contains(",")) {
									colorsplit = allfields[4].split(",");
									if(colorsplit.length > 2) {
										if(colorsplit.length == 3) {
											addColor = new Color(Integer.parseInt(colorsplit[0]), Integer.parseInt(colorsplit[1]), Integer.parseInt(colorsplit[2]));
											if(track.getColors().containsKey(addColor.getRGB())) {
												addNode.color = addColor.getRGB();
											}
											else {
												
												track.getColors().put(addColor.getRGB(),addColor);
												addNode.color = addColor.getRGB();
											}
										}
										else {
											addColor = new Color(Integer.parseInt(colorsplit[0]), Integer.parseInt(colorsplit[1]), Integer.parseInt(colorsplit[2]), Integer.parseInt(colorsplit[3]));
											if(track.getColors().containsKey(addColor.getRGB())) {
												addNode.color = addColor.getRGB();
											}
											else {
												
												track.getColors().put(addColor.getRGB(),addColor);
												addNode.color = addColor.getRGB();
											}
										}
										
									}
								}
							}
						}
					}
					
				}			
		
			
	    	}
	    	catch(Exception e) {
	  //  		tabixReader.close();
	    		e.printStackTrace();
	    		for(int i = 0 ; i<allfields.length; i++) {
	    			System.out.println(allfields[i] +"\t");
	    		}
	    	}
	    
	     }     
	   if(track.minvalue < 0) {
		   track.negatives = true;
	   }
	   if(track.getBBfileReader() == null && ((track.minvalue == 0 && track.maxvalue == 0) || track.minvalue == Double.MAX_VALUE)) {
		   track.hasvalues = false;
	   }
	   else {
		   track.hasvalues = true;
		   track.scale = Math.max(Math.abs(track.maxvalue), Math.abs(track.minvalue));
	   }
	   
	     addNode = null;
	     track.cleared = false;
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}
void iterateGFF(TabixReaderMod.Iterator iterator, BedTrack track) {
	try {
		if(iterator == null) {
			return;
		}
		String line;
		String[] split; 
	 	BedNode addNode = track.getHead();	 	
	 	addNode.putNext(null);
	 	track.maxvalue = 0;
		track.minvalue = Double.MAX_VALUE;
		track.setCurrent(track.getHead());
		boolean first = true;
		track.maxvalue = 0;
		track.minvalue = Double.MAX_VALUE;	
		
	   while((line = iterator.next()) != null) {
		   if(!Main.drawCanvas.loading) {
			   track.used = false;
				track.intersect = false;
				removeBedhits(track);
				
				repaint();
	    		break;
	    	}
	    	try {
				split = line.split("\t");	    			    		
	 			
				if(split.length > 4) {
					if(track.limitValue != (double)Integer.MIN_VALUE) {	    			
	    				if(!Double.isNaN(Double.parseDouble(split[5]))) {
		    				if(Double.parseDouble(split[5]) < track.limitValue) {
		    					continue;
		    				}
	    				}    			
		    		}
					//System.out.println(track.iszerobased);
					if(Integer.parseInt(split[4])-(Integer.parseInt(split[3])) == 1) {
						addNode.putNext(new BedNode(split[0],(Integer.parseInt(split[3]))-track.iszerobased,Integer.parseInt(split[4])-(Integer.parseInt(split[3])), track));
						
					}
					else {
						addNode.putNext(new BedNode(split[0],(Integer.parseInt(split[3]))-track.iszerobased,Integer.parseInt(split[4])-(Integer.parseInt(split[3]))+track.iszerobased, track));
						
					}
					
					addNode.getNext().putPrev(addNode);
					
					addNode = addNode.getNext();
					
					if(first) {				
					
						if(Main.SELEXhash.containsKey(split[2].replace(".pfm", ""))) {						
							track.selex = true;
						}						
					}
					if(track.selex) {
						addNode.id = split[2].replace(".pfm", "");
						addNode.name = Main.factorNames.get(addNode.id);
						
					}
					else {
						
						addNode.name = split[2];
					}
					if(split.length > 5) {						
						try {
							addNode.value = Double.parseDouble(split[5]);
						}
						catch(Exception e) {
							e.printStackTrace();
							addNode.value = null;
						}
						if(addNode.value != null && !Double.isNaN(addNode.value)) {
							if(track.maxvalue < addNode.value) {
								track.maxvalue = addNode.value;
								
							}
							if(track.minvalue > addNode.value) {
								track.minvalue = addNode.value;
							}
						}
						
						if(split.length > 6) {						
							addNode.forward = split[6].contains("-") ? false : true;	
							
							if(split.length > 8) {
								addNode.info = split[8];
								
								String[] infoFields = split[8].split(";");
								for(int i = 0 ; i<infoFields.length; i++) {
									if(infoFields[i].startsWith("Name=")) {
										if(addNode.id == null) {
											addNode.secondaryName = split[2];
											addNode.name = infoFields[i].substring(5);
										}
										
									}
									else if(infoFields[i].startsWith("pfm=")) {
										if(Main.factorNames.get(infoFields[i].substring(4).replace(".pfm", "")) != null) {
											addNode.secondaryName = split[2];
											addNode.id = infoFields[i].substring(4).replace(".pfm", "");
											addNode.name = Main.factorNames.get(addNode.id);
										}									
									}
								}
								
							}
						}
						
						
					}				
				}
				else {
					break;
				}			
			
	    	}
	    	catch(Exception e) {
	 //   		tabixReader.close();
	    		e.printStackTrace();
	    	}
	    	first = false;
	     }     
	   if(track.minvalue < 0) {
		   track.negatives = true;
	   }
	   if(track.getBBfileReader() == null && ((track.minvalue == 0 && track.maxvalue == 0) || track.minvalue == Double.MAX_VALUE)) {
		   track.hasvalues = false;
	   }
	   else {
		   track.hasvalues = true;
		   track.scale = Math.max(Math.abs(track.maxvalue), Math.abs(track.minvalue));
	   }
	  
	     addNode = null;
	    
	     track.cleared = false;
	    
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}
void iterateTSV(TabixReaderMod.Iterator iterator, BedTrack track) {
	try {
		if(iterator == null) {
			return;
		}
		String line;
		String[] split; 
	 	BedNode addNode = track.getHead();	 	
	 	addNode.putNext(null);
	 	track.maxvalue = 0;
		track.minvalue = Double.MAX_VALUE;
		track.setCurrent(track.getHead());
		boolean first = true;
		track.maxvalue = 0;
		track.minvalue = Double.MAX_VALUE;	
		
	   while((line = iterator.next()) != null) {
		   if(!Main.drawCanvas.loading) {
			   track.used = false;
				track.intersect = false;
				removeBedhits(track);
				
				repaint();
	    		break;
	    	}
	    	try {
				split = line.split("\t"); 			
				
					if(track.limitValue != (double)Integer.MIN_VALUE) {	    
						if(track.valuecolumn != null) {
							if(!Double.isNaN(Double.parseDouble(split[track.valuecolumn]))) {
			    				if(Double.parseDouble(split[track.valuecolumn]) < track.limitValue) {
			    					continue;
			    				}
		    				} 
						}	    				   			
		    		}
					if(track.endcolumn != null) {
						addNode.putNext(new BedNode(split[track.chromcolumn],(Integer.parseInt(split[track.startcolumn]))-track.iszerobased,Integer.parseInt(split[track.endcolumn])-(Integer.parseInt(split[track.startcolumn])), track));
					}
					else {
						addNode.putNext(new BedNode(split[track.chromcolumn],Integer.parseInt(split[track.startcolumn])-track.iszerobased,1, track));
						
					}
					addNode.getNext().putPrev(addNode);
					
					addNode = addNode.getNext();
					
					if(first) {				
						
						if(Main.SELEXhash.containsKey(split[2].replace(".pfm", ""))) {						
							track.selex = true;
						}						
					}
				
						if(track.namecolumn != null) {
							addNode.name = split[track.namecolumn];
						}
						if(track.strandcolumn != null) {
							addNode.forward = split[track.strandcolumn].contains("-") ? false : true;
						}
						if(track.basecolumn != null) {
							addNode.name = split[track.basecolumn];
						}
									
						try {
							if(track.valuecolumn != null) {
								addNode.value = Double.parseDouble(split[track.valuecolumn]);
							}
						}
						catch(Exception e) {
							e.printStackTrace();
							addNode.value = null;
						}
						if(addNode.value != null && !Double.isNaN(addNode.value)) {
							if(track.maxvalue < addNode.value) {
								track.maxvalue = addNode.value;
								
							}
							if(track.minvalue > addNode.value) {
								track.minvalue = addNode.value;
							}
						}
						/*
						if(split.length > 6) {						
							addNode.forward = split[6].contains("-") ? false : true;	
							
							if(split.length > 8) {
								addNode.info = split[8];
							}
						
					}	*/			
						
			
	    	}
	    	catch(Exception e) {
	 //   		tabixReader.close();
	    		e.printStackTrace();
	    	}
	    	first = false;
	     }     
	   if(track.minvalue < 0) {
		   track.negatives = true;
	   }
	   if(track.getBBfileReader() == null && ((track.minvalue == 0 && track.maxvalue == 0) || track.minvalue == Double.MAX_VALUE)) {
		   track.hasvalues = false;
	   }
	   else {
		   track.hasvalues = true;
		   track.scale = Math.max(Math.abs(track.maxvalue), Math.abs(track.minvalue));
	   }	  
	   addNode = null;	    
	   track.cleared = false;
	    
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}
void iterateWig(BigWigIterator iterator, BedTrack track) {
	
	try {
		
		if(iterator == null) {
			return;
		}
		BedNode addNode = track.getHead();
		track.maxvalue = 0;
		track.minvalue = Double.MAX_VALUE;
		track.setCurrent(track.getHead());		
		addNode.putNext(null);
		
		WigItem features;
		
	   while(iterator.hasNext()) {
		   if(!Main.drawCanvas.loading) {
			   track.used = false;
				track.intersect = false;
				removeBedhits(track);
				 addNode = null;
				//repaint();
			   return;
		   }
	    	try {
			
	    		features = iterator.next();
	    		try {
		    		if(features == null || features.getStartBase()-track.iszerobased > track.bedend) {
		    			break;
		    		}
	    		}
	    		catch(Exception e) {
	    			e.printStackTrace();
	    		}
	    		/*if(testfield == 778291 ) {
	    			System.out.println(features.getStartBase());
	    		}*/
	    		if(track.limitValue != null && track.limitValue != (double)Integer.MIN_VALUE) {	   
	    			
    				if((double)features.getWigValue() < track.limitValue) {    					
    					continue;
    				}    				
	    		}
	    		
	   			addNode.putNext(new BedNode(features.getChromosome(),features.getStartBase()-track.iszerobased, (features.getEndBase()-(features.getStartBase())), track));
				addNode.getNext().putPrev(addNode);
				addNode = addNode.getNext();			
				addNode.value = (double)features.getWigValue();		
				//addNode.name = addNode.getChrom() +":" +features.getStartBase();
				if(!Double.isNaN(addNode.value)) {
					if(track.maxvalue < addNode.value) {
						track.maxvalue = addNode.value;		
						track.scale = Math.max(Math.abs(track.maxvalue), Math.abs(track.minvalue));
					}
					if(track.minvalue > addNode.value) {
						track.minvalue = addNode.value;	
						track.scale = Math.max(Math.abs(track.maxvalue), Math.abs(track.minvalue));
					}
					 
				}					
	    	}
	    	catch(Exception e) {
	  //  		tabixReader.close();
	    		e.printStackTrace();
	    		
	    		break;
	    	}
	    	
	     }     
	   if(track.minvalue < 0) {
		   track.negatives = true;
	   }
	   if(track.getBBfileReader() == null && ((track.minvalue == 0 && track.maxvalue == 0) || track.minvalue == Double.MAX_VALUE)) {
		   track.hasvalues = false;
	   }
	   else {
		   track.hasvalues = true;
		   track.scale = Math.max(Math.abs(track.maxvalue), Math.abs(track.minvalue));
	   }	  
	  
	  
	     addNode = null;
	     track.cleared = false;
	    
	}
	catch(Exception e) {
		
		e.printStackTrace();
		
	}
}
void iterateZoom(ZoomLevelIterator iterator, BedTrack track) {
	try {
	
		BedNode addNode = track.getHead();
		track.maxvalue = 0;
		track.minvalue = Double.MAX_VALUE;
		track.setCurrent(track.getHead());		
		addNode.putNext(null);
		
		ZoomDataRecord features;
		
	   while(iterator.hasNext()) {
		   if(!Main.drawCanvas.loading) {
			   track.used = false;
				track.intersect = false;
				removeBedhits(track);
				
				repaint();
			   break;
		   }
	    	try {
	    		
	    		features = iterator.next();
	    		if(features.getChromStart()-track.iszerobased > track.bedend) {
	    			
	    			break;
	    		}
	    		
	    		if(track.limitValue != (double)Integer.MIN_VALUE) {   		
	    			
    				if((double)features.getMeanVal() < track.limitValue) {
    					
    					continue;
    				}    
    				
	    		}
	    		
	   			addNode.putNext(new BedNode(features.getChromName(),features.getChromStart()-track.iszerobased, (features.getChromEnd()-(features.getChromStart())), track));
				addNode.getNext().putPrev(addNode);
				addNode = addNode.getNext();			
				addNode.value = (double)features.getMeanVal();		
				//addNode.name = addNode.getChrom() +":" +features.getStartBase();
				if(!Double.isNaN(addNode.value)) {
					if(track.maxvalue < addNode.value) {
						track.maxvalue = addNode.value;		
						track.scale = Math.max(Math.abs(track.maxvalue), Math.abs(track.minvalue));
					}
					if(track.minvalue > addNode.value) {
						track.minvalue = addNode.value;	
						track.scale = Math.max(Math.abs(track.maxvalue), Math.abs(track.minvalue));
					}
					 
				}					
	    	}
	    	catch(Exception e) {
	    		e.printStackTrace();	   
	    	}
	    	
	     }     
	   if(track.minvalue < 0) {
		   track.negatives = true;
	   }
     if(track.getBBfileReader() == null && ((track.minvalue == 0 && track.maxvalue == 0) || track.minvalue == Double.MAX_VALUE)) {
		   track.hasvalues = false;
	   }
	   else {
		   track.hasvalues = true;
		   track.scale = Math.max(Math.abs(track.maxvalue), Math.abs(track.minvalue));
	   }	  
	  
	  	
	     addNode = null;
	     track.cleared = false;
	}
	catch(Exception e) {
		
		e.printStackTrace();
		
	}
}
void getBEDfeatures(BedTrack track, int start, int end) {
	
	if(track.url != null) {
		track.used = false;	
	}
	
	track.bedstart = start;
	track.bedend = end;
	 
	
	BedReader reader = new BedReader(track, start, end);	
	
	if(!Main.drawCanvas.loading) {
		
		reader.execute();
	}
	else {
	//if(Main.nothread) {
		reader.getNodes();
	}
	//}
	/*else {
		
		reader.execute();
	}*/
	
}
void removeBeds(BedTrack track) {
	
	track.nulled = true;
	BedNode addNode = track.getHead();
	try {
		while(addNode.getNext()!= null) {
			addNode = addNode.getNext();
			if(addNode == null) {
				break;
			}
			if(addNode.getPrev() != null) {
				addNode.getPrev().putNext(null);
				addNode.getPrev().putPrev(null);
			}			
		}
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	addNode = null;
	track.setCurrent(null);
	track.getHead().putNext(null);
	track.getBedLevels().clear();
	track.prepixel = 0;	
	track.cleared = true;
	Boolean ison = false;
	
	for(int i = 0 ; i<bedTrack.size();i++) {
		if(bedTrack.get(i).intersect) {
			ison = true;
			break;
		}
	}
	if(!ison) {
		Main.bedCanvas.bedOn = false;
	}
}

public void paint(Graphics g) {
	try {	
		
		drawScreen(g);
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}
public void pressIntersect(BedTrack track) {
	BedCanvas.annoTrack = track;
	
	if(track.intersect) {
		track.intersect = false;
		
		boolean bedon = false;
		for(int i = 0 ; i<bedTrack.size(); i++) {
			if(bedTrack.get(i).intersect) {
				bedon = true;
				break;
			}
		}
		if(!bedon) {
			bedOn = false;
		}
		
		//if(bedTrack.get(selectedPlay).small || (!bedTrack.get(selectedPlay).small && bedTrack.get(selectedPlay).getHead().getNext() != null && Main.drawCanvas.splits.get(0).viewLength < 1000000)) {
			annotate(track.getHead(), FileRead.head.getNext());
			annotator = false;
			BedCanvas.annoTrack = null;
		/*}
		else {
			
		}*/
		if(VariantHandler.commonSlider.getValue() > 1 && !VariantHandler.clusterBox.getText().equals("0") && !VariantHandler.clusterBox.getText().equals("")) {
	    	Main.drawCanvas.clusterId = 1;
			Main.drawCanvas.calcClusters(FileRead.head);
		}
		
	}
	else {
		
		track.intersect = true;
		/*if((track.small && track.getZoomlevel() == null)|| (!track.small && track.getHead().getNext() != null && Main.drawCanvas.splits.get(0).viewLength < Settings.windowSize)) {
			Annotator annotator = new Annotator(track);
			annotator.execute();
			/*Main.drawCanvas.loading("Annotating variants with " +track.file.getName());
			try {
				
				annotate(track.getHead());
				BedCanvas.annoTrack = null;
				annotator = false;
				//bedTrack.get(selectedPlay).used = true;
				
				if(VariantHandler.commonSlider.getValue() > 1 && !VariantHandler.clusterBox.getText().equals("0") && !VariantHandler.clusterBox.getText().equals("")) {
			    	Main.drawCanvas.clusterId = 1;
					Main.drawCanvas.calcClusters(FileRead.head);
				}
				
				VariantHandler.clusterTable.repaint();
				Draw.updatevars = true;
				Draw.calculateVars = true;
				Main.drawCanvas.ready("Annotating variants with " +track.file.getName());
				Main.drawCanvas.repaint();
			}
			catch(Exception e) {
				Main.drawCanvas.ready("Annotating variants with " +track.file.getName());
				e.printStackTrace();
			}*/
		/*}
		else {
			//TODO*/
			
			if(!track.used) {
				if(!bedOn) {
					VarNode current = FileRead.head.getNext();
					while(current != null) {	
						
						current.bedhit = true;	
						current = current.getNext();
					}
				}
				//System.out.println(Main.drawCanvas.loading +" " +Main.drawCanvas.loadingtext);
				bedOn = true;
				Annotator annotator = new Annotator(track);
				annotator.execute();
				
			}
			else {
				
				bedOn = true;
				annotate(track.getHead(), FileRead.head.getNext());
				BedCanvas.annoTrack = null;
				
			}
			
	//	}
	}
	
}
public void mouseClicked(MouseEvent event) {
	switch(event.getModifiers()) {	
		case InputEvent.BUTTON1_MASK: {	
			if(selectedPlay > -1 && bedTrack.get(selectedPlay).playbox.intersects(sideMouseRect)) {
				if(!Main.drawCanvas.loading) {
					
					pressIntersect(bedTrack.get(selectedPlay));			
				}
			}
			else if(selectedPlay > -1 && bedTrack.get(selectedPlay).graphBox.intersects(sideMouseRect)) {
				pressGraph(bedTrack.get(selectedPlay));
			}
			else if(removeTrack > -1) {
				removeTrack(removeTrack);
			}	
			break;
		}
		case InputEvent.BUTTON3_MASK: {	
			
			if(!sidebar) {
				this.bedTrack.get(hoverIndex).getPopup().show(this, mouseX, mouseY);
			}
		}
	}
	Draw.updatevars = true;
	repaint();	
	
}	

public void pressGraph(BedTrack track) {
	
	if(track.graph) {
		track.graph =false;
		track.getCollapseBox().setText("Auto collapse");
		repaint();
	}
	else {
		track.graph =true;
		track.getCollapseBox().setText("Auto scale");
		repaint();
	}
}
void removeBedhits(BedTrack remtrack) {
	boolean bedon = false;
	
	for(int i = 0 ; i<bedTrack.size(); i++) {
		if(bedTrack.get(i).intersect && !bedTrack.get(i).equals(remtrack)) {
			bedon = true;
		}
	}
	
	if(!bedon) {
		bedOn = false;
	}
	VarNode current = FileRead.head.getNext();
	while(current != null) {
		if(current.getBedHits() != null) {
			for(int i = 0; i<current.getBedHits().size(); i++) {
				if(current.getBedHits().get(i).getTrack().equals(remtrack)) {
					current.getBedHits().remove(i);
					i--;
				}
			}	
			if(current.getBedHits().size() == 0) {
				current.remBedhits();
				current.setBedhit(false);
			}
		}
		
		current = current.getNext();
	}
}

void removeTrack(int removeTrack) {
	if(removeTrack >=this.bedTrack.size() ) {
		removeTrack = this.bedTrack.size()-1;
	}
	BedTrack remtrack = this.bedTrack.get(removeTrack);
	remtrack.intersect = false;
	MethodLibrary.removeHeaderColumns(remtrack);
	removeBedhits(remtrack);
	remtrack.setDrawNode(null);
	boolean bedon = false;
	for(int i = 0 ; i<bedTrack.size(); i++) {
		if(bedTrack.get(i).intersect) {
			bedon = true;
		}
	}
	if(!bedon) {
		bedOn = false;
	}
	annotate(remtrack.getHead(), FileRead.head.getNext());
	if(VariantHandler.commonSlider.getValue() > 1 && !VariantHandler.clusterBox.getText().equals("0") && !VariantHandler.clusterBox.getText().equals("")) {
    	Main.drawCanvas.clusterId = 1;
		Main.drawCanvas.calcClusters(FileRead.head);
	}
	FileRead.removeTable(remtrack);
	this.bedTrack.remove(removeTrack);
	this.trackDivider.remove(removeTrack);

	if(this.bedTrack.size() == 0) {
		
		track = null;
		infoNode = null;
		preInfoNode = null;
		Main.bedScroll.setVisible(false);
		Main.trackPane.setDividerLocation(0);
		Main.trackPane.setDividerSize(0);
		if(!Main.controlScroll.isVisible()) {
			Main.varpane.setDividerSize(0);
			Main.trackPane.setVisible(false);
			Main.trackPane.setDividerSize(0);
			Main.varpane.setResizeWeight(0.0);   
		}
	}
	Boolean ison = false;
	
	for(int i = 0 ; i<bedTrack.size();i++) {
		bedTrack.get(i).trackIndex = i;
		if(bedTrack.get(i).intersect) {
			ison = true;
			
		}
	}
	if(!ison) {
		Main.bedCanvas.bedOn = false;
	}
	remtrack = null;
}
/*
void annotate(BedTrack track) {
	
	int varlength = 0;
	BedNode currentBed = track.getHead().getNext();
	Object[] t1 = new Object[Main.bedCanvas.bedTrack.size()];
	Object[] t2 = new Object[Main.bedCanvas.bedTrack.size()];
	for(int i = 0; i<t1.length; i++) {
		if(Main.bedCanvas.bedTrack.get(i).intersect) {
			t1[i] = 1;
		}
		else {
			t1[i] = null;
		}		
	}	

	if(currentBed == null) {
		return;
	}
	VarNode current = FileRead.head.getNext();	 
	int position = 0;	
	
	try {
	
	if(current != null) {
	
		while(currentBed != null) {
			
			if(current.getNext() == null && currentBed.getPosition() > current.getPosition()) {				
				break;
			}
			
			if(current.getPrev() == null && currentBed.getPosition()+currentBed.getLength() < current.getPosition()) {	
				currentBed = currentBed.getNext();
				continue;
			}
			
			if(current != null && current.getPrev() != null) {
				while(current.getPrev().getPosition() >= currentBed.getPosition()) {					
					if(current.getPrev() != null) {			
						
						current = current.getPrev();						
					}						
				}
			}
			
			if(current.indel) {
				position = current.getPosition()+1;	
				varlength = MethodLibrary.getBaseLength(current.vars, 1);
				
			}
			else {		
				position = current.getPosition();
				varlength = 0;
			}
		
			while(position < currentBed.getPosition()+currentBed.getLength()) {
				
				try {
					
				if(position+varlength >= currentBed.getPosition()) {
					
					if(!currentBed.getTrack().used) {	
						if(current.getBedHits() == null) {
							current.setBedhits();
						}
						
						current.getBedHits().add(currentBed);			
						current.setBedhit(true);
						
						if(currentBed.varnodes == null) {
							currentBed.varnodes = new ArrayList<VarNode>();
						}
						
						if(!currentBed.varnodes.contains(current)) {
							currentBed.varnodes.add(current);				
						}
					}					
				}
											
				if(current != null) {
					
					current.bedhit = checkIntersect(current, t1, t2);					
					current = current.getNext();
					
					if(current == null) {
						break;
					}
					position = current.getPosition();	
					if(current.indel) {
						position = current.getPosition()+1;	
						varlength = MethodLibrary.getBaseLength(current.vars, 1);
						
					}
					else {		
						position = current.getPosition();
						varlength = 0;
					}
					
				}
				else {
					break;
				}		
			}
			catch(Exception e) {
				
				e.printStackTrace();
				break;
			}
			}	
			if(current == null) {
				break;
			}
			currentBed = currentBed.getNext();
		}			
	}
	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
		
	}
	
	track.used = true;
	current = null;
	currentBed = null;
}
*/
boolean annotate(BedNode head, VarNode node) {
	
	if(FileRead.head.getNext() == null && head.getTrack().small) {
		head.getTrack().used = false;
		
		return false; 
	}
	//if(current.getPosition() == 36852366) {}
	
	int varlength = 0;
	BedNode currentBed = head.getNext();
	Object[] t2 = new Object[Main.bedCanvas.bedTrack.size()];
	Object[] t1 = new Object[Main.bedCanvas.bedTrack.size()];
	
	for(int i = 0; i<t1.length; i++) {
		if(!Main.bedCanvas.bedTrack.get(i).intersect) {
			continue;
		}
		if(Main.bedCanvas.bedTrack.get(i).getIntersectBox().isSelected() && !Main.bedCanvas.bedTrack.get(i).annotator) {			
			t1[i] = 1;		
		}				
	}	
	if(head.getTrack().used) {
		VarNode current = FileRead.head.getNext();	 
		while(current != null) {	
			
			current.bedhit = checkIntersect(current, t1, t2);	
			current = current.getNext();
		}
		
		Draw.updatevars = true;
		Draw.calculateVars = true;
		Main.drawCanvas.repaint();
		return true;
	}
	
	if(currentBed == null) {
		return false;
	}
	
	VarNode current = node;	 
	node = null;
	int position = 0;	
	
	try {
		
	if(current != null) {
		
		while(currentBed != null) {
			if(!Main.drawCanvas.loading) {
				break;
			}
			if(current.getNext() == null) {				
				if(currentBed.getPosition() > current.getPosition()+1) {
					break;
				}				
			}
			
			if(current.getPrev().getPrev() == null && currentBed.getPosition()+currentBed.getLength() < current.getPosition()) {	
				currentBed = currentBed.getNext();
				continue;
			}
			
			if(current != null && current.getPrev() != null) {
				while(current.getPrev().getPosition() >= currentBed.getPosition()) {					
					if(current.getPrev() != null) {						
						current = current.getPrev();						
					}						
					else  {
						
						current = current.getNext();						
						break;
					}
				}
			}			
			
			position = current.getPosition();	
			
			if(current.indel) {
				
				varlength = MethodLibrary.getBaseLength(current.vars, 1);				
			}
			else {		
				
				varlength = 0;
			}
			
			/*if(position == 34286255) {
				
				if(position+varlength >= currentBed.getPosition()) {
					System.out.println(position +" " +(currentBed.getPosition()+currentBed.getLength()) +" " +currentBed.getPosition());
				}
				
			}*/
			while(position < currentBed.getPosition()+currentBed.getLength()) {
				
				if(!Main.drawCanvas.loading) {
					break;
				}
				
				try {
					
					
					
					if(position+varlength >= currentBed.getPosition()) {
						
						if(!currentBed.getTrack().used) {				
							
							if(current.getBedHits() == null) {
								current.setBedhits();
							}
							found = true;
							if(currentBed.id != null || currentBed.getTrack().selex) {
							
								if(currentBed.getTrack().getAffinityBox().isSelected()) {
									if(currentBed.getTrack().limitValue != (double)Integer.MIN_VALUE && currentBed.getTrack().limitValue != null) {
										found = false;
										for(int i=0;i<current.vars.size();i++) {
											if(Math.abs(MethodLibrary.calcAffiniyChange(current, current.vars.get(i).getKey(), currentBed)) > Math.abs(currentBed.getTrack().limitValue)) {
												found = true;
												break;
											}
										}
										
									}
								}
							}
							if(found) {
								if(currentBed.getTrack().basecolumn != null) {
									for(int i = 0 ; i< current.vars.size(); i++) {										
										if(currentBed.name.equals(current.vars.get(i).getKey())) {
											current.getBedHits().add(currentBed);			
											current.setBedhit(true);
											
											if(currentBed.varnodes == null) {
												currentBed.varnodes = new ArrayList<VarNode>();
											}											
											if(!currentBed.varnodes.contains(current)) {												
												currentBed.varnodes.add(current);				
											}
										}
									}
								}
								else {
									
									current.getBedHits().add(currentBed);			
									current.setBedhit(true);
									
									if(currentBed.varnodes == null) {
										currentBed.varnodes = new ArrayList<VarNode>();
									}									
									if(!currentBed.varnodes.contains(current)) {										
										currentBed.varnodes.add(current);				
									}
								}
							}
						}							
					}
					
					if(current != null) {
						
						if(!head.getTrack().annotator) {
							
							current.bedhit = checkIntersect(current, t1, t2);		
						}
						if(current.getNext() != null) {
							current = current.getNext();
						}
						else {
							break;
						}
						
						position = current.getPosition();	
						//System.out.println(position +" " +(currentBed.getPosition()+currentBed.getLength()));
						if(current.indel) {
						//	position = current.getPosition()+1;	
							varlength = MethodLibrary.getBaseLength(current.vars, 1);
							
						}
						else {		
						//	position = current.getPosition();
							varlength = 0;
						}
						
					}
					else {
						
						break;
					}		
				}
				catch(Exception e) {
					
					e.printStackTrace();
					break;
				}
			}	
					
			currentBed = currentBed.getNext();
			
		}
		
	}
	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
		
	}
	
	if(!head.getTrack().annotator && !FileRead.bigcalc) {
		head.getTrack().used = true;
		while(current != null) {	
			
			current.bedhit = checkIntersect(current, t1, t2);	
			current = current.getNext();
		}
	}
	
	
	current = null;
	currentBed = null;
	
	return true;
}

public void checkIntersectAll(VarNode node, Object[] t1) {
	
	Object[] t2 = new Object[Main.bedCanvas.bedTrack.size()];
	
	for(int i = 0; i<t1.length; i++) {
		if(Main.bedCanvas.bedTrack.get(i).intersect) { 
			continue;			
		}
		if(Main.bedCanvas.bedTrack.get(i).getIntersectBox().isSelected()) {
			t1[i] = 1;
			continue;
		}
		if(Main.bedCanvas.bedTrack.get(i).getSubtracttBox().isSelected()) { 
			t1[i] = null;
			continue;
		}
	
	}	
	while(node != null) {
		
		node.bedhit = checkIntersect(node, t1, t2);
		
		node = node.getNext();
	}
	node = null;
}

boolean checkIntersect(VarNode current, Object[] t1, Object[] t2) {
	t2 = new Object[t1.length];
	
	if(current == null || current.getBedHits() == null) {
		for(int i = 0 ;i<this.bedTrack.size(); i++) {
			if(this.bedTrack.get(i).intersect && this.bedTrack.get(i).getIntersectBox().isSelected()) {
				return false;
			}				
		}
		return true;
	}

	for(int i = 0; i<current.getBedHits().size(); i++) {	
		if(!current.getBedHits().get(i).getTrack().intersect) {			
			continue;
		}
		if(current.getBedHits().get(i).getTrack().getSubtracttBox().isSelected()) {
			t2[current.getBedHits().get(i).getTrack().trackIndex] = 1;			
		}
		else if(current.getBedHits().get(i).getTrack().getIntersectBox().isSelected()) {					
			t2[current.getBedHits().get(i).getTrack().trackIndex] = 1;					
		}		
	}
	
	if(Arrays.deepEquals(t1, t2)) {		
		return true;	
	}
	else {		
		return false;
	}
}
public void mouseEntered(MouseEvent arg0) {}	
public void mouseExited(MouseEvent arg0) {}	
@SuppressWarnings("unchecked")
public void mousePressed(MouseEvent event) {
	resize =false;
	pressY = event.getY();
	mouseDrag = true;
	
	this.requestFocus();
	Main.drawCanvas.pressX = event.getX();
	
	Main.drawCanvas.tempDrag = Main.drawCanvas.pressX;
	switch(event.getModifiers()) {	
	
	case InputEvent.BUTTON1_MASK: {	
		
	//	if(this.settingsbutton > -1) {
			if(sidebar && this.hoverIndex < bedTrack.size() && this.sideMouseRect.intersects(bedTrack.get(this.hoverIndex).settingsButton)) {
				
				this.bedTrack.get(this.hoverIndex).getPopup().show(this, mouseX, mouseY);
				//settingsbutton = -1;
				break;
			}
		
	
		//}
		if(resizer) {
			preresize = mouseY;
			preresizer = this.trackDivider.get(this.resizeDivider)*this.getHeight();
			tempDivider = (ArrayList<Double>)trackDivider.clone();			
		}		
		
		break;
	}
	case InputEvent.BUTTON3_MASK: {
		
		this.zoomDrag = false;
		if(sidebar) {
			this.bedTrack.get(hoverIndex).getPopup().show(this, mouseX, mouseY);
		}
		
		
	}
}
	
}	

void calcScale(BedTrack track) {
	
	if(track.getCollapseBox().isSelected()) {
		Double min = (double)Integer.MAX_VALUE, max = (double)Integer.MIN_VALUE;
		BedNode node = null;
		if(track.getCurrent() != null) {
			node = track.getCurrent();
			
		}
		else {
			
			return;
		}
		
		while(node != null && node.getPosition() < Main.drawCanvas.splits.get(0).end) {
			try {
				if(node.value < min) {
					min = node.value;				
				}
				if(node.value > max) {
					max = node.value;
				}
			}
			catch(Exception e) {
				node = node.getNext();
				continue;
			}
			node = node.getNext();		
		}
		//track.minvalue = min;
		//track.maxvalue = max;
		if(track.minvalue < 0) {
			   track.negatives = true;
			   track.scale = Math.max(Math.abs(max), Math.abs(min));
			   
	    }
	    else {
	    	   track.scale = max;
	    }
		
		node = null;
	}
	else {
		if(track.minvalue < 0) {
			   track.negatives = true;
			   track.scale = Math.max(Math.abs(track.maxvalue), Math.abs(track.minvalue));
	    }
	    else {
	    	   track.scale = track.maxvalue;
	    }
	}
	
	repaint();
}

public void mouseReleased(MouseEvent event) {
	
	if(!sidebar) {
		if(zoomDrag) {			
			if(mouseX-Main.drawCanvas.pressX > 0) {			
				Main.drawCanvas.gotoPos(Main.drawCanvas.selectedSplit.start+(Main.drawCanvas.pressX-Main.drawCanvas.selectedSplit.offset)/Main.drawCanvas.selectedSplit.pixel, Main.drawCanvas.selectedSplit.start+(mouseX-Main.drawCanvas.selectedSplit.offset)/Main.drawCanvas.selectedSplit.pixel);
			}
		}
		bedFeatureFetcher fetch = new bedFeatureFetcher();
		fetch.execute();
		Main.chromDraw.updateExons = true;
		repaint();
		Main.drawCanvas.mouseDrag = false;
		Main.chromDraw.repaint();
		Draw.updatevars = true;
		Main.drawCanvas.repaint();
	}
	
	mouseDrag = false; 
	zoomDrag = false;	
	Main.drawCanvas.lineZoomer = false;
	lineZoomer = false;
	
	/*bedFeatureFetcher fetch = new bedFeatureFetcher();
	fetch.execute();*/
}	


public void getMoreBeds(BedTrack track) {
	try {
		
	if(Main.drawCanvas.splits.get(0).viewLength < Settings.windowSize && !track.small) {			
		
		if(!Main.drawCanvas.lineZoomer && !Main.drawCanvas.zoomDrag && !Main.drawCanvas.mouseDrag) {
			
			if(!FileRead.cancelfileread && ((Main.drawCanvas.splits.get(0).start < track.bedstart || Main.drawCanvas.splits.get(0).end > track.bedend) || track.nulled)) {		
				
				track.nulled = false;
				track.bedstart = (int)Main.drawCanvas.splits.get(0).start-1000;
				track.bedend = (int)Main.drawCanvas.splits.get(0).end+1000;
				
				getBEDfeatures(track, (int)track.bedstart, (int)track.bedend);
				track.cleared = false;
				if(track.getHead().getNext() != null) {
					
				track.setCurrent(track.getHead().getNext());
				track.setDrawNode(track.getCurrent());
				
				}
			}
		}
	}
	else if(track.getZoomlevel() != null && track.getBBfileReader() != null) {
	
		int zoom = 1;
		for(int i =2;i<track.getBBfileReader().getZoomLevels().getZoomHeaderCount();i++) {
			if(track.getBBfileReader().getZoomLevels().getZoomLevelHeader(i).getReductionLevel() < (Main.drawCanvas.splits.get(0).viewLength/(Main.drawCanvas.splits.get(0).pixel*Main.drawCanvas.splits.get(0).viewLength))) {
				zoom = i;
			}
			else {
				break;
			}
		}
		
		if(zoom != track.getZoomlevel() || ((Main.drawCanvas.splits.get(0).start < track.bedstart || Main.drawCanvas.splits.get(0).end > track.bedend) || track.nulled)) {
			
			track.setZoomlevel(zoom);
			if(!Main.drawCanvas.lineZoomer && !Main.drawCanvas.zoomDrag && !Main.drawCanvas.mouseDrag) {
				
			//	if(!FileRead.cancelfileread && ((Main.drawCanvas.splits.get(0).start < track.bedstart || Main.drawCanvas.splits.get(0).end > track.bedend) || track.nulled)) {		
				
					track.nulled = false;
					track.bedstart = (int)Main.drawCanvas.splits.get(0).start-1000;
					track.bedend = (int)Main.drawCanvas.splits.get(0).end+1000;					
					getBEDfeatures(track, (int)Main.drawCanvas.splits.get(0).start-1000, (int)Main.drawCanvas.splits.get(0).end+1000);
					track.cleared = false;
					if(track.getHead().getNext() != null) {						
						track.setCurrent(track.getHead().getNext());
						track.setDrawNode(track.getCurrent());
					
				//	}
				}
			}
		}
	}
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}

@SuppressWarnings("unchecked")
public void mouseDragged(MouseEvent event) {
	switch(event.getModifiers()) {	
	
	case InputEvent.BUTTON1_MASK: {	
		if(sidebar) {
			return;
		}
		this.mouseX = event.getX();
		this.mouseY = event.getY();		
		
		if(resizer && this.getHeight() > this.trackDivider.size()*20) {				
				
				this.trackDivider.set(this.resizeDivider, mouseY/(double)this.getHeight());
			
				if((positivelock || negative) && this.trackDivider.get(this.resizeDivider)*this.getHeight()-preresizer >= 0) {	
					positivelock = false;
					negative = false;
					preresize = mouseY;					
					tempDivider = (ArrayList<Double>)trackDivider.clone();
					
				}
				else if ((!negative || negativelock ) && this.trackDivider.get(this.resizeDivider)*this.getHeight()-preresizer < 0) {					
					negativelock = false;
					preresize = mouseY;					
					tempDivider = (ArrayList<Double>)trackDivider.clone();
					negative = true;
				}
				
				if(negativelock) { // || (this.trackDivider.get(this.resizeDivider+1)*this.getHeight()) -(this.trackDivider.get(this.resizeDivider)*this.getHeight()) < 20) {
					this.trackDivider.set(this.resizeDivider, (this.trackDivider.get(this.resizeDivider+1)*this.getHeight()-19)/(double)this.getHeight());
				}
				if(positivelock) { // || (this.trackDivider.get(this.resizeDivider+1)*this.getHeight()) -(this.trackDivider.get(this.resizeDivider)*this.getHeight()) < 20) {
					this.trackDivider.set(this.resizeDivider, 19/(double)this.getHeight());
				}
				if(this.trackDivider.get(this.resizeDivider)*this.getHeight()-preresizer < 0) {
					
					negative = true;
					positivelock = true;
					if(this.resizeDivider > 0) {
						
						for(int i =1 ; i<this.resizeDivider+1; i++) {								
							if((this.trackDivider.get(i)*this.getHeight())-(this.trackDivider.get(i-1)*this.getHeight()) < 20) {
								
								this.trackDivider.set(i,((this.trackDivider.get(i-1)*this.getHeight())+19)/(double)this.getHeight());	
							}
							else {
								positivelock = false;
								if(i != this.resizeDivider) {
									this.trackDivider.set(i,(this.tempDivider.get(i)/preresize)*mouseY);			
								}
							}
						}						
					}
					if((this.trackDivider.get(0)*this.getHeight()) >= 20) {
						this.trackDivider.set(0,(this.tempDivider.get(0)/preresize)*mouseY);			
						positivelock = false;
					}
					else {
						this.trackDivider.set(0, 19/(double)this.getHeight());
					}
					
				}
				else {
					
					negative = false;
					negativelock = true;
					if(this.resizeDivider < this.trackDivider.size()-1) {						
						
						for(int i =this.resizeDivider ; i<this.trackDivider.size()-1; i++) {
							
							if((this.trackDivider.get(i+1)*this.getHeight())-(this.trackDivider.get(i)*this.getHeight()) < 20) {
								
								this.trackDivider.set(i,((this.trackDivider.get(i+1)*this.getHeight())-19)/(double)this.getHeight());
							}
							else {
								
								negativelock = false;
								if(i !=this.resizeDivider) {
									try {
										this.trackDivider.set(i,1-((1-this.tempDivider.get(i))/(this.getHeight()-preresize))*(this.getHeight()-mouseY));
									}
									catch(Exception e) {
										
									//	e.printStackTrace();
									}
								}
							}							
						}		
						if(this.getHeight()-(this.trackDivider.get(this.trackDivider.size()-2)*this.getHeight()) >= 20) {
							negativelock = false;
						}
						else {
							this.trackDivider.set(this.trackDivider.size()-2, (this.getHeight()-19)/(double)this.getHeight());
						}
						
					}
				}				
				
				preresizer = this.trackDivider.get(this.resizeDivider)*this.getHeight();
				repaint();
				
			
		}
		else if(lineZoomer) {	
			
			 if (Main.drawCanvas.selectedSplit.start > 1 || Main.drawCanvas.selectedSplit.end < Main.drawCanvas.selectedSplit.chromEnd ) {					
				 Main.drawCanvas.gotoPos(Main.drawCanvas.selectedSplit.start-(Main.drawCanvas.tempDrag-mouseX)/Main.drawCanvas.selectedSplit.pixel*2, Main.drawCanvas.selectedSplit.end+(Main.drawCanvas.tempDrag-mouseX)/Main.drawCanvas.selectedSplit.pixel*2);
			 }
			
			Main.drawCanvas.tempDrag = mouseX;
			Main.chromDraw.updateExons = true;
			repaint();
			Main.chromDraw.repaint();
			Draw.updatevars = true;
			Main.drawCanvas.repaint();
			
		}
		else {				
			if(getCursor().getType() != Cursor.N_RESIZE_CURSOR) {
			zoomDrag = true;				
			
			repaint();
			}
			return;
		}
		break;
	}
	case InputEvent.BUTTON3_MASK: {	
		if(sidebar) {
			return;
		}
		if((int)Main.drawCanvas.selectedSplit.start == 1 && (int)Main.drawCanvas.selectedSplit.end == Main.drawCanvas.selectedSplit.chromEnd) {
			break;
		}			
	
		Main.drawCanvas.mouseDrag = true;
		Main.drawCanvas.moveX = event.getX();
		Main.drawCanvas.drag(Main.drawCanvas.moveX);
		break;
	}
	case 17: {
		if(sidebar) {
			return;
		}
		if((int)Main.drawCanvas.selectedSplit.start == 1 && (int)Main.drawCanvas.selectedSplit.end == Main.drawCanvas.selectedSplit.chromEnd) {
			break;
		}			
	
		Main.drawCanvas.mouseDrag = true;
		Main.drawCanvas.moveX = event.getX();
		Main.drawCanvas.drag(Main.drawCanvas.moveX);
		break;
	}
	}
}	
public void mouseMoved(MouseEvent event) {
	mouseY = event.getY();
	this.mouseX = event.getX();
	this.mouseRect.setBounds(event.getX()-Main.sidebarWidth, event.getY(), 1, 1);
	this.sideMouseRect.setBounds(event.getX(), event.getY(), 1, 1);
	
	if(!sidebar && event.getX() < Main.sidebarWidth) {
		sidebar = true;
	}
	else {
		if(sidebar && event.getX() >= Main.sidebarWidth) {
			sidebar = false;
		}
	}
	
	for(int i = 0; i<this.trackDivider.size(); i++) {
		if(this.trackDivider.get(i)*this.getHeight() < mouseY) {
			continue;
		}
		if(this.trackDivider.get(i)*this.getHeight() >= mouseY) {			
			if(hoverIndex != i) {
				hoverIndex = i;				
			}
			break;
		}	
	}
	repaint();
	
}


@Override
public void keyTyped(KeyEvent e) {
	// TODO Auto-generated method stub
	
}


@Override
public void keyPressed(KeyEvent e) {
	// TODO Auto-generated method stub
	int keyCode = e.getKeyCode();
	
	
	if(keyCode == KeyEvent.VK_PLUS || keyCode == 107) {
		heightchanged = true;
		for(int i = 0 ; i< bedTrack.size(); i++) {
			if(bedTrack.get(i).nodeHeight < 20) {		
				bedTrack.get(i).nodeHeight++;
				
			}
		}
		
		repaint();
	}
	if(keyCode == KeyEvent.VK_MINUS || keyCode == 109) {
		
		heightchanged = true;
		for(int i = 0 ; i< bedTrack.size(); i++) {
			if(bedTrack.get(i).nodeHeight > 1) {
				bedTrack.get(i).nodeHeight--;
			}
		}
		repaint();
	}
	if(keyCode == KeyEvent.VK_PAGE_UP) {
		
		if(this.forwardColor.getAlpha() < 255) {
			this.forwardColor = new Color(this.forwardColor.getRed(), this.forwardColor.getGreen(), this.forwardColor.getBlue(), this.forwardColor.getAlpha()+15);
			this.reverseColor = new Color(this.reverseColor.getRed(), this.reverseColor.getGreen(), this.reverseColor.getBlue(), this.reverseColor.getAlpha()+15);
			for(int i = 0; i<bedTrack.size(); i++) {
				
				if(bedTrack.get(i).getColors().size() > 0) {
					 Iterator<Map.Entry<Integer, Color>> it;
					 Map.Entry<Integer, Color> pair;
					Color color;
					it = bedTrack.get(i).getColors().entrySet().iterator();
				    while (it.hasNext()) {
				    	pair = (Map.Entry<Integer, Color>)it.next();
				    	color = pair.getValue();
				    	color = new Color(color.getRed(), color.getGreen(), color.getBlue(), color.getAlpha()+15);
				    	bedTrack.get(i).getColors().put(pair.getKey(), color);
				    }					
				}
			}
		
		}
		
		repaint();
	}
	if(keyCode == KeyEvent.VK_PAGE_DOWN) {
		
		if(this.forwardColor.getAlpha() > 0) {
			this.forwardColor = new Color(this.forwardColor.getRed(), this.forwardColor.getGreen(), this.forwardColor.getBlue(), this.forwardColor.getAlpha()-15);
			this.reverseColor = new Color(this.reverseColor.getRed(), this.reverseColor.getGreen(), this.reverseColor.getBlue(), this.reverseColor.getAlpha()-15);
			
			for(int i = 0; i<bedTrack.size(); i++) {
				if(bedTrack.get(i).getColors().size() > 0) {					
					Iterator<Map.Entry<Integer, Color>> it;
					Map.Entry<Integer, Color> pair;
					Color color;
					it = bedTrack.get(i).getColors().entrySet().iterator();
					
				    while (it.hasNext()) {
				    	
				    	pair = (Map.Entry<Integer, Color>)it.next();
				    	color = pair.getValue();
				    	
				    	color = new Color(color.getRed(), color.getGreen(), color.getBlue(), color.getAlpha()-15);
				    	bedTrack.get(i).getColors().put(pair.getKey(), color);
				    	
				    }					
				}
			}
		}
		
		repaint();
	}
}


@Override
public void keyReleased(KeyEvent e) {
	// TODO Auto-generated method stub
	
}
	
	
	
	
}
