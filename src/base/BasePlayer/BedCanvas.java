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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDCodec;
import base.BBfile.BBFileReader;


import base.BBfile.BedFeature;
import base.BBfile.BigBedIterator;
import base.BBfile.BigWigIterator;
import base.BBfile.WigItem;
import htsjdk.tribble.bed.BEDFeature;


import javax.swing.JPanel;
import javax.swing.SwingWorker;

public class BedCanvas extends JPanel implements MouseMotionListener, MouseListener, KeyListener {

private static final long serialVersionUID = 1L;
	Color forwardColor = new Color(171,194,171, 255), reverseColor = new Color(194,171,171, 255);
	BufferedImage bufImage, nodeImage;
	Graphics2D buf, nodebuf;
	int width, height, bedwidth;
//	private TabixReader tabixReader;
	Color zoomColor = new Color(115,115,100,100);
	//ArrayList<BedNode> bednodes = new ArrayList<BedNode>();
	ArrayList<String> infolist = new ArrayList<String>();
	Composite backupComposite;
	ArrayList<BedTrack> bedTrack = new ArrayList<BedTrack>();
	TabixReader.Iterator iterator = null;
	TabixReader.Iterator bedIterator = null;
	//BBFileReader bbfileReader = null;
	ArrayList<Double> trackDivider = new ArrayList<Double>();
	ArrayList<Double> tempDivider;
	BedNode drawNode;
	BedTrack track;
	Rectangle remoBox = new Rectangle();
	private int bedEndPos;
	private boolean foundlevel;
//	BEDCodec bedcodec = null;
	private int[][] matrix;
	private int sum;
	private int prevY;
	private Color bedcolor;
	private int bedcounter;
	private int trackstart;
	private FontMetrics fm;
	private Rectangle2D textWidth;
	private int lettercount;
	private int mouseY;
	private int hoverIndex;
	boolean bedOn = false;
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
	public boolean graph = false;
	private int bedheight=0;
	private int mouseX =0;
	private boolean zoomDrag=false;
	boolean lineZoomer = false;
	private String zoomtext= "";
	private int zoompostemp=0;
	private int pressY = 0;
	int nodeHeight = 10;

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
	
	BedCanvas(int width, int height) {	
		
		this.width = width;
		this.height = height;
		bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)Main.screenSize.getWidth(), (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
		buf = (Graphics2D)bufImage.getGraphics();
		buf.setRenderingHints(Draw.rh);
		nodeImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)Main.screenSize.getWidth(), (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));
		nodebuf = (Graphics2D)nodeImage.getGraphics();
		nodebuf.setRenderingHints(Draw.rh);
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
	
	if(this.trackDivider.get(this.trackDivider.size()-1) != 1.0) {	
		 for(int i = 0 ; i<Main.bedCanvas.trackDivider.size(); i++) {
			 Main.bedCanvas.trackDivider.set(i, ((i+1)*(this.getHeight()/(double)trackDivider.size())/this.getHeight()));
		 }		
	}
	
	drawSidebar();
	buf.setStroke(Draw.doubleStroke);
	
	buf.drawImage(Draw.image, Main.sidebarWidth, 0, this.getWidth(),(int)Main.screenSize.getHeight(), this);
	buf.setColor(Draw.backColor);
	buf.fillRect(Main.sidebarWidth, 0, this.getWidth(), nodeImage.getHeight());	
	buf.setColor(Color.gray);
	
	if(!zoomDrag && !resize) {		
		drawNodes();		
	}	
	
	buf.setColor(Color.black);
	buf.drawLine(Main.sidebarWidth-2, 0, Main.sidebarWidth-2, this.getHeight());	
	
	if(resize) {		
		buf.drawImage(nodeImage, Main.sidebarWidth, 0, nodeImage.getWidth(),(int)(Main.vardivider*Main.varPaneDivider.getY()), null);
	}
	else {		
		buf.drawImage(nodeImage, Main.sidebarWidth, 0, null);
	}
	if(resizer && !mouseDrag) {
		resizer = false;
	}
	for(int i = 0 ; i<bedTrack.size(); i++) {
		if(i <bedTrack.size()-1) {
			buf.drawLine(0, (int)(trackDivider.get(i)*this.getHeight()), this.getWidth(), (int)(trackDivider.get(i)*this.getHeight()));
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
			buf.drawString("Scale [" +MethodLibrary.round(bedTrack.get(i).minvalue, 2) +", " +MethodLibrary.round(bedTrack.get(i).maxvalue,2) +"]", Main.sidebarWidth+10, (int)(trackDivider.get(i)*this.getHeight())-5);
			buf.setColor(Color.black);
		}
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
		buf.drawLine((Main.drawCanvas.getDrawWidth())/2+Main.sidebarWidth+2, 0, ((Main.drawCanvas.getDrawWidth()))/2+Main.sidebarWidth+2, Main.bedScroll.getViewport().getHeight());
		buf.drawLine((int)((Main.drawCanvas.getDrawWidth())/2+Main.drawCanvas.splits.get(0).pixel+Main.sidebarWidth+2), 0, (int)(((Main.drawCanvas.getDrawWidth()))/2+Main.drawCanvas.splits.get(0).pixel+Main.sidebarWidth+2), Main.bedScroll.getViewport().getHeight());
		buf.setStroke(Draw.basicStroke);
	}	
	drawZoom();
	g.drawImage(bufImage, 0, 0, null);	
}	

void drawSidebar() {
	
	buf.setColor(Draw.sidecolor);
	buf.fillRect(0, 0, Main.sidebarWidth, this.getHeight());
	buf.setColor(Draw.softColor);
	buf.fillRect(0, 0, Main.sidebarWidth, this.getHeight());
	buf.setColor(Color.black);
	
	if(bedTrack.size() > 0) {
	//	buf.setStroke(Draw.doubleStroke);
		overlapping = false;
		for(int i = 0; i<bedTrack.size(); i++) {
			if(!buf.getColor().equals(Color.black)) {
				buf.setColor(Color.black);
			}
			if(i == 0) {
				trackstart = 10;
				trackheight = (int)(trackDivider.get(i)*this.getHeight());
			}
			else {
				
				trackstart = 10+(int)(trackDivider.get(i-1)*this.getHeight());
				trackheight = (int)(trackDivider.get(i)*this.getHeight())-(int)(trackDivider.get(i-1)*this.getHeight());
			}
			
			if(bedTrack.get(i).file != null) {
				buf.drawString(bedTrack.get(i).file.getName(), 10, trackstart);
			}
			else {
				buf.drawString(bedTrack.get(i).url.getFile(), 10, trackstart);
			}
			
			if(trackheight > 40) {
				if((int)bedTrack.get(i).playbox.getBounds().getY() != trackstart+10) {
					
					bedTrack.get(i).playbox.setBounds(10, trackstart+10, 20, 20);
					bedTrack.get(i).playTriangle.reset();
					
					bedTrack.get(i).playTriangle.addPoint(14, trackstart+14);
					bedTrack.get(i).playTriangle.addPoint(14, trackstart+26);
					bedTrack.get(i).playTriangle.addPoint(26, trackstart+20);						
					bedTrack.get(i).graphBox.setBounds(bedTrack.get(i).playbox.x+30, trackstart+10, 20, 20);
				}
				if(Main.varsamples > 0) {
					if(bedTrack.get(i).getCurrent() != null) {
						buf.setColor(Color.lightGray);
						buf.fillRect(bedTrack.get(i).playbox.getBounds().x,bedTrack.get(i).playbox.getBounds().y,bedTrack.get(i).playbox.getBounds().width,bedTrack.get(i).playbox.getBounds().height );
						
						if(bedTrack.get(i).intersect) {
							buf.setColor(Color.green);
							buf.fillRect(bedTrack.get(i).playTriangle.getBounds().x, bedTrack.get(i).playTriangle.getBounds().y, 12,12);
						}
						else {
							buf.setColor(Color.red);
							buf.fillPolygon(bedTrack.get(i).playTriangle);
						}
						if(sideMouseRect.intersects(bedTrack.get(i).playbox)) {
							
							overlapping = true;
							selectedPlay = i;
							if(getCursor().getType() != Cursor.HAND_CURSOR) {
								
								setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
							}	
							buf.setColor(Color.white);
							buf.drawRect(bedTrack.get(i).playbox.getBounds().x-1,bedTrack.get(i).playbox.getBounds().y-1,bedTrack.get(i).playbox.getBounds().width+1, bedTrack.get(i).playbox.getBounds().height+1);
							
						}
					}
					else {
						bedTrack.get(i).intersect = false;
					}					
				}
				
				if(bedTrack.get(i).maxvalue != 0 || bedTrack.get(i).minvalue != 0) {
					if(bedTrack.get(i).getCurrent() != null) {
					buf.setColor(Color.lightGray);
					buf.fillRect(bedTrack.get(i).graphBox.getBounds().x,bedTrack.get(i).graphBox.getBounds().y,bedTrack.get(i).graphBox.getBounds().width,bedTrack.get(i).graphBox.getBounds().height );
					
					if(bedTrack.get(i).graph) {
						buf.setColor(Color.white);
						buf.drawLine(bedTrack.get(i).graphBox.getBounds().x, bedTrack.get(i).graphBox.getBounds().y+20, bedTrack.get(i).graphBox.getBounds().x+4, bedTrack.get(i).graphBox.getBounds().y+10);
						buf.drawLine(bedTrack.get(i).graphBox.getBounds().x+4, bedTrack.get(i).graphBox.getBounds().y+10, bedTrack.get(i).graphBox.getBounds().x+10, bedTrack.get(i).graphBox.getBounds().y+16);
						buf.drawLine(bedTrack.get(i).graphBox.getBounds().x+10, bedTrack.get(i).graphBox.getBounds().y+16, bedTrack.get(i).graphBox.getBounds().x+16, bedTrack.get(i).graphBox.getBounds().y+6);
						buf.drawLine(bedTrack.get(i).graphBox.getBounds().x+16, bedTrack.get(i).graphBox.getBounds().y+6, bedTrack.get(i).graphBox.getBounds().x+20, bedTrack.get(i).graphBox.getBounds().y+12);
						
					}
					else {
						buf.setColor(Color.black);
						buf.drawLine(bedTrack.get(i).graphBox.getBounds().x, bedTrack.get(i).graphBox.getBounds().y+20, bedTrack.get(i).graphBox.getBounds().x+4, bedTrack.get(i).graphBox.getBounds().y+10);
						buf.drawLine(bedTrack.get(i).graphBox.getBounds().x+4, bedTrack.get(i).graphBox.getBounds().y+10, bedTrack.get(i).graphBox.getBounds().x+10, bedTrack.get(i).graphBox.getBounds().y+16);
						buf.drawLine(bedTrack.get(i).graphBox.getBounds().x+10, bedTrack.get(i).graphBox.getBounds().y+16, bedTrack.get(i).graphBox.getBounds().x+16, bedTrack.get(i).graphBox.getBounds().y+6);
						buf.drawLine(bedTrack.get(i).graphBox.getBounds().x+16, bedTrack.get(i).graphBox.getBounds().y+6, bedTrack.get(i).graphBox.getBounds().x+20, bedTrack.get(i).graphBox.getBounds().y+12);
						
					
					}
					
					if(sideMouseRect.intersects(bedTrack.get(i).graphBox)) {
						
						overlapping = true;
						selectedPlay = i;
						if(getCursor().getType() != Cursor.HAND_CURSOR) {
							
							setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
						}	
						buf.setColor(Color.white);
						buf.drawRect(bedTrack.get(i).graphBox.getBounds().x-1,bedTrack.get(i).graphBox.getBounds().y-1,bedTrack.get(i).graphBox.getBounds().width+1, bedTrack.get(i).graphBox.getBounds().height+1);
						
					}
					}
				}
				
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
			buf.setStroke(Draw.basicStroke);		
		
			if(sidebar && sideMouseRect.intersects(this.remoBox)) {				
			//	buf.setStroke(Draw.doubleStroke);
				removeTrack  = hoverIndex;
				buf.setColor(Color.white);
				buf.fillRect(this.remoBox.x-2, this.remoBox.y-1, this.remoBox.width+2, this.remoBox.height+2);
			}
			else {				
				removeTrack = -1;
				if(this.remoBox.getBounds().y != (int)(trackDivider.get(hoverIndex)*this.getHeight())-11) {
					this.remoBox.setBounds(Main.sidebarWidth-11, (int)(trackDivider.get(hoverIndex)*this.getHeight())-11, 8, 8);
				}
			}
		
			buf.setColor(Color.black);
			buf.drawRect(this.remoBox.x-1, this.remoBox.y, this.remoBox.width, this.remoBox.height);
			buf.drawLine(this.remoBox.x-1,this.remoBox.y,this.remoBox.x+7, this.remoBox.y +8);
			buf.drawLine(this.remoBox.x-1,this.remoBox.y+8,this.remoBox.x+7, this.remoBox.y);
		
		}
	}
	
}

void drawZoom() {	
	
	if(lineZoomer ) {
		
		buf.setColor(Color.black);
		buf.setStroke(Draw.dashed);
		
		buf.drawLine(Main.drawCanvas.pressX, pressY, mouseX, mouseY);
	
		
	}
	else if(zoomDrag) {
		
		
		buf.setStroke(Draw.dashed);
		buf.setColor(Color.black);
		buf.setFont(ChromDraw.seqFont);
		if(this.mouseX-Main.drawCanvas.pressX >= 0) {
			buf.drawRect(Main.drawCanvas.pressX, 0, this.mouseX-Main.drawCanvas.pressX, this.getHeight());	
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
	buf.setStroke(Draw.doubleStroke);
}

void drawLogo(BedTrack track, BedNode node, int trackstart) {	
	
	if(Main.SELEXhash.containsKey(node.id)) {	
		matrix = Main.SELEXhash.get(node.id);
		
		if((trackstart +this.trackheight) - (trackstart+15+((node.level-1)*(selexheight+15))+(track.mouseWheel*((selexheight+15)*2))) < selexheight) {
			selexheight = (trackstart +this.trackheight) - (trackstart+15+((node.level-1)*(selexheight+15))+(track.mouseWheel*((selexheight+15)*2)));
		}		
		
		
		if(matrix != null && (trackstart +this.trackheight) > trackstart+15+((node.level-1)*(40+15))+(track.mouseWheel*((40+15)*2))) {
			
			sum = 0;
			if(node.forward == null || node.forward) {
				for(int j = 0; j<matrix[0].length; j++) {
					prevY = trackstart+15+((node.level-1)*(40+15))+(track.mouseWheel*((40+15)*2));
					sum = matrix[0][j] + matrix[1][j] +matrix[2][j] +matrix[3][j];
				//	System.out.print((100*(matrix[0][j]/(double)sum)) +"\t" +(100*(matrix[1][j]/(double)sum)) +"\t" +(100*(matrix[2][j]/(double)sum)) +"\t" +(100*(matrix[3][j]/(double)sum)) +"\n");
					nodebuf.drawImage(Main.A, (int)(((node.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+(j*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[0][j]/(double)sum)), null);
					prevY += (int)(selexheight*(matrix[0][j]/(double)sum));
					nodebuf.drawImage(Main.C, (int)(((node.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+(j*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[1][j]/(double)sum)), null);
					prevY += (int)(selexheight*(matrix[1][j]/(double)sum));
					nodebuf.drawImage(Main.G, (int)(((node.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+(j*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[2][j]/(double)sum)), null);
					prevY += (int)(selexheight*(matrix[2][j]/(double)sum));
					nodebuf.drawImage(Main.T, (int)(((node.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+(j*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[3][j]/(double)sum)), null);
				}
			}
			else {
				bedcounter = 0;
				for(int j = matrix[0].length-1; j>=0; j--) {
					prevY = trackstart+15+((node.level-1)*(40+15))+(track.mouseWheel*((40+15)*2));
					sum = matrix[0][j] + matrix[1][j] +matrix[2][j] +matrix[3][j];
				//	System.out.print((100*(matrix[0][j]/(double)sum)) +"\t" +(100*(matrix[1][j]/(double)sum)) +"\t" +(100*(matrix[2][j]/(double)sum)) +"\t" +(100*(matrix[3][j]/(double)sum)) +"\n");
					nodebuf.drawImage(Main.A, (int)(((node.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+(bedcounter*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[3][j]/(double)sum)), null);
					prevY += (int)(selexheight*(matrix[3][j]/(double)sum));
					nodebuf.drawImage(Main.C, (int)(((node.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+(bedcounter*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[2][j]/(double)sum)), null);
					prevY += (int)(selexheight*(matrix[2][j]/(double)sum));
					nodebuf.drawImage(Main.G, (int)(((node.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+(bedcounter*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[1][j]/(double)sum)), null);
					prevY += (int)(selexheight*(matrix[1][j]/(double)sum));
					nodebuf.drawImage(Main.T, (int)(((node.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+(bedcounter*Main.drawCanvas.splits.get(0).pixel)), prevY, (int)Main.drawCanvas.splits.get(0).pixel, (int)(selexheight*(matrix[0][j]/(double)sum)), null);
					bedcounter++;
				}
			}
		}
	}	
	if(selexheight != 40) {
		selexheight = 40;
	}
	
}

void drawNodes() {
	
	//nodebuf.setColor(Draw.backColor);
	nodebuf.setComposite( Main.drawCanvas.composite);		
	nodebuf.fillRect(0,0, bufImage.getWidth(),this.getHeight());	
	nodebuf.setComposite(this.backupComposite);	
	
	for(int i = 0; i<bedTrack.size(); i++) {
		track = bedTrack.get(i);
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
		
//		nodebuf.drawImage(Draw.image, 0, trackstart-5, bufImage.getWidth(),Main.bedScroll.getViewport().getHeight()/bedTrack.size()+1, this);
//		nodebuf.setColor(Draw.loadColor);
//		nodebuf.fillRect(0, trackstart-5, bufImage.getWidth(), Main.bedScroll.getViewport().getHeight()/bedTrack.size()+1);
		
		if(Main.drawCanvas.splits.get(0).viewLength < 1000000 && !track.small) {
			if(!Main.drawCanvas.lineZoomer && !Main.drawCanvas.zoomDrag && !Main.drawCanvas.mouseDrag) {
				
				if(Main.drawCanvas.splits.get(0).start < track.bedstart || Main.drawCanvas.splits.get(0).end > track.bedend || track.nulled) {		
					track.nulled = false;
					track.bedstart = (int)Main.drawCanvas.splits.get(0).start-1000;
					track.bedend = (int)Main.drawCanvas.splits.get(0).end+1000;
					getBEDfeatures(track, (int)Main.drawCanvas.splits.get(0).start-1000, (int)Main.drawCanvas.splits.get(0).end+1000);
					track.cleared = false;
					if(track.getHead().getNext() == null) {
						continue;
					}
				}
			}
		}
		else if(!track.small) {
			if(!track.cleared) {
				removeBeds(track);
			}
			nodebuf.setColor(zoomColor);
			nodebuf.setFont(new Font("SansSerif", Font.BOLD, 50));
			nodebuf.drawString("ZOOM IN < 1Mbp", 100, trackstart+50);
			nodebuf.setFont(Draw.defaultFont);
			continue;
		}
		
		if(track.getHead().getNext() == null) {
			continue;
		}
		if(track.getCurrent().getPrev() != null && track.getCurrent().getPosition() > Main.drawCanvas.splits.get(0).start) {
			while(track.getCurrent().getPosition() > Main.drawCanvas.splits.get(0).start && !track.getCurrent().getPrev().equals(track.getHead())) {
				track.setCurrent(track.getCurrent().getPrev());
			}
		}
		else {
			while(track.getCurrent() != null && track.getCurrent().getPosition()+track.getCurrent().getLength() < Main.drawCanvas.splits.get(0).start && track.getCurrent().getNext() != null && track.getCurrent().getNext().getPosition()+track.getCurrent().getLength() < Main.drawCanvas.splits.get(0).start ) {
				track.setCurrent(track.getCurrent().getNext());
			}
		}
		drawNode = track.getCurrent();
		overlap = false;
		while(drawNode != null && drawNode.getPosition() < Main.drawCanvas.splits.get(0).end) {
			bedwidth = (int)(drawNode.getLength()*Main.drawCanvas.selectedSplit.pixel);
			if(bedwidth < 1) {
				bedwidth = 1;
			}
			if(!track.getColors().isEmpty()) {
				bedcolor = track.getColors().get(drawNode.color);
			}
			else if(drawNode.forward == null) {
				bedcolor = forwardColor;
			}
			else if(drawNode.forward) {
				bedcolor = forwardColor;
				
			}
			else {
				bedcolor= reverseColor;
			}
			nodebuf.setColor(bedcolor);
			if(track.graph) {
				if(!track.negatives) {
					if(drawNode.value == null || Double.isNaN(drawNode.value)) {
						track.prepixel = (int)((drawNode.getPosition()-Main.drawCanvas.selectedSplit.start)*Main.drawCanvas.selectedSplit.pixel);
						testRect.setBounds(track.prepixel,trackstart,bedwidth,nodeHeight);
						
						nodebuf.fillRect(track.prepixel, trackstart,bedwidth,nodeHeight);	
					}
					else {
						track.prepixel = (int)((drawNode.getPosition()-Main.drawCanvas.selectedSplit.start)*Main.drawCanvas.selectedSplit.pixel);
						bedheight = (int)((drawNode.value/(track.scale))*this.trackheight-10);
						testRect.setBounds(track.prepixel,(int)(trackDivider.get(i)*this.getHeight())-bedheight,bedwidth,bedheight);
						nodebuf.fillRect(track.prepixel,testRect.y,bedwidth,bedheight);	
						
					}
				}
				else if(drawNode.value!=null) {
					track.prepixel = (int)((drawNode.getPosition()-Main.drawCanvas.selectedSplit.start)*Main.drawCanvas.selectedSplit.pixel);
					if(!Double.isNaN(drawNode.value) && drawNode.value >= 0) {
						bedheight = (int)((drawNode.value/(track.scale*2))*this.trackheight-10);
						
						if(drawNode.value == 0) {
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
						
						bedheight = (int)((Math.abs(drawNode.value)/(track.scale*2))*this.trackheight-10);
						
						if(drawNode.value == 0) {
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
				if(mouseRect.intersects(testRect)) {					
					nodebuf.setColor(Color.white);
					nodebuf.fillRect(testRect.getBounds().x, testRect.getBounds().y,(int)testRect.getWidth(), (int)testRect.getHeight());
					if(infoNode != drawNode) {
						infoNode = drawNode;
					}
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
						bedEndPos = drawNode.getPosition()+drawNode.getLength();
																
						if(track.getBedLevels().isEmpty()) {
						
							track.getBedLevels().add(bedEndPos);
							drawNode.level = 1;
					//		level = 1;
						}
						else {
							
							foundlevel = false;
							for(int c=0;c<track.getBedLevels().size(); c++) {
								
								if(track.getBedLevels().get(c) < drawNode.getPosition()) {
									
								//	level = c+1;
									drawNode.level = c+1;
									foundlevel = true;				
									track.getBedLevels().set(c, bedEndPos);
									break;
								}
							
							}
						if(!foundlevel) {			
							track.getBedLevels().add(bedEndPos);
							//	level = track.bedLevelMatrix.size();
								drawNode.level = track.getBedLevels().size();
							}
						}					
						
					if(Main.drawCanvas.splits.get(0).viewLength <= 300 && track.selex) {
						
						if((drawNode.level-1)*(selexheight+15)+(track.mouseWheel*((selexheight+15)*2))+5 >= 5 && (drawNode.level-1)*(selexheight+15)+(track.mouseWheel*((selexheight+15)*2))+15 < this.trackheight) {
							testRect.setBounds((int)((drawNode.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel),trackstart+((drawNode.level-1)*(selexheight+15))+(track.mouseWheel*((selexheight+15)*2)),bedwidth,10);
							if(testRect.intersects(mouseRect)) {
								nodebuf.setColor(Color.white);
								nodebuf.fillRect(testRect.x,testRect.y,bedwidth,10);
								if(infoNode != drawNode) {
									infoNode = drawNode;
								}
								overlap = true;
							}
							else {
								nodebuf.fillRect(testRect.x,testRect.y,bedwidth,10);
							}
							
							
							if(drawNode.name != null) {
								fm = nodebuf.getFontMetrics();
								textWidth = fm.getStringBounds(drawNode.name, nodebuf);							
								nodebuf.setColor(Color.black);
								
								if(bedwidth > textWidth.getWidth()) {
									nodebuf.drawString(drawNode.name, (int)((drawNode.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel),trackstart+(drawNode.level*(selexheight+15))-(selexheight+7)+(track.mouseWheel*((selexheight+15)*2))+1);
								}
								else {
									lettercount = (int)(bedwidth/fm.getStringBounds("R", nodebuf).getWidth());
									nodebuf.drawString(drawNode.name.substring(0, lettercount), (int)((drawNode.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel),trackstart+(drawNode.level*(selexheight+15))-(selexheight+7)+(track.mouseWheel*((selexheight+15)*2))+1);
									
								}						
							}
							
							drawLogo(track,drawNode, trackstart);
							
						}
					}
					else {
						if(((drawNode.level-1)*(nodeHeight+(nodeHeight/2))+(track.mouseWheel*75)+(nodeHeight*2) > this.trackheight)) {
							drawNode = drawNode.getNext();
							continue;
						}
						
						if((drawNode.level-1)*(nodeHeight+(nodeHeight/2))+(track.mouseWheel*75)+5 >= 5) {
							testRect.setBounds((int)((drawNode.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel), trackstart+((drawNode.level-1)*(nodeHeight+(nodeHeight/2)))+(track.mouseWheel*75),bedwidth, nodeHeight);
							if(testRect.intersects(mouseRect)) {
								nodebuf.setColor(Color.white);
								nodebuf.fillRect(testRect.x, testRect.y,bedwidth, nodeHeight);	
								if(infoNode != drawNode) {
									infoNode = drawNode;
								}
								overlap = true;
								
								
							}
							else {
								nodebuf.fillRect(testRect.x, testRect.y,bedwidth, nodeHeight);
							}
						
			//				nodebuf.fillRect((int)((drawNode.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel), trackstart+((drawNode.level-1)*15)+(track.mouseWheel*75),bedwidth, 10);
						
						if(drawNode.name != null && testRect.width > 10) {
							fm = nodebuf.getFontMetrics();
							textWidth = fm.getStringBounds(drawNode.name, nodebuf);							
							nodebuf.setColor(Color.black);
							if((int)((drawNode.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel) < 0 && testRect.x+testRect.width > 0) {
								if(bedwidth > textWidth.getWidth()) {
									nodebuf.drawString(drawNode.name, 0,testRect.getBounds().y+(nodeHeight)-1);
								}
								else {
									lettercount = (int)(bedwidth/fm.getStringBounds("R", nodebuf).getWidth());
									nodebuf.drawString(drawNode.name.substring(0, lettercount),0,testRect.getBounds().y+(nodeHeight)-1);
									
								}
							}
							else {
								if(bedwidth > textWidth.getWidth()) {
									nodebuf.drawString(drawNode.name, (int)((drawNode.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel),testRect.getBounds().y+(nodeHeight)-1);
								}
								else {
									lettercount = (int)(bedwidth/fm.getStringBounds("R", nodebuf).getWidth());
									nodebuf.drawString(drawNode.name.substring(0, lettercount), (int)((drawNode.getPosition()-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel),testRect.getBounds().y+(nodeHeight)-1);
									
								}
							}
						}
					}
				}
	//		}
			}
			drawNode = drawNode.getNext();
		}
		
		if(overlap ) {
			drawInfo();
		}
		drawNode = null;
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

void drawInfo() {
	
//	if(boxmetrics.width != 0) {
	
		if(!infoNode.equals(preInfoNode)) {	
			
			preInfoNode = infoNode;
			infolist.clear();
			int maxlength = 0;
			if(infoNode.name != null) {
				infolist.add(infoNode.name);
				if(infoNode.name.length() > infolist.get(maxlength).length()) {
					maxlength = infolist.size()-1;
				}
			}
			infolist.add("Position: " +Main.drawCanvas.selectedSplit.chrom +":" +infoNode.getPosition());
			if(infolist.get(infolist.size()-1).length() > infolist.get(maxlength).length()) {
				maxlength = infolist.size()-1;
			}
			infolist.add("Length: " +infoNode.getLength());
			if(infolist.get(infolist.size()-1).length() > infolist.get(maxlength).length()) {
				maxlength = infolist.size()-1;
			}
			if(infoNode.value != null) {
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
			if(fm== null) {
				fm = nodebuf.getFontMetrics();
			}
			boxmetrics = fm.getStringBounds(infolist.get(maxlength), nodebuf).getBounds();				
			
		}
		
		nodebuf.setColor(Draw.loadColor);
		nodebuf.fillRect(10, trackstart+10, boxmetrics.width+8, boxmetrics.height*infolist.size());
		nodebuf.setColor(Color.black);
		for(int i = 0; i<infolist.size(); i++) {
			nodebuf.drawString(infolist.get(i), 14, (trackstart+6)+(boxmetrics.height*(i+1)));
		}	
	//}
}

public class BedReader extends SwingWorker<String, Object> {
	int start, end;
	BedTrack track;
	
	public BedReader(BedTrack track, int start, int end) {
		this.start = start;
		if(this.start < 0) {
			this.start = 0;
		}
		this.end = end;
		this.track = track;
		track.loading = true;
	}

	protected String doInBackground() {
		
		Main.drawCanvas.loading("Reading BED-file");
		try {		
			SeekableStream stream = null;
			TabixReader tabixReader = null;
			BBFileReader bbfileReader = null;
			if(track.file != null) {
				if(track.file.getName().toLowerCase().endsWith(".bb") || track.file.getName().toLowerCase().endsWith(".bigwig") || track.file.getName().toLowerCase().endsWith(".bigbed")|| track.file.getName().toLowerCase().endsWith(".bw")) {
					bbfileReader = new BBFileReader(track.file.getCanonicalPath(), track);
					 if(track.first) {
						 if(bbfileReader.getChromosomeName(0).contains("chr")) {
							 track.chr = "chr";
						 }
					
						 track.first = false;
					 }
						
				}
				else {
					
					tabixReader = new TabixReader(track.file.getCanonicalPath());  
				}
			}
			else {
				
				stream = SeekableStreamFactory.getInstance().getStreamFor(track.url);
				
				if(track.index != null) {
					tabixReader = new TabixReader(track.url.toString(), track.index.toString(), stream);
				}
				else {
					try {
						
						bbfileReader = new BBFileReader(track.url.toString(), stream, track);
						 if(track.first) {
							 if(bbfileReader.getChromosomeName(0).contains("chr")) {
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
			
			TabixReader.Iterator bedIterator = null;
			BigWigIterator wigiter = null;
			BigBedIterator bediter = null;
			 
				 if(tabixReader == null) {
					
					 if(bbfileReader.isBigWigFile()) {
						
						 wigiter = bbfileReader.getBigWigIterator(track.chr+Main.chromosomeDropdown.getSelectedItem().toString() , start, track.chr+Main.chromosomeDropdown.getSelectedItem().toString(), end, false);
					 }
					 else {
						 bediter  = bbfileReader.getBigBedIterator(track.chr+Main.chromosomeDropdown.getSelectedItem().toString() , start, track.chr+Main.chromosomeDropdown.getSelectedItem().toString(), end, false);
					 }
					
					
					 
				 }
				 else {
					 if(track.first) {
						 if(track.chr == null) {
							 track.chr = "";
						 }
						 String chr = tabixReader.getChromosomes().iterator().next();
						 if(chr.contains("chr")) {
							
							 track.chr = "chr";
						 }
						 
						 track.first = false;
						
					 }					
					try {
					 bedIterator = tabixReader.query(track.chr +Main.chromosomeDropdown.getSelectedItem().toString() +":" +start +"-" +end);
					}
					catch(Exception e) {
						System.out.println(track.chr +Main.chromosomeDropdown.getSelectedItem().toString() +":" +start +"-" +end);
						e.printStackTrace();
					}
				 }
				/*
				  if((wigiter != null && !wigiter.hasNext()) || (bedIterator == null) || bedIterator.next() == null) {
					  try {
						  	if(tabixReader == null) {
						  		if(bbfileReader.isBigWigFile()) {
						  			
						  			wigiter = bbfileReader.getBigWigIterator("chr" +Main.chromosomeDropdown.getSelectedItem().toString() , start, "chr" +Main.chromosomeDropdown.getSelectedItem().toString(), end, true);
						  			
						  		}
						  		else {
						  			
						  			
						  			 bediter  = bbfileReader.getBigBedIterator("chr" +Main.chromosomeDropdown.getSelectedItem().toString() , start,"chr" + Main.chromosomeDropdown.getSelectedItem().toString(), end, true);
										
						  		}
							 }
							 else {
								 bedIterator = tabixReader.query("chr" +Main.chromosomeDropdown.getSelectedItem().toString()+":" +start +"-" +end );
							 }
					  }
					  catch(Exception ex) {
						  ex.printStackTrace();
						  ErrorLog.addError(ex.getStackTrace());
						
					  }
				  }
			 */
			  
			if(track.file != null) {  
				if(track.file.getName().toLowerCase().endsWith("bed.gz") || track.file.getName().toLowerCase().endsWith("bedgraph.gz") || track.file.getName().toLowerCase().endsWith("bedgraph")) {
				  iterateBED(bedIterator, track);
				}
				else if (track.file.getName().toLowerCase().endsWith("gff.gz")){
					iterateGFF(bedIterator, track);
				}
				else if(track.file.getName().toLowerCase().endsWith(".bigwig") || track.file.getName().toLowerCase().endsWith(".bw")) {
					iterateWig(wigiter, track);
				}
				else if(track.file.getName().toLowerCase().endsWith(".bb") || track.file.getName().toLowerCase().endsWith(".bigbed")) {
				
					iterateBigBed(bediter, track);
				}
			}
			else {
				if(track.url.toString().toLowerCase().endsWith("bed.gz") || track.url.toString().toLowerCase().endsWith("bedgraph.gz")) {
					  iterateBED(bedIterator, track);
					}
					else if (track.url.toString().toLowerCase().endsWith("gff.gz")){
						iterateGFF(bedIterator, track);
					}
					else if(bbfileReader.isBigBedFile()) {
						iterateBigBed(bediter, track);
					}
					else if(bbfileReader.isBigWigFile()) {
						iterateWig(wigiter, track);
					}
			}
			track.loading = false;
			  Main.drawCanvas.ready("Reading BED-file");
			  if(tabixReader != null) {
				  tabixReader.close();
			  }
	//     track.tail = addNode;
	     if(stream != null) {
	    	 stream.close();
	     }
	   
		}
		catch(Exception e) {
			ErrorLog.addError(e.getStackTrace());
			e.printStackTrace();
	//		tabixReader.close();
			Main.drawCanvas.ready("Reading BED-file");
		}
		track.loading = false;
		repaint();
		if(track.waiting) {
			track.waiting = false;
			Main.bedCanvas.annotate(track.getHead());
			Draw.updatevars = true;
		}
		return "";
	}
	
}

void iterateBED(TabixReader.Iterator iterator, BedTrack track) {
	if(iterator == null) {
		return;
	}
	try {
	String line;
	 BedNode addNode = track.getHead();
	 track.maxvalue = 0;
		track.minvalue = Double.MAX_VALUE;
		while(addNode.getNext()!= null) {
			addNode = addNode.getNext();
			addNode.getPrev().putNext(null);
			addNode.getPrev().putPrev(null);			
		}
		
		addNode = track.getHead();
		
		 track.setCurrent(track.getHead());
		 boolean novalue = false, nostrand = false, first = true, nocolor = true, bedgraph = false;
		 if(track.file.getCanonicalPath().endsWith("bedgraph.gz") || track.file.getCanonicalPath().endsWith("bedgraph")) {
			 bedgraph = true;
		 }
		 BEDFeature features;
		 BEDCodec bedcodec = new BEDCodec();
	   while((line = iterator.next()) != null) {
	    	try {
	    		features = bedcodec.decode(line);
	    		addNode.putNext(new BedNode(features.getContig(),features.getStart()-1, (features.getEnd()-features.getStart())+2, track));	    		
	    		addNode.putNext(new BedNode(features.getContig(),features.getStart()-1, (features.getEnd()-features.getStart())+2, track));	    		
				addNode.getNext().putPrev(addNode);
				addNode = addNode.getNext();
				
				if(features.getName() != null) {
					
					/*
					if(!track.getNames().contains(features.getName())) {
						track.getNames().put(features.getName(), features.getName());
						addNode.name = track.getNames().get(features.getName());
					}
					else {
						addNode.name = track.getNames().get(features.getName());		
					}*/
					if(first) {						
						
						if(Main.SELEXhash.containsKey(features.getName())) {
						
							track.selex = true;
						}
						
					}
					if(track.selex) {
						addNode.id = features.getName();
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
					if(track.getColors().containsKey(features.getColor().getRGB())) {
						addNode.color = features.getColor().getRGB();
					}
					else {
						
						track.getColors().put(features.getColor().getRGB(),features.getColor());
						addNode.color = features.getColor().getRGB();
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
	 track.maxvalue = 0;
		track.minvalue = Double.MAX_VALUE;
		 track.setCurrent(track.getHead());
		 boolean novalue = false, nostrand = false, first = true, nocolor = true, bedgraph = false;
		
		 BedFeature features;
		String[] allfields;
	
	   while(iterator.hasNext()) {
	    	try {
				//split = line.split("\t");
	    			    		
	    		features = iterator.next();
	    		allfields = features.getRestOfFields();
	    //		position = features.getStart();
				addNode.putNext(new BedNode(features.getChromosome(),features.getStartBase()-1, (features.getEndBase()-features.getStartBase()), track));
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
	    	}
	    	first = false;
	     }     
	   if(track.minvalue < 0) {
		   track.negatives = true;
	   }
	   track.scale = Math.max(Math.abs(track.maxvalue), Math.abs(track.minvalue));
	     addNode = null;
	     track.cleared = false;
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}
void iterateGFF(TabixReader.Iterator iterator, BedTrack track) {
	try {
		
		String line;
		String[] split; 
	 	BedNode addNode = track.getHead();
	 	track.maxvalue = 0;
		track.minvalue = Double.MAX_VALUE;
		 track.setCurrent(track.getHead());
		 boolean first = true;
			
		
	   while((line = iterator.next()) != null) {
	    	try {
				split = line.split("\t");	    			    		
	 			
				if(split.length > 4) {
					addNode.putNext(new BedNode(split[0],(Integer.parseInt(split[3])),Integer.parseInt(split[4])-(Integer.parseInt(split[3]))+1, track));
					addNode.getNext().putPrev(addNode);
					addNode = addNode.getNext();
			
					if(first) {				
						
						if(Main.SELEXhash.containsKey(split[2])) {						
							track.selex = true;
						}						
					}
					if(track.selex) {
						addNode.id = split[2];
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
	   track.scale = Math.max(Math.abs(track.maxvalue), Math.abs(track.minvalue));
	     addNode = null;
	     track.cleared = false;
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}
void iterateWig(BigWigIterator iterator, BedTrack track) {
	try {
	
	 BedNode addNode = track.getHead();
	 track.maxvalue = 0;
		track.minvalue = Double.MAX_VALUE;
		 track.setCurrent(track.getHead());
		
			
		WigItem features;
		
	   while(iterator.hasNext()) {
	    	try {
			
	    		features = iterator.next();
	   			addNode.putNext(new BedNode(features.getChromosome(),features.getStartBase()+1, (features.getEndBase()-(features.getStartBase())), track));
				addNode.getNext().putPrev(addNode);
				addNode = addNode.getNext();			
				addNode.value = (double)features.getWigValue();		
				
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
	    		System.out.println(e);
	    		
	    		break;
	    	}
	    	
	     }     
	   if(track.minvalue < 0) {
		   track.negatives = true;
	   }
	   track.scale = Math.max(Math.abs(track.maxvalue), Math.abs(track.minvalue));
	  
	     addNode = null;
	     track.cleared = false;
	}
	catch(Exception e) {
		
		e.printStackTrace();
		
	}
}
void getBEDfeatures(BedTrack track, int start, int end) {

	track.used = false;
	BedReader reader = new BedReader(track, start, end);
	reader.execute();
}
void removeBeds(BedTrack track) {
	
	track.nulled = true;
	drawNode = null;
	BedNode addNode = track.getHead();
	while(addNode.getNext()!= null) {
		addNode = addNode.getNext();
		addNode.getPrev().putNext(null);
		addNode.getPrev().putPrev(null);
	}
	addNode = null;
//	track.getHead().putNext(null);
	track.setCurrent(null);
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

public void mouseClicked(MouseEvent event) {
	
	if(selectedPlay > -1 && bedTrack.get(selectedPlay).playbox.intersects(sideMouseRect)) {
		if(bedTrack.get(selectedPlay).intersect) {
			bedTrack.get(selectedPlay).intersect = false;
			
			boolean bedon = false;
			for(int i = 0 ; i<bedTrack.size(); i++) {
				if(bedTrack.get(i).intersect) {
					bedon = true;
				}
			}
			if(!bedon) {
				bedOn = false;
			}
			annotate(bedTrack.get(selectedPlay).getHead());
			if(VariantHandler.commonSlider.getValue() > 1 && !VariantHandler.clusterBox.getText().equals("0") && !VariantHandler.clusterBox.getText().equals("")) {
		    	Main.drawCanvas.clusterId = 1;
				Main.drawCanvas.calcClusters(FileRead.head);
			}
		}
		else {
			bedOn = true;
			bedTrack.get(selectedPlay).intersect = true;
			annotate(bedTrack.get(selectedPlay).getHead());
		}
		if(VariantHandler.commonSlider.getValue() > 1 && !VariantHandler.clusterBox.getText().equals("0") && !VariantHandler.clusterBox.getText().equals("")) {
	    	Main.drawCanvas.clusterId = 1;
			Main.drawCanvas.calcClusters(FileRead.head);
		}
		Draw.updatevars = true;
		Draw.calculateVars = true;
		Main.drawCanvas.repaint();
	}
	else if(selectedPlay > -1 && bedTrack.get(selectedPlay).graphBox.intersects(sideMouseRect)) {
		if(bedTrack.get(selectedPlay).graph) {
			bedTrack.get(selectedPlay).graph =false;
			repaint();
		}
		else {
			bedTrack.get(selectedPlay).graph =true;
			repaint();
		}
	}
	else if(removeTrack > -1) {
		FileRead.removeTable(this.bedTrack.get(removeTrack));
		this.bedTrack.remove(removeTrack);
		this.trackDivider.remove(removeTrack);
		if(this.bedTrack.size() == 0) {
			drawNode = null;
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
	}
	repaint();
	
	
}	
void annotate(BedNode head) {
	int varlength = 0;
	BedNode currentBed = head.getNext();
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
				if(position+varlength > currentBed.getPosition()) {
					
					if(!currentBed.getTrack().used) {	
						if(current.getBedHits() == null) {
							current.setBedhits();
						}
						
						current.getBedHits().add(currentBed);			
						current.setBedhit();
						
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
	
	head.getTrack().used = true;
	current = null;
	currentBed = null;
}

boolean checkIntersect(VarNode current, Object[] t1, Object[] t2) {
	t2 = new Object[t1.length];
	if(current.getBedHits() == null) {
		return false;
	}
	for(int i = 0; i<current.getBedHits().size(); i++) {	
		if(current.getBedHits().get(i).getTrack().intersect) {
			if(t2[current.getBedHits().get(i).getTrack().trackIndex] == null) {
				t2[current.getBedHits().get(i).getTrack().trackIndex] = 1;
			}	
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
public void mousePressed(MouseEvent event) {
	resize =false;
	pressY = event.getY();
	mouseDrag = true;
	
	this.requestFocus();
	Main.drawCanvas.pressX = event.getX();
	
	Main.drawCanvas.tempDrag = Main.drawCanvas.pressX;
	switch(event.getModifiers()) {	
	
	case InputEvent.BUTTON1_MASK: {	
		if(resizer) {
			preresize = mouseY;
			preresizer = this.trackDivider.get(this.resizeDivider)*this.getHeight();
			tempDivider = (ArrayList<Double>)trackDivider.clone();	
			
		}
		
		
		break;
	}
	case InputEvent.BUTTON3_MASK: {
		
		this.zoomDrag = false;
		
	//	repaint();
	}
}
	
}	
public void mouseReleased(MouseEvent event) {
	
	if(!sidebar) {
		if(zoomDrag) {			
			if(mouseX-Main.drawCanvas.pressX > 0) {			
				Main.drawCanvas.gotoPos(Main.drawCanvas.selectedSplit.start+(Main.drawCanvas.pressX-Main.drawCanvas.selectedSplit.offset)/Main.drawCanvas.selectedSplit.pixel, Main.drawCanvas.selectedSplit.start+(mouseX-Main.drawCanvas.selectedSplit.offset)/Main.drawCanvas.selectedSplit.pixel);
			}
		}
	
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
}	
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
										System.out.println(this.tempDivider.get(i));
										e.printStackTrace();
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
			zoomDrag = true;				
			
			repaint();
			return;
		}
		break;
	}
	case InputEvent.BUTTON3_MASK: {	
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
//	this.mouseX = event.getX();

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
		nodeHeight++;
		
		repaint();
	}
	if(keyCode == KeyEvent.VK_MINUS || keyCode == 109) {
		nodeHeight--;
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
