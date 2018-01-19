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
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.io.Serializable;
import java.net.URL;
import java.util.ArrayList;
import java.util.Hashtable;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JPopupMenu;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;

import base.BBfile.BBFileHeader;


public class BedTrack implements Serializable, ActionListener, KeyListener, MouseListener, PopupMenuListener, DocumentListener, ChangeListener {
	
	private static final long serialVersionUID = 1L;
	final File file;
	final URL url, index;
	Integer chromcolumn, startcolumn, endcolumn, valuecolumn, namecolumn, strandcolumn, basecolumn;  
	private transient BedNode current, head = new BedNode("",0,0, this);
	int prepixel = 0, mouseWheel = 0, bedstart=0, bedend=0, trackIndex;
	private transient ArrayList<Integer> bedLevelMatrix= new ArrayList<Integer>();
	private transient Hashtable<Integer, Color> colors = new Hashtable<Integer, Color>();
	boolean selex = false, small = false, cleared = false, intersect = false, graph = false;
	Rectangle playbox = new Rectangle(), graphBox = new Rectangle(), settingsButton = new Rectangle();	
	Polygon playTriangle = new Polygon();
	private transient BedTable table;
	double maxvalue = 0, minvalue = Double.MAX_VALUE, scale = 0;
	boolean negatives = false, first = true, used = false, loading = false, waiting = false, hasvalues = false, updated = false;
	String chr = "";
	
	ArrayList<JCheckBox> menuBoxes;
	private transient BedNode drawNode = null;
	short iszerobased = 0, nodeHeight = 10;
	private transient BBFileHeader bbfileheader = null;
	public boolean nulled = false;
	private transient JPopupMenu menu;
	private transient JCheckBox zerobased, logscale, intersectBox, subtractBox, collapseBox, affinityChange, varcalc;
	private transient JTextField limitField;
	private transient JButton columnSelector;
	Double limitValue = (double)Integer.MIN_VALUE;
	private transient ColumnSelector selector = null;
	private transient JSlider sizeSlider;
	private transient base.BBfile.BBFileHeader fileheader;
	private transient base.BBfile.BBZoomLevels zoomlevels;
	private transient BBFileReader bbfileReader;
	private transient Integer zoomLevel = null;	
	
	public Integer getZoomlevel() {
		return this.zoomLevel;
	}
	public void setZoomlevel(Integer i) {
		this.zoomLevel = i; 
	}
	public BedNode getDrawNode() {
		return this.drawNode;
	}
	public void setDrawNode(BedNode node) {
		this.drawNode = node;
	}
	
	public BBFileReader getBBfileReader() {
		return this.bbfileReader;
	}
	public void setBBfileReader(BBFileReader reader) {
		this.bbfileReader = reader;
	}
	public base.BBfile.BBFileHeader getFileHeader() {
		return this.fileheader;
	}
	public void setFileHeader(base.BBfile.BBFileHeader header) {
		this.fileheader = header;
	}
	public base.BBfile.BBZoomLevels getZoomlevels() {
		return this.zoomlevels;
	}
	public void setZoomlevels(base.BBfile.BBZoomLevels levels) {
		this.zoomlevels = levels;
	}
	public ColumnSelector getSelector() {
		return this.selector;
	}
	public JPopupMenu getPopup() {
		return this.menu;
	}
	public JCheckBox getIntersectBox() {
		return this.intersectBox;
	}
	public JCheckBox getSubtracttBox() {
		return this.subtractBox;
	}
	public JTextField getLimitField() {
		return this.limitField;
	}
	public JCheckBox getCollapseBox() {
		return this.collapseBox;
	}
	public JCheckBox getAffinityBox() {
		return this.affinityChange;
	}
	public JCheckBox getLogscale() {
		return this.logscale;
	}
	public JCheckBox getZerobased() {
		return this.zerobased;
	}
	public JCheckBox getVarcalc() {
		return this.varcalc;
	}
	public BedTrack(File file, int indexnro) {		
		this.file = file;
		url = null;
		index = null;
		this.trackIndex = indexnro;
		selector = new ColumnSelector(this);
		setmenu();
	}
	
	public BedTrack(URL url, URL index, int indexnro) {
		this.file = null;
		this.url = url;
		this.index = index;
		this.trackIndex = indexnro;
		setmenu();
	}
	public void setCollapsebox() {
		collapseBox = new JCheckBox("Auto collapse");
	}
	
	public void setmenu() {
		try {
			settingsButton = new Rectangle();
			if(menuBoxes == null || menuBoxes.size() != 7) {
				intersectBox = new JCheckBox("Intersect");			
				subtractBox = new JCheckBox("Subtract");	
				zerobased = new JCheckBox("Zero Based");
				varcalc = new JCheckBox("Apply in annotation");
				logscale = new JCheckBox("Log Scale");				
				setCollapsebox();
				affinityChange = new JCheckBox("Report affinity changes");
				menuBoxes = new ArrayList<JCheckBox>();				
				menuBoxes.add(intersectBox);
				menuBoxes.add(subtractBox);
				menuBoxes.add(zerobased);
				menuBoxes.add(logscale);
				menuBoxes.add(collapseBox);
				menuBoxes.add(affinityChange);
				menuBoxes.add(varcalc);
				zerobased.setSelected(true);
				intersectBox.setSelected(true);
				affinityChange.setVisible(false);				
			}
			else {					
				intersectBox = menuBoxes.get(0);
				subtractBox = menuBoxes.get(1);
				zerobased = menuBoxes.get(2);
				logscale = menuBoxes.get(3);
				collapseBox = menuBoxes.get(4);
				affinityChange = menuBoxes.get(5);				
				varcalc = menuBoxes.get(6);						
				
			}
			
			columnSelector = new JButton("File format");
			
			menu = new JPopupMenu("Options");
			menu.setLayout(new GridBagLayout());
			GridBagConstraints con = new GridBagConstraints();
			
			limitField = new JTextField("Value limit");
			
			sizeSlider = new JSlider(JSlider.VERTICAL,1,20,nodeHeight);
			sizeSlider.setInverted(true);
			nodeHeight = (short)sizeSlider.getValue();
			sizeSlider.addChangeListener(this);
			sizeSlider.setValue(nodeHeight);
			con.anchor = GridBagConstraints.NORTHWEST;
			con.gridheight = 8;
			menu.add(sizeSlider,con);
			con.gridx = 1;
			con.gridy++;
			con.gridheight = 1 ;
			menu.add(intersectBox,con);
			con.gridy++;
			menu.add(subtractBox,con);
			con.gridy++;
			menu.add(zerobased,con);
			con.gridy++;
			menu.add(logscale,con);
			con.gridy++;
			menu.add(collapseBox,con);
			con.gridy++;
			menu.add(affinityChange,con);	
			con.gridy++;
			menu.add(varcalc, con);
			con.gridy++;
			menu.add(limitField,con);
			con.gridy++;
			menu.add(columnSelector,con);
			
			limitField.setPreferredSize(new Dimension(menu.getFontMetrics(Main.menuFont).stringWidth("__Value limit__"), Main.defaultFontSize+6));
			limitField.setMinimumSize(new Dimension(menu.getFontMetrics(Main.menuFont).stringWidth("__Value limit__"), Main.defaultFontSize+6));
			limitField.setToolTipText("Value limit");
			if(limitValue != (double)Integer.MIN_VALUE) {				
				limitField.setText(""+limitValue);
			}
			
			menu.addPopupMenuListener(this);
			collapseBox.addActionListener(this);
			collapseBox.setSelected(true);
			limitField.addMouseListener(this);
			intersectBox.addActionListener(this);
			subtractBox.addActionListener(this);
			zerobased.addActionListener(this);
			logscale.addActionListener(this);
			limitField.addKeyListener(this);
			limitField.getDocument().addDocumentListener(this);
			columnSelector.addActionListener(this);
			for(int c = 0 ; c<getPopup().getComponentCount(); c++) {
				getPopup().getComponent(c).setFont(Main.menuFont);
			}
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
	}
	void setTable(BedTable table) {
		this.table = table;
	}
	void setSelector() {
		this.selector = new ColumnSelector(this);	
	}
	
	public JButton getSelectorButton() {
		return this.columnSelector;
	}
	BedTable getTable() {
		return this.table;
	}
	void setHead() {
		this.head = new BedNode("",0,0, this);
	}
	
	void setColors() {
		this.colors = new Hashtable<Integer, Color>();
	}
	public BBFileHeader getBBheader() {
		return this.bbfileheader;
	}
	public void setBBheader(BBFileHeader header) {
		this.bbfileheader = header;
	}
	BedNode getCurrent() {
		return current;
	}
	void setCurrent(BedNode current) {
		this.current = current;
	}
	BedNode getHead() {
		return this.head;
	}
	ArrayList<Integer> getBedLevels() {
		return this.bedLevelMatrix;
	}
	void setBedLevels() {
		bedLevelMatrix= new ArrayList<Integer>();
	}
	
	Hashtable<Integer, Color> getColors() {
		return colors;
	}
	@Override
	public void actionPerformed(ActionEvent e) {
		updated = true;
		if(e.getSource() == zerobased) {
			updated = false;
			if(zerobased.isSelected()) {
				iszerobased = 0;
				
			}
			else {
				iszerobased = 1;
			}
			if(!small) {
				if(Main.drawCanvas.splits.get(0).viewLength < Settings.windowSize) {
					int start =  (int)Main.drawCanvas.splits.get(0).start+(int)((Main.drawCanvas.splits.get(0).end-Main.drawCanvas.splits.get(0).start)/2)- Settings.windowSize/2;
					Main.bedCanvas.getBEDfeatures(this,start, start+ Settings.windowSize);
				}
			}
			else {
				Main.bedCanvas.getBEDfeatures(this,1, Main.drawCanvas.splits.get(0).chromEnd);
			}
		}
		else if(e.getSource() == intersectBox) {
			subtractBox.setSelected(false);
			Object[] t1 = new Object[Main.bedCanvas.bedTrack.size()];
			
			for(int i = 0; i<t1.length; i++) {
				if(!Main.bedCanvas.bedTrack.get(i).intersect) {
					continue;
				}
				if(Main.bedCanvas.bedTrack.get(i).getIntersectBox().isSelected()) {			
					t1[i] = 1;		
				}				
			}	
			
			Main.bedCanvas.checkIntersectAll(FileRead.head.getNext(),t1);			
			Draw.updatevars = true;
			Main.drawCanvas.repaint();
		}
		else if(e.getSource() == subtractBox) {
			intersectBox.setSelected(false);
			Object[] t1 = new Object[Main.bedCanvas.bedTrack.size()];
			
			for(int i = 0; i<t1.length; i++) {
				if(!Main.bedCanvas.bedTrack.get(i).intersect) {
					continue;
				}
				if(Main.bedCanvas.bedTrack.get(i).getIntersectBox().isSelected()) {			
					t1[i] = 1;		
				}				
			}	
			
			Main.bedCanvas.checkIntersectAll(FileRead.head.getNext(),t1);
			Draw.updatevars = true;
			Main.drawCanvas.repaint();
		}
		else if(e.getSource() == columnSelector) {
			if(selector != null) {
				selector.frame.setVisible(true);
			}
			
		}
		
		Main.bedCanvas.repaint();
		
	}
	@Override
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void keyPressed(KeyEvent e) {
		if(e.getSource() == limitField) {
			if(e.getKeyCode() == KeyEvent.VK_ENTER) {
				try {
					limitValue = Double.parseDouble(limitField.getText());
					
					if(!small || (file !=null && file.getName().endsWith("bw")) || (url != null && url.getFile().endsWith("bw"))) {
						
						if(Main.drawCanvas.splits.get(0).viewLength <  Settings.windowSize) {
							int start =  (int)Main.drawCanvas.splits.get(0).start- Settings.windowSize;
							Main.bedCanvas.getBEDfeatures(this,start, start+ Settings.windowSize*2);
						}
						else {
							if(this.intersect) {
								
								Main.bedCanvas.removeBedhits(this);
								BedCanvas.Annotator annotator = Main.bedCanvas.new Annotator(this);
								annotator.execute();	
							}
						}
					}
					else {
						Main.bedCanvas.removeBedhits(this);
						Main.bedCanvas.getBEDfeatures(this,1, Main.drawCanvas.splits.get(0).chromEnd);
						
					}
					updated = true;
				}
				catch(Exception ex ) {
					limitValue = (double)Integer.MIN_VALUE;
				}
			}
		}
		
	}
	@Override
	public void keyReleased(KeyEvent e) {
		
	}

	@Override
	public void mouseClicked(MouseEvent e) {
	
	}

	@Override
	public void mousePressed(MouseEvent e) {
		
		if(e.getSource() == limitField) {
			
			if(limitField.isEditable()) {
				if(limitField.getText().contains("limit")) {
					limitField.setText("");
				}
			}
		}
	}

	@Override
	public void mouseReleased(MouseEvent e) {
	
		
	}

	@Override
	public void mouseEntered(MouseEvent e) {
		
		
	}

	@Override
	public void mouseExited(MouseEvent e) {
	
		
	}
	

	@Override
	public void popupMenuWillBecomeVisible(PopupMenuEvent e) {
		// TODO Auto-generated method stub
		updated = true;
	}

	@Override
	public void popupMenuWillBecomeInvisible(PopupMenuEvent e) {
	
		if(updated) {
			return;
		}
		try {
		
			limitValue = Double.parseDouble(limitField.getText());
			
			if(!small) {
				
				if(Main.drawCanvas.splits.get(0).viewLength <  Settings.windowSize) {
					int start =  (int)Main.drawCanvas.splits.get(0).start+(int)((Main.drawCanvas.splits.get(0).end-Main.drawCanvas.splits.get(0).start)/2)- Settings.windowSize/2;
					Main.bedCanvas.getBEDfeatures(this,start, start+ Settings.windowSize);
				}
				else {
					if(this.intersect) {
						
						Main.bedCanvas.removeBedhits(this);
						BedCanvas.Annotator annotator = Main.bedCanvas.new Annotator(this);
						annotator.execute();	
					}
				}
			}
			else {
				Main.bedCanvas.removeBedhits(this);
				Main.bedCanvas.getBEDfeatures(this,1, Main.drawCanvas.splits.get(0).chromEnd);
				
			}
			
			
		}
		catch(Exception ex) {
			
			limitValue = (double)Integer.MIN_VALUE;
		
			limitField.setForeground(Color.red);
			limitField.revalidate();
			if(!small) {
				
				if(Main.drawCanvas.splits.get(0).viewLength <  Settings.windowSize) {
					int start =  (int)Main.drawCanvas.splits.get(0).start+(int)((Main.drawCanvas.splits.get(0).end-Main.drawCanvas.splits.get(0).start)/2)- Settings.windowSize/2;
					Main.bedCanvas.getBEDfeatures(this,start, start+ Settings.windowSize);
				}
			}
			else {
				Main.bedCanvas.removeBedhits(this);
				Main.bedCanvas.getBEDfeatures(this,1, Main.drawCanvas.splits.get(0).chromEnd);
				
			}
		}
		
	}

	@Override
	public void popupMenuCanceled(PopupMenuEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void insertUpdate(DocumentEvent e) {
		updated = true;
		
		if(limitField.getText().length() == 0) {
			limitValue = (double)Integer.MIN_VALUE;
			
			updated = false;
			/*if(!small) {
				
				if(Main.drawCanvas.splits.get(0).viewLength < 1000000) {
					int start =  (int)Main.drawCanvas.splits.get(0).start+(int)((Main.drawCanvas.splits.get(0).end-Main.drawCanvas.splits.get(0).start)/2)-500000;
					Main.bedCanvas.getBEDfeatures(this,start, start+1000000);
				}
			}
			else {
				Main.bedCanvas.removeBedhits(this);
				Main.bedCanvas.getBEDfeatures(this,1, Main.drawCanvas.splits.get(0).chromEnd);
				
			}*/
			return;
		}
		
		
		
			try {
				limitValue = Double.parseDouble(limitField.getText());
				limitField.setForeground(Color.black);
				updated = false;
				limitField.revalidate();
			/*	if(!small) {
					
					if(Main.drawCanvas.splits.get(0).viewLength < 1000000) {
						int start =  (int)Main.drawCanvas.splits.get(0).start+(int)((Main.drawCanvas.splits.get(0).end-Main.drawCanvas.splits.get(0).start)/2)-500000;
						Main.bedCanvas.getBEDfeatures(this,start, start+1000000);
					}
				}
				else {
					Main.bedCanvas.removeBedhits(this);
					Main.bedCanvas.getBEDfeatures(this,1, Main.drawCanvas.splits.get(0).chromEnd);
					
				}
				*/
				
			}
			catch(Exception ex) {
				limitValue = (double)Integer.MIN_VALUE;
				updated = false;
				limitField.setForeground(Color.red);
				limitField.revalidate();
			}
			
		}		
		
	

	@Override
	public void removeUpdate(DocumentEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void changedUpdate(DocumentEvent e) {
		
		
		
	
	}
	@Override
	public void stateChanged(ChangeEvent e) {
		if(e.getSource() == sizeSlider) {
			Main.bedCanvas.heightchanged = true;
			nodeHeight = (short)sizeSlider.getValue();
			Main.bedCanvas.updateTrack = this;
			Main.bedCanvas.repaint();
			Main.bedCanvas.updateTrack = null;
		}
		
	}
}
