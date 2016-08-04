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
import java.awt.Polygon;
import java.awt.Rectangle;
import java.io.File;
import java.io.Serializable;
import java.net.URL;
import java.util.ArrayList;
import java.util.Hashtable;

import base.BBfile.BBFileHeader;

public class BedTrack implements Serializable {
	
	private static final long serialVersionUID = 1L;
	final File file;
	final URL url, index;
	
	private transient BedNode current, head = new BedNode("",0,0, this);
	int prepixel = 0, mouseWheel = 0, bedstart=0, bedend=0, trackIndex;
	private transient ArrayList<Integer> bedLevelMatrix= new ArrayList<Integer>();
//	private transient Hashtable<String, String> names = new Hashtable<String, String>();
	private transient Hashtable<Integer, Color> colors = new Hashtable<Integer, Color>();
	boolean selex = false, small = false, cleared = false, intersect = false, graph = false;
	Rectangle playbox = new Rectangle(), graphBox = new Rectangle(), collapseBox = new Rectangle();	
	Polygon playTriangle = new Polygon();
	private transient BedTable table;
	double maxvalue = 0, minvalue = Double.MAX_VALUE, scale = 0;
	boolean negatives = false, first = true, used = false, loading = false, waiting = false;
	String chr = "";
	private transient BBFileHeader bbfileheader = null;
	public boolean nulled = false;
	
	public BedTrack(File file, int indexnro) {
		this.file = file;
		url = null;
		index = null;
		this.trackIndex = indexnro;
	}
	public BedTrack(URL url, URL index, int indexnro) {
		this.file = null;
		this.url = url;
		this.index = index;
		this.trackIndex = indexnro;
	}
	void setTable(BedTable table) {
		this.table = table;
	}
	BedTable getTable() {
		return this.table;
	}
	void setHead() {
		this.head = new BedNode("",0,0, this);
	}
	/*void setNames() {
		this.names = new Hashtable<String, String>();
	}*/
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
	/*Hashtable<String, String> getNames() {
		return this.names;
	}*/
	Hashtable<Integer, Color> getColors() {
		return colors;
	}
}
