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
import java.awt.Font;
import java.awt.Rectangle;
import java.io.Serializable;
import java.util.ArrayList;


public class Reads implements Serializable{
	
	private static final long serialVersionUID = 1L;

	Sample sample;
	private transient double maxcoverage = 0.0;
	private transient double[][] coverages; // = new int[(int)Main.screenSize.getWidth()];
	private transient ArrayList<ReadNode> reads = new ArrayList<ReadNode>();
	private transient ArrayList<ReadNode[]> headAndTail = new ArrayList<ReadNode[]>();

	private transient ReadNode firstRead, lastRead;
	private transient int maxReadSize = 0;
	int startpos, endpos;
	int readwheel = 0;
	private transient Rectangle scrollbar = new Rectangle(), scroller = new Rectangle();
	private transient int readstart=0, readend=0, coveragestart = 0, coverageend = 0;
	private transient boolean readScroll = false;
	public boolean complete = false, loading = false, nodraw = false;
	Font readfont =new Font("SansSerif", Font.BOLD, 10);
	int readHeight = 10;
	public ReferenceSeq reference;
	
	Rectangle getScrollBar() {
		return this.scrollbar;
	}
	Rectangle getScroller() {
		return this.scroller;
	}
	boolean isReadScroll() {
		return this.readScroll;
	}
	void setReadScroll(boolean value) {
		this.readScroll = value;
	}
	void setMaxcoverage(double coverage) {
		this.maxcoverage = coverage;
	}
	double getMaxcoverage() {
		return this.maxcoverage;
	}
	int getReadSize() {
		return this.maxReadSize;
	}
	void setReadSize(int max) {
		this.maxReadSize = max;
	}
	double[][] getCoverages() {
		return this.coverages;
	}
	void setCoverages(double[][] coverages) {
		this.coverages = coverages;
	}
	int getReadStart() {
		return readstart;
	}
	int getCoverageStart() {
		return this.coveragestart;
	}
	int getCoverageEnd() {
		return this.coverageend;
	}
	void setCoverageStart(int start) {
		this.coveragestart = start;		
	}
	void setCoverageEnd(int end) {
		this.coverageend = end;		
	}
	int getReadEnd() {
		return readend;
	}
	void setReadStart(int value) {
		this.readstart = value;
	}
	void setReadEnd(int value) {
		this.readend = value;
	}
	ReadNode getFirstRead() {
		return this.firstRead;
	}
	ReadNode getLastRead() {
		return this.lastRead;
	}
	void setFirstRead(ReadNode read) {
		this.firstRead = read;
	}
	void setLastRead(ReadNode read) {
		this.lastRead = read;		
	}
	void resetReads() {
		this.reads = new ArrayList<ReadNode>();
		this.headAndTail = new ArrayList<ReadNode[]>();
		this.loading = false;
	}
	ArrayList<ReadNode> getReads() {
		return this.reads;
	}
	ArrayList<ReadNode[]> getHeadAndTail() {
		return this.headAndTail;
	}
}


