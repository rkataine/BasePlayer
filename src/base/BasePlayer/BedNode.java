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

import java.util.ArrayList;

public class BedNode {
	
	private static final long serialVersionUID = 1L;
	private final int position;	
	private final int length;
	int level;
	Double value;
	Boolean forward, inVarlist = false;
	int mutations;
	private final String chrom;
	private BedNode next, prev;
	String name;
	String id;
	String info;
	int color;
	private final BedTrack track;
	ArrayList<VarNode> varnodes;
	
	public BedNode(String chr, int position, int length, BedTrack track) {
		this.position = position;
		if(length < 1) {
			this.length = 1;
		}
		else {
			this.length = length;	
		}
		this.track = track;
		this.chrom = chr;
	}
	
	public String getChrom() {
		return this.chrom;
	}
	public BedTrack getTrack() {
		return this.track;
	}
	public void removeNode() {
		if(this.getPrev() != null) {
			this.getPrev().putNext(this.getNext());
		}
		if(this.getNext() != null) {
			this.getNext().putPrev(this.getPrev());
		}
		this.putNext(null);
		this.putPrev(null);
	}
	public int getPosition() {
		return this.position;
	}
	public int getDrawPosition() {
		return this.position+1;
	}
	public int getLength() {
		return this.length;
	}
	public BedNode getNext() {
		return this.next;
	}
	public BedNode getPrev() {
		return this.prev;
	}
	public void putNext(BedNode node) {
		this.next = node;
	}
	public void putPrev(BedNode node) {
		this.prev = node;
	}
}
