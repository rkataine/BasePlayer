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
import java.util.AbstractMap;
import java.util.ArrayList;

//import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class VarNode {
	
	private final int position;
	
	private VarNode prev, next;
	private final byte refBase;
	private ArrayList<Transcript.Exon> exons;
	private ArrayList<Transcript> transcripts;
	private boolean inGene;
	boolean inVarList, bedhit;
	Integer clusterId;
	private ArrayList<BedNode> bedHits;
	String rscode;
	ClusterNode clusterNode;
	ArrayList<Map.Entry<String, ArrayList<SampleNode>>> vars = new ArrayList<Map.Entry<String, ArrayList<SampleNode>>>();	
	boolean incluster, indel, found;
	public String codon;
	boolean controlled;
	private Entry<String, ArrayList<SampleNode>> entry;

	public boolean coding = false;
	
	public VarNode(int position, byte ref, String variation, short coverage, short calls, boolean genotype, short quality, String rscode, Sample sample, VarNode prev, VarNode next) {
		this.position = position;
		this.rscode = rscode;
		this.refBase = ref; //Main.bases.get(ref);
		if (variation.length() > 1 ) indel = true;
		if(vars.size() == 0) {
			Map.Entry<String, ArrayList<SampleNode>> entry = new AbstractMap.SimpleEntry<String, ArrayList<SampleNode>>(variation, new ArrayList<SampleNode>());
			
			entry.getValue().add(new SampleNode(coverage, calls, genotype, quality, sample));
			vars.add(entry);
		}
		else {
			found = false;
			for(Map.Entry<String, ArrayList<SampleNode>> entry : vars) {
				if(entry.getKey().equals(variation)) {
					entry.getValue().add(new SampleNode(coverage, calls, genotype, quality, sample));
					found = true;
					break;
				}
				
			}
			if(!found) {
				Map.Entry<String, ArrayList<SampleNode>> entry = new AbstractMap.SimpleEntry<String, ArrayList<SampleNode>>(variation, new ArrayList<SampleNode>());
				entry.getValue().add(new SampleNode(coverage, calls, genotype, quality, sample));
				
				vars.add(entry);
			}
		}
		
	//	samples.add(new SampleNode(variation, coverage, calls, genotype, quality, sample));
		this.prev = prev;
		this.next = next;		
	}
	
	public void setTranscripts() {
		this.transcripts = new ArrayList<Transcript>();
	}
	public void setExons() {
		this.exons = new ArrayList<Transcript.Exon>();
	}
	public void setBedhits() {
		this.bedHits = new ArrayList<BedNode>();
	}
	public void remBedhits() {
		this.bedHits = null;
	}
/*	public VarNode clone() throws CloneNotSupportedException {
        return (VarNode) super.clone();
	}*/	
	public void putNext(VarNode node) {
		this.next = node;
	}
	public void putPrev(VarNode node) {
		this.prev = node;
	}		
	public byte getRefBase() {
		return this.refBase;
	}
	public int getPosition() {
		return this.position;
	}
	public VarNode getPrev() {
		return this.prev;
	}
	
	public VarNode getNextVisible(VarNode node) {
		node = node.getNext();
		while(Main.drawCanvas.hideNodeTotal(node)) {
			if(node != null) {
				node = node.next;
			}
			else {
				break;
			}
		}
		return node;	
		
	}
	
	private VarNode getNextVisibleCluster(VarNode node) {
		while(Main.drawCanvas.hideNodeCluster(node)) {
			if(node != null) {
				node = node.next;
			}
			else {
				break;
			}
		}
		return node;
	}
	public String getChrom() {
		if(getTranscripts() != null) {
			return getTranscripts().get(0).getChrom();
		}
		else if(getExons() != null) {
			return getExons().get(0).getTranscript().getChrom();
		}
		else {
			return null;
		}
		
	}
	public VarNode getNext(VarNode node ) {
		if(Main.drawCanvas.hideNodeCluster(this.next)) {
			if(this.next != null) {
				return getNextVisibleCluster(this.next.next);
			}
			else {
				return null;
			}
		}
		else {
			
			return this.next;
		}
	}
	public VarNode getNext() {
		return this.next;
		
	}
	public boolean isInGene() {
		return this.inGene;
	}
	public void setInGene() {
		this.inGene = true;
	}
	public void setBedhit(boolean value) {
		this.bedhit = value;
	}
	public boolean isBedhit() {
		return this.bedhit;		
	}
	public ArrayList<BedNode> getBedHits() {
		return this.bedHits;
	}
	public void setCodon(String codon) {
		this.codon = ChromDraw.codons.get(codon);
	}
	public void setRscode(String value) {
		this.rscode = value;
	}
	public String isRscode() {
		return this.rscode;
	}
	public String getCodon() {
		return this.codon;
	}
	public java.util.List<Transcript.Exon> getExons() {
		return this.exons;
	}
	public java.util.List<Transcript> getTranscripts() {
		return this.transcripts;
	}
	public List<SampleNode> getSamples(String variation) {
		for(Map.Entry<String, ArrayList<SampleNode>> entry : vars) {
			if(entry.getKey().equals(variation)) {
				return entry.getValue();
			}			
		}
		return null;
	}

	public void addSample(String variation, short coverage, short calls, boolean genotype, short quality, Sample sample) {	
		found = false;	
		if (variation.length() > 1 ) indel = true;
		for(Map.Entry<String, ArrayList<SampleNode>> entry : vars) {
			if(entry.getKey().equals(variation)) {
				entry.getValue().add(new SampleNode(coverage, calls, genotype, quality, sample));
				found = true;
				break;
			}			
		}
		if(!found) {
			Map.Entry<String, ArrayList<SampleNode>> entry = new AbstractMap.SimpleEntry<String, ArrayList<SampleNode>>(variation, new ArrayList<SampleNode>());
			entry.getValue().add(new SampleNode(coverage, calls, genotype, quality, sample));
			vars.add(entry);
		}

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
	public void removeSample(Sample sample) {
				
		for(int v = 0; v<vars.size();v++) {			
			entry = vars.get(v);
			for(int i = 0 ;i<entry.getValue().size(); i++) {
				if(entry.getValue().get(i).getSample() == null) {
					break;
				}
				if(entry.getValue().get(i).getSample().equals(sample)) {
					entry.getValue().remove(i);
					
					if(entry.getValue().size() == 0) {
						vars.remove(entry);						
					}
					break;					
				}
			}			
		}
		if(vars.size() == 0) {
			if(this.getNext() != null) {
				this.getNext().putPrev(this.getPrev());
			}
			this.getPrev().putNext(this.getNext());
			
		}
	/*	SampleNode node; // TODO
		for(int i = 0;i<samples.size(); i++) {
			node = samples.get(i);
			if(node.getSample().getIndex() > sample.getIndex()) {
				break;
			}
			if(node.getSample().equals(sample)) {
				samples.remove(node);				
				i--;
			}
		}
		if(samples.size() == 0) {
			if(this.getNext() != null) {
				this.getNext().putPrev(this.getPrev());
			}
			this.getPrev().putNext(this.getNext());
			Main.datasize--;
		}
		node = null;*/
	}
}
