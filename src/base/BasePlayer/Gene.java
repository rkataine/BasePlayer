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
import java.util.HashMap;

public class Gene {
	
	private boolean strand = false, showIsoforms = false;
	private String description = "-", name, ID, chrom;
	private Transcript canonical;
	private int start, end;
	private ArrayList<Transcript> transcripts = new ArrayList<Transcript>();
	
	
	public Gene(String[] line) {
		this.chrom = line[0];
		this.name = line[3];
		this.ID = line[6];
		this.start = Integer.parseInt(line[1]);
		this.end = Integer.parseInt(line[2]);
		this.description = line[16];
		if(line[5].equals("+")) {
			this.strand = true;			
		}			
	}
	public Gene(String chrom, HashMap<String, String> gffhash) {
		if(chrom.equals("-1")) {
			this.chrom = gffhash.get("seqid");
		}
		else {
			this.chrom = chrom;
		}
		
		this.name = FileRead.getInfoValue(gffhash,"name");
		if(gffhash.containsKey("dbxref")) {
			this.ID = gffhash.get("dbxref");
		}
		else {
			this.ID = FileRead.getInfoValue(gffhash,"id");
		}
		this.start = Integer.parseInt(gffhash.get("start"));
		
		this.end = Integer.parseInt(gffhash.get("end"));
		if(gffhash.containsKey("description")) {
			this.description = FileRead.getInfoValue(gffhash,"description").replace("%20", " ");
		}
		else if (gffhash.containsKey("note")) {
			this.description = FileRead.getInfoValue(gffhash,"note").replace("%20", " ");
		}
		if(FileRead.getInfoValue(gffhash,"strand").equals("+")) {
			this.strand = true;			
		}	
	}
	public boolean getStrand() {
		return this.strand;
	}
	public boolean showIsoforms() {
		return this.showIsoforms;
	}
	
	public String getChrom() { return this.chrom; }
	public int getStart() { return this.start; }
	public int getEnd() {return this.end;}
	public String getDescription() {return this.description;}
	public String getName() {return this.name;}
	public String getID() {return this.ID;}
	public Transcript getCanonical() {return this.canonical;}
	public ArrayList<Transcript> getTranscripts() {	return this.transcripts;}
	public void addTranscript(Transcript trans) {
		this.transcripts.add(trans);
		if (this.start > trans.getStart()) {
			this.start = trans.getStart();
		}
		if(this.end < trans.getEnd()) {
			this.end = trans.getEnd();
		}
		if(trans.isCanonical()) {
			this.canonical = trans;
		}
	
	}
	public void setShowIsoforms(boolean value) {
		this.showIsoforms = value;
	}
	
	public void setStart(int start) { this.start = start; }
	public void setEnd(int end) { this.end = end;}
	public void setDescription(String descr) { this.description = descr;}
	public void setCanonical(Transcript trans) { this.canonical = trans;}
	
	
	
	
}
