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
import java.awt.Rectangle;

import java.util.ArrayList;
import java.util.HashMap;

public class Transcript {
	
	
	private boolean canonical = true;
	private final String biotype,uniprot,ID;

	private int start, end, codingStart, codingEnd, exonstart, exonend, length;
	private Exon[] exons;
	ArrayList<Exon> exonArray;
	private String[] exonstarts, exonends, exonphases;
	private Gene gene;
	
	short startphase, endphase, endphasetemp;
	int firstamino = 0, ypos = 0;
/*	ArrayList<VarNode> varnodes = new ArrayList<VarNode>();
	int mutations = 0, nonsense=0, missense=0, synonymous=0, intronic = 0, utr =0;
	ArrayList<Sample> samples = new ArrayList<Sample>();
	*/
	
	class Exon {
		
		
		private int start, end;
		short startPhase = -1,endPhase =-1, nro;
		private int firstamino;
		
		private Rectangle rectangle;
		Transcript transcript;		
		
		public Exon(int start, int end, short startphase, short endphase, short nro, int firstAmino, Transcript transcript) {
			this.start = start;
			this.end = end;
			
			this.startPhase = startphase;
			this.endPhase = endphase;
			this.rectangle = new Rectangle();
			this.nro = nro;
			this.transcript = transcript;
			this.firstamino = firstAmino;
		}
		public Exon(int start, int end, short nro, Transcript transcript) {
			this.start = start;
			this.end = end;			
			this.nro = nro;
			this.transcript = transcript;
			
		}
		
		public int getFirstAmino() {
			return this.firstamino;
		}
		public int getStart() {
			return this.start;
		}
		public int getEnd() {
			return this.end;
		}
		public short getStartPhase() {
			return this.startPhase;
		}
		public short getEndPhase() {
			return this.endPhase;
		}
		public Rectangle getRectangle() {
			return this.rectangle;
		}
		public short getNro() {
			return this.nro;
		}	
		public Transcript getTranscript() {
			return this.transcript;
		}
	}
	
	public Transcript(String[] line) {
	
		this.biotype = line[10];
		this.uniprot = line[8];
		
		if(line[9].equals("-")) {
			this.canonical = false;			
		}	
		if(line[7].contains(":")) {
			this.ID = line[7].split(":")[1];
		}
		else {
			this.ID = line[7];
		}
		
		
		this.start = Integer.parseInt(line[1]);
		this.end = Integer.parseInt(line[2]);
		this.length = this.end-this.start;
		this.codingStart = Integer.parseInt(line[11]);
		this.codingEnd = Integer.parseInt(line[12]);
		
		Boolean strand = false;
		if(line[5].equals("+")) {
			strand = true;			
		}	
		
		this.exons = new Exon[Integer.parseInt(line[4])];		
		exonstarts = line[13].split(",");
		exonends = line[14].split(",");
		exonphases = line[15].split(",");
		//Forward strand
		if(strand) {
			firstamino = 0;
			for(short i = 0; i<this.exons.length; i++) {
				exonstart = Integer.parseInt(exonstarts[i]);
				exonend = Integer.parseInt(exonends[i]);
				startphase = Short.parseShort(exonphases[i]);
				
				if(startphase > 0) {
					startphase = (short)(3-startphase);
				}
				
				if(startphase > -1) {
					endphasetemp = 0;
					if(this.codingStart >= exonstart && this.codingStart < exonend) {
						firstamino = 1;
						endphase = (short)((exonend-this.codingStart)%3);
						if(endphase > 0) {
							endphasetemp = (short)(3-endphase);
						}
						if(this.codingEnd > exonstart && this.codingEnd < exonend) {
							endphase = 0;
						}
						this.exons[i] = new Exon(exonstart, exonend, startphase, endphase, (short)(i+1), firstamino, this);
						firstamino += ((exonend-this.codingStart+endphasetemp)/3);
					}
					else {
						if(this.codingEnd > exonstart && this.codingEnd < exonend) {
							endphase = 0;
						}
						else {
							endphase = (short)((exonend-(exonstart+startphase))%3);		
							if(endphase > 0) {
								endphasetemp = (short)(3-endphase);
							}
						}
						this.exons[i] = new Exon(exonstart, exonend, startphase, endphase, (short)(i+1), firstamino, this);
						firstamino += ((exonend-exonstart-startphase+endphasetemp)/3);
					}
				}
				else {
					this.exons[i] = new Exon(exonstart, exonend, startphase, (short)-1, (short)(i+1), (short)0, this);
				}
			}
		}
		//Reverse strand
		else {
			
			firstamino = (short)0;
			for(int i = this.exons.length-1; i>=0; i--) {	
				exonstart = Integer.parseInt(exonstarts[i]);
				exonend = Integer.parseInt(exonends[i]);
				startphase = Short.parseShort(exonphases[i]);
				
				if(startphase > 0)  {
					startphase = (short)(3-startphase);
				}
				
				if(startphase > -1) {
					endphasetemp = 0;
					if(this.codingEnd > exonstart && this.codingEnd <= exonend) {
						
						firstamino = 1;
						endphase = (short)((this.codingEnd-exonstart)%3);
						if(endphase > 0) {
							endphasetemp = (short)(3-endphase);
						}
						if(this.codingStart > exonstart && this.codingStart < exonend) {
							endphase = 0;
						}
						this.exons[i] = new Exon(exonstart, exonend, startphase, endphase, (short)(this.exons.length-i), firstamino, this);
						firstamino += ((this.codingEnd-exonstart+endphasetemp)/3);
					}
					else {
						if(this.codingStart > exonstart && this.codingStart < exonend) {
							endphase = 0;
						}
						else {
							endphase = (short)(((exonend-startphase)-exonstart)%3);
							if(endphase > 0) {
								endphasetemp = (short)(3-endphase);
							}
						}
						this.exons[i] = new Exon(exonstart, exonend, startphase, endphase, (short)(this.exons.length-i), firstamino, this);
						firstamino += ((exonend-exonstart-startphase+endphasetemp)/3);
					}
				}
				else {
					this.exons[i] = new Exon(exonstart, exonend, startphase, (short)-1, (short)(this.exons.length-i), (short)0, this);
				}
			}
		}
		exonstarts = null;
		exonends = null;
		exonphases = null;
	}
	public Transcript(HashMap<String,String> gffLine, Gene gene) {
		exonArray = new ArrayList<Exon>();
		this.uniprot = FileRead.getInfoValue(gffLine, "uniprot");
		this.biotype = FileRead.getInfoValue(gffLine, "biotype");
		if(gffLine.containsKey("transcript_id")) {
			this.ID = "TranscriptID:" +gffLine.get("transcript_id");
		}
		else {
			this.ID = FileRead.getInfoValue(gffLine, "id");
		}
		
		this.start = Integer.parseInt(FileRead.getInfoValue(gffLine, "start"));
		this.end = Integer.parseInt(FileRead.getInfoValue(gffLine, "end"));
		this.length = this.end -this.start;
		this.gene = gene;
		this.gene.addTranscript(this);
		this.codingStart = Integer.MAX_VALUE;
		this.codingEnd = -1;
	}
	
	
	public void addExon(HashMap<String, String> gffHash, Transcript trans) {
		//ADDAA exoneita arrayhyn... älä heti luo uuutta... kato jos overlappaa edellisten kanssa ja jos löytyy CDS
		if(trans.start > Integer.parseInt(gffHash.get("start"))) {
			trans.start = Integer.parseInt(gffHash.get("start"));
		}
		if(trans.end < Integer.parseInt(gffHash.get("end"))) {
			trans.end = Integer.parseInt(gffHash.get("end"));
		}
		if(exonArray.size() == 0) {
			Exon addExon = new Exon(Integer.parseInt(gffHash.get("start")),Integer.parseInt(gffHash.get("end")), (short)(exonArray.size()+1), trans );
			
			exonArray.add(addExon);
			if(!gffHash.get("phase").equals(".")) {
				if(trans.getStrand()) {
					if(trans.codingStart > Integer.parseInt(gffHash.get("start"))) {
						trans.codingStart = Integer.parseInt(gffHash.get("start"));
					}
					if(trans.codingEnd < Integer.parseInt(gffHash.get("end"))) {
						trans.codingEnd = Integer.parseInt(gffHash.get("end"));
					}
					
				}
				else {
					if(trans.codingStart > Integer.parseInt(gffHash.get("start"))) {
						trans.codingStart = Integer.parseInt(gffHash.get("start"));
					}
					if(trans.codingEnd < Integer.parseInt(gffHash.get("end"))) {
						trans.codingEnd = Integer.parseInt(gffHash.get("end"));
					}
				}				
				
				addExon.startPhase = Short.parseShort(gffHash.get("phase"));
			}
		}
		else if(exonArray.get(0).getStart() > Integer.parseInt(gffHash.get("end"))) {
			Exon addExon = new Exon(Integer.parseInt(gffHash.get("start")),Integer.parseInt(gffHash.get("end")), (short)(exonArray.size()+1), trans );
			exonArray.add(0,addExon);
			
			if(!gffHash.get("phase").equals(".")) {
				if(trans.getStrand()) {
					if(trans.codingStart > Integer.parseInt(gffHash.get("start"))) {
						trans.codingStart = Integer.parseInt(gffHash.get("start"));
					}
					if(trans.codingEnd < Integer.parseInt(gffHash.get("end"))) {
						trans.codingEnd = Integer.parseInt(gffHash.get("end"));
					}
				}
				else {
					if(trans.codingStart > Integer.parseInt(gffHash.get("start"))) {
						trans.codingStart = Integer.parseInt(gffHash.get("start"));
					}
					if(trans.codingEnd < Integer.parseInt(gffHash.get("end"))) {
						trans.codingEnd = Integer.parseInt(gffHash.get("end"));
					}
				}				
				
				addExon.startPhase = Short.parseShort(gffHash.get("phase"));
			}
		}
		else if(exonArray.get(exonArray.size()-1).getEnd() < Integer.parseInt(gffHash.get("start"))) {
			Exon addExon = new Exon(Integer.parseInt(gffHash.get("start")),Integer.parseInt(gffHash.get("end")), (short)(exonArray.size()+1), trans );
			exonArray.add(addExon);
			
			if(!gffHash.get("phase").equals(".")) {
				
					if(trans.getStrand()) {
						if(trans.codingStart > Integer.parseInt(gffHash.get("start"))) {
							trans.codingStart = Integer.parseInt(gffHash.get("start"));
						}
						if(trans.codingEnd < Integer.parseInt(gffHash.get("end"))) {
							trans.codingEnd = Integer.parseInt(gffHash.get("end"));
						}
					}
					else {
						if(trans.codingStart > Integer.parseInt(gffHash.get("start"))) {
							trans.codingStart = Integer.parseInt(gffHash.get("start"));
						}
						if(trans.codingEnd < Integer.parseInt(gffHash.get("end"))) {
							trans.codingEnd = Integer.parseInt(gffHash.get("end"));
						}
					}					
					
				addExon.startPhase = Short.parseShort(gffHash.get("phase"));
			}
		}
		else {
			
			int start = Integer.parseInt(gffHash.get("start"));
			
			for(int i = 0; i<exonArray.size(); i++) {
				
				if(start >= exonArray.get(i).getStart() && start <exonArray.get(i).getEnd()) {
					
					if(!gffHash.get("phase").equals(".")) {
						
						if(trans.getStrand()) {
							if(trans.codingStart > Integer.parseInt(gffHash.get("start"))) {
								trans.codingStart = Integer.parseInt(gffHash.get("start"));
							}
							if(trans.codingEnd < Integer.parseInt(gffHash.get("end"))) {
								trans.codingEnd = Integer.parseInt(gffHash.get("end"));
							}
						}
						else {
							
							if(trans.codingStart > Integer.parseInt(gffHash.get("start"))) {
								trans.codingStart = Integer.parseInt(gffHash.get("start"));
							}
							if(trans.codingEnd < Integer.parseInt(gffHash.get("end"))) {
								trans.codingEnd = Integer.parseInt(gffHash.get("end"));
							}
						}					
						if(!gffHash.get("type").contains("codon")) {
							exonArray.get(i).startPhase = Short.parseShort(gffHash.get("phase"));
						}
					}
					
					
					if(Integer.parseInt(gffHash.get("end")) > exonArray.get(i).getEnd()) {
						
						exonArray.get(i).end = Integer.parseInt(gffHash.get("end"));
						
					}
					break;
				}
			}
			
		}
		
		//Exon addExon = new Exon(Integer.parseInt(gffHash.get("start")),Integer.parseInt(gffHash.get("end")), (short)(exonArray.size()+1), trans );
	
	}
	
	public Gene getGene() {
		return this.gene;
	}
	public void setGene(Gene gene) {
		this.gene = gene;
	}
	public int getStart() {
		return start;
	}
	public int getEnd() {
		return end;
	}
	public int getLength() {
		return this.length;
	}
	public Integer getCodingStart() {
		return codingStart;
	}
	public Integer getCodingEnd() {
		return codingEnd;
	}			
	public boolean isCanonical() {
		return canonical;
	}
	public boolean getStrand() {
		return gene.getStrand();
	}
	public String getChrom() {
		return this.gene.getChrom();
	}
	public String getBiotype() {
		return biotype;
	}
	public String getUniprot() {
		return uniprot;
	}
	public String getGenename() {
		return this.gene.getName();
	}
	public String getDescription() {
		return this.gene.getDescription();
	}
	
	public String getENSG() {
		return this.gene.getID();
	}
	public String getENST() {
		return ID;
	}
	
	public Exon[] getExons() {
		return this.exons;
	}
	
	
}
