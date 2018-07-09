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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

//import java.util.ArrayList;

public class ReadNode {

		private final int position, startoffset;
		private final int matepos;
		private int insertsize;
		private final Integer readWidth;
		private final short quality;
		java.awt.Rectangle rect;
		private final String readName;
		private final boolean forward;
		private final boolean isDiscordant;
		private final boolean primary;
		private final boolean mateforward;
		private String mateChrom;
		//boolean moved = false;
		SplitClass split;
		ReadNode prev, next;
		int yposition = 0;
		ArrayList<ReadNode> mates; //  = new ArrayList<ReadNode>();
		final String SA;
//		private CigarElement element;
		private final java.util.ArrayList<java.util.Map.Entry<Integer,Byte>> mismatches;
		private final Cigar cigar;
//		private boolean foundsame;
//		private int mispos;
		public ReadNode() {
			this.primary = true;
			this.position = -1;			
			this.rect = null;
			this.readWidth = 0;
			this.readName = null;			
			this.quality = 0;
			this.insertsize = 0;
			this.forward = false;
			SA = null;
			this.matepos = -1;
			this.mateChrom = null;
			this.mateforward = false;
			this.cigar = null;
			this.isDiscordant = false;
			this.mismatches = null;
			this.startoffset = 0;
		}
	/*	public ReadNode(int pos, boolean forward, Cigar cigar, String readname, StringBuffer sequence, int seqStartPos) {
		}
		*/
		public ReadNode(SAMRecord record, boolean cg, String chrom, Sample sample, SplitClass split, Reads readClass,java.util.ArrayList<java.util.Map.Entry<Integer,Byte>> mismatches) {
			
			this.rect = new java.awt.Rectangle();			
			this.primary = !record.getSupplementaryAlignmentFlag();			
			this.readName = record.getReadName();
			
			this.quality = (short)record.getMappingQuality();
			this.mismatches = mismatches;
			if(record.getInferredInsertSize() == 0) {
				this.insertsize = record.getReadLength();
			}
			else {
				this.insertsize = record.getInferredInsertSize();
			}
			this.forward = !record.getReadNegativeStrandFlag();
			SA = record.getStringAttribute("SA");
		
			isDiscordant = MethodLibrary.isDiscordant(record, cg);
			
				
			
			this.split = split;
			if(record.getCigarString().contains("H") || (SA != null/* && !Settings.softClips*/)) {
				
				this.position = record.getAlignmentStart();		
				
				this.readWidth = record.getAlignmentEnd()-this.position+1;
				this.startoffset = record.getAlignmentStart()-record.getUnclippedStart();
				
			}
			else {
				
				this.position = record.getUnclippedStart();	
				this.readWidth = record.getUnclippedEnd()- this.position+1;
				this.startoffset = 0;
			}
			
		//	System.out.println(sample.mates.size());
			
		 //TODO MT & M
			
			if(record.getReadPairedFlag() && !record.getMateUnmappedFlag()) {
				if(sample.hasMates == false) {
					sample.hasMates = true;
				}
				if(record.getMateReferenceName().contains("chr")) {
					if(record.getMateReferenceName().contains("M") ) {
						this.mateChrom = record.getMateReferenceName().replace("chr", "") +"T";
					}
					else {
						this.mateChrom = record.getMateReferenceName().replace("chr", "");
					}
				}
				else {
					this.mateChrom = record.getMateReferenceName();
				}
				this.mateforward = !record.getMateNegativeStrandFlag();
				this.matepos = record.getMateAlignmentStart();
								
				try {
					if(this.isDiscordant && sample.getMates() != null && sample.getMates().size() > 0) {
					
						this.split = split;
					
					if(sample.getMates().containsKey(record.getReadName())) {
						
						this.mates = sample.getMates().get(record.getReadName());
					
							this.mates.add(this);
						
					
					}
				
					}
				}
				catch(Exception e) {
					ErrorLog.addError(e.getStackTrace());
					e.printStackTrace();
				}
			}
			else if(SA != null) {
				
					
				this.mateChrom = null;
				this.mateforward = false;
				this.matepos = -1;
				
				if(sample.getMates() != null) {
					if(sample.getMates().containsKey(record.getReadName())) {							
							this.split = split;	
							this.mates = sample.getMates().get(record.getReadName());
							this.mates.add(this);												
					}
						
				}
				}				
			
			else {
					
				this.mateChrom = null;
				this.mateforward = false;
				this.matepos = -1;
			}
			
			if(record.getCigarLength() > 1) {
				this.cigar = record.getCigar();
				
		
			}
			else {
				cigar = null;
			
			}
			
		}			
	
		void setRectBounds(int x, int y, int width, int height) {
			this.rect.setBounds(x,y,width,height);
		}
		java.awt.Rectangle getRect() {return this.rect;}
		
		int getPosition() {	return this.position;	}
		int getStartOffset() {	return this.startoffset;	}
		int getWidth() {return this.readWidth;	}
		String getName() {	return this.readName;	}
		int getMatePos() {	return this.matepos;	}
		String getMateChrom() {	return this.mateChrom;	}
		
		int getInsertSize() {	return this.insertsize;	}
		boolean isForward() {	return this.forward;	}
		boolean isMateForward() {	return this.mateforward;	}
		ReadNode getPrev() {	return this.prev;	}
		ReadNode getNext() {	return this.next;	}
		ArrayList<ReadNode> getMates() {	return this.mates;	}
		boolean getPrimary() { return this.primary; }
		void addMate(ReadNode mate) {
		
		//	if(mate.getInsertSize() != this.insertsize) {
				this.mates.add(mate);
		//	}
		}
		void setPrev(ReadNode node) {
			this.prev = node;
		}
		void setNext(ReadNode node) {
			this.next = node;
		}	
		boolean isDiscordant() {
			return this.isDiscordant;
		}
	/*	void setDiscordant() {
			this.isDiscordant = true;
		}
		*/
		Cigar getCigar() {
			return this.cigar;
		}
		java.util.ArrayList<java.util.Map.Entry<Integer,Byte>> getMismatches() {
			return this.mismatches;
		}
		Short getMappingQuality() {
			return this.quality;
		}
}
