package base.BasePlayer;

import java.io.RandomAccessFile;

public class ReferenceSeq {

	byte[] seq;
	int startpos = 0, endpos = 0;
	String chrom;
	RandomAccessFile ref;
	
	public ReferenceSeq(String chrom, int start, int end, RandomAccessFile file) {
		this.seq = getSeq(chrom, start, end, file);		
		this.startpos = start;
		this.endpos = end;
		this.ref = file;
		this.chrom = chrom;		
	}
	
	public byte[] getSeq() {
		return this.seq;
	}
	public int getStartPos() {
		return this.startpos;
	}
	public int getEndPos() {
		return this.endpos;
	}
	void appendToStart(int pos) {
		byte[] first = getSeq(chrom, pos, startpos, ref);		
		byte[] combined = new byte[first.length+seq.length];
		System.arraycopy(first,0,combined,0,first.length);
		System.arraycopy(seq,0,combined,first.length,seq.length);
		seq = combined;		
		this.startpos = pos;
	}
	public void append(int pos) {
		byte[] second = getSeq(chrom, this.endpos, pos, ref);		
		byte[] combined = new byte[second.length+seq.length];
		System.arraycopy(seq,0,combined,0,seq.length);
		System.arraycopy(second,0,combined,seq.length,second.length);
		seq = combined;				
		this.endpos = pos;
		
	}
	public byte[] getSeq(String chrom, int start, int end, RandomAccessFile seqchrom) {
		try {
			
		byte[] seqresult = new byte[(end-start+1)+((end-start)/(Main.chromIndex.get(Main.refchrom +chrom)[2].intValue()-1))];
		byte[] temp = new byte[end-start];
		
		seqchrom.seek((Main.chromIndex.get(Main.refchrom +chrom)[0]+(start)+((start)/Main.chromIndex.get(Main.refchrom +chrom)[2].intValue())));
		seqchrom.readFully(seqresult);
		
		if(seqresult[0] == 10) {
			
			seqchrom.seek((Main.chromIndex.get(Main.refchrom +chrom)[0]+(start-1)+((start)/Main.chromIndex.get(Main.refchrom +chrom)[2].intValue())));
			seqchrom.readFully(seqresult);			
		}	
		int pointer =0;
		for(int i = 0 ; i<seqresult.length;i++) {
			if(pointer > temp.length-1) {				
				break;
			}
			if(seqresult[i] != 10) {
				temp[pointer] = seqresult[i];
				pointer++;
			}			
		}
		seqresult = null;
		return temp;
		}
		catch(Exception e) {
			
			e.printStackTrace();
		}
		return null;
	}
	
}
