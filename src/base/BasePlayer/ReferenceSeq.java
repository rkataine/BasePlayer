package base.BasePlayer;

import java.io.RandomAccessFile;

public class ReferenceSeq {

	byte[] seq;
	int startpos = 0, endpos = 0;
	String chrom;
	RandomAccessFile ref;
	static boolean wait = false;
	public ReferenceSeq(String chrom, int start, int end, RandomAccessFile file) {
		
		wait = true;
		this.startpos = start;
		this.endpos = end;
		this.seq = getSeq(chrom, start, end, file);		
	
		
		this.ref = file;
		this.chrom = chrom;
		wait = false;
	}
	public ReferenceSeq() {
		
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
		try {
			while(wait) {			
				Thread.sleep(100);			
			}
				Thread.sleep(0);
		}
		catch(Exception ex) {
			ex.printStackTrace();
		}
		if(pos >= this.startpos) {
			return;
		}
		
		wait = true;
		byte[] first = getSeq(chrom, pos, startpos, ref);		
		byte[] combined = new byte[first.length+seq.length];
		System.arraycopy(first,0,combined,0,first.length);
		System.arraycopy(seq,0,combined,first.length,seq.length);
		seq = combined;		
		this.startpos = pos;
		wait = false;
	}
	public void append(int pos) {
		try {
			while(wait) {			
				Thread.sleep(100);			
			}
				Thread.sleep(0);
		}
		catch(Exception ex) {
			ex.printStackTrace();
		}
		if(pos <= this.endpos) {
			
			return;
		}		
	
		wait = true;
		byte[] second = getSeq(chrom, this.endpos, pos, ref);		
		byte[] combined = new byte[second.length+seq.length];
		System.arraycopy(seq,0,combined,0,seq.length);
		System.arraycopy(second,0,combined,seq.length,second.length);
		seq = combined;				
		this.endpos = pos;
		wait = false;
	}
	public byte[] getSeq(String chrom, int start, int end, RandomAccessFile seqchrom) {
		try {
			if(chrom == null || end-start < 0) {
				return null;
			}
			if(!Main.chromIndex.containsKey(Main.refchrom +chrom)) {
				if(!Main.chromIndex.containsKey(chrom.replace("chr", ""))) {
					return null;
				}
				else {
					chrom = chrom.replace("chr", "");
				}			
			}		
			if(start < 0) {
				start = 0;
				this.startpos = 0;
			}		
			byte[] seqresult = null;
			try {
				seqresult = new byte[(end-start+1)+((end-start)/(Main.chromIndex.get(Main.refchrom +chrom)[2].intValue()-1))];
			}
			catch(Exception e) {
				System.out.println((end-start) +" " +Main.chromIndex.get(Main.refchrom +chrom)[2].intValue());
				e.printStackTrace();
				seqresult = new byte[(end-start+1)+((end-start)/40)];
			}
			byte[] temp = new byte[end-start];
			
			seqchrom.seek((Main.chromIndex.get(Main.refchrom +chrom)[0]+(start)+((start)/Main.chromIndex.get(Main.refchrom +chrom)[2].intValue())));
			if(seqchrom.getFilePointer() + seqresult.length > seqchrom.length()) {
				seqresult = new byte[(int)((seqchrom.length())-(seqchrom.getFilePointer()))];
			}
			seqchrom.readFully(seqresult);
			if(seqresult.length == 0) {
				byte[] result = {0};
				return result;
			}
			if(seqresult[0] == 10) {
				
				seqchrom.seek((Main.chromIndex.get(Main.refchrom +chrom)[0]+(start+1)+((start)/Main.chromIndex.get(Main.refchrom +chrom)[2].intValue())));
				if(seqchrom.getFilePointer() + seqresult.length > seqchrom.length()-1) {
					seqresult = new byte[(int)((seqchrom.length())-(seqchrom.getFilePointer()))];
				}
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
