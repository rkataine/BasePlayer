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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.Feature;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.tribble.util.TabixUtils;
import java.awt.Graphics2D;
import java.awt.GraphicsConfiguration;
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Locale;
import java.util.Map;
import java.util.zip.GZIPInputStream;

public class MethodLibrary {

	static String formatNumber(int number) {
		
		return NumberFormat.getNumberInstance(Locale.US).format(number);	
	}

	public static class ReadListSorter implements Comparator<Object[]> {

		public int compare(Object[] o1, Object[] o2) {  
	        Object[] f1 = o1;  
	        Object[] f2 = o2;  
	     
	        if ( (Integer)f1[0] < (Integer)(f2[0]) ) {  
	                return -1;  
	        } 
	        else if((Integer)f1[0] < (Integer)(f2[0])) {  
	        		return 1;  
	        }
	        else {
	        	return 0;
	        }        
	       
	}  
	}
	public static class ChromSorter implements Comparator<String> {

		public int compare(String o1, String o2) {  
			Integer n1= null, n2=null;
			if(o1.matches("^-?\\d+$")) { n1 = Integer.parseInt(o1); }
			else if(o1.equals("X")) { n1 = 23; }			
			else if(o1.equals("Y")) { n1 = 24; }			
			else if(o1.equals("M")) { n1 = 25; }			
			else if(o1.equals("MT")) { n1 = 25; }		
			
			if(o2.matches("^-?\\d+$")) { n2 = Integer.parseInt(o2); }
			else if(o2.equals("X")) { n2 = 23; }
			else if(o2.equals("Y")) { n2 = 24; }
			else if(o2.equals("M")) { n2 = 25; }
			else if(o2.equals("MT")) { n2 = 25; }
			
			if(n1 != null && n2 != null) {
			
		        if ( n1 < n2 ) {  
		                return -1;  
		        } 
		        else if(n1 > n2) {  
		        		return 1;  
		        }
		        else {
		        	return 0;
		        }        
			}
			else if (n1 == null && n2 == null) {
				if ( o1.compareTo(o2) < 0 ) {  
	                return -1;  
	        } 
	        else if(o1.compareTo(o2) > 0 ) {  
	        		return 1;  
	        }
	        else {
	        	return 0;
	        }   
				
			}
			else {
				if(n1 != null) {
					return -1;
				}
				else {
					return 1;
				}
			}
	}  
	}
	public static class mateListSorter implements Comparator<ReadNode> {
		
		public int compare(ReadNode o1, ReadNode o2) {  
		
	            
	        if ( o1.split.offset + (o1.getPosition()-o1.split.start)*o1.split.pixel <  o2.split.offset +(o2.getPosition()-o2.split.start)*o2.split.pixel ) {  
	                return -1;  
	        } 
	        else  {  
	        		return 1;  
	        }              
	       
	}  
	}
	public int getLongestDel(VarNode node) {
		int maxlength = 0;
		for(int i = 0 ; i<node.vars.size(); i++) {
			if(node.vars.get(i).getKey().length() == 1) {
				continue;
			}
			if(node.vars.get(i).getKey().startsWith("i")) {
				continue;
			}
			if(node.vars.get(i).getKey().length() > maxlength) {
				maxlength = node.vars.get(i).getKey().length();
			}
		}		
		return 0;
	}
	public static class BEDSorter implements Comparator<String[]> {

		public int compare(String[] o1, String[] o2) {  
	        String[] f1 = o1;  
	        String[] f2 = o2;  
	    //    String[] f3 = o3;  
	        int retStatus=0;  
	        if ( Integer.parseInt(f1[0]) < Integer.parseInt(f2[0]) ) {  
	                retStatus = -1;  
	        } else if ( Integer.parseInt(f1[0]) > Integer.parseInt(f2[0]) ) {  
	                retStatus = 1;  
	        } else if ( Integer.parseInt(f1[1]) < Integer.parseInt(f2[1]) ) {  
	                        retStatus = -1;  
	                } else if ( Integer.parseInt(f1[1]) > Integer.parseInt(f2[1]) ) {  
	                        retStatus = 1;  
	                } else if ( Integer.parseInt(f1[2]) < Integer.parseInt(f2[2])) {  
	                                retStatus = -1;  
	                        } else if( Integer.parseInt(f1[2]) > Integer.parseInt(f2[2]) ){  
	                                retStatus = 1;  
	                        } else {  
	                                retStatus = 0;  
	                        }  
	                  
	         
	        return retStatus;  
	}  
	}
	static BufferedImage toCompatibleImage(BufferedImage image)
	{
		// obtain the current system graphical settings
		GraphicsConfiguration gfx_config = GraphicsEnvironment.
			getLocalGraphicsEnvironment().getDefaultScreenDevice().
			getDefaultConfiguration();

		/*
		 * if image is already compatible and optimized for current system 
		 * settings, simply return it
		 */
		if (image.getColorModel().equals(gfx_config.getColorModel()))
			return image;

		// image is not optimized, so create a new image that is
		BufferedImage new_image = gfx_config.createCompatibleImage(
				image.getWidth(), image.getHeight(), image.getTransparency());

		// get the graphics context of the new image to draw the old image on
		Graphics2D g2d = (Graphics2D) new_image.getGraphics();

		// actually draw the image and dispose of context no longer needed
		g2d.drawImage(image, 0, 0, null);
		g2d.dispose();

		// return the new optimized image
		return new_image; 
	}
	
	public static class controlsorter implements Comparator<SampleNode> {

		@Override
		public int compare(SampleNode s1, SampleNode s2) {

			if(s1.getControlSample().getIndex() < s2.getControlSample().getIndex()) {
				return -1;
			}
			else {
				return 1;
			}
			
			
		}
		
	}
	static void unzip(File infile, File outfile) {
		
			try {
				GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(infile));
				BufferedReader reader = new BufferedReader(new InputStreamReader(gzip));
			    BufferedWriter fastaWriter = new BufferedWriter(new FileWriter(outfile));				
				String line;
				
				while((line = reader.readLine()) != null) {
					if(fastaWriter != null) {
						fastaWriter.write(line +"\n");
					}
				}
				reader.close();
				fastaWriter.close();
			}
			catch(Exception e ) {
				e.printStackTrace();
			}
		
	}
	
	static void blockCompressAndIndex(ArrayList<String[]> in, String bgzfOut, boolean deleteOnExit,SAMSequenceDictionary dict) throws IOException {

	 
	    
	    File outFile= new File(bgzfOut);
	    TabixIndexCreator indexCreator =null;
	   if(dict.getSequences().size() > 300) {
		
		  indexCreator = new TabixIndexCreator(dict,TabixFormat.BED);	  
	   }
	   else {
		  indexCreator = new TabixIndexCreator(TabixFormat.BED);
	   }
	    
	    BlockCompressedOutputStream writer = new BlockCompressedOutputStream(outFile);
	    String header = "#Chrom\tGeneStart\tGeneEnd\tName\tExonCount\tStrand\tENSG\tENST\tUniProt\tCanonical\tBiotype\tCodingStart\tCodingEnd\tExonStarts\tExonEnds\tStartPhases\tDescription\n";
		writer.write(header.getBytes());
	    long filePosition= writer.getFilePointer();
	    OWNCodec bedCodec= new OWNCodec();
	   for(int i = 0 ; i<in.size(); i++) {
	        String[] line = in.get(i);
	       
	        Feature bed = bedCodec.decode(line);
	        if(bed==null) continue;
	        for(int j = 0 ; j<line.length; j++) {
	        	if(line[j] == null) {
	        		System.out.println(j);
	        		continue;
	        	}
	    		writer.write(line[j].getBytes());
	    		if(j < line.length-1) {
	    			writer.write('\t');
	    		}
	    	}
	    	if(i < in.size()-1) {
	    		writer.write('\n');
	    	}
	    	
	        indexCreator.addFeature(bed, filePosition);
	    	
	        filePosition = writer.getFilePointer();
	    }
	 
	    writer.flush();

	    System.err.print("Indexing... ");

	    File tbi= new File(bgzfOut + TabixUtils.STANDARD_INDEX_EXTENSION);
	    if(tbi.exists() && tbi.isFile()){
	        System.err.println("Index file exists: " + tbi);
	        tbi.delete();
	    }
	    Index index = indexCreator.finalizeIndex(writer.getFilePointer());
	    index.writeBasedOnFeatureFile(outFile);
	    
	    writer.close();

	    System.err.println("Done");

	    if(deleteOnExit){
	        outFile.deleteOnExit();
	        File idx= new File(outFile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
	        idx.deleteOnExit();
	    }
	}
	public static boolean isDiscordant(SAMRecord record, boolean cg) {
		
		if(record.getReadPairedFlag() && !record.getMateUnmappedFlag()) {
			if(record.getReferenceIndex() != record.getMateReferenceIndex()) {
				return true;
			}
			else if(record.getMateNegativeStrandFlag() == record.getReadNegativeStrandFlag() && !cg) {
				return true;
			}
			else if (cg && record.getMateNegativeStrandFlag() != record.getReadNegativeStrandFlag()) {
				return true;
			}
			else if (Math.abs(record.getInferredInsertSize()) > Settings.getMaxInsertSize()) {
				return true;			
			}
		}
		
		return false;
	}
	
	public static int getRegion(int position, SplitClass split, int pointer) {
		
		Gene gene = split.getGenes().get(0);	
		
		Boolean inGene = false;
		while(position > gene.getStart()) {
			try {
				pointer++;
			if(position >= gene.getStart() && position <= gene.getEnd()) {	
				
				gene = null;
				return pointer;
						
			}
			if(pointer > split.getGenes().size()-1) {
				return -1;
			}
			gene = split.getGenes().get(pointer);	
			
		}
		catch(Exception e) {
			System.out.println(position);
			e.printStackTrace();
			break;
		}
		}	
		if(!inGene) {		
			gene = null;
			return -1;			
		}
	
		gene = null;
		return -1;
	}
	
	public static int getBaseLength(ArrayList<Map.Entry<String, ArrayList<SampleNode>>> vars, int baselength ) {
		int len = 1;
		for(Map.Entry<String, ArrayList<SampleNode>> entry : vars) {			
			if(entry.getKey().length() < 2 || entry.getKey().startsWith("i")) {
				continue;
			}	
			if(Character.isDigit(entry.getKey().charAt(3))) {
				len = Integer.parseInt(entry.getKey().substring(3));
				if(baselength < len) {					
					baselength = len;
					
				}	
			}
			else {
				continue;
			}
			//f(entry.getKey().substring(3).matches("\\d+")) {
						
			//}
		/*	else {
				if(baselength < entry.getKey().substring(3).length()) {
					baselength = entry.getKey().substring(3).length();
				}
			}*/				
		}
		
		return baselength;	
	}
	public static double round(double value, int places) {
	    if (places < 0) throw new IllegalArgumentException();
	   
	    if(Double.isNaN(value) || Double.isInfinite(value)) {
	    	return 0;
	    }
	    
	    BigDecimal bd = new BigDecimal(value);
	//    System.out.println(bd +" " +(int)-Math.log10(bd.doubleValue()));
	    if(bd.setScale(places, RoundingMode.HALF_UP).doubleValue() == 0.0) {
	    	
	    	return bd.setScale((int)-Math.log10(bd.doubleValue())+places, RoundingMode.HALF_UP).doubleValue();
	    	
	    }
	    else {
	    	return bd.setScale(places, RoundingMode.HALF_UP).doubleValue();
	    }
	    	
	}
	
	public static String reverseComplement(String string) {
		
		StringBuffer newstringBuffer = new StringBuffer(); 
		
		for(int i =string.length()-1; i >= 0; i--) {
			if(string.charAt(i) == 'A') {
				newstringBuffer.append('T');
			//	newString += 'T';
			}
			else if(string.charAt(i) == 'C') {
				//	newString += 'G';
				newstringBuffer.append('G');
			}
			else if(string.charAt(i) == 'G') {
				//	newString += 'C';
				newstringBuffer.append('C');
			}
			else if(string.charAt(i) == 'T') {
				//	newString += 'A';
				newstringBuffer.append('A');
			}
			else if(string.charAt(i) == 'N') {
				//	newString += 'N';
				newstringBuffer.append('N');
			}	
		}	
		return newstringBuffer.toString();
	}
	
	public static String getStrand(boolean value) {
		if(value) {
			return "+";
		}
		else {
			return "-";
		}
	}
	
	public static String aminoEffect(String amino) {
		if(amino.length() < 4) {
			
			return "";
		}
	
		if(!amino.contains(";")) {
			if(amino.contains("fs") || amino.contains("Stop") || (amino.length() == 7 && amino.startsWith("Met1") && !amino.substring(4).equals("Met")) || amino.contains("spl")) {
				return "nonsense";
			}
			if(amino.length() > 6 && amino.substring(0,3).equals(amino.substring(amino.length()-3))) {
				return "synonymous";
			}
			if(amino.contains("UTR")) {
				return "UTR";
			}
			if(!amino.contains("UTR") && !amino.contains("intro") && (amino.contains("if") || !amino.substring(0, 3).equals(amino.substring(amino.length()-3)))) {
				return "missense";
			}			
				return "intronic";
			
		}
		else {
			
			if(amino.contains("fs") || amino.contains("Stop") ||amino.contains("spl") ||  amino.matches(".*Met1\\w+") && !amino.matches(".*Met1Met")) {
				return "nonsense";
			}
			String[] aminoTable = amino.split(";");
			boolean syn = false, utr = amino.contains("UTR");
			for(int i = 0; i< aminoTable.length; i++) {
				
				if(!syn && aminoTable[i].substring(0,3).equals(aminoTable[i].substring(aminoTable[i].length()-3))) {
					
					syn = true;
					continue;
				}
				
				if(!aminoTable[i].contains("UTR") && !aminoTable[i].contains("intro") && (aminoTable[i].contains("if") || !aminoTable[i].substring(0, 3).equals(aminoTable[i].substring(aminoTable[i].length()-3)))) {
					
					return "missense";
				}				
			}
			if(syn) {
				return "synonymous";
			}
			if(utr) {
				return "UTR";
			}
			return "intronic";
			
		}
				
		
	}
	
/*
public static String getSpliceAmino(String chrom, Transcript.Exon exon, boolean phase) {
	try {
		
			Transcript transcript = exon.getTranscript();
			//Forward strand
			if(transcript.getStrand()) {
				if (phase) {
						Transcript.Exon preExon = transcript.getExons()[exon.getNro()-2];
						return getAminoAcid(Main.chromDraw.getSeq(chrom,preExon.getEnd()-preExon.getEndPhase(), preExon.getEnd(), Main.referenceFile) +Main.chromDraw.sequence.substring(exon.getStart()- Main.chromDraw.seqStartPos,(exon.getStart()+exon.getStartPhase()- Main.chromDraw.seqStartPos)));
									
				}
				else {
						Transcript.Exon nextExon = transcript.getExons()[exon.getNro()];
						return getAminoAcid(Main.chromDraw.sequence.substring(exon.getEnd()-exon.getEndPhase()-Main.chromDraw.seqStartPos,exon.getEnd()- Main.chromDraw.seqStartPos) +Main.chromDraw.getSeq(chrom,nextExon.getStart(), nextExon.getStart()+nextExon.getStartPhase(), Main.referenceFile));
									
				}
			}
			//Reverse strand
			else {
				if (phase) {
					//TODO korjaa
					Transcript.Exon preExon = transcript.getExons()[exon.getNro()-2];
					
					return getAminoAcid(reverseComplement(Main.chromDraw.sequence.substring(exon.getEnd()-exon.getStartPhase()-Main.chromDraw.seqStartPos,exon.getEnd()- Main.chromDraw.seqStartPos) +Main.chromDraw.getSeq(chrom,preExon.getStart(), preExon.getStart()+preExon.getEndPhase(), Main.referenceFile)));
								
			}
			else {
					Transcript.Exon nextExon = transcript.getExons()[exon.getNro()];
					return getAminoAcid(reverseComplement(Main.chromDraw.getSeq(chrom,nextExon.getEnd()-nextExon.getStartPhase(), nextExon.getEnd(), Main.referenceFile) +Main.chromDraw.sequence.substring(exon.getStart()- Main.chromDraw.seqStartPos,(exon.getStart()+exon.getEndPhase()- Main.chromDraw.seqStartPos))));
								
			}
			}
		
		
		
		
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	return null;
	
}
	*/
static String getAminoAcid(String codon) {
	return ChromDraw.aminoacids.get(codon);
}

}
