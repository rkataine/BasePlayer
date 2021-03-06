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
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.Feature;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;

import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.GraphicsConfiguration;
import java.awt.GraphicsEnvironment;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

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
	public static class varSorter implements Comparator<SampleNode> {

		public int compare(SampleNode o1, SampleNode o2) {			
			if(o1.getSample() == null || o2.getSample() == null) {
				return 1;
			}
			if(o1.getSample().getIndex() < o2.getSample().getIndex()) {
				return -1;
				
			}
			else if(o1.getSample().getIndex() > o2.getSample().getIndex()) {
				return 1;
			}
			else {
				return 0;
			}			 
		}
	}
	public static class ChromSorter implements Comparator<String> {

		public int compare(String o1, String o2) {  
			Long n1= null, n2=null;
			
			if(o1.matches("^-?\\d+$")) { n1 = Long.parseLong(o1); }
			
		
			if(o2.matches("^-?\\d+$")) { n2 = Long.parseLong(o2); }
			
			if(n1 != null && n2 != null) {
			
		        if ( n1 < n2 ) {  
		                return -1;  
		        } 
		        else if(n1 > n2) {  
		        		return 1;  
		        }
		             
			}
			else if (n1 == null && n2 == null) {
				if(o1.length() != o2.length()) {
					if(o1.length() < o2.length()) {
						return -1;
					}
					else if(o1.length() > o2.length()) {
						return 1;
					}
				}
				else {
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
			}
			else {
				if(n1 != null) {
					return -1;
				}
				else {
					return 1;
				}
			}
			return 0;
	}  
	}
	static HashMap<String, Integer> mapChrnameToIndex(SAMFileHeader header) {
		HashMap<String, Integer> chrresult = new HashMap<String, Integer>();
		
		SAMSequenceDictionary dict = header.getSequenceDictionary();
		
		List<SAMSequenceRecord> chrs = dict.getSequences();
		String chr = "";
		for(int i = 0 ; i<chrs.size(); i++) {			
			chr = chrs.get(i).getSequenceName();
			if(Main.chromModel.getIndexOf(chr.replace("chr", "")) > -1) {
				chrresult.put(chr, Main.chromModel.getIndexOf(chr.replace("chr", "")));
				
			}
			else {
				if(chr.replace("chr", "").equals("M")) {
					if(Main.chromModel.getIndexOf("MT") > -1) {
						chrresult.put(chr, Main.chromModel.getIndexOf("MT"));
						
					}
				}
				else if (chr.replace("chr", "").equals("MT")) {
					if(Main.chromModel.getIndexOf("M") > -1) {
						
						chrresult.put(chr, Main.chromModel.getIndexOf("M"));
					}
				}
			}
			
		}
		
		return chrresult;
	}
	public static void showVariantMenu(Component component, final VarNode varOverLap, SampleNode sampleOverLap, int x, int y, final String altbase) {
	String line = "";
	ArrayList<String> annolines = new ArrayList<String>();
	if(sampleOverLap != null) {
		line = FileRead.getVCFLine(varOverLap.getChrom(), varOverLap.getPosition(), varOverLap.getPosition()+1, sampleOverLap.getSample());
	}
	if(Main.drawCanvas.annotationOn) {
		for(int i = 0;i<Main.drawCanvas.sampleList.size(); i++) {
			if(Main.drawCanvas.sampleList.get(i).annotation) {
				String annoline = FileRead.getVCFLine(varOverLap.getChrom(), varOverLap.getPosition(), varOverLap.getPosition()+1, Main.drawCanvas.sampleList.get(i));
				
				if(annoline != null) {
					annolines.add(Main.drawCanvas.sampleList.get(i).getName() +"XX" +annoline);
				}
			}
		}
	}
	if(line != null && line.length() > 1) {
		
		JPopupMenu menu = new JPopupMenu();
		JTextArea area = new JTextArea();
		JScrollPane menuscroll = new JScrollPane();
		JButton hg19 = new JButton("Show variant in VarSome (hg19)");
		
		area.setFont(Main.menuFont);
		menu.add(hg19);
		menu.add(menuscroll);		
		//menu.setPreferredSize(new Dimension(300, 300));
		//area.setMaximumSize(new Dimension(300, 600));
		area.setLineWrap(true);
		area.setWrapStyleWord(true);
		final String[] split = line.split("\t");
		hg19.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent arg0) {	
				Main.gotoURL("https://varsome.com/variant/hg19/" +split[0] +"-" +split[1] +"-" +split[3] +"-" +split[4]);
				
			}
			
		});
		boolean first = true;
		if(varOverLap.getBedHits() != null) {
			
			for(int i = 0; i<varOverLap.getBedHits().size(); i++) {
				if(varOverLap.getBedHits().get(i).getTrack().selex) {
					if(varOverLap.getBedHits().get(i).getTrack().getAffinityBox().isSelected()) {
						if(first) {
							area.append("Affinity Changes:\n");
							first = false;
						}
						Main.drawCanvas.baseHover = altbase;
						Double value = MethodLibrary.calcAffiniyChange(varOverLap,Main.drawCanvas.baseHover,varOverLap.getBedHits().get(i));
						area.append(MethodLibrary.shortName(varOverLap.getBedHits().get(i).name, 7, true) +"=" +MethodLibrary.round(varOverLap.getBedHits().get(i).value,3) +" (" +MethodLibrary.round(value,3) +")\n");
					}
				}
			}
		}
		if(!first) {
			area.append("\n");
		}
		int longest = 0;
		String longestText = "";
		try {
			
			area.append("VCF info: " +sampleOverLap.getSample().getName() +"\n\n");
			area.append("POS: " +split[0] +":" +split[1]+"\nID: "+split[2]+"\n");
			area.append("ALT: " +split[3] +" > " +split[4] +"\n");
			area.append("QUAL: " +split[5] +"\n");
			area.append("FILTER: " +split[6] +"\n");
			area.append("INFO:\n");
			if(split.length > 7) {
				String[] infosplit = split[7].split(";");
				for(int i = 0 ; i< infosplit.length; i++) {
					if(infosplit[i].length() > 30) {
						infosplit[i] = infosplit[i].substring(0, 30) +"...";
					}
					if(longest < infosplit[i].length()) {
						longest = infosplit[i].length();
						longestText = infosplit[i];
					}
					area.append(infosplit[i] +"\n");
				}
			}
			
			if(annolines.size() != 0) {
				for(int i = 0;i<annolines.size();i++) {
					String[] temp = annolines.get(i).split("XX");
					String[] annosplit = temp[1].split("\t");
					if(annosplit.length > 7) {
						area.append("\n" +temp[0] +"\n");
						String[] infosplit = annosplit[7].split(";");
						
						for(int j = 0 ; j< infosplit.length; j++) {
							if(infosplit[j].length() > 30) {
								infosplit[j] = infosplit[j].substring(0, 30) +"...";
							}
							if(longest < infosplit[j].length()) {
								longest = infosplit[j].length();
								longestText = infosplit[j];
							}
							area.append(infosplit[j] +"\n");
						}
					}
				}
			}
			
			if(split.length > 8) {
				area.append("\nFORMAT: " +split[8] +"\n");
			}
			
		if(split.length > 9) {
			for(int i = 9 ; i<split.length; i++) {
				area.append(split[i] +"\n");
			}		
		}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
		menu.setMaximumSize(new Dimension((int)(Main.width*0.6), 300));
		menu.setPreferredSize(new Dimension(Main.bedCanvas.buf.getFontMetrics().stringWidth(longestText)+100, 300));
		
		area.setCaretPosition(0);
		area.revalidate();
		menuscroll.getViewport().add(area);
		menu.pack();
		menu.show(component, x,y);		
		
	}
	else {
		JPopupMenu menu = new JPopupMenu();
		JTextArea area = new JTextArea();
		//JScrollPane menuscroll = new JScrollPane();
		JButton hg19 = new JButton("Show variant in VarSome (hg19)");
		
		area.setFont(Main.menuFont);
		menu.add(hg19);
		
		
		hg19.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent arg0) {	
				Main.gotoURL("https://varsome.com/variant/hg19/" +varOverLap.getChrom() +"-" +(varOverLap.getPosition()+1) +"-" +Main.getBase.get(varOverLap.getRefBase()) +"-" +altbase);
				
			}
			
		});
		
		menu.pack();
		menu.show(component, x,y);		
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
	static void blockCompressAndIndex(String[] in,BlockCompressedOutputStream writer) throws IOException {		
		       
		        for(int j = 0 ; j<in.length; j++) {
		        	if(in[j] == null) {
		        		System.out.println(j);
		        		continue;
		        	}
		    		writer.write(in[j].getBytes());
		    		if(j < in.length-1) {
		    			writer.write('\t');
		    		}
		    	}
		    	
		    	writer.write('\n');	    
		
	}
	
	static void blockCompressAndIndex(File infile, File outFile, boolean deleteOnExit) throws IOException {
		GZIPInputStream gzip = null;
		
		  
	    gzip = new GZIPInputStream(new FileInputStream(infile));
	    BufferedReader reader = new BufferedReader(new InputStreamReader(gzip));	
			  
		  
		
	    TabixIndexCreator indexCreator =null;
	    indexCreator = new TabixIndexCreator(TabixFormat.BED);
	    BlockCompressedOutputStream writer = new BlockCompressedOutputStream(outFile);
	    //String header = "#Chrom\tGeneStart\tGeneEnd\tName\tExonCount\tStrand\tENSG\tENST\tUniProt\tCanonical\tBiotype\tCodingStart\tCodingEnd\tExonStarts\tExonEnds\tStartPhases\tDescription\n";
		String header = "#chrom\tstart\tend\tname\tscore\tstrand\n";
	    writer.write(header.getBytes());
		String bedline;
	    long filePosition= writer.getFilePointer();
	    OWNCodec bedCodec= new OWNCodec();
	   //for(int i = 0 ; i<in.size(); i++) {
	    while((bedline = reader.readLine()) != null) { 
	    	String[] line = bedline.split("\t");
	       
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
	    	//if(i < in.size()-1) {
	    	writer.write('\n');
	    	//}
	    	
	        indexCreator.addFeature(bed, filePosition);
	    	
	        filePosition = writer.getFilePointer();
	    }
	    reader.close();
	    writer.flush();

	    System.err.print("Indexing... ");

	    File tbi= new File(outFile.getCanonicalPath() + TabixUtils.STANDARD_INDEX_EXTENSION);
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
	static void blockCompressAndIndexVCF(File infile, File outFile, boolean deleteOnExit) throws IOException {
		GZIPInputStream gzip = null;
		
		  
	    gzip = new GZIPInputStream(new FileInputStream(infile));
	    BufferedReader reader = new BufferedReader(new InputStreamReader(gzip));	
			  
		  
		
	    TabixIndexCreator indexCreator =null;
	    indexCreator = new TabixIndexCreator(TabixFormat.VCF);
	    BlockCompressedOutputStream writer = new BlockCompressedOutputStream(outFile);
	    //String header = "#Chrom\tGeneStart\tGeneEnd\tName\tExonCount\tStrand\tENSG\tENST\tUniProt\tCanonical\tBiotype\tCodingStart\tCodingEnd\tExonStarts\tExonEnds\tStartPhases\tDescription\n";
		//String header = "#chrom\tstart\tend\tname\tscore\tstrand\n";
	    //writer.write(header.getBytes());
		String bedline;
	    long filePosition= writer.getFilePointer();
	    VCFHeader vcfheader = new VCFHeader();
		VCFHeaderLine headerline = new VCFHeaderLine("format","##fileformat=VCFv4.1");
		vcfheader.addMetaDataLine(headerline);
		VariantHandler.vcfCodec.setVCFHeader(vcfheader, VCFHeaderVersion.VCF4_1);
	   // OWNCodec bedCodec= new OWNCodec();
	   //for(int i = 0 ; i<in.size(); i++) {
	    while((bedline = reader.readLine()) != null) { 
	    	String[] line = bedline.split("\t");
	       
	    	Feature vcf = VariantHandler.vcfCodec.decode(bedline);			
	        if(vcf==null) continue;
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
	    	//if(i < in.size()-1) {
	    	writer.write('\n');
	    	//}
	    	
	        indexCreator.addFeature(vcf, filePosition);
	    	
	        filePosition = writer.getFilePointer();
	    }
	    reader.close();
	    writer.flush();

	    System.err.print("Indexing... ");

	    File tbi= new File(outFile.getCanonicalPath() + TabixUtils.STANDARD_INDEX_EXTENSION);
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
	static void blockCompressAndIndex(ArrayList<String[]> in, String bgzfOut, boolean deleteOnExit,SAMSequenceDictionary dict) throws IOException {

	 
	    
	    File outFile= new File(bgzfOut);
	    TabixIndexCreator indexCreator =null;
	   if(dict != null && dict.getSequences().size() > 300) {
		
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
			else if (Math.abs(record.getInferredInsertSize()) > Settings.insertSize) {
				return true;			
			}
		}
		
		return false;
	}
	
	public static String getOverlappingGenes(int start, int end, SplitClass split) {
		
		
		if(split.getGenes().size() == 0) {
			return null;
		}
		
		int pointer = 0;
		Gene gene = split.getGenes().get(pointer);		
		Boolean inGene = false;
		ArrayList<String> genes = new ArrayList<String>();
		while(gene.getEnd() < start) {
			if(pointer > split.getGenes().size()-2) {
				pointer--;
				break;
			}
			pointer++;
			
			gene = split.getGenes().get(pointer);	
		}
		//pointer--;
		//gene = split.getGenes().get(pointer);	
		
		while(end > gene.getStart()) {
			//System.out.println(end +" " +gene.getStart() +" " +gene.getName());
			try {
				
			if((start >= gene.getStart() && start <= gene.getEnd()) || (end >= gene.getStart() && end <= gene.getEnd())) {		
				inGene = true;
			}
			if(!genes.contains(gene.getName())) {				
				genes.add(gene.getName());		
			}						
			if(pointer > split.getGenes().size()-2) {
				break;
			}
			pointer++;
			gene = split.getGenes().get(pointer);			
			
		}
		catch(Exception e) {
			System.out.println(start);
			e.printStackTrace();
			break;
		}
		}	
		if(!inGene) {	
			if(pointer > 0) {
				pointer--;
				gene = split.getGenes().get(pointer);	
				genes.add(gene.getName());
				pointer++;
			}			
			if(pointer < split.getGenes().size()) {
				gene = split.getGenes().get(pointer);	
				genes.add(gene.getName());
			}
			
		}		
		StringBuffer gens = new StringBuffer("");
		if(inGene) {
			for(int i = 0 ;i<genes.size(); i++) {			
				if(i > 0) {
					gens.append(",");
				}
				gens.append(genes.get(i));			
			}
		}
		else {
			if(genes.size() == 2) {
				gens.append(genes.get(0) +"..." +genes.get(1));
			}
			else if (pointer == 0) {
				gens.append("..." +genes.get(0));
			}
			else {
				gens.append(genes.get(0) +"...");
			}
		}
		return gens.toString();
	}
	
	public static int getRegion(int position, SplitClass split, int pointer) {
		if(split.getGenes().size() == 0) {
			return -1;
		}
		Gene gene = split.getGenes().get(pointer);		
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
	
	public static int getControlBaseLength(String ref, String alt, int baselength) {
		if(ref.length() == 1) {
			return 0;
		}
		else if(!alt.contains(",")) {
			return ref.length() -alt.length();
		}
		else {
			String[] altsplit = alt.split(",");
			int len = 0;
			String altbase;
			for(int i = 0; i<altsplit.length; i++) {
				altbase = FileRead.getVariant(ref,altsplit[i]);
				if(altbase.length() < 2 || altbase.startsWith("i")) {
					continue;
				}	
				if(Character.isDigit(altbase.charAt(3))) {
					len = Integer.parseInt(altbase.substring(3));
					if(baselength < len) {					
						baselength = len;
						
					}	
				}
			}
			return baselength;
		}
		
		
	}
	public static int map(int value, int low1, int high1, int low2, int high2) {
			if((high1 - low1) == 0) {
				return 0;
			}
			else {
				 return low2 + (high2 - low2) * (value - low1) / (high1 - low1);
			}
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
			
		}
		
		return baselength;	
	}
	
	public static double round(double value, int places) {
	    if (places < 0) throw new IllegalArgumentException();
	   
	    if(Double.isNaN(value) || Double.isInfinite(value)) {
	    	return 0;
	    }
	    
	    BigDecimal bd = new BigDecimal(value);

	    if(bd.setScale(places, RoundingMode.HALF_UP).doubleValue() == 0.0) {
	    	
	    	return bd.setScale((int)-Math.log10(bd.doubleValue())+places, RoundingMode.HALF_UP).doubleValue();
	    	
	    }
	    else {
	    	return bd.setScale(places, RoundingMode.HALF_UP).doubleValue();
	    }	    	
	}
	public static void removeHeaderColumns(Object column) {		
		
		for(int i = VariantHandler.table.geneheader.size()-1; i>0; i--) {
			
			if(VariantHandler.table.geneheader.get(i)[0].equals(column)) {
				if(VariantHandler.table.geneheader.get(i)[0] instanceof ControlFile) {
					VariantHandler.table.geneheader.remove(i);
					VariantHandler.table.geneheader.remove(i);
				}
				else {
					VariantHandler.table.geneheader.remove(i);
				}
				for(int j = i;j<VariantHandler.table.geneheader.size(); j++) {
					VariantHandler.table.geneheader.get(j)[1]=(int)VariantHandler.table.geneheader.get(j-1)[1]+(int)VariantHandler.table.geneheader.get(j-1)[2];
				}
				break;
			}
		}
		VariantHandler.table.repaint();
		for(int i = VariantHandler.clusterTable.geneheader.size()-1; i>0; i--) {
			if(VariantHandler.clusterTable.geneheader.get(i)[0].equals(column)) {
				if(VariantHandler.clusterTable.geneheader.get(i)[0] instanceof ControlFile) {
					VariantHandler.clusterTable.geneheader.remove(i);
					VariantHandler.clusterTable.geneheader.remove(i);
				}
				else {
					VariantHandler.clusterTable.geneheader.remove(i);
				}
				for(int j = i;j<VariantHandler.clusterTable.geneheader.size(); j++) {
					VariantHandler.clusterTable.geneheader.get(j)[1]=(int)VariantHandler.clusterTable.geneheader.get(j-1)[1]+(int)VariantHandler.clusterTable.geneheader.get(j-1)[2];
				}
				break;
			}
		}
		for(int i = VariantHandler.clusterTable.header.size()-1; i>0; i--) {
			if(VariantHandler.clusterTable.header.get(i)[0].equals(column)) {				
					VariantHandler.clusterTable.header.remove(i);				
					for(int j = i;j<VariantHandler.clusterTable.header.size(); j++) {
						VariantHandler.clusterTable.header.get(j)[1]=(int)VariantHandler.clusterTable.header.get(j-1)[1]+(int)VariantHandler.clusterTable.header.get(j-1)[2];
					}
					break;
			}
		}
		VariantHandler.clusterTable.repaint();
		for(int j = 0 ; j<VariantHandler.tables.size(); j++) {
			
			for(int i = VariantHandler.tables.get(j).geneheader.size()-1; i>0; i--) {
				
				if(VariantHandler.tables.get(j).geneheader.get(i)[0].equals(column)) {
					if(VariantHandler.tables.get(j).geneheader.get(i)[0] instanceof ControlFile) {
						VariantHandler.tables.get(j).geneheader.remove(i);
						VariantHandler.tables.get(j).geneheader.remove(i);
					}
					else {
						VariantHandler.tables.get(j).geneheader.remove(i);
					}
					for(int k = i;k<VariantHandler.tables.get(j).geneheader.size(); k++) {
						VariantHandler.tables.get(j).geneheader.get(k)[1]=(int)VariantHandler.tables.get(j).geneheader.get(k-1)[1]+(int)VariantHandler.tables.get(j).geneheader.get(k-1)[2];
					}
					break;
				}
			}
			VariantHandler.tables.get(j).revalidate();
			VariantHandler.tables.get(j).repaint();
		}	
	}
	
public static StringBuffer[] makeTrackArray(ArrayList<VarNode> nodes) {
	StringBuffer[] bedarray = new StringBuffer[Main.bedCanvas.bedTrack.size()];
	VarNode node = null;
		for(int j = 0; j<nodes.size(); j++) {
			node = nodes.get(j);
		if(node.getBedHits() == null) {
			
			return null;
		}
			makeTrackArray(node, null, true);
		}
		node = null;
		return bedarray;	
	}
	
	static String shortString(String stringi, int length) {
		if(stringi == null) {
			return "";
		}
		if(stringi.length() > length) {
			return stringi.substring(0, length) +"...";
		}
		else {
			return stringi;
		}
		
	}
	
	public static StringBuffer[] makeTrackArray(VarNode node, String base, boolean shorten) {
	
		if(node.getBedHits() == null) {
			return null;
		}
	
		StringBuffer[] bedarray = new StringBuffer[Main.bedCanvas.bedTrack.size()];
		for(int v = 0; v < node.getBedHits().size(); v++) {
			
			if(bedarray[node.getBedHits().get(v).getTrack().trackIndex] == null) {
				
				if(base != null && node.getBedHits().get(v).getTrack().basecolumn != null) {
					//System.out.println(node.getBedHits().get(v).getPosition() +":" +node.getBedHits().get(v).name +" " +node.getPosition()+":"+base);
					if(node.getBedHits().get(v).name.equals(base)) {						
						bedarray[node.getBedHits().get(v).getTrack().trackIndex] = new StringBuffer(""+MethodLibrary.round(node.getBedHits().get(v).value,2));
					}					
				}
				else {
					if(node.getBedHits().get(v).name != null && node.getBedHits().get(v).name.length()> 0) {						
						
						bedarray[node.getBedHits().get(v).getTrack().trackIndex] = new StringBuffer(shortName(node.getBedHits().get(v).name,10,shorten));						
						
						if(node.getBedHits().get(v).getTrack().hasvalues) {
							bedarray[node.getBedHits().get(v).getTrack().trackIndex].append("=" +MethodLibrary.round(node.getBedHits().get(v).value,2));
						}					
					}
					else if(node.getBedHits().get(v).getTrack().hasvalues) {
						
						bedarray[node.getBedHits().get(v).getTrack().trackIndex] = new StringBuffer(""+MethodLibrary.round(node.getBedHits().get(v).value,2));
						
					}
					else {
						bedarray[node.getBedHits().get(v).getTrack().trackIndex] = new StringBuffer("hit");
					}
					if(node.getBedHits().get(v).getTrack().selex) {
						if(node.getBedHits().get(v).getTrack().getAffinityBox().isSelected()) {
							Double value = calcAffiniyChange(node, base, node.getBedHits().get(v));
							bedarray[node.getBedHits().get(v).getTrack().trackIndex].append(">"+round(value,3));
						}
					}
				}
			}
			else {
				
				if(base != null && node.getBedHits().get(v).getTrack().basecolumn != null) {
					if(bedarray[node.getBedHits().get(v).getTrack().trackIndex] == null) {
						bedarray[node.getBedHits().get(v).getTrack().trackIndex] = new StringBuffer("-");
					}
					continue;
				}
				else {
					bedarray[node.getBedHits().get(v).getTrack().trackIndex].append(";");
					if(node.getBedHits().get(v).name != null) {
						
						bedarray[node.getBedHits().get(v).getTrack().trackIndex].append(shortName(node.getBedHits().get(v).name,10, shorten));
						//bedarray[node.getBedHits().get(v).getTrack().trackIndex].append(node.getBedHits().get(v).name);
						
						if(node.getBedHits().get(v).getTrack().hasvalues) {
							bedarray[node.getBedHits().get(v).getTrack().trackIndex].append("=" +MethodLibrary.round(node.getBedHits().get(v).value,2));
						}					
					}
					else if(node.getBedHits().get(v).getTrack().hasvalues) {
						bedarray[node.getBedHits().get(v).getTrack().trackIndex].append(""+MethodLibrary.round(node.getBedHits().get(v).value,2));
					}
					else {
						bedarray[node.getBedHits().get(v).getTrack().trackIndex].append("hit");
					}	
					if(node.getBedHits().get(v).getTrack().selex) {
						if(node.getBedHits().get(v).getTrack().getAffinityBox().isSelected()) {
							Double value = calcAffiniyChange(node, base, node.getBedHits().get(v));
							bedarray[node.getBedHits().get(v).getTrack().trackIndex].append(">"+round(value,3));
						}
					}
				}			
			}
		}		
		return bedarray;		
	}
	public static String[] makeMultiAlt(String chrom, int pos, String ref, VarNode node) {
		int longestdel = 0;
		String[] result = new String[2];
		for(int i = 0; i<node.vars.size(); i++) {
			if(node.vars.get(i).getKey().startsWith("del")) {
				if(node.vars.get(i).getKey().matches("del\\d+")) {
					int len = Integer.parseInt(node.vars.get(i).getKey().substring(3));
					if(longestdel < len) {
						longestdel = len;
					}
				}
				else {
					if(longestdel == 0) {
						longestdel = 1;
					}	
				}			
			}
		}
		StringBuffer buffer = new StringBuffer("");
		if(longestdel > 0) {
			if(Main.drawCanvas.splits.get(0).getReference() == null ) {
				Main.drawCanvas.splits.get(0).setReference(new ReferenceSeq());
			}
			result[0] = new String(Main.drawCanvas.splits.get(0).getReference().getSeq(chrom, pos, pos+longestdel+1, Main.referenceFile));
			for(int i = 0; i<node.vars.size(); i++) {
				if(node.vars.get(i).getKey().length() == 1) {
					buffer.append(node.vars.get(i).getKey() +result[0].substring(1) +",");
				}
				else {
					if(node.vars.get(i).getKey().startsWith("ins")) {
						buffer.append(result[0].charAt(0) +"" +node.vars.get(i).getKey().substring(3) +result[0].substring(1) +",");
					}
					else {
						if(node.vars.get(i).getKey().matches("del\\d+")) {
							int len = Integer.parseInt(node.vars.get(i).getKey().substring(3));
							buffer.append(result[0].substring(0, result[0].length()-len)+",");	
						}
						else {
							
							buffer.append(result[0].substring(0, result[0].length()-1)+",");	
						}
					}
				}
			}
		}
		else {
			result[0] = ref;
			
			for(int i = 0; i<node.vars.size(); i++) {
				if(node.vars.get(i).getKey().length() == 1) {
					buffer.append(node.vars.get(i).getKey() +",");
				}
				else {
					buffer.append(ref +node.vars.get(i).getKey().substring(3) +",");
				}
			}
			
		}
		
		result[1] = buffer.toString();
	
		return result;
	}
	
	public static String[] makeIndelColumns(String chrom, int pos, String ref, String indel) {
		String[] result = new String[2];
		
		if(indel.contains("del")) {
			if(Main.drawCanvas.splits.get(0).getReference() == null) {
				Main.drawCanvas.splits.get(0).setReference(new ReferenceSeq());
			}
			if(indel.length() == 4 && !indel.matches(".*\\d+")) {
				result[0] = ref +indel.substring(3);
				result[1] = ref;
			}
			else {
				
				result[0] = new String(Main.drawCanvas.splits.get(0).getReference().getSeq(chrom, pos, pos+Integer.parseInt(indel.substring(3))+1, Main.referenceFile));
				result[1] = ref;
			}
		}
		else {
			if(indel.length() == 4 && !indel.matches(".*\\d+")) {
				result[0] = ref;
				result[1] = ref +indel.substring(3);
			}
			else {
				result[0] = ref;
				
				result[1] = ref +indel.substring(3);
				
			}
		}
		
		return result;	

	}
	
	static void createVCFIndex2(File file) {
		try {
		 BufferedReader reader = new BufferedReader(new FileReader(file));
    	
	   	  String line;
	   	FileReader read = new FileReader(file);
    	int addbyte = 1;
    	int max = 10000, count = 0;
    	while(read.read() != 10) {
    		if(read.read() == 13) {
    			addbyte = 2;
    			break;
    		}
    		if(count > max) {
    			break;
    		}
    		count++;
    	}
    	
    	read.close();
		  TabixIndexCreator indexCreator = null;	
		  SAMSequenceDictionary dict = null;
		
		VCFHeader vcfheader = new VCFHeader();
		VCFHeaderLine headerline = new VCFHeaderLine("format","##fileformat=VCFv4.1");
		vcfheader.addMetaDataLine(headerline);
		VariantHandler.vcfCodec.setVCFHeader(vcfheader, VCFHeaderVersion.VCF4_1);
		long filepointer = 0;
		
		boolean cancelled = false, first = true;
	    while((line = reader.readLine()) != null) {
	    	
	    	if(line.startsWith("#") || line.startsWith("\"") || line.length() < 10) {
	    		filepointer += line.length()+addbyte;
	    		continue;
	    	}
	    	if(first) {
	    		first = false;
	    		if(line.startsWith("chr")) {
	    			dict = AddGenome.ReadDictChr(Main.ref);
	    		}
	    		else {
	    			dict = AddGenome.ReadDict(Main.ref);
	    		}
	    		indexCreator = new TabixIndexCreator(dict, TabixFormat.VCF);	
	    	}
	    	Feature vcf = VariantHandler.vcfCodec.decode(line);			
	    	if(!Main.drawCanvas.loading) {
	    		cancelled = true;
	    		break;
	    	}
	    	
			indexCreator.addFeature(vcf, filepointer);
			filepointer += line.length()+addbyte;
	    }
	   
	   
	    if(!cancelled) {
	    	 Index index = indexCreator.finalizeIndex(filepointer);
	    	 index.writeBasedOnFeatureFile(file);
	    }  
	
	    reader.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}	
	}
	
	
	static File createVCFIndex(File file) {
		try {
			FileInputStream testfile = new FileInputStream(file);
			BufferedInputStream testbuf = new BufferedInputStream(testfile);
			boolean ok = BlockCompressedInputStream.isValidFile(testbuf);
			testbuf.close();
			testfile.close();
			if(!ok) {
				if (JOptionPane.showConfirmDialog(Main.drawScroll, "File is not bgzipped, want to bgzip it now?", "Bgzipped?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE) == JOptionPane.YES_OPTION){
					 Main.drawCanvas.loadingtext = "Bgzipping and creating index for " +file.getName();
					blockCompressAndIndexVCF(file, new File(file.getCanonicalPath().replace(".vcf.gz", "_bgzip.vcf.gz")), false);
					return new File(file.getCanonicalPath().replace(".vcf.gz", "_bgzip.vcf.gz"));
				
				}
				else {
					return null;
				}
				
				
			}
			BlockCompressedInputStream reader = new BlockCompressedInputStream(file);
			
			
		   	String line;
			
			TabixIndexCreator indexCreator;	
			SAMSequenceDictionary dict = AddGenome.ReadDict(Main.ref);
			indexCreator = new TabixIndexCreator(dict, TabixFormat.VCF);	
			VCFHeader vcfheader = new VCFHeader();
			VCFHeaderLine headerline = new VCFHeaderLine("format","##fileformat=VCFv4.1");
			vcfheader.addMetaDataLine(headerline);
			VariantHandler.vcfCodec.setVCFHeader(vcfheader, VCFHeaderVersion.VCF4_1);
			long filepointer = 0;
			
			boolean cancelled = false;
			    while((line = reader.readLine()) != null) {
			    	if(line.startsWith("#")) {
			    	filepointer =  reader.getFilePointer();
			    		continue;
			    	}
			    	Feature vcf = VariantHandler.vcfCodec.decode(line);			
			    	if(!Main.drawCanvas.loading) {
			    		cancelled = true;
			    		break;
			    	}
			    	
					indexCreator.addFeature(vcf, filepointer);
					filepointer = reader.getFilePointer();
			    }	   
		   
		    if(!cancelled) {
		    	 Index index = indexCreator.finalizeIndex(filepointer);
		    	 index.writeBasedOnFeatureFile(file);
		    }	
		    
		    	reader.close();
		    	return file;
		}
		catch(Exception e) {
			e.printStackTrace();
			Main.showError(e.getMessage(), "Error");
			return null;
		}	
		
	}
	static void createBEDIndex(File file) {
		try {
		 BlockCompressedInputStream reader = new BlockCompressedInputStream(file);
		
	   	 String line;
		
		  TabixIndexCreator indexCreator;	
		  SAMSequenceDictionary dict = AddGenome.ReadDict(Main.ref);
		  indexCreator = new TabixIndexCreator(dict, TabixFormat.BED);	
		
	
		BEDCodecMod bedCodec= new BEDCodecMod();		      
		long filepointer = 0;		
		boolean cancelled = false;
		
	    while((line = reader.readLine()) != null) {
	    	try {
	    	if(line.startsWith("#")) {
	    	filepointer =  reader.getFilePointer();
	    		continue;
	    	}
	    	Feature bed = bedCodec.decode(line);	    	
			indexCreator.addFeature(bed, filepointer);
			filepointer = reader.getFilePointer();
	    	}
		    catch(Exception ex) {
		    	ex.printStackTrace();
		    }
	    }	   
	   
	    if(!cancelled) {
	    	 Index index = indexCreator.finalizeIndex(filepointer);
	    	 index.writeBasedOnFeatureFile(file);
	    }  
	
	    reader.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}	
	}
	static void createBEDIndex2(File file) {
		try {
		 BufferedReader reader = new BufferedReader(new FileReader(file));
    	 FileReader read = new FileReader(file);
    	int addbyte = 1;
    	int max = 10000, count = 0;
    	while(read.read() != 10) {
    		if(read.read() == 13) {
    			addbyte = 2;
    			break;
    		}
    		if(count > max) {
    			break;
    		}
    		count++;
    	}
    	
    	read.close();
	   	  String line;
		
		  TabixIndexCreator indexCreator;	
		  SAMSequenceDictionary dict = AddGenome.ReadDict(Main.ref);
		  indexCreator = new TabixIndexCreator(dict, TabixFormat.BED);	
		
	
		BEDCodecMod bedCodec= new BEDCodecMod();		       
	  
		long filepointer = 0;
		
		boolean cancelled = false;
		boolean isChr = false, first = true;
	    while((line = reader.readLine()) != null) {
	    	try {
	    	
	    	if(line.startsWith("#") || line.startsWith("track")) {
	    		filepointer +=  line.length()+addbyte;
	    		continue;
	    	}
	    	if(first) {
	    		if(line.startsWith("chr")) {
	    			isChr = true;
	    		}
	    		first = false;
	    	}
	    	if(isChr) {
	    		try {
	    			Feature bed = bedCodec.decode(line.substring(3));
	    			indexCreator.addFeature(bed, filepointer);
	    		}
	    		catch(Exception e) {
	    			
	    			e.printStackTrace();
	    		}
	    	}
	    	else {
	    		try {
	    			Feature bed = bedCodec.decode(line);
	    			indexCreator.addFeature(bed, filepointer);
	    		}
	    		catch(Exception e) {
	    			filepointer +=  line.length()+addbyte;
	    			continue;
	    		}
	    	}
	    	
	    	
			
			filepointer +=  line.length()+addbyte;
			
	    	}
		    catch(Exception ex) {
		    	
		    	ex.printStackTrace();
		    	break;
		    }
	    }
	   
	   
	    if(!cancelled) {
	    	 Index index = indexCreator.finalizeIndex(filepointer);
	    	 index.writeBasedOnFeatureFile(file);
	    }  
	
	    reader.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}	
	}
	public static String shortName(String line, int maxlength, boolean shorten) {
		if(line.length() > maxlength && shorten) {
			return line.substring(0, maxlength) +"...";
		}
		else {
			return line;
		}
	}
	public static void addHeaderColumns(Object column) {		
		
		if(Main.bedCanvas.bedTrack.size() == 0 || column instanceof BedTrack) {
			Object[] obj = new Object[3]; 
			obj[0] = column; 	
			obj[1] = (int)VariantHandler.table.geneheader.get(VariantHandler.table.geneheader.size()-1)[1] + (int)VariantHandler.table.geneheader.get(VariantHandler.table.geneheader.size()-1)[2];
			obj[2] = 100;
			VariantHandler.table.geneheader.add(obj);
			obj = new Object[3]; 
			obj[0] = column; 	
			obj[1] = (int)VariantHandler.clusterTable.geneheader.get(VariantHandler.clusterTable.geneheader.size()-1)[1] + (int)VariantHandler.clusterTable.geneheader.get(VariantHandler.clusterTable.geneheader.size()-1)[2];
			obj[2] = 100;
			VariantHandler.clusterTable.geneheader.add(obj);
			for(int i = 0 ;i<VariantHandler.tables.size(); i++) {
				if(VariantHandler.tables.get(i).bedtrack.equals(column)) {
					
					for(int j = 0 ;j<VariantHandler.tables.size()-1; j++) {
						obj = new Object[3]; 
						obj[0] = VariantHandler.tables.get(j).bedtrack; 	
						obj[1] = (int)VariantHandler.tables.get(i).geneheader.get(VariantHandler.tables.get(i).geneheader.size()-1)[1] + (int)VariantHandler.tables.get(i).geneheader.get(VariantHandler.tables.get(i).geneheader.size()-1)[2];
						obj[2] = 100;
						VariantHandler.tables.get(i).geneheader.add(obj);
						
					}
					continue;
				}
				obj = new Object[3]; 
				obj[0] = column; 	
				obj[1] = (int)VariantHandler.tables.get(i).geneheader.get(VariantHandler.tables.get(i).geneheader.size()-1)[1] + (int)VariantHandler.tables.get(i).geneheader.get(VariantHandler.tables.get(i).geneheader.size()-1)[2];
				obj[2] = 100;
				VariantHandler.tables.get(i).geneheader.add(obj);
				if(column instanceof ControlFile) {
					obj = new Object[3]; 
					obj[0] = "OR"; 	
					obj[1] = (int)VariantHandler.tables.get(i).geneheader.get(VariantHandler.tables.get(i).geneheader.size()-1)[1] + (int)VariantHandler.tables.get(i).geneheader.get(VariantHandler.tables.get(i).geneheader.size()-1)[2];
					obj[2] = 100;
					VariantHandler.tables.get(i).geneheader.add(obj);
				}
			}
			if(column instanceof ControlFile) {
				obj = new Object[3]; 
				obj[0] = "OR"; 	
				obj[1] = (int)VariantHandler.table.geneheader.get(VariantHandler.table.geneheader.size()-1)[1] + (int)VariantHandler.table.geneheader.get(VariantHandler.table.geneheader.size()-1)[2];
				obj[2] = 100;
				VariantHandler.table.geneheader.add(obj);
				obj = new Object[3]; 
				obj[0] = "OR"; 	
				obj[1] = (int)VariantHandler.clusterTable.geneheader.get(VariantHandler.clusterTable.geneheader.size()-1)[1] + (int)VariantHandler.clusterTable.geneheader.get(VariantHandler.clusterTable.geneheader.size()-1)[2];
				obj[2] = 100;
				VariantHandler.clusterTable.geneheader.add(obj);
				
			}
			else {
				obj = new Object[3]; 
				obj[0] = column; 	
				obj[1] = (int)VariantHandler.clusterTable.header.get(VariantHandler.clusterTable.header.size()-1)[1] + (int)VariantHandler.clusterTable.header.get(VariantHandler.clusterTable.header.size()-1)[2];
				obj[2] = 100;
				VariantHandler.clusterTable.header.add(obj);
				
				if((int)obj[1]+(int)obj[2] > VariantHandler.clusterTable.getWidth()) {
					if(VariantHandler.clusterTable.bufImage.getWidth() < (int)obj[1]+(int)obj[2]) {
						VariantHandler.clusterTable.bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)VariantHandler.clusterTable.width*2, (int)VariantHandler.clusterTable.height, BufferedImage.TYPE_INT_ARGB));	
						VariantHandler.clusterTable.buf = (Graphics2D)VariantHandler.clusterTable.bufImage.getGraphics();
					}					
				}
				VariantHandler.clusterTable.setPreferredSize(new Dimension((int)obj[1]+(int)obj[2], VariantHandler.clusterTable.getHeight()));
				VariantHandler.clusterTable.revalidate();
			}
		}
		else {
			for(int i = VariantHandler.table.geneheader.size()-2; i> 0; i--) {
				if(VariantHandler.table.geneheader.get(i)[0] instanceof BedTrack == false) {
					Object[] obj = new Object[3]; 
					obj[0] = column; 	
					obj[1] = (int)VariantHandler.table.geneheader.get(i)[1] + (int)VariantHandler.table.geneheader.get(i)[2];
					obj[2] = 100;
					VariantHandler.table.geneheader.add(i+1,obj);
					
					obj = new Object[3]; 
					obj[0] = "OR"; 	
					obj[1] = (int)VariantHandler.table.geneheader.get(VariantHandler.table.geneheader.size()-1)[1] + (int)VariantHandler.table.geneheader.get(VariantHandler.table.geneheader.size()-1)[2];
					obj[2] = 100;
					VariantHandler.table.geneheader.add(i+2,obj);
					
					for(int j = i;j<VariantHandler.table.geneheader.size(); j++) {
						VariantHandler.table.geneheader.get(j)[1]=(int)VariantHandler.table.geneheader.get(j-1)[1]+(int)VariantHandler.table.geneheader.get(j-1)[2];
					}
					break;
				}
			}			
			for(int i = VariantHandler.clusterTable.geneheader.size()-2; i> 0; i--) {
				if(VariantHandler.clusterTable.geneheader.get(i)[0] instanceof BedTrack == false) {
					Object[] obj = new Object[3]; 
					obj[0] = column; 	
					obj[1] = (int)VariantHandler.clusterTable.geneheader.get(i)[1] + (int)VariantHandler.clusterTable.geneheader.get(i)[2];
					obj[2] = 100;
					VariantHandler.clusterTable.geneheader.add(i+1,obj);
					
					obj = new Object[3]; 
					obj[0] = "OR"; 	
					obj[1] = (int)VariantHandler.clusterTable.geneheader.get(VariantHandler.clusterTable.geneheader.size()-1)[1] + (int)VariantHandler.clusterTable.geneheader.get(VariantHandler.clusterTable.geneheader.size()-1)[2];
					obj[2] = 100;
					VariantHandler.clusterTable.geneheader.add(i+2,obj);
					
					for(int j = i;j<VariantHandler.clusterTable.geneheader.size(); j++) {
						VariantHandler.clusterTable.geneheader.get(j)[1]=(int)VariantHandler.clusterTable.geneheader.get(j-1)[1]+(int)VariantHandler.clusterTable.geneheader.get(j-1)[2];
					}
					break;
				}
			}	
			for(int t=0;t<VariantHandler.tables.size();t++) {
				for(int i = VariantHandler.tables.get(t).geneheader.size()-2; i> 0; i--) {
					if(VariantHandler.tables.get(t).geneheader.get(i)[0] instanceof BedTrack == false) {
						Object[] obj = new Object[3]; 
						obj[0] = column; 	
						obj[1] = (int)VariantHandler.tables.get(t).geneheader.get(i)[1] + (int)VariantHandler.tables.get(t).geneheader.get(i)[2];
						obj[2] = 100;
						VariantHandler.tables.get(t).geneheader.add(i+1,obj);
						
						obj = new Object[3]; 
						obj[0] = "OR"; 	
						obj[1] = (int)VariantHandler.tables.get(t).geneheader.get(i+1)[1] + (int)VariantHandler.tables.get(t).geneheader.get(i+1)[2];
						obj[2] = 100;
						VariantHandler.tables.get(t).geneheader.add(i+2,obj);
						
						for(int j = i;j<VariantHandler.tables.get(t).geneheader.size(); j++) {
							VariantHandler.tables.get(t).geneheader.get(j)[1]=(int)VariantHandler.tables.get(t).geneheader.get(j-1)[1]+(int)VariantHandler.tables.get(t).geneheader.get(j-1)[2];
						}
						if((int)VariantHandler.tables.get(t).geneheader.get(VariantHandler.tables.get(t).geneheader.size()-1)[1]+(int)VariantHandler.tables.get(t).geneheader.get(VariantHandler.tables.get(t).geneheader.size()-1)[2] > VariantHandler.tables.get(t).getWidth()) {
							if(VariantHandler.tables.get(t).bufImage.getWidth() < (int)VariantHandler.tables.get(t).geneheader.get(VariantHandler.tables.get(t).geneheader.size()-1)[1]+(int)VariantHandler.tables.get(t).geneheader.get(VariantHandler.tables.get(t).geneheader.size()-1)[2] ) {
								VariantHandler.tables.get(t).bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)VariantHandler.tables.get(t).width*2, (int)VariantHandler.tables.get(t).height, BufferedImage.TYPE_INT_ARGB));	
								VariantHandler.tables.get(t).buf = (Graphics2D)VariantHandler.tables.get(t).bufImage.getGraphics();
							}
							VariantHandler.tables.get(t).setPreferredSize(new Dimension((int)VariantHandler.tables.get(t).geneheader.get(VariantHandler.tables.get(t).geneheader.size()-1)[1]+(int)VariantHandler.tables.get(t).geneheader.get(VariantHandler.tables.get(t).geneheader.size()-1)[2] , VariantHandler.tables.get(t).getHeight()));
							
						}
						VariantHandler.tables.get(t).revalidate();
						VariantHandler.tables.get(t).repaint();
						break;
					}
				}		
			}
		}
				
		if((int)VariantHandler.table.geneheader.get(VariantHandler.table.geneheader.size()-1)[1]+(int)VariantHandler.table.geneheader.get(VariantHandler.table.geneheader.size()-1)[2] > VariantHandler.table.getWidth()) {
			if(VariantHandler.table.bufImage.getWidth() < (int)VariantHandler.table.geneheader.get(VariantHandler.table.geneheader.size()-1)[1]+(int)VariantHandler.table.geneheader.get(VariantHandler.table.geneheader.size()-1)[2] ) {
				VariantHandler.table.bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)VariantHandler.table.width*2, (int)VariantHandler.table.height, BufferedImage.TYPE_INT_ARGB));	
				VariantHandler.table.buf = (Graphics2D)VariantHandler.table.bufImage.getGraphics();
			}
			VariantHandler.table.setPreferredSize(new Dimension((int)VariantHandler.table.geneheader.get(VariantHandler.table.geneheader.size()-1)[1]+(int)VariantHandler.table.geneheader.get(VariantHandler.table.geneheader.size()-1)[2] , VariantHandler.table.getHeight()));
			
		}
		if((int)VariantHandler.clusterTable.geneheader.get(VariantHandler.clusterTable.geneheader.size()-1)[1]+(int)VariantHandler.clusterTable.geneheader.get(VariantHandler.clusterTable.geneheader.size()-1)[2] > VariantHandler.clusterTable.getWidth()) {
			if(VariantHandler.clusterTable.bufImage.getWidth() < (int)VariantHandler.clusterTable.geneheader.get(VariantHandler.clusterTable.geneheader.size()-1)[1]+(int)VariantHandler.clusterTable.geneheader.get(VariantHandler.clusterTable.geneheader.size()-1)[2] ) {
				VariantHandler.clusterTable.bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)VariantHandler.clusterTable.width*2, (int)VariantHandler.clusterTable.height, BufferedImage.TYPE_INT_ARGB));	
				VariantHandler.clusterTable.buf = (Graphics2D)VariantHandler.clusterTable.bufImage.getGraphics();
			}
			VariantHandler.clusterTable.setPreferredSize(new Dimension((int)VariantHandler.clusterTable.geneheader.get(VariantHandler.clusterTable.geneheader.size()-1)[1]+(int)VariantHandler.clusterTable.geneheader.get(VariantHandler.clusterTable.geneheader.size()-1)[2] , VariantHandler.clusterTable.getHeight()));
			
		}		
		
		VariantHandler.table.revalidate();
		VariantHandler.table.repaint();
		VariantHandler.clusterTable.revalidate();
		VariantHandler.clusterTable.repaint();
		
	}
	public static int[][] reverseMatrix(int[][] matrix) {
		int[][] newMatrix = new int[4][matrix[0].length]; 
		int pointer = 0;
	
		for(int i = matrix[0].length-1; i>=0; i--) {
			
			newMatrix[3][pointer] = matrix[0][i];
			newMatrix[2][pointer] = matrix[1][i];
			newMatrix[1][pointer] = matrix[2][i];
			newMatrix[0][pointer] = matrix[3][i];
			pointer++;
		}
		
		return newMatrix;
		
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
			
			if(!amino.contains("UTR") && (amino.contains("fs") || (amino.length() == 7 && amino.startsWith("Met1") && !amino.substring(4).equals("Met")) || amino.contains("spl"))) {
				
				return "nonsense";
			}
			if(amino.length() > 6 && amino.substring(0,3).equals(amino.substring(amino.length()-3))) {
				return "synonymous";
			}
			if(amino.contains("UTR")) {
				return "UTR";
			}
			if(amino.contains("Stop") ) {
				if(amino.startsWith("Stop") && amino.endsWith("Stop")) {
					return "synonymous";
				}
				else {
					return "nonsense";					
				}
			
			}
			if(!amino.contains("UTR") && !amino.contains("intro") && (amino.contains("if") || !amino.substring(0, 3).equals(amino.substring(amino.length()-3)))) {
				return "missense";
			}			
			
				return "intronic";
			
		}
		else {
			
			if((amino.contains("fs") || amino.contains("spl") ||  (amino.matches(".*Met1\\D+.*") && !amino.matches(".*Met1Met.*")))) {
				return "nonsense";
			}
			String[] aminoTable = amino.split(";");
			boolean syn = false, utr = amino.contains("UTR"), nonsense = false;
			for(int i = 0; i< aminoTable.length; i++) {
				
				if(!syn && aminoTable[i].substring(0,3).equals(aminoTable[i].substring(aminoTable[i].length()-3))) {
					
					syn = true;
					continue;
				}
				if( amino.contains("Stop")) {
					if(amino.startsWith("Stop") && amino.endsWith("Stop")) {
						syn = true;
					}
					else {
						return "nonsense";
					}
					
				}
				try {
				if(!aminoTable[i].contains("UTR") && !aminoTable[i].contains("intro") && (aminoTable[i].contains("if") || !aminoTable[i].substring(0, 3).equals(aminoTable[i].substring(aminoTable[i].length()-3)))) {
					
					return "missense";
				}	
				}
				catch(Exception e) {
					System.out.println(aminoTable[i]);
					System.out.println(e);
				}
				
			}
			if(nonsense) {
				return "nonsense";
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
	static Double calcAffiniyChange(VarNode node, String alt, BedNode bednode) {
		if(alt == null || alt.length() > 1) {
			return 0.0;
		}			
		if(bednode.getTrack().selex) {
			
			int index = node.getPosition()-bednode.getPosition();
			
			if(index < 0) {
				
				return 0.0;
			}
			int[][] matrix=null;
			if(bednode.forward) {
				matrix = Main.SELEXhash.get(bednode.id);
			}
			else {
				matrix = MethodLibrary.reverseMatrix(Main.SELEXhash.get(bednode.id));
			}					
			
			double sum;										
			Double value;					
			double mutatedvalue;
			sum = matrix[0][index] + matrix[1][index] + matrix[2][index] + matrix[3][index];
		//	System.out.println(matrix[0][index] +" " +matrix[1][index] +" " +matrix[2][index] +" " +matrix[3][index]);
			value = matrix[Main.baseMap.get(node.getRefBase())-1][index]/(double)sum;
			
			mutatedvalue = matrix[Main.baseMap.get((byte)alt.charAt(0))-1][index]/(double)sum;
			
			value = value*Math.log(value/Main.background.get(node.getRefBase()));
			if(mutatedvalue != 0) {
				mutatedvalue = mutatedvalue*Math.log(mutatedvalue/Main.background.get((byte)alt.charAt(0)));
			}
			
			
			
			return mutatedvalue-value;
			//System.out.println(mutatedvalue-value);
								
		}
		else {
			return 0.0;
		}
	
}
static String getAminoAcid(String codon) {
	return ChromDraw.aminoacids.get(codon);
}

}
