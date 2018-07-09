package base.BasePlayer;
import htsjdk.samtools.CRAMFileReader;
import htsjdk.samtools.QueryInterval;
//import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.ref.ReferenceSource;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;

public class CoveragePattern {
	static SamReader samFileReader;	
	static CRAMFileReader CRAMReader = null;
	static File indexFile;
	static ReferenceSource reference;
	static RandomAccessFile referenceFile;
	static boolean CRAM = false;
	static Iterator<SAMRecord> bamIterator;
	static File readfile;
	private static SAMRecord samRecord;
	public static void main(String[] args) {
		try {
			if(args.length < 3) {
				System.out.println("Start the script: java -jar CoveragePattern.jar bamfile.bam bedfile.bed.gz reference.fa");
				System.exit(0);
			}
			else {
				  reference = new ReferenceSource(new File(args[2]));
				  referenceFile = new RandomAccessFile(new File(args[2]), "r");
				  BufferedReader reader = null;
				  GZIPInputStream gzip = null;
				  String line;
				  readfile = new File(args[0]);
				  if(args[0].endsWith(".cram")) {
					  CRAM = true;
				  }
				  if(args[1].endsWith(".gz")) {					  
					  gzip = new GZIPInputStream(new FileInputStream(args[1]));
					  reader = new BufferedReader(new InputStreamReader(gzip));							
				  }
				  else {
					  reader = new BufferedReader(new FileReader(args[1]));		
				  }
				  String[] split;
				  reader.readLine();
				  int[] coverages = new int[2039];
				  int start, end, index;
				  while((line = reader.readLine())!=null) {
					  if(line.startsWith("2")) {
						  break;
					  }					  
					  split = line.split("\\t");
					  start = Integer.parseInt(split[1]);
					  end = Integer.parseInt(split[2]);
					  if(split[5].equals("-")) {
						  continue;
					  }
					  bamIterator = getBamIterator(readfile,split[0],start-1000,end+1000);
					  
					  while(bamIterator != null && bamIterator.hasNext()) {	
							
							try {						
							
								samRecord = bamIterator.next(); 							
							
								if(samRecord.getReadUnmappedFlag()) {									
									continue;
								}
								if(samRecord.getReadNegativeStrandFlag()) {
									continue;
								}
								index = samRecord.getUnclippedStart()-start+1000;
								
								if(index < 0) {									
									continue;
								}
								if(index > 2038) {
									
									break;
								}
								coverages[index]++;
								
							}
							catch(Exception e) {
								e.printStackTrace();
							}
					  }
					 
					  
				  }			     
				  for(int i=0; i<coverages.length; i++) {
					  System.out.print(coverages[i] +"\t");
				  }
				  System.out.println();
			}
		  }
		  catch(Exception e) {
			  e.printStackTrace();
		  }		
	}
	static Iterator<SAMRecord> getBamIterator(File readfile, String chrom, int startpos, int endpos) {
		try {	
			if(samFileReader !=null) {
				samFileReader.close();
			}
			if(CRAM) {					
				 if(indexFile == null) {
					 indexFile = new File(readfile +".crai");
				 }
				 CRAMReader = new CRAMFileReader(readfile, indexFile,reference,ValidationStringency.SILENT);
				 if(CRAMReader != null && !CRAMReader.hasIndex()) {			
						
						return null;
					}
			}
			else {
				try {
				samFileReader = SamReaderFactory.make().open(readfile); //new SAMFileReader(readfile);
				if(samFileReader != null && !samFileReader.hasIndex()) {	
					
					return null;
				}
				}
				catch(Exception e) {
					
					e.printStackTrace();
				}
			}	
			
			
			if(CRAM) {
				
				QueryInterval[] interval = { new QueryInterval(CRAMReader.getFileHeader().getSequence(chrom).getSequenceIndex(), startpos, endpos) };				
				Iterator<SAMRecord> value = CRAMReader.query(interval, false);
				
				return value;
			}
			else {
				Iterator<SAMRecord> value = null;
			
				try {
					//SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
					value = samFileReader.queryOverlapping(chrom, startpos, endpos);	
				}
				catch(htsjdk.samtools.SAMFormatException e) {
					e.printStackTrace();
				}
								
				return value;
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		return null;
	}
}
