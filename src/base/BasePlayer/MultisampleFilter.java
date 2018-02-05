package base.BasePlayer;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;


@SuppressWarnings({ "unused", "deprecation" })

public class MultisampleFilter {

	static Hashtable<String, Long[]> chromIndex = new Hashtable<String, Long[]>();
	
	private static byte[] seqresult;
	
	private static StringBuffer seqBuffer;
	
	private static RandomAccessFile chromo;
	private static String resString;
	private static int headerLen;
	private static RandomAccessFile infile;
	private static boolean flankfound = false;

	public static void main(String[] args) {
		
		
		if(args.length < 2) {
			System.out.println("Give bam directory and reference file");
			System.out.println("Example: zcat samples/sample.vcf.gz | java -jar MultisampleFilter.jar all ref/hs37d5.fa > results/sample.vcf");
		}
		else {
		String bamdir = args[0];
		String reffile = args[1];
		
		setChromDrop(reffile);
		heterogeneity(new File("X:/cg7/projects/Spomyoma/MY5014_heterogeneity/filtered/Somatic_list.vcf.gz"));
		somaticOrNot(bamdir);		
		
		}
	}
	
	public static void somaticOrNot(String bamDir) {
		try {
			
			List<File> bamList = Collections.synchronizedList(new ArrayList<File>()); 
			List<String[]> clusterList = Collections.synchronizedList(new ArrayList<String[]>()); 
			List<String> clusterMutations = Collections.synchronizedList(new ArrayList<String>()); 
			List<String> sampleList = Collections.synchronizedList(new ArrayList<String>()); 
			File bam = new File(bamDir);
			
			File[] bams = bam.listFiles();
			String[] locus = {};
			String samples = "";
					
			for(File file: bams) {
				if(file.getAbsolutePath().endsWith(".bam")) {
					samples += file.getName().substring(0,10)+"\t";
				//	if(file.getName().startsWith("c289") || file.getName().startsWith("c291") || file.getName().startsWith("s1161") || file.getName().startsWith("c115") || file.getName().startsWith("c222") || file.getName().startsWith("c232") || file.getName().startsWith("c286") || file.getName().startsWith("c299") || file.getName().startsWith("c594") || file.getName().startsWith("s1113") ) {
						bamList.add(file);
				//	}
				}
			}
			
			String chrom ="";
	//		BufferedReader bufReader = new BufferedReader(new FileReader("C:/HY-Data/RKATAINE/c182_6780_T_LP6005135-DNA_H01.somatic.vcf"));
			String line, base = "", basetemp="", header=""; //, info = "";
			String[] split= null;
			int postemp = 0, pos = 0, samplecount = 0, windowLength = 0, chromcounter = 0;
			double mutrate = 0;
			boolean found = false;
//			System.out.println("#Chrom\tPosition\tRef\tAlt\tTumor\t" +samples);
			
			InputStreamReader isReader = new InputStreamReader(System.in);
			BufferedReader bufReader = new BufferedReader(isReader);
			boolean setheader = true;
			while(true) {
				line = bufReader.readLine();
				if(line == null) {
					break;
				}					
		               
		        
				if(line.startsWith("#")) {
					if(line.startsWith("##INFO") && setheader) {
						header += "##INFO=<ID=FILTER,Number=.,Type=String,Description=\"20bp flanks around variation and variant reads/coverage in 10 blood samples\">\n";
						setheader = false;
					}
					header += line +"\n";
					continue;
				}
				if(header.length() > 1) {
					System.out.print(header);
					header = "";
				}
				split = line.split("\\t");
				if(split[0].equals("23")) {
					split[0] = "X";
				}
				else if(split[0].equals("24")) {
					split[0] = "Y";
				}
				else if(split[0].contains("25")) {
					split[0] = "MT";
				}				
				
				if(!split[0].equals(chrom)) {
					chrom = split[0];
					
					if(chrom.equals("23")) {
						chrom = "X";
					}
					else if(chrom.equals("24")) {
						chrom = "Y";
					}
					else if(chrom.contains("25")) {
						chrom = "MT";
					}					
					chromcounter = 0;
				}
				
				
				//info= split[10];
				pos = Integer.parseInt(split[1]);
				
				if(split[3].length() > 1) {
					base = "DEL" +(split[3].length()-1);
					pos = Integer.parseInt(split[1])+1;
				}
				else if(split[4].length() > 1) {
					base = "INS" +(split[4].length()-1);
					pos = Integer.parseInt(split[1])+1;
				}
				else {					
					base = ""+split[4];					
				}				
			//	checkReads(chrom, pos, base, bamList, base, "", split);
			}		
		}
		catch(Exception e) {
			e.printStackTrace();			
		}
	}
	public static void heterogeneity(File multivcf) {
		try {
			
			List<File> bamList = Collections.synchronizedList(new ArrayList<File>());			
			String[] locus = {};
			String samples = "";
			File[] addDir =multivcf.getParentFile().listFiles(new FilenameFilter() {
 	    	     public boolean accept(File dir, String name) {
 	    	        return name.toLowerCase().endsWith(".bam");
 	    	     }
 	    	});		
			
			String outfile = multivcf.getCanonicalPath().replace(".vcf.gz", "") +".tsv";
			BufferedWriter writer = new BufferedWriter(new FileWriter(outfile));
			String chrom ="";

			String line, base = "", basetemp=""; //, info = "";
			StringBuffer header = new StringBuffer(""), row = new StringBuffer("");
			String[] split= null;
			int postemp = 0, pos = 0, samplecount = 0, windowLength = 0, chromcounter = 0;
			
			GZIPInputStream gzip = null;
			BufferedReader reader = null;
		 
			  try {
				  gzip = new GZIPInputStream(new FileInputStream(multivcf));
				  reader = new BufferedReader(new InputStreamReader(gzip));		
			  }
			  catch(Exception e) {
				e.printStackTrace();					
			  }
			boolean setheader = true;
			while((line = reader.readLine()) !=null) {				
				if(setheader) {
					if(line.startsWith("##")) {											
						continue;
					}
					if(line.startsWith("#CHROM")) {
						split = line.split("\t");
						header.append(split[0] +"\t" +split[1] +"\t" +split[2] +"\t" +split[3] +"\t" +split[4] +"\t" +split[5]);
						for(int i = 0 ; i<split.length-9; i++) {
							header.append("\t"+split[9+i]);
							for(File file: addDir) {
								if(split[9+i].contains(file.getName().subSequence(0,file.getName().indexOf(".")))) {
									bamList.add(file);		
								}										
							}
						}
						writer.write(header.toString() +"\n");
						//System.out.println(header.toString());
						setheader = false;
						continue;
					}
				}				
				
				split = line.split("\t");
								
				pos = Integer.parseInt(split[1]);
				chrom = split[0];
				if(split[3].length() > 1) {
					base = "DEL" +(split[3].length()-1);
					pos = Integer.parseInt(split[1])+1;
				}
				else if(split[4].length() > 1) {
					base = "INS" +(split[4].length()-1);
					pos = Integer.parseInt(split[1])+1;
				}
				else {					
					base = ""+split[4];					
				}
				row = new StringBuffer(split[0] +"\t" +split[1] +"\t" +split[2] +"\t" +split[3] +"\t" +split[4] +"\t" +split[5]);
				for(int i = 0 ; i<split.length-9; i++) {
					if(split[9+i].equals("0/0")) {
						String result = checkReads(chrom, pos, base, bamList.get(i), base, "", split);						
						row.append("\t"+result);
					}
					else {
						row.append("\t" +split[9+i].split(":")[2]);						
					}					
				}
				writer.write(row.toString() +"\n");
				System.out.println(row.toString());			
			}		
			writer.close();
		}
		catch(Exception e) {
			e.printStackTrace();		
			
		}
		
	}
	public static String checkReads(String chromString, int pos, String variation, File bamfile, String ref, String info, String[] split) {
		CigarElement cigar = null;
		SAMRecord samRecord = null;
		
		try {
		
			boolean indelfound = false, flankfound=false, first = true;
			String chrom = chromString;
			
			SAMFileReader inputSam = new SAMFileReader(bamfile);
			samRecord = new SAMRecord(inputSam.getFileHeader());	 
			int start = pos-200, end = pos +200;
			Iterator<SAMRecord> ite;
			
			if(chrom.equals("X")) {
				chrom = "23";
			}
			else if(chrom.equals("Y")) {
				chrom = "24";
			}
			else if(chrom.contains("M")) {
				chrom = "25";
			}
			ite = inputSam.iterator(inputSam.getIndex().getSpanOverlapping(Integer.parseInt(chrom)-1, start, end));	
			int poscount = 0;
			int readPos = 0, position = 0;
			String base = variation;
			
			int altcount = 0, refcount = 0;
			
			while(ite.hasNext() && ite != null) {
				
				samRecord =ite.next(); 	
				
				if(samRecord.getUnclippedEnd() < pos) {
					continue;
				}
				if(samRecord.getUnclippedStart() > pos) {
					break;
				}
		/*		if(!flankfound && (pos - samRecord.getUnclippedStart() > 11 && samRecord.getUnclippedEnd() -pos > 11)) {
					readseq = samRecord.getReadString().substring((pos-samRecord.getUnclippedStart())-10, (pos-samRecord.getUnclippedStart())+11);
					readseq = readseq.substring(0,10).toLowerCase() +readseq.substring(10, 11) +readseq.substring(11,21).toLowerCase();
					flankfound = true;
					
				}
			*/
				position = 0;
				poscount = 0;
				indelfound = false;
				
				if(samRecord.getCigar().numCigarElements() > 1) {		
					if(base.length() > 1) {
						if(samRecord.getCigarString().contains("I") || samRecord.getCigarString().contains("D")) {
							altcount++;
						}
						else {
							refcount++;
						}
					}
					else {
					cigar = null;
					for(int k = 0; k<samRecord.getCigar().numCigarElements(); k++) {
						
						cigar = samRecord.getCigar().getCigarElement(k);
						if ((samRecord.getUnclippedStart() + position) > pos) {								
							break;
						}
						if(cigar.getOperator().compareTo(CigarOperator.MATCH_OR_MISMATCH)== 0) {
							position += cigar.getLength();
						}
						else if(cigar.getOperator().compareTo(CigarOperator.DELETION) == 0) {	
							
				//			System.out.println(bamfile.getName().substring(0, 5) +" " +samRecord.getUnclippedStart() +" " +samRecord.getCigarString() +" " +pos +" " +(samRecord.getUnclippedStart() + position));												
						/*	if(base.length() > 1 && pos == samRecord.getUnclippedStart() + position) {
								indelfound = true;
								altcount++;
								break;
							}*/
							position += cigar.getLength();
							poscount -=	cigar.getLength();	
						}
						else if(cigar.getOperator().compareTo(CigarOperator.INSERTION) == 0) {		
						/*	if(base.length() > 1 && samRecord.getUnclippedStart() + position == pos) {
								altcount++;
								indelfound = true;
								break;
							}*/
							poscount += cigar.getLength();
						}
						else if(cigar.getOperator().compareTo(CigarOperator.SOFT_CLIP) == 0) {
							position+=cigar.getLength()+1;
						}
					
					}
					}
				}
			/*	if(base.length() > 1) {
					if(!indelfound) {
						refcount++;
					}
				}*/
				if(base.length() == 1) {
					readPos = (pos - samRecord.getUnclippedStart()) + poscount;							
					
								
					try {
						base = Character.toString(samRecord.getReadString().charAt(readPos));
					}
					catch(Exception e) {
						
						continue;
					}
					
					if(base.equals(variation)) {
						altcount++;				
					}
					else {
						refcount++;
					}
				}
				/*
				if(count > 0) { 
					samplecount++;
					count = 0;
					if(samplecount == 1) {
						inputSam.close();
						return true;
					}
					break;
				}
				*/
			}
			inputSam.close();
			
			if(refcount == 0 && altcount == 0) {
				return "0,0";
				
			}
			else {
				return refcount +"," +altcount;
			}
			
			
		
			
		}
		catch(Exception e) {
			
			e.printStackTrace();
		}
		
		return "";
	}	

	static void setChromDrop(String dir) {
		try {
		File chromindex = null;
		
			infile = new RandomAccessFile(new File(dir), "r");
			java.util.List<String> chromnamevector = Collections.synchronizedList(new ArrayList<String>());
			chromindex = new File(dir +".fai");
			BufferedReader reader = new BufferedReader(new FileReader(chromindex));
			String line;
			String[] split;
			 
			while((line = reader.readLine()) != null) {
				split = line.split("\t");
			//	chromnamevector.add(split[0]);
				Long[] add = {Long.parseLong(split[2]), Long.parseLong(split[1]),  Long.parseLong(split[3])};
				
				chromIndex.put(split[0], add);
			}
			
			
			reader.close();
		
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
}
