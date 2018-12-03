package base.BasePlayer;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
//import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

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



public class MultisampleFilter {

	static Hashtable<String, Long[]> chromIndex = new Hashtable<String, Long[]>();
	
	private static byte[] seqresult;
	
	private static StringBuffer seqBuffer;
	
	private static RandomAccessFile chromo;
	private static String resString;
	private static long headerLen;
	private static RandomAccessFile infile;


	public static void main(String[] args) {
		
		
		if(args.length < 2) {
			System.out.println("Give bam file/directory and reference file");
			System.out.println("Example: zcat samples/sample.vcf.gz | java -jar MultisampleFilter.jar bamfile/directory ref/hs37d5.fa > results/sample.vcf");
		}
		else {
			String bamdir = args[0];
			String reffile = args[1];
			//String reffile = "C:/LocalData/rkataine/BasePlayer_dev/BasePlayer/genomes/Homo_sapiens.GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa";
			//String bamdir = "X:/cg8/scripts/SomaticFiltering/GermlineToSomatic/normalpool";
			//String file = "X:/cg8/scripts/SomaticFiltering/GermlineToSomatic/samples/test.vcf";
			setChromDrop(reffile);
			//heterogeneity(new File("X:/cg7/projects/Spomyoma/MY5014_heterogeneity/filtered/Somatic_list.vcf.gz"));
			somaticOrNot(bamdir);		
		
		}
	}
	
	public static void somaticOrNot(String bamDir) {
		try {
			
			List<File> bamList = Collections.synchronizedList(new ArrayList<File>()); 
			/*List<String[]> clusterList = Collections.synchronizedList(new ArrayList<String[]>()); 
			List<String> clusterMutations = Collections.synchronizedList(new ArrayList<String>()); 
			List<String> sampleList = Collections.synchronizedList(new ArrayList<String>()); 
			*/
			File bam = new File(bamDir);
			File[] bams = null;
			if(bam.isDirectory()) {
				bams = bam.listFiles();
			}
			else {
				bams = new File[1];
				bams[0] = bam;
			}
			
			//String[] locus = {};
			//String samples = "";
					
			for(File file: bams) {
				if(file.getAbsolutePath().endsWith(".bam")) {
					//samples += file.getName().substring(0,10)+"\t";
				//	if(file.getName().startsWith("c289") || file.getName().startsWith("c291") || file.getName().startsWith("s1161") || file.getName().startsWith("c115") || file.getName().startsWith("c222") || file.getName().startsWith("c232") || file.getName().startsWith("c286") || file.getName().startsWith("c299") || file.getName().startsWith("c594") || file.getName().startsWith("s1113") ) {
						bamList.add(file);
				//	}
				}
			}
			
			String chrom ="";
	//		BufferedReader bufReader = new BufferedReader(new FileReader("C:/HY-Data/RKATAINE/c182_6780_T_LP6005135-DNA_H01.somatic.vcf"));
			String line, base = "";
			String[] split= null;
			int pos = 0;
			
			int bloodCount = 0;
			boolean infoFound = false, setheader=true;
//			System.out.println("#Chrom\tPosition\tRef\tAlt\tTumor\t" +samples);
			ArrayList<String> headerList = new ArrayList<String>();
			InputStreamReader isReader = new InputStreamReader(System.in);
			BufferedReader bufReader = new BufferedReader(isReader);
			//BufferedReader bufReader = new BufferedReader(new FileReader(vcffile));
			while(true) {
				line = bufReader.readLine();
				if(line == null) {
					break;
				}		
					        
				if(line.startsWith("#")) {
					if(line.startsWith("##INFO")) {
						
						infoFound = true;
					}					
					
					if(line.contains("flanks")) {
						String[] blood = line.split("\\s+");
						bloodCount = Integer.parseInt(blood[blood.length-3]);
						
					}
					headerList.add(line.trim());
					//header += line +"\n";
					continue;
				}
				if(headerList.size() > 0) {
					for (int i = 0 ; i< headerList.size(); i++) {
						if(infoFound && headerList.get(i).startsWith("##INFO") && setheader) {
							System.out.println("##INFO=<ID=FILTER,Number=.,Type=String,Description=\"20bp flanks around variation and variant reads/coverage in " +(bamList.size()+bloodCount) +" blood samples\">");						
							setheader = false;
							infoFound = true;
						}
						if(!infoFound) {
							if(headerList.size() == 1) {
								System.out.println("##INFO=<ID=FILTER,Number=.,Type=String,Description=\"20bp flanks around variation and variant reads/coverage in " +(bamList.size()+bloodCount) +" blood samples\">");						
								infoFound = true;
							}
							else if(i == 1) {
								System.out.println("##INFO=<ID=FILTER,Number=.,Type=String,Description=\"20bp flanks around variation and variant reads/coverage in " +(bamList.size()+bloodCount) +" blood samples\">");						
								infoFound = true;
							}
						}
						System.out.println(headerList.get(i));
					}
					
					headerList.clear();
				}
				split = line.split("\\t");
				
				if(!split[0].equals(chrom)) {
					chrom = split[0];
					
					headerLen = chromIndex.get(chrom)[0];
					
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
				
				if(split[7].contains("FILTER=")) {
					String filterfield = "";
					String infofield = "";
					String[] infoSplit = split[7].split(";");
					for(int i = 0; i < infoSplit.length; i++) {
						if(infoSplit[i].startsWith("FILTER=")) {
							filterfield = infoSplit[i];
						}
						else {
							infofield += infoSplit[i] +";";
						}						
					}
					System.out.print(split[0] +"\t" +split[1] +"\t" +split[2] +"\t" +split[3] +"\t" +split[4] +"\t" +split[5] +"\t" +split[6] +"\t" +infofield +filterfield);
				}
				else {	
					String seq = getSeq(pos-10, pos+11, infile);
					String error = "", foxog = "";
					if(seq.substring(9, 12).equals("CTC") && base.equals("C") || seq.substring(9, 12).equals("GAG") && base.equals("G")) {
						error = "errorCCC";
					}
					if(seq.substring(9, 12).equals("CCG") && base.equals("A") || seq.substring(9, 12).equals("CGG") && base.equals("T")) {
						foxog = "FOXOG";
					}
					String errors = "";
					if(error.length() > 0) {
						errors = ";errorCCC";
					}
					else if(foxog.length() > 0) {
						errors = ";FOXOG";
					}
					seq = seq.substring(0,10).toLowerCase() +seq.substring(10, 11) +seq.substring(11,21).toLowerCase();
				
					System.out.print(split[0] +"\t" +split[1] +"\t" +split[2] +"\t" +split[3] +"\t" +split[4] +"\t" +split[5] +"\t" +split[6] +"\t" +split[7]+errors +";FILTER=" +seq);
				
				}
				for(int i = 0 ; i< bamList.size(); i++) {
					String result = checkReads(chrom, pos, base, bamList.get(i));
					System.out.print(","+result);
				}
				for(int i = 8; i<split.length; i++) {
					System.out.print("\t" +split[i]);
				}
				System.out.println();
			}		
			
			
		}
		catch(Exception e) {
			e.printStackTrace();			
		}
	}

	public static String checkReads(String chromString, int pos, String variation, File bamfile) {
		CigarElement cigar = null;
		SAMRecord samRecord = null;
		
		try {
		
		
			String chrom = chromString;
			SamReader inputSam = SamReaderFactory.make().open(bamfile);
			//SAMFileReader inputSam = new SAMFileReader(bamfile);
			samRecord = new SAMRecord(inputSam.getFileHeader());	 
			int start = pos-200, end = pos +200;
			Iterator<SAMRecord> ite;
			
			/*if(chrom.equals("X")) {
				chrom = "23";
			}
			else if(chrom.equals("Y")) {
				chrom = "24";
			}
			else if(chrom.contains("M")) {
				chrom = "25";
			}*/
			ite = inputSam.queryOverlapping(chrom, start, end);
			//ite = inputSam.iterator(inputSam.getIndex().getSpanOverlapping(Integer.parseInt(chrom)-1, start, end));	
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
				else {
					if(base.length() > 1) {						
						refcount++;					
					}
				}
				/*if(base.length() > 1) {
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
				return "0/0";
				
			}
			else {
				return altcount +"/" +(altcount+refcount);
			}
			
			
		
			
		}
		catch(Exception e) {
			
			e.printStackTrace();
		}
		
		return "";
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
						String result = checkReads(chrom, pos, base, bamList.get(i));						
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
public static String getSeq(int start, int end, RandomAccessFile seqchrom) {
		
		try {
			seqresult = new byte[end-start+100];
			if(seqresult.length > 30000) {
				return "";
			}
		}
		catch(Exception e) {
			e.printStackTrace();
			return "";
		}
		
//		resString ="";
		seqBuffer = new StringBuffer();
		chromo = seqchrom;
		
		try {
		
			if((headerLen+(start)+((start)/60))-1 >= (seqchrom.length() -seqresult.length)) {
				
				seqresult = new byte[(int)(seqchrom.length() - (headerLen+(start)+((start)/60)-1))-1];
			}
			
			chromo.seek((headerLen+(start)+((start)/60))-1);
			chromo.readFully(seqresult);
			
			if(seqresult[0] == 10) {
				chromo.seek((headerLen+(start-1)+((start)/60))-1);
				seqchrom.readFully(seqresult);
			}			
					
			for(int i= 0; i< seqresult.length; i++){
				
				if(seqresult[i] != 10) {
					seqBuffer.append((char)seqresult[i]);
			//		resString+=(char)seqresult[i];					
				}
				
				resString = seqBuffer.toString().toUpperCase();
				if(resString.length() >= end-start) {
					break;
				}			
			}
		}
		catch (Exception e) {
			e.printStackTrace();
		//	System.out.println(e +" getseq");
		}	
		
		if(resString.length() < end-start) {
			return "";
		}
		
		return resString.substring(0, end-start);			
	}
	
}
