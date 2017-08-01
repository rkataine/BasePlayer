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
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.zip.GZIPInputStream;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.SwingWorker;
import javax.swing.UIManager;
import javax.swing.table.DefaultTableModel;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;

public class AddGenome  extends JPanel implements ActionListener, MouseListener {

	private static final long serialVersionUID = 1L;
	static JFrame frame = new JFrame("Add new genome");
	static boolean annotation = false;
	JLabel genomeFileText;
	JLabel annotationFileText;
	static JTextField genomeName;	
	JButton openRef, openAnno, add, checkUpdates;
	static JButton download, remove;
	static int longestName = 0;
	static String userDir = new File(Main.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParent().replace("%20", " ");
	File genomeFile, annotationFile;
	static boolean downloading = false;
	static HashMap<String, URL[]> genomeHash = new HashMap<String, URL[]>();
	static ArrayList<String> removables = new ArrayList<String>();
	static Object[] headers = {"Download genomes"}, remheaders = {"Installed genomes"}; 
    static String[][] data = {}, remdata = {}; 
    static HashMap<String, Integer[]> sizeHash = new HashMap<String, Integer[]>();
    static JPanel panel = new JPanel(new GridBagLayout());
	
	static DefaultTableModel model = new DefaultTableModel(data, headers) {
	 
		private static final long serialVersionUID = 1L;

		@Override
		
	    public boolean isCellEditable(int row, int column) {
	      
	       return false;
	    }
	};
	static DefaultTableModel remmodel = new DefaultTableModel(remdata, remheaders) {
		 
		private static final long serialVersionUID = 1L;

		@Override
		
	    public boolean isCellEditable(int row, int column) {
	      
	       return false;
	    }
	};
	static JTable genometable = new JTable(model);
	static JTable remtable = new JTable(remmodel);
	
	
	
	static void checkGenomes() {
		
		File checkdir = new File(userDir +"/genomes/");
		File[] addDir;
		for(int i = model.getRowCount() - 1; i >=0; i--) {
		   model.removeRow(i); 
		}
		for(int i = remmodel.getRowCount() - 1; i >=0; i--) {
			remmodel.removeRow(i); 
		}	
		
		removables.clear();
		
		int currentlen = 0, length = 0;
  		addDir = checkdir.listFiles();
  		for(int f = 0; f<addDir.length; f++) { 	
  			currentlen = genometable.getFontMetrics(genometable.getFont()).stringWidth(addDir[f].getName());
			if(currentlen > length) {
				length = currentlen;				
			}	
			AddGenome.removables.add(addDir[f].getName());
			Object[] row = {addDir[f].getName()};	
  			remmodel.addRow(row);
  		}	  		
	  	
	  	
	  	Iterator<Entry<String, URL[]>> iterator = genomeHash.entrySet().iterator();
		ArrayList<String> organisms = new ArrayList<String>();
		
		while(iterator.hasNext()) {
			organisms.add(iterator.next().getKey());			
		}
		Collections.sort(organisms);
		for(int i = 0 ; i< organisms.size(); i++) {
			if(!AddGenome.removables.contains(organisms.get(i))) {
				Object[] row = {organisms.get(i)};			
				model.addRow(row);			
				currentlen = genometable.getFontMetrics(genometable.getFont()).stringWidth(organisms.get(i));
				if(currentlen > length) {
					length = currentlen;				
				}	
			}
		}
		
		AddGenome.longestName = length;
		
		if(genometable.getRowCount() > 15) {
			genometable.setPreferredScrollableViewportSize((new Dimension(AddGenome.longestName+20,genometable.getRowHeight()*15)));
		}
		else {			
			genometable.setPreferredScrollableViewportSize((new Dimension(AddGenome.longestName+20,genometable.getRowHeight()*(genometable.getRowCount()+1))));
		}
		if(remtable.getRowCount() > 15) {
			remtable.setPreferredScrollableViewportSize((new Dimension(AddGenome.longestName+20,remtable.getRowHeight()*15)));
		}
		else {			
			remtable.setPreferredScrollableViewportSize((new Dimension(AddGenome.longestName+20,remtable.getRowHeight()*(remtable.getRowCount()+1))));
		}
		frame.pack();
	}
	
	
	public AddGenome() {
		super(new BorderLayout());	
		makeGenomes();
		
		checkGenomes();
		genomeFileText = new JLabel("Select reference fasta-file");
		annotationFileText = new JLabel("Select annotation gff3-file");
		genomeName = new JTextField("Give name of the genome");	
		openRef = new JButton("Browse");	
		openAnno = new JButton("Browse");
		add = new JButton("Add");
		download = new JButton("Download");
		remove = new JButton("Remove");
		checkUpdates = new JButton("Check updates");
		download.setEnabled(false);
		remove.setEnabled(false);
		download.addActionListener(this);
		remove.addActionListener(this);
		panel.setBackground(Draw.sidecolor);
		checkUpdates.addActionListener(this);
		this.setBackground(Draw.sidecolor);
		frame.getContentPane().setBackground(Draw.sidecolor);
		GridBagConstraints c = new GridBagConstraints();
		c.anchor = GridBagConstraints.NORTHWEST;	
		c.fill = GridBagConstraints.BOTH;
		c.gridx = 0;
		c.gridy = 0;
		c.insets = new Insets(2,4,2,4);		
		c.gridwidth = 2;
		genometable.setSelectionMode(0);
		genometable.setShowGrid(false);		
		remtable.setSelectionMode(0);
		remtable.setShowGrid(false);	
		JScrollPane scroll = new JScrollPane();		
		scroll.getViewport().setBackground(Color.white);
		scroll.getViewport().add(genometable);
		JScrollPane remscroll = new JScrollPane(remtable);
		
		remscroll.getViewport().setBackground(Color.white);
		
		genometable.addMouseListener(this);
		remtable.addMouseListener(this);
		panel.add(new JLabel("Manage genomes"), c);
		c.gridy++;		
		//c.fill = GridBagConstraints.NONE;
		panel.add(scroll,c);		
		c.gridy++;
		c.fill = GridBagConstraints.NONE;
		panel.add(download,c);	
		c.gridy++;
		c.fill = GridBagConstraints.BOTH;
		panel.add(remscroll,c);
		c.gridy++;
		
		c.fill = GridBagConstraints.NONE;
		c.gridwidth = 1;
		panel.add(remove, c);
		c.gridx = 1;
		panel.add(checkUpdates, c);
		checkUpdates.setEnabled(false);
		c.gridwidth = 2;
		c.gridx = 0;
		c.fill = GridBagConstraints.BOTH;
		c.gridy++;
		panel.add(new JLabel("Add genome manually"),c);
		c.gridy++;
		
	
		c.gridwidth = 2;
		
		panel.add(new JSeparator(),c);
		c.gridwidth = 1;
		c.gridy++;
		panel.add(genomeFileText, c);
		c.gridx = 1;
		panel.add(openRef, c);		
		
		c.gridx = 0;
		openRef.addActionListener(this);
		c.gridy++;
		panel.add(annotationFileText,c);
		c.gridx=1;
		panel.add(openAnno, c);
		c.gridy++;
		panel.add(add,c);
		openAnno.addActionListener(this);
		add.addActionListener(this);
		add.setEnabled(false);
		add(panel,BorderLayout.NORTH);
		
		
	/*	html.append("<a href=http:Homo_sapiens_GRCh37:Ensembl_genes> Homo sapiens GRCh37 with Ensembl</a> or <a href=http:Homo_sapiens_GRCh37:RefSeq_genes>RefSeq</a> gene annotations<br>");
		html.append("<a href=http:Homo_sapiens_GRCh38:Ensembl_genes> Homo sapiens GRCh38 with Ensembl</a> or <a href=http:Homo_sapiens_GRCh38:RefSeq_genes>RefSeq</a> gene annotations<br><br>");
		
		html.append("<a href=http:Mus_musculus_GRCm38:Ensembl_genes> Mus musculus GRCm38 with Ensembl</a> or <a href=http:Mus_musculus_GRCm38:RefSeq_genes>RefSeq</a> gene annotations<br>");
		html.append("<a href=http:Rattus_norvegicus:Ensembl_genes> Rattus norvegicus with Ensembl gene annotations</a><br>");
		html.append("<a href=http:Saccharomyces_cerevisiae:Ensembl_genes> Saccharomyces cerevisiae with Ensembl gene annotation</a><br>");					
		html.append("<a href=http:Ciona_intestinalis:Ensembl_genes> Ciona intestinalis with Ensembl gene annotation</a><br>");
		Object[] row = {"Homo_sapiens_GRCh37"};
		Object[] row = {"Homo_sapiens_GRCh38"};
		
		model.addRow(row);
		/*	genomeName.setPreferredSize(new Dimension(300,20));
		this.add(genomeName);
		this.add(new JSeparator());
		
		this.add(openRef);
		openRef.addActionListener(this);
		this.add(genomeFileText);
		this.add(openAnno);
		openAnno.addActionListener(this);
		this.add(annotationFileText);
		this.add(add);
		add.addActionListener(this);
		if(annotation) {
			openRef.setVisible(false);
			genomeFileText.setVisible(false);
			genomeName.setEditable(false);
		}
		genomeFileText.setEditable(false);
		annotationFileText.setEditable(false);*/
	}
	void makeGenomes() {
		try {
			URL[] urls = {new URL("ftp://ftp.ensembl.org/pub/grch37/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"),
						  new URL("ftp://ftp.ensembl.org/pub/grch37/update/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.gff3.gz"),
						  new URL("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz")};
			Integer[] sizes = new Integer[2];
			sizes[0] = 900;
			sizes[1] = 40;
			sizeHash.put("Homo_sapiens_GRCh37", sizes);
			genomeHash.put("Homo_sapiens_GRCh37", urls);
			
		/*	URL[] urls2 = {new URL("ftp://ftp.ensembl.org/pub/grch37/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"),
					  new URL("ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.37.3/GFF/ref_GRCh37.p5_top_level.gff3.gz"),
					  new URL("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz") };
			sizes = new Integer[2];
			sizes[0] = 900;
			sizes[1] = 35;
			sizeHash.put("Homo_sapiens_GRCh37:RefSeq_genes", sizes);
			genomeHash.put("Homo_sapiens_GRCh37:RefSeq_genes", urls2);
			*/
			URL[] urls3 = { new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/ciona_intestinalis/dna/Ciona_intestinalis.KH.dna.toplevel.fa.gz"),
						new URL("ftp://ftp.ensembl.org/pub/release-85/gff3/ciona_intestinalis/Ciona_intestinalis.KH.85.gff3.gz")
			};
			sizes = new Integer[2];
			sizes[0] = 34;
			sizes[1] = 5;
			sizeHash.put("Ciona_intestinalis", sizes);
			genomeHash.put("Ciona_intestinalis", urls3);
			
			URL[] urls4 = { new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"),
					new URL("ftp://ftp.ensembl.org/pub/release-85/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.85.gff3.gz")
			};
			sizes = new Integer[2];
			sizes[0] = 4;
			sizes[1] = 1;
			sizeHash.put("Saccharomyces_cerevisiae", sizes);
			genomeHash.put("Saccharomyces_cerevisiae", urls4);
			
			URL[] urls5 = {new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz"),
					  new URL("ftp://ftp.ensembl.org/pub/release-85/gff3/mus_musculus/Mus_musculus.GRCm38.85.gff3.gz"),
					  new URL("http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/cytoBand.txt.gz")};
			
			sizes = new Integer[2];
			sizes[0] = 801;
			sizes[1] = 58;
			sizeHash.put("Mus_musculus_GRCm38", sizes);
			genomeHash.put("Mus_musculus_GRCm38", urls5);
			
		/*	URL[] urls6 = {new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz"),
					  new URL("ftp://ftp.ncbi.nlm.nih.gov/genomes/M_musculus/GFF/ref_GRCm38.p4_top_level.gff3.gz") };
			sizes = new Integer[2];
			sizes[0] = 801;
			sizes[1] = 28;
			sizeHash.put("Mus_musculus_GRCm38:RefSeq_genes", sizes);
			genomeHash.put("Mus_musculus_GRCm38:RefSeq_genes", urls6);
			*/
			URL[] urls7 = {new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz"),
					  new URL("ftp://ftp.ensembl.org/pub/release-85/gff3/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.85.gff3.gz") };
			sizes = new Integer[2];
			sizes[0] = 812;
			sizes[1] = 15;
			sizeHash.put("Rattus_norvegicus", sizes);
			genomeHash.put("Rattus_norvegicus", urls7);
			
			
			URL[] urls8 = {new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
					  new URL("ftp://ftp.ensembl.org/pub/release-85/gff3/homo_sapiens/Homo_sapiens.GRCh38.85.gff3.gz"),
					  new URL("http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz")};
			sizes = new Integer[2];
			sizes[0] = 860;
			sizes[1] = 38;
			sizeHash.put("Homo_sapiens_GRCh38", sizes);
			genomeHash.put("Homo_sapiens_GRCh38", urls8);
		
			/*URL[] urls9 = {new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
					  new URL("ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/GFF/ref_GRCh38.p7_top_level.gff3.gz"),
					  new URL("http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz")};
			sizes = new Integer[2];
			sizes[0] = 860;
			sizes[1] = 41;
			sizeHash.put("Homo_sapiens_GRCh38:RefSeq_genes", sizes);
			genomeHash.put("Homo_sapiens_GRCh38:RefSeq_genes", urls9);
		*/
		}
		catch(Exception ex) {
			ex.printStackTrace();
		}
	}
	public static class OutputRunner extends SwingWorker<String, Object> {
		File genomefile, annotationFile;
		String genomeName, urls;
		
		public OutputRunner(String genomeName, File genomefile, File annotationFile) {
			this.genomefile = genomefile;
			this.genomeName = genomeName;
			this.annotationFile = annotationFile;
			
		}
		public OutputRunner(String genomeName, File genomefile, File annotationFile, String urls) {
			this.genomefile = genomefile;
			this.genomeName = genomeName;
			this.annotationFile = annotationFile;
			this.urls = urls;
		}
		protected String doInBackground() {
			if(Main.drawCanvas != null) {
				Main.drawCanvas.ready("all");
				Main.drawCanvas.loading("Processing files...");
			}
			
			createGenome(genomeName, genomefile, annotationFile);
			if(Main.drawCanvas != null) {
				Main.drawCanvas.ready("Processing files...");
			}
			return "";
		}
		void createGenome(String genomeName, File genomeFile, File annotationFile) {
			
			try {				
				
			File annofile = null;
			if(Main.drawCanvas != null) {
				Main.drawCanvas.splits.get(0).getGenes().clear();
			}
			Boolean ok = false;
			String annotationUrl = "";		
			File fastatest = new File("test");
			
			if(this.urls != null) {			
				URL[] urlList = AddGenome.genomeHash.get(urls);
				URL fastafile= urlList[0];
				String targetDir = userDir +"/genomes/" +urls +"/";
				fastatest = new File(targetDir +FilenameUtils.getName(fastafile.getFile()).substring(0,FilenameUtils.getName(fastafile.getFile()).indexOf(".gz") ));
				
				if(!new File(targetDir).exists()) {
					File fasta = new File(targetDir +FilenameUtils.getName(fastafile.getFile()));
				//	fastatest = new File(targetDir +FilenameUtils.getName(fastafile.getFile()).substring(0,FilenameUtils.getName(fastafile.getFile()).indexOf(".gz") ));
					new File(targetDir).mkdir();
					if(!fasta.exists() && !fastatest.exists()) {
						Main.drawCanvas.loadingtext ="Downloading " +genomeFile.getName();
						Main.downloadFile(fastafile, targetDir, AddGenome.sizeHash.get(urls)[0]);				
					}
				}
				annotationUrl = annotationFile.getName();
				URL gfffile= urlList[1];
				targetDir = userDir +"/genomes/" +urls +"/annotation/";
				Main.drawCanvas.loadingtext ="Downloading " +gfffile.getFile();
				String filetest = Main.downloadFile(gfffile, targetDir, AddGenome.sizeHash.get(urls)[1]);				
				if(!filetest.equals(FilenameUtils.getName(gfffile.getFile()))) {
					annotationUrl = filetest;
					annotationFile = new File( userDir +"/genomes/" +urls +"/annotation/" +filetest +"/" +filetest);
				}				
				if(Main.drawCanvas != null) {
					Main.drawCanvas.loadingtext ="Downloading " +filetest;
				}
				else {
					System.out.println("Downloading " +filetest);
				}
				
				if(urlList.length == 3) {
					URL bandFile = urlList[2];
					targetDir = userDir +"/genomes/" +urls +"/";
					Main.downloadFile(bandFile, targetDir, 0);
					File infile = new File(targetDir+FilenameUtils.getName(bandFile.getFile()));
					File outfile = new File(targetDir+"bands.txt");
					MethodLibrary.unzip(infile, outfile);
				}
				
			}
			
			if(urls == null) {
			
				ok = new File(userDir +"/genomes/" +genomeName +"/").mkdir();
				
			}
			else {
				ok = true;
			}
			
			if(ok) {
				String genomedir = userDir +"/genomes/" +genomeName +"/";
				if(!fastatest.exists()) {
					if(Main.drawCanvas != null) {
						Main.drawCanvas.loadingtext = "Unpacking and indexing: " +genomeFile.getName();			
					}
					else {
						System.out.println("Unpacking and indexing: " +genomeFile.getName());
					}
							
					indexFasta(genomeFile, genomedir);
					Main.addGenomeFile(genomeName);
				}
			
				if(annotationFile != null) {
					if(Main.drawCanvas != null) {
						Main.drawCanvas.loadingtext = "Unpacking and indexing: " +annotationFile.getName();
					}
					else {
						System.out.println("Unpacking and indexing: " +annotationFile.getName());
					}
					
					if(urls == null) {
						ok = new File(genomedir +"annotation/" +annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gff"))+"/").mkdirs();
					}
					else {
						ok = true;
					}
					
					if(urls != null) {
						annofile = new File(genomedir +"annotation/" +annotationUrl +"/"+annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gff")) +".bed.gz");
						
					}
					else {
						annofile = new File(genomedir +"annotation/"+annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gff")) +"/" +annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gff")) +".bed.gz");
						
					}
					if(ok) {
					
						
						parseGFF(annotationFile, annofile, genomeName, genomedir +genomeFile.getName());
					
						Main.addAnnotationFile(genomeName, annofile);
					}
					else {
						
						Main.addAnnotationFile(genomeName, annofile);
					}
					
				}
				if(urls != null) {
					genomeFile.delete();
					annotationFile.delete();
					
				}
				Main.defaultGenome = genomeName;
				Main.defaultAnnotation = annofile.getName();	
				if(Main.drawCanvas != null) {
					Main.setChromDrop(genomeName);
					Main.getBands();
					Main.geneDropdown.setSelectedItem(Main.defaultAnnotation);
					Main.getExons();
					Main.chromosomeDropdown.setSelectedIndex(0);
					AddGenome.downloading = false;
				}			
				else {
					System.out.println("Genome ready: " +genomeFile.getName() +" : " +Main.defaultAnnotation);
				}
			   
			    checkGenomes();
				}
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		
	}

static String downloadAnnotation(URL fileurl, String genome) {
	try {
		
		String annotationUrl = fileurl.getFile();	
		String targetDir = userDir +"/genomes/" +genome +"/annotation/";
		File dir = new File(userDir +"/genomes/" +genome);
		File[] files = dir.listFiles();
		String genomestring = "";
		for(int i = 0 ; i < files.length; i++) {
			if(files[i].getName().endsWith(".fa") || files[i].getName().endsWith(".fasta")) {
				genomestring = files[i].getCanonicalPath();
				break;
			}
		}
		if(Main.drawCanvas != null) {
			Main.drawCanvas.loadingtext = "Downloading " +annotationUrl;
		}
		else {
			System.out.println("Downloading " +annotationUrl);
		}
		String filetest = Main.downloadFile(fileurl, targetDir, AddGenome.sizeHash.get(genome)[1]);	
		File annotationFile = new File( userDir +"/genomes/" +genome +"/annotation/" +filetest +"/" +filetest);
		String genomedir = userDir +"/genomes/" +genome+"/";
		System.out.println(AddGenome.genomeHash.get(genome)[0]);
		File annofile = new File(genomedir +"annotation/" +annotationUrl +"/"+annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gff")) +".bed.gz");
		System.out.println(genome +" " +genomestring);
	//	parseGFF(annotationFile, annofile, genome, genomestring);		
		
		
		Main.addAnnotationFile(genome, annofile);
		checkGenomes();
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	
	return "";
}
	
static SAMSequenceDictionary ReadDict(File fastafile) {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		
		try {
			BufferedReader reader = null;
			if(fastafile.getName().endsWith(".gz")) {
				reader = new BufferedReader(new FileReader(fastafile.getCanonicalPath().substring(0,fastafile.getCanonicalPath().lastIndexOf(".gz") ) +".fai"));
			}
			else {
				reader = new BufferedReader(new FileReader(fastafile.getCanonicalPath() +".fai"));
			}
			
		
			String line;
			String[] split;
			
			while((line = reader.readLine())!=null) {				
				split = line.split("\t");
				SAMSequenceRecord record = new SAMSequenceRecord(split[0], Integer.parseInt(split[1]));
				dict.addSequence(record);				
			}
			reader.close();
		}
		catch(Exception e){
			
			e.printStackTrace();
		}
 		return dict;
	}
	
	static void parseGFF(File gffFile, File outfile, String genome, String genomeFile) {		
		
			try {
				
				SAMSequenceDictionary dict = ReadDict(new File(genomeFile));				
				
				FileRead.readGFF(gffFile, outfile.getCanonicalPath(), dict);
			}
			catch(Exception e) {
				e.printStackTrace();
			}
			/*	Main.addGenome
				Main.genomehash.get(genome).add(outfile);	
				JMenuItem additem = new JMenuItem(outfile.getName().substring(0,outfile.getName().indexOf(".gff")));
				additem.setName(outfile.getName().substring(0,outfile.getName().indexOf(".gff")));
				additem.addMouseListener(Main.frame);
				Main.addMenu.add(additem);
				
				break;
			*/
			
		
	}
	
	static void indexFasta(File fastafile, String targetDir) {
		try {
			
			ArrayList<String[]> faiList = new ArrayList<String[]>();
			
			BufferedReader reader;
			Boolean zip = false;
			GZIPInputStream gzip = null;
			BufferedWriter fastaWriter = null;
			BufferedWriter indexWriter = null;
		//	BufferedWriter dictWriter;
			if(fastafile.getName().endsWith(".gz")) {
				zip = true;
				gzip = new GZIPInputStream(new FileInputStream(fastafile));
			    reader = new BufferedReader(new InputStreamReader(gzip));
			    if(!new File (targetDir +fastafile.getName().substring(0, fastafile.getName().indexOf(".gz"))).exists()) {
			    	fastaWriter = new BufferedWriter(new FileWriter(targetDir +fastafile.getName().substring(0, fastafile.getName().indexOf(".gz"))));
			    }
			    if(!new File(targetDir+fastafile.getName().substring(0, fastafile.getName().indexOf(".gz")) +".fai").exists()) {
					indexWriter = new BufferedWriter(new FileWriter(targetDir+fastafile.getName().substring(0, fastafile.getName().indexOf(".gz")) +".fai"));
				}
			}
			else {
				reader = new BufferedReader(new FileReader(fastafile));	
				if(!new File(targetDir+"/" +fastafile.getName()).exists()) {
					
					fastaWriter = new BufferedWriter(new FileWriter(targetDir +fastafile.getName()));
				}
				if(!new File(targetDir+"/" +fastafile.getName() +".fai").exists()) {
					indexWriter = new BufferedWriter(new FileWriter(targetDir+fastafile.getName() +".fai"));
				}
				
			}
//			dictWriter = new BufferedWriter(new FileWriter(targetDir+"/" +fastafile.getName().substring(0, fastafile.getName().indexOf(".fa")) +".dict"));
//			dictWriter.write(header);
			String line;
			
			long counter = 0, chromlength = 0, pointer = 0;
			String[] split;
			
			while((line = reader.readLine()) != null) {
				if(fastaWriter != null) {
					fastaWriter.write(line +"\n");
				}
				counter += line.length()+1;
				if(line.startsWith(">")) {
					if(faiList.size() > 0) {
						faiList.get(faiList.size()-1)[1] = ""+chromlength;
			//			dictWriter.write("@SQ\tSN:" +faiList.get(faiList.size()-1)[0] +"\tLN:" +chromlength +"\tUR:file:" +fastafile.getCanonicalPath() +"\tM5:none\n");
					}					
					
					String[] row = new String[5];
					faiList.add(row);
					split = line.split("\\s+");
					pointer = counter;
					row[0] = split[0].substring(1);
					
					row[2] = ""+pointer;						
					line = reader.readLine();
					if(fastaWriter != null) {
						fastaWriter.write(line +"\n");
					}
					chromlength = line.length();
					row[3] = ""+line.length();
					row[4] = ""+(line.length()+1);
					
					counter += line.length()+1;						
					
				}
				else {
					chromlength += line.length();
				}					
			}
			faiList.get(faiList.size()-1)[1] = ""+chromlength;
	//		dictWriter.write("@SQ\tSN:" +faiList.get(faiList.size()-1)[0] +"\tLN:" +chromlength +"\tUR:file:" +fastafile.getCanonicalPath() +"\tM5:none");
			if(fastaWriter != null) {
				fastaWriter.close();
			}
	//		dictWriter.close();
			reader.close();
		
			if(zip) {
				gzip.close();
			}
			if(indexWriter != null) {
				for(int i = 0 ; i<faiList.size(); i++) {
					for(int j = 0; j<4; j++) {
						indexWriter.write(faiList.get(i)[j] +"\t");
					}
					indexWriter.write(faiList.get(i)[4] +"\n");
				}
				indexWriter.close();
			
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private static void createAndShowGUI() {	
	
		if(Main.drawCanvas != null) {
		 	frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE); 
		 	frame.setVisible(false); 
		}
		else {
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
		 	frame.setVisible(true); 
		}
		 	 frame.setResizable(false);    
		frame.setTitle("Genome select");
		    JComponent newContentPane = new AddGenome();
		//    newContentPane.setOpaque(false); 
		   
		    frame.setContentPane(newContentPane);
		    frame.pack();
	   
	}
	public static void main(String[] args) {
		try {
			//UIManager.setLookAndFeel("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel"); 
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		 javax.swing.SwingUtilities.invokeLater(new Runnable() {
	          public void run() {	          
	          		createAndShowGUI();
	          }
	      });	
	}
	public void actionPerformed(ActionEvent event) {
		if(event.getSource() == download) {
			if(!downloading) {
	    		downloading = true;
	    		Main.downloadGenome(genometable.getValueAt(genometable.getSelectedRow(), 0).toString());	    		
	    	}			
		}
		else if(event.getSource() == checkUpdates) {
			//TODO kattoo onko päivityksiä annotaatioon
				String ref = remtable.getValueAt(remtable.getSelectedRow(), 0).toString();
				if(AddGenome.genomeHash.get(ref) != null) {
					ArrayList<String> testfiles = new ArrayList<String>();
					if(Main.drawCanvas != null) {
						for(int i = 0 ; i<Main.genomehash.get(ref).size();i++) {
							testfiles.add(Main.genomehash.get(ref).get(i).getName().replace(".bed.gz",""));
							
						}
					}
					URL testfile = AddGenome.genomeHash.get(ref)[1];
					try {
						String result = Main.checkFile(testfile, testfiles);
						if(result.length() == 0) {
							
							Main.showError("You have newest annotation file.", "Note");
						}
						else {
							int n = JOptionPane.showConfirmDialog(Main.drawCanvas, "New annotation file found: " +result +"\nDownload it now?", "Note", JOptionPane.YES_NO_OPTION);
							if(n == JOptionPane.YES_OPTION) {
								downloadAnnotation(new URL(testfile.getProtocol() +"://" +testfile.getHost() +testfile.getPath().substring(0,testfile.getPath().lastIndexOf("/")+1)+result), ref);
							
							}
							
						}
						
					}
					catch(Exception e) {
						e.printStackTrace();
					}
				}
				
		}
		else if(event.getSource() == remove) {		
		
				String removeref = remtable.getValueAt(remtable.getSelectedRow(), 0).toString();
			
				
					
				try {
					if(Main.drawCanvas != null ) {
						if(removeref.equals(Main.refDropdown.getSelectedItem().toString())) {
							Main.referenceFile.close();
							if(ChromDraw.exonReader != null) {
								ChromDraw.exonReader.close();
							}							
						}
					}
					FileUtils.deleteDirectory(new File(userDir +"/genomes/" +removeref));
					if(Main.drawCanvas != null ) {
						Main.refDropdown.removeItem(removeref);
					}
					remmodel.removeRow(remtable.getSelectedRow());
					
					
					for(int i = 0 ; i<Main.genome.getItemCount(); i++) {
						if(Main.genome.getItem(i).getName() != null) {
							
							if(Main.genome.getItem(i).getName().equals(removeref)) {
								Main.genome.remove(Main.genome.getItem(i));							
								break;
							}
						}
					}
					checkGenomes();
					
				}
				catch(Exception e) {
					Main.showError("Could not delete genome folder.\nYou can do it manually by deleting folder " +userDir +"/genomes/"+removeref, "Note");
				}
		}
		else if(event.getSource() == add) {
		
			if(genomeFile == null) {
				if(new File(genomeFileText.getText()).exists()) {
					genomeFile = new File(genomeFileText.getText());		
					
				}	
				else {
					genomeFileText.setText("Select reference genome fasta-file.");
					genomeFileText.setForeground(Color.red);
					return;
				}
			}
			if(annotationFile == null) {
				if(new File(annotationFileText.getText()).exists()) {
					annotationFile = new File(annotationFileText.getText());
				}
			}
			
			/*if(genomeName.getText().contains("Give name") || genomeName.getText().length() == 0) {
				genomeName.setText("Give name of the genome");
				genomeName.setForeground(Color.red);
				genomeName.revalidate();
				
			}
			else if(!annotation && new File(Main.userDir +"/genomes/"+genomeName.getText().trim().replace("\\s+", "_")).exists()) {
				genomeName.setText("This genome exists already.");
				genomeName.setForeground(Color.red);
				genomeName.revalidate();
			}
			else */
		
			if((genomeFileText.getText().length() == 0 || genomeFileText.getText().startsWith("Select reference"))) {				
				genomeFileText.setText("Select reference genome fasta-file.");
				genomeFileText.setForeground(Color.red);
				genomeFileText.revalidate();
			}
			
			else {	
				
				OutputRunner runner = new OutputRunner(genomeFile.getName(), genomeFile, annotationFile);
				runner.execute();
			}
			
		}	
		else if(event.getSource() == openRef) {
			try {
				
			 JFileChooser chooser = new JFileChooser(Main.downloadDir);	 
	    	  chooser.setMultiSelectionEnabled(false);
	    	  chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
	    	  chooser.setAcceptAllFileFilterUsed(false);	    	
	    	  MyFilterFasta fastaFilter = new MyFilterFasta();	    
	    	 
	    	  chooser.addChoosableFileFilter(fastaFilter); 	    	
	    	  chooser.setDialogTitle("Select reference fasta-file");
	    	  if(Main.screenSize != null) {
	    		  chooser.setPreferredSize(new Dimension((int)Main.screenSize.getWidth()/3, (int)Main.screenSize.getHeight()/3));
	    	  }
	   	  
	          int returnVal = chooser.showOpenDialog((Component)this.getParent());	         
	          
	         if (returnVal == JFileChooser.APPROVE_OPTION) {
	        	 
	        	 genomeFile = chooser.getSelectedFile(); 
	        	 Main.writeToConfig("DownloadDir=" +genomeFile.getParent());
	        	 genomeFileText.setText(genomeFile.getName());
	        	 genomeFileText.revalidate();  	
	        	 frame.pack();
	          }
		 }
		 catch(Exception ex) {
			 ex.printStackTrace();
		 }		
		}
		else if(event.getSource() == openAnno) {
			try {
				 JFileChooser chooser = new JFileChooser(Main.downloadDir);	 
		    	  chooser.setMultiSelectionEnabled(false);
		    	  chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
		    	  chooser.setAcceptAllFileFilterUsed(false);	    	
		    	  MyFilterGFF gffFilter = new MyFilterGFF();	    
		    	 
		    	  chooser.addChoosableFileFilter(gffFilter); 	    	
		    	  chooser.setDialogTitle("Select annotation gff3-file");
		    	  if(Main.screenSize != null) {
		    		  chooser.setPreferredSize(new Dimension((int)Main.screenSize.getWidth()/3, (int)Main.screenSize.getHeight()/3));
		    	  }
		          int returnVal = chooser.showOpenDialog((Component)this.getParent());	         
		          
		         if (returnVal == JFileChooser.APPROVE_OPTION) {
		        	 if(genomeFile == null) {
		        		 genomeFile = Main.fastahash.get(Main.hoverGenome);
		        	 }
		        	 annotationFile = chooser.getSelectedFile(); 
		        	 Main.writeToConfig("DownloadDir=" +annotationFile.getParent());
		        	 annotationFileText.setText(annotationFile.getName());
		        	 annotationFileText.revalidate();  	
		          }
			 }
			 catch(Exception ex) {
				 ex.printStackTrace();
			 }		
		}
	}
	
	static class MyFilterFasta extends javax.swing.filechooser.FileFilter {
		public boolean accept(File file) { 
				if (file.isDirectory()) {
					return true;
			    }	
				if (file.getName().endsWith(".fa")) {
					return true;
				}					
				if (file.getName().endsWith(".fa.gz")) {
					return true;
				}
				if (file.getName().endsWith(".fasta")) {
					return true;
				}
				if (file.getName().endsWith(".fasta.gz")) {
					return true;
				}
				else {
					return false;
				}	
		} 
		
		public String getDescription() { return "*.fa, *.fasta"; }
}
	static class MyFilterGFF extends javax.swing.filechooser.FileFilter {
		public boolean accept(File file) { 
				if (file.isDirectory()) {
					return true;
			    }	
				if (file.getName().endsWith(".gff3")) {
					return true;
				}
				if(file.getName().endsWith(".gff3.gz")) {
					return true;
				}
				else {
					return false;
				}	
		} 
		
		public String getDescription() { return "*.gff3, *.gff3.gz"; }
}

	@Override
	public void mouseClicked(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}


	@Override
	public void mousePressed(MouseEvent e) {
		
		if(e.getSource() == remtable) {
			remove.setEnabled(true);
			genometable.clearSelection();
			download.setEnabled(false);
			checkUpdates.setEnabled(true);
			
		}
		if(e.getSource() == genometable) {
			download.setEnabled(true);
			remtable.clearSelection();
			remove.setEnabled(false);
			checkUpdates.setEnabled(false);
		}
		
	}


	@Override
	public void mouseReleased(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}


	@Override
	public void mouseEntered(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}


	@Override
	public void mouseExited(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
}
