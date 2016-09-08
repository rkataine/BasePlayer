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

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.zip.GZIPInputStream;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.SwingWorker;

import org.apache.commons.io.FilenameUtils;

public class AddGenome  extends JPanel implements ActionListener {

	private static final long serialVersionUID = 1L;
	static JFrame frame = new JFrame("Add new genome");
	static boolean annotation = false;
	static JTextField genomeName = new JTextField("Give name of the genome");
	JTextField genomeFileText = new JTextField("Select reference genome fasta-file");
	JButton openRef = new JButton("Open");
	JTextField annotationFileText = new JTextField("Select annotation gff3-file");
	JButton openAnno = new JButton("Open");
	JButton add = new JButton("Add");
	File genomeFile, annotationFile;
	
	
	public AddGenome() {
		super(new GridLayout(6,2));		
		this.setBackground(Color.black);
		genomeName.setPreferredSize(new Dimension(300,20));
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
			Main.drawCanvas.ready("all");
			Main.drawCanvas.loading("Processing files...");
			createGenome(genomeName, genomefile, annotationFile);
			Main.drawCanvas.ready("Processing files...");
			return "";
		}
		void createGenome(String genomeName, File genomeFile, File annotationFile) {
			try {
				
				
			File annofile = null;
			Main.drawCanvas.splits.get(0).getGenes().clear();
			Boolean ok = false;
			String[] urlsplit = null;
			File fastatest = new File("test");
			if(this.urls != null) {
				urlsplit = urls.split(":");
				URL[] urlList = WelcomeScreen.genomeHash.get(urls);
				URL fastafile= urlList[0];
				String targetDir = Main.userDir +"/genomes/" +urlsplit[0] +"/";
				File fasta = new File(targetDir +FilenameUtils.getName(fastafile.getFile()));
				fastatest = new File(targetDir +FilenameUtils.getName(fastafile.getFile()).substring(0,FilenameUtils.getName(fastafile.getFile()).indexOf(".gz") ));
				new File(targetDir).mkdir();
				if(!fasta.exists() && !fastatest.exists()) {
					Main.drawCanvas.loadingtext ="Downloading " +genomeFile.getName();
					Main.downloadFile(fastafile, targetDir, WelcomeScreen.sizeHash.get(urls)[0]);				
				}
				
				URL gfffile= urlList[1];
				targetDir = Main.userDir +"/genomes/" +urlsplit[0] +"/annotation/" +urlsplit[1] +"/";
				new File(targetDir).mkdirs();
				File gff = new File(targetDir +FilenameUtils.getName(gfffile.getFile()));				
				
				Main.drawCanvas.loadingtext ="Downloading " +gff.getName();
				Main.downloadFile(gfffile, targetDir, WelcomeScreen.sizeHash.get(urls)[1]);				
				
				
				if(urlList.length == 3) {
					URL bandFile = urlList[2];
					targetDir = Main.userDir +"/genomes/" +urlsplit[0] +"/";
					Main.downloadFile(bandFile, targetDir, 0);
					File infile = new File(targetDir+FilenameUtils.getName(bandFile.getFile()));
					File outfile = new File(targetDir+"bands.txt");
					MethodLibrary.unzip(infile, outfile);
				}
				
			}
			if(urls == null && !annotation) {
				ok = new File(Main.userDir +"/genomes/" +genomeName +"/").mkdir();
			}
			else {
				ok = true;
			}
			
			if(ok) {
				String genomedir = Main.userDir +"/genomes/" +genomeName +"/";
				if(!annotation && !fastatest.exists()) {
					Main.drawCanvas.loadingtext = "Unpacking and indexing: " +genomeFile.getName();
					
					indexFasta(genomeFile, genomedir);
					Main.addGenomeFile(genomeName);
				}
			
				if(annotationFile != null) {
					
					Main.drawCanvas.loadingtext = "Unpacking and indexing: " +annotationFile.getName();
					if(urls == null) {
						ok = new File(genomedir +"annotation/" +annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gff"))+"/").mkdirs();
					}
					else {
						ok = true;
					}
					
					if(urls != null) {
						annofile = new File(genomedir +"annotation/" +urlsplit[1] +"/"+annotationFile.getName().substring(0,annotationFile.getName().indexOf(".gff")) +".bed.gz");
						
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
					Main.defaultGenome = genomeName;
					Main.setChromDrop(genomeName);
					Main.getBands();
					Main.defaultAnnotation = annofile.getName();	
				    Main.getExons();
				    Main.chromosomeDropdown.setSelectedIndex(0);
				}
			}
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		
	}
	
static SAMSequenceDictionary ReadDict(File fastafile) {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		
		try {
			BufferedReader reader = null;
			if(fastafile.getName().endsWith(".gz")) {
				reader = new BufferedReader(new FileReader(fastafile.getCanonicalPath().substring(0,fastafile.getCanonicalPath().indexOf(".gz") ) +".fai"));
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
				System.out.println(dict);
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
		//	ArrayList<String> dictList = new ArrayList<String>();
		//	String header = "@HD\tVN:1.0\tSO:unsorted\n";
		//	dictList.add(header);
			
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
	
		
		 	frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE); 
			 frame.setVisible(true); 
		
		 	 frame.setResizable(true);    
		
		    JComponent newContentPane = new AddGenome();
		    newContentPane.setOpaque(false); 
		   
		    frame.setContentPane(newContentPane);
		    frame.pack();
	   
	}
	public static void main(String[] args) {
		
		 javax.swing.SwingUtilities.invokeLater(new Runnable() {
	          public void run() {	          
	          		createAndShowGUI();
	          }
	      });	
	}
	public void actionPerformed(ActionEvent event) {
		
		if(event.getSource() == add) {
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
			
			if(genomeName.getText().contains("Give name") || genomeName.getText().length() == 0) {
				genomeName.setText("Give name of the genome");
				genomeName.setForeground(Color.red);
				genomeName.revalidate();
				
			}
			else if(!annotation && new File(Main.userDir +"/genomes/"+genomeName.getText().trim().replace("\\s+", "_")).exists()) {
				genomeName.setText("This genome exists already. Give new name.");
				genomeName.setForeground(Color.red);
				genomeName.revalidate();
			}
			else if(!annotation && (genomeFileText.getText().length() == 0 || genomeFileText.getText().startsWith("Select reference"))) {
				
				genomeFileText.setText("Select reference genome fasta-file.");
				genomeFileText.setForeground(Color.red);
				genomeFileText.revalidate();
			}
			
			else {	
				OutputRunner runner = new OutputRunner(genomeName.getText().trim().replace("\\s+", "_"), genomeFile, annotationFile);
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
	    	  chooser.setPreferredSize(new Dimension((int)Main.screenSize.getWidth()/3, (int)Main.screenSize.getHeight()/3));
	          int returnVal = chooser.showOpenDialog((Component)this.getParent());	         
	          
	         if (returnVal == JFileChooser.APPROVE_OPTION) {
	        	 
	        	 genomeFile = chooser.getSelectedFile(); 
	        	 Main.writeToConfig("DownloadDir=" +genomeFile.getParent());
	        	 genomeFileText.setText(genomeFile.getName());
	        	 genomeFileText.revalidate();  	
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
		    	  chooser.setPreferredSize(new Dimension((int)Main.screenSize.getWidth()/3, (int)Main.screenSize.getHeight()/3));
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
}
