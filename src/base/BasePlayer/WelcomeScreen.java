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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridLayout;
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.net.URL;
import java.util.HashMap;

import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;
import javax.swing.text.html.HTMLEditorKit;
import javax.swing.text.html.StyleSheet;

public class WelcomeScreen  extends JPanel implements ActionListener {
	static boolean downloading = false;
	private static final long serialVersionUID = 1L;
	static JFrame frame = new JFrame("Welcome screen");	
	JEditorPane htmlPage;
	boolean editorial = false;
	static Image image;
	static HashMap<String, URL[]> genomeHash = new HashMap<String, URL[]>();
	static HashMap<String, Integer[]> sizeHash = new HashMap<String, Integer[]>();
	static JCheckBox dontShow = new JCheckBox("Don't show this message again.");
	public WelcomeScreen() {
		super(new BorderLayout());		
		try {
		makeGenomes();
		this.setBackground(Color.black);
				
		StringBuffer html = new StringBuffer("<html><body style='color:black;background:rgb(228,228,218);'><h1>Welcome to BasePlayer</h1>");
				if(editorial) {
					html.append("<p>This is editors/reviewers version of BasePlayer for testing purposes.<br>"
							+ "Check http://baseplayer.fi for instructions<br>"
							+ "This package includes following items:"
							+ "<ul><li> human reference sequence and Ensembl gene annotation for chromosome 20.</li>"
							+ "<li>6 whole-genome samples from 1000 Genomes Project, including variants (vcf.gz) and read sequences (cram). File->Open samples</li>"
							+ "<li>Variant data from ExAC for allele frequency filtering. File->Add controls</li>"
							+ "<li>Additional tracks for ENCODE regulatory regions, transcription factor binding sites from Ensembl Biomart and replication timing. File->Add tracks</li>"
							+ "</ul>"
							+ "You can download whole reference genome and full annotation from links below.");					
				}
		
		
				if(Main.genomehash == null || Main.genomehash.size() == 0) {
					html.append("<p>Before we start, download your favorite genome using following links or add new genome by hand in File->Change/add genome->Add new genome</p><br>");					
				}
				else {
					html.append("<p>Download your favorite genome using following links or add new genome by hand in File->Change/add genome->Add new genome</p><br>");			
				}
				
					html.append("<a href=http:Homo_sapiens_GRCh37:Ensembl_genes> Homo sapiens GRCh37 with Ensembl</a> or <a href=http:Homo_sapiens_GRCh37:RefSeq_genes>RefSeq</a> gene annotations<br>");
					html.append("<a href=http:Homo_sapiens_GRCh38:Ensembl_genes> Homo sapiens GRCh38 with Ensembl</a> or <a href=http:Homo_sapiens_GRCh38:RefSeq_genes>RefSeq</a> gene annotations<br><br>");
					
					html.append("<a href=http:Mus_musculus_GRCm38:Ensembl_genes> Mus musculus GRCm38 with Ensembl</a> or <a href=http:Mus_musculus_GRCm38:RefSeq_genes>RefSeq</a> gene annotations<br>");
					html.append("<a href=http:Rattus_norvegicus:Ensembl_genes> Rattus norvegicus with Ensembl gene annotations</a><br>");
					html.append("<a href=http:Saccharomyces_cerevisiae:Ensembl_genes> Saccharomyces cerevisiae with Ensembl gene annotation</a><br>");					
					html.append("<a href=http:Ciona_intestinalis:Ensembl_genes> Ciona intestinalis with Ensembl gene annotation</a><br>");
				
					html.append("</body> </html>");
					
			htmlPage = new JEditorPane();
			
			//htmlPage.setBackground(Draw.sidecolor);
			
		//	htmlPage.setOpaque(false);	
			htmlPage.setEditable(false);
			htmlPage.setEditorKit(JEditorPane.createEditorKitForContentType("text/html"));
			htmlPage.setText(html.toString());
			
			htmlPage.setPreferredSize(new Dimension(400, 500));
			htmlPage.addHyperlinkListener(new HyperlinkListener() {
			public void hyperlinkUpdate(HyperlinkEvent hyperlinkEvent) {
			    HyperlinkEvent.EventType type = hyperlinkEvent.getEventType();
			    final URL url = hyperlinkEvent.getURL();
			    if (type == HyperlinkEvent.EventType.ACTIVATED) {
			    	if(!downloading) {
			    		downloading = true;
			    		Main.downloadGenome(url.toString().substring(5));
			    		downloading = false;
			    	}
			    }
			  }
		});
		htmlPage.setCaretPosition(0);
		JScrollPane scrollpane = new JScrollPane(htmlPage);
		
		scrollpane.setPreferredSize(new Dimension(400, 600));
		
		add(scrollpane);
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	void makeGenomes() {
		try {
			URL[] urls = {new URL("ftp://ftp.ensembl.org/pub/grch37/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"),
						  new URL("ftp://ftp.ensembl.org/pub/grch37/update/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz"),
						  new URL("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz")};
			Integer[] sizes = new Integer[2];
			sizes[0] = 900;
			sizes[1] = 40;
			sizeHash.put("Homo_sapiens_GRCh37:Ensembl_genes", sizes);
			genomeHash.put("Homo_sapiens_GRCh37:Ensembl_genes", urls);
			
			URL[] urls2 = {new URL("ftp://ftp.ensembl.org/pub/grch37/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"),
					  new URL("ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.37.3/GFF/ref_GRCh37.p5_top_level.gff3.gz"),
					  new URL("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz") };
			sizes = new Integer[2];
			sizes[0] = 900;
			sizes[1] = 35;
			sizeHash.put("Homo_sapiens_GRCh37:RefSeq_genes", sizes);
			genomeHash.put("Homo_sapiens_GRCh37:RefSeq_genes", urls2);
			
			URL[] urls3 = { new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/ciona_intestinalis/dna/Ciona_intestinalis.KH.dna.toplevel.fa.gz"),
						new URL("ftp://ftp.ensembl.org/pub/release-85/gff3/ciona_intestinalis/Ciona_intestinalis.KH.85.gff3.gz")
			};
			sizes = new Integer[2];
			sizes[0] = 34;
			sizes[1] = 5;
			sizeHash.put("Ciona_intestinalis:Ensembl_genes", sizes);
			genomeHash.put("Ciona_intestinalis:Ensembl_genes", urls3);
			
			URL[] urls4 = { new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"),
					new URL("ftp://ftp.ensembl.org/pub/release-85/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.85.gff3.gz")
			};
			sizes = new Integer[2];
			sizes[0] = 4;
			sizes[1] = 1;
			sizeHash.put("Saccharomyces_cerevisiae:Ensembl_genes", sizes);
			genomeHash.put("Saccharomyces_cerevisiae:Ensembl_genes", urls4);
			
			URL[] urls5 = {new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz"),
					  new URL("ftp://ftp.ensembl.org/pub/release-85/gff3/mus_musculus/Mus_musculus.GRCm38.85.gff3.gz"),
					  new URL("http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/cytoBand.txt.gz")};
			
			sizes = new Integer[2];
			sizes[0] = 801;
			sizes[1] = 58;
			sizeHash.put("Mus_musculus_GRCm38:Ensembl_genes", sizes);
			genomeHash.put("Mus_musculus_GRCm38:Ensembl_genes", urls5);
			
			URL[] urls6 = {new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz"),
					  new URL("ftp://ftp.ncbi.nlm.nih.gov/genomes/M_musculus/GFF/ref_GRCm38.p4_top_level.gff3.gz") };
			sizes = new Integer[2];
			sizes[0] = 801;
			sizes[1] = 28;
			sizeHash.put("Mus_musculus_GRCm38:RefSeq_genes", sizes);
			genomeHash.put("Mus_musculus_GRCm38:RefSeq_genes", urls6);
			
			URL[] urls7 = {new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz"),
					  new URL("ftp://ftp.ensembl.org/pub/release-85/gff3/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.85.gff3.gz") };
			sizes = new Integer[2];
			sizes[0] = 812;
			sizes[1] = 15;
			sizeHash.put("Rattus_norvegicus:Ensembl_genes", sizes);
			genomeHash.put("Rattus_norvegicus:Ensembl_genes", urls7);
			
			
			URL[] urls8 = {new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
					  new URL("ftp://ftp.ensembl.org/pub/release-85/gff3/homo_sapiens/Homo_sapiens.GRCh38.85.gff3.gz"),
					  new URL("http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz")};
			sizes = new Integer[2];
			sizes[0] = 860;
			sizes[1] = 38;
			sizeHash.put("Homo_sapiens_GRCh38:Ensembl_genes", sizes);
			genomeHash.put("Homo_sapiens_GRCh38:Ensembl_genes", urls8);
			URL[] urls9 = {new URL("ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
					  new URL("ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/GFF/ref_GRCh38.p7_top_level.gff3.gz"),
					  new URL("http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz")};
			sizes = new Integer[2];
			sizes[0] = 860;
			sizes[1] = 41;
			sizeHash.put("Homo_sapiens_GRCh38:RefSeq_genes", sizes);
			genomeHash.put("Homo_sapiens_GRCh38:RefSeq_genes", urls9);
		}
		catch(Exception ex) {
			ex.printStackTrace();
		}
	}
	private static void createAndShowGUI() {	
	
			if(Main.frame != null) {
				frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE); 
			}
			else {
				frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
			}
			frame.setAlwaysOnTop(true);
		 	frame.setVisible(true); 
		
		 	frame.setResizable(true);    
		
		    JComponent newContentPane = new WelcomeScreen();
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
	@Override
	public void actionPerformed(ActionEvent e) {
		if(e.getSource() == dontShow) {

			if(Launcher.firstStart) {
				Main.writeToConfig("FirstStart=false");
			}
		}
		
	}


}
