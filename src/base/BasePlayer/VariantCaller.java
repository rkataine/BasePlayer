package base.BasePlayer;

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

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.SwingWorker;

import base.BasePlayer.BedCanvas.Annotator;

public class VariantCaller extends JPanel implements ActionListener {
	
	private static final long serialVersionUID = 1L;
	static Boolean loading = false;
	static JFrame frame = new JFrame("Variant Caller");
	JLabel testing = new JLabel("Test version");
	static JTextField minreads = new JTextField("4");
	JLabel minreadlabel = new JLabel("Min. alt-read count");
	static JTextField mincoverage = new JTextField("10");
	JLabel mincoveragelabel = new JLabel("Min. coverage");
	static JTextField minbaseq = new JTextField("5");
	JLabel minbaseqlabel = new JLabel("Min. base quality");
	static JTextField minmappingq = new JTextField("10");
	JLabel minmappingqlabel = new JLabel("Min. mapping quality");
	static JCheckBox bothstrand = new JCheckBox("Require both strands");
	static JCheckBox bothruns = new JCheckBox("Require multiple runs");
	
	JCheckBox onlySel = new JCheckBox("Calc only for selected");
	static JCheckBox inanno = new JCheckBox("Before annotation");
	JButton execute = new JButton("Execute");
	public VariantCaller(boolean caller) {
		
	}
	public VariantCaller() {
		super(new GridLayout(8,2));		
		try {
			testing.setForeground(Color.red);
			execute.addActionListener(this);
			inanno.addActionListener(this);
			add(testing);
			add(new JLabel());
			add(minreads);
			add(minreadlabel);
			add(mincoverage);
			add(mincoveragelabel);
			add(minmappingq);
			add(minmappingqlabel);
			add(minbaseq);
			add(minbaseqlabel);
			add(bothstrand);
			add(bothruns);
			add(onlySel);
			add(inanno);
			//add(new JSeparator());
			add(execute);			
			bothstrand.setToolTipText("Variant call is discared if mismatches are present only in the other strand.");
			bothruns.setToolTipText("If bam file is merged from multiple runs, variant call must be present in reads from at least two separate runs.");
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private static void createAndShowGUI() {	
	
			if(Main.frame != null) {
				frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE); 
				frame.setVisible(false);		
			}
			else {
				frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
				frame.setVisible(true);		
			}
			frame.setAlwaysOnTop(true);
		 	
		    JComponent newContentPane = new VariantCaller();
		    newContentPane.setOpaque(false);		   
		    frame.setContentPane(newContentPane);
		    setFonts(Main.menuFont);
		    frame.pack();
		   
		    for(int i =0; i<frame.getContentPane().getComponentCount(); i++) {
		    	frame.getContentPane().getComponent(i).setMinimumSize(frame.getContentPane().getComponent(i).getPreferredSize());
		    }
		    minreads.setMinimumSize(new Dimension(minreads.getWidth(),minreads.getHeight()));	
		    mincoverage.setMinimumSize(new Dimension(mincoverage.getWidth(),mincoverage.getHeight()));
			minbaseq.setMinimumSize(new Dimension(minbaseq.getWidth(),minbaseq.getHeight()));
			minmappingq.setMinimumSize(new Dimension(minmappingq.getWidth(),minmappingq.getHeight()));
			
	}
	public static void main(String[] args) {
		
		 javax.swing.SwingUtilities.invokeLater(new Runnable() {
	          public void run() {	          
	          		createAndShowGUI();
	          }
	      });	
	}
	
	static void setFonts(Font menuFont) {
		for(int i = 0 ; i<VariantCaller.frame.getContentPane().getComponentCount(); i++) {
			VariantCaller.frame.getContentPane().getComponent(i).setFont(menuFont);
		}
		VariantCaller.frame.pack();
	}
	@Override
	public void actionPerformed(ActionEvent event) {		
		if(event.getSource() == execute) {
			try {				
				VarCaller caller = new VarCaller();
				caller.execute();												
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		else if(event.getSource() == inanno) {
			if(inanno.isSelected()) {				
				FileRead.caller = true;
				FileRead.checkSamples();
			}
			else {
				FileRead.caller = false;
				FileRead.checkSamples();
			}
		}
	}

public class VarCaller extends SwingWorker<String, Object> {
		int readlevel,minquality, minreadquality, mincoverage; 
		boolean onlysel;
		public VarCaller(int readlevel, int minquality, int minreadquality, int mincoverage, boolean onlysel) {
			this.readlevel = readlevel;
			this.minquality = minquality;
			this.minreadquality = minreadquality;
			this.mincoverage = mincoverage;
			this.onlysel = onlysel;
		}
		public VarCaller() {
			int readlevel,minquality, minreadquality, mincoverage; 
			try {
				readlevel = Integer.parseInt(VariantCaller.minreads.getText());
			}
			catch(Exception e) {
				VariantCaller.minreads.setForeground(Color.red);
				return;
			}
			try {
				minquality = Integer.parseInt(VariantCaller.minbaseq.getText());
			}
			catch(Exception e) {
				VariantCaller.minbaseq.setForeground(Color.red);
				return;
			}
			try {
				minreadquality = Integer.parseInt(VariantCaller.minmappingq.getText());
			}
			catch(Exception e) {
				VariantCaller.minmappingq.setForeground(Color.red);
				return;
			}
			try {
				mincoverage = Integer.parseInt(VariantCaller.mincoverage.getText());
			}
			catch(Exception e) {
				VariantCaller.mincoverage.setForeground(Color.red);
				return;
			}
			this.readlevel = readlevel;
			this.minquality = minquality;
			this.minreadquality = minreadquality;
			this.mincoverage = mincoverage;
			this.onlysel = onlySel.isSelected();
		}
	
		void callVariants() {

			if(onlysel && (Main.drawCanvas.selectedSample == null || Main.drawCanvas.selectedSample.getTabixFile() != null || Main.drawCanvas.selectedSample.samFile == null)) {
				VariantCaller.loading = false;
				return;			
			}
		
			int startpos=(int)Main.drawCanvas.splits.get(0).start, endpos=(int)Main.drawCanvas.splits.get(0).end;
			FileRead reader = new FileRead();
			FileRead.head.putNext(null);
			VariantCaller.loading = true;
			
			if(!Main.drawCanvas.loading) {
				Main.drawCanvas.loading("Calling variants");
			}
			else {
				Main.drawCanvas.loadingtext = "Calling variants";
			}
			
			
			Main.bedCanvas.bedOn = false;
				for(int s = 0; s<Main.samples;s++) {
					
					try {	
						if(!Main.drawCanvas.loading) {
							break;
						}
						Main.drawCanvas.loadBarSample = 0;  
						Main.drawCanvas.loadbarAll  = s/Main.samples*100;
						Sample sample = Main.drawCanvas.sampleList.get(s);
						Reads readClass = sample.getreadHash().get(Main.drawCanvas.splits.get(0));
						sample.basequalsum = 0;
						sample.basequals = 0;
						if(onlysel && !sample.equals(Main.drawCanvas.selectedSample)) {
							continue;
						}
						if(sample.getTabixFile() != null || sample.samFile == null) {
							continue;
						}
						
						VarNode node = FileRead.head.getNext();						
						
						while(node != null) {							
							node.removeSample(sample);				
							node = node.getNext();					
						}
						
						node = null;
						Reads reads = null;
						try {
							reads = (Reads)sample.getreadHash().get(Main.drawCanvas.splits.get(0)).clone();
							
						}
						catch(Exception e) {
							e.printStackTrace();
						}
						
						if(!sample.calledvariants) {
							Main.varsamples++;
							VariantHandler.commonSlider.setMaximum(Main.varsamples);
							VariantHandler.commonSlider.setUpperValue(Main.varsamples);
							VariantHandler.geneSlider.setMaximum(Main.varsamples);
						}
						sample.maxCoveragecaller = 0;
						sample.calledvariants = true;
						String[] coveragebases = { "", "A", "C", "G", "T", "N"};
						double homlevel = 0.95;						
						
						
						int interval;
						
						StringBuffer genotypes = new StringBuffer("");
						int calls = 0, refs = 0;
									
						String refbase = "";
						
						//boolean genotype = false;
						reader.current = FileRead.head;
						if(endpos-startpos < Settings.windowSize) {
							interval = endpos-startpos;
						}
						else {
							interval = Settings.windowSize;
						}
						Draw.updatevars = true;
						Main.drawCanvas.repaint();
						
						String[] line = new String[10];
						line[2] = ".";
						line[5] = "99";
						line[6] = "PASS";
						line[7] = "INFO";
						line[8] = "GT:AD:DP";
					//	BedNode bednode = null;
						int start, end;
					/*	if(Main.bedCanvas.bedTrack.size() > 0) {
							if(Main.bedCanvas.bedTrack.get(0).getHead().getNext() == null) {
								break;
							}
							else {
								bednode = Main.bedCanvas.bedTrack.get(0).getHead().getNext();
								start =	bednode.getPosition();
								end =  bednode.getPosition() + bednode.getLength()+1;
							//	bednode = bednode.getNext();
							//	while(bednode.getPosition() <= end) {
							//		end = bednode.getPosition() + bednode.getLength()+1;
							//		bednode = bednode.getNext();
							//	}
								
							}
						}
						else {
						*/
							start = startpos;
							end = startpos + interval;
						//}
						
						Main.drawCanvas.variantsStart = (int)start;
						while(start < endpos) {
							
							if(!Main.drawCanvas.loading) {
								Main.drawCanvas.ready("all");
								
								break;
							}
							
							if(end > endpos) {
								end = endpos;
							}
							Main.drawCanvas.variantsEnd = end;
							ReferenceSeq reference = new ReferenceSeq(Main.drawCanvas.splits.get(0).chrom, start-300, end+300, Main.referenceFile);
							VariantCall[][] coverages = reader.variantCaller(Main.drawCanvas.splits.get(0).chrom, start,end, reads, minquality,minreadquality, reference);
							
							for(int i = 0 ; i<coverages.length; i++) {		
								if(!Main.drawCanvas.loading) {
									break;
								}
								if(coverages[i][0] == null) {
									continue;
								}
								if(coverages[i][0].calls > sample.maxCoveragecaller) {
									sample.maxCoveragecaller = coverages[i][0].calls;
								}
								
								for(int j = 1; j< coverages[i].length-2; j++) {
									if(coverages[i][j] == null) {
										continue;
									}
									
									if(coverages[i][j].calls >= readlevel && coverages[i][0].calls >= mincoverage /*&& coverages[i][j].runs.size() > 1 && coverages[i][j].strands.size() > 1*/) {
										if(VariantCaller.bothstrand.isSelected() && coverages[i][j].strands.size() < 2) {
											continue;
										}
										if(VariantCaller.bothruns.isSelected() && readClass.getRuns() > 1 && coverages[i][j].runs.size() < 2) {
											continue;
										}
										genotypes = new StringBuffer("");
										calls = coverages[i][j].calls;
										refs = coverages[i][0].calls;
										if(calls/(double)refs >= homlevel) {
											genotypes.append("1/1:");
										//	genotype = true;
										}
										else {
											genotypes.append("0/1:");
										//	genotype = false;
										}
										
										genotypes.append((refs-calls) +"," +calls +":" +refs);
										refbase = Main.getBase.get(reference.getSeq()[(reads.getCoverageStart() + i -1) -reference.getStartPos()]);								
										line[0] = Main.drawCanvas.splits.get(0).chrom;										
										line[1] = ""+(reads.getCoverageStart() + i);							
										line[3] = refbase;										
										line[4] = coveragebases[j];
										line[9] = genotypes.toString();										
										reader.readLine(line, sample);
										
									}									
								}
							}
							/*if(bednode != null) {
								bednode = bednode.getNext();
								
								if(bednode != null) {
									start = bednode.getPosition();
									end = bednode.getPosition()+bednode.getLength()+1;
								}
								else {
									break;
								}
							}
							else {*/
								start = end;						
								end += interval;
						//	}
							Draw.updatevars = true;
							Main.drawCanvas.repaint();
							
						}
						
						node = null;
						reads = null;
						readClass = null;
						System.out.println(sample.getName() +": " +(sample.basequalsum/(double)sample.basequals));
					}
					catch(Exception e) {
						e.printStackTrace();
						break;
					}					
					
				}			
			
			Main.drawCanvas.variantsEnd = endpos;
			reader.current = FileRead.head;
			reader = null;
			FileRead.annotate();	
			if(Control.controlData.controlsOn) {		    	
			    Control.applyControl();
			}						
			
			Main.drawCanvas.calcClusters(FileRead.head,1);
				
			Main.bedCanvas.bedOn = true;
			boolean ison = false;
			for(int i = 0 ; i<Main.bedCanvas.bedTrack.size();i++) {
				if(Main.bedCanvas.bedTrack.get(i).intersect) {
				
					ison = true;
					break;
				}
			}
			if(!ison) {
				Main.bedCanvas.bedOn = false;
			}
		
			if(Main.bedCanvas.bedOn) {
			
				Main.drawCanvas.loadingtext = "Annotating variants";
				for(int i = 0 ; i<Main.bedCanvas.bedTrack.size(); i++) {
					if(Main.bedCanvas.bedTrack.get(i).small && Main.bedCanvas.bedTrack.get(i).getBBfileReader() == null) {						
						if( Main.bedCanvas.bedTrack.get(i).intersect && !Main.bedCanvas.bedTrack.get(i).loading) {
							Main.bedCanvas.annotate(Main.bedCanvas.bedTrack.get(i).getHead(), FileRead.head.getNext());	
							Main.bedCanvas.intersected = true;
							
						}							
						else if(Main.bedCanvas.bedTrack.get(i).intersect && Main.bedCanvas.bedTrack.get(i).loading) {
							Main.bedCanvas.bedTrack.get(i).waiting = true;							
						}
					}
					else if(Main.bedCanvas.bedTrack.get(i).intersect) {		
						
						BedCanvas.Annotator annotator = Main.bedCanvas.new Annotator(Main.bedCanvas.bedTrack.get(i));
						annotator.annotateVars();
						
						Main.bedCanvas.intersected = true;
						
					}
				}	
				
			}
			
			Main.drawCanvas.ready("Calling variants");			
					
			Draw.updatevars = true;
			Main.drawCanvas.repaint();
			VariantCaller.loading = false;
			
		}
		@Override
		protected String doInBackground() throws Exception {
			callVariants();
			return null;
		}
	}
}
