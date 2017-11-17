package base.BasePlayer;

import java.awt.Color;
import java.awt.Graphics;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import javax.swing.JPanel;

public class Plotter extends JPanel {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	int[][] array = null;
	int width, barwidth;
	Color lightgray = new Color(220,220,220);
	int samplewidth = 0, scale = 0, valuecount = 0, chromsize;
	int bottom = 310, highest = 0, windowsize = 3000000;
	int[] simplearray;
	double[] chromNormalize = {
			1.02702702702703,
			1.08108108108108,
			1.10810810810811,
			1.13513513513514,
			1.10810810810811,
			1.24324324324324,
			0.972972972972973,
			0.972972972972973,
			0.891891891891892,
			0.945945945945946,
			0.945945945945946,
			1,
			1.02702702702703,
			1.02702702702703,
			0.972972972972973,
			0.783783783783784,
			0.837837837837838,
			1.08108108108108,
			0.810810810810811,
			0.837837837837838,
			0.891891891891892,
			0.783783783783784
	};	
	
	public Plotter(int[][] array, int width) {
		this.width = width;
		samplewidth = this.width; ///Main.samples;
		this.valuecount = array[0].length;
		barwidth = samplewidth/valuecount;
		this.array = array;
		int sum = 0, countsum = 0;
		for(int i = 0; i<array.length; i++) {			
			sum = 0;
			countsum = 0;
			for(int j = 0; j<array[0].length; j++) {
				if(array[i][j] > highest) {
					highest = array[i][j];					
				}
				sum+=array[i][j]*(j+1);
				countsum += (j+1);
				System.out.print(array[i][j] +" ");
			}
			System.out.println();
			System.out.println(sum/(double)countsum);
		}		
	}
	public Plotter(int width) {
		try {
			boolean merit = true;
			String control = "";
			if(VariantHandler.synonymous.isSelected()) {
				control = "_control";
			}
			
			String path = "/mnt/cg8/Riku/Ohutsuoli/bafseg" +control +"/";
			if(!merit) {
				path = "X:/cg8/Riku/Ohutsuoli/bafseg" +control +"/";
			}
			//BufferedReader coverages = new BufferedReader(new FileReader(path +"/cg8/Riku/Ohutsuoli/MSS_coverages.txt"));
			BufferedReader coverages = new BufferedReader(new FileReader(path +"coverages.txt"));
			String line;
			String[] split;
			HashMap<String, Integer> covs = new HashMap<String, Integer>();
						
			while((line = coverages.readLine()) != null) {
				split = line.split("\\s+");
				covs.put(split[0], Integer.parseInt(split[3]));
			}
			coverages.close();
			this.width = width;			
			chromsize = Main.drawCanvas.splits.get(0).chromEnd;
			simplearray = new int[chromsize/windowsize];
			barwidth =width/simplearray.length;			
			Map.Entry<String, ArrayList<SampleNode>> entry;
			HashMap<String, BufferedWriter> samplewriter = new HashMap<String, BufferedWriter>();
			
			//BufferedWriter namewriter = new BufferedWriter(new FileWriter(path +"extracted/sample_names.txt"));
			//namewriter.write("Assay\tFilename\tIGV_index\n");
			for(int i = 0 ; i<Main.drawCanvas.sampleList.size(); i++) {
				if(Main.drawCanvas.sampleList.get(i).multiVCF || Main.drawCanvas.sampleList.get(i).getTabixFile() == null) {
					continue;
				}
				BufferedWriter writer = new BufferedWriter(new FileWriter(path +"extracted/"+Main.drawCanvas.sampleList.get(i).getName() +".tsv"));
				writer.write("Name\tChr\tPosition\tGenotype\tB Allele Frequency\tLog R Ratio\tamplified\tdepleted\n");
				samplewriter.put(Main.drawCanvas.sampleList.get(i).getName(), writer);
			//	namewriter.write(Main.drawCanvas.sampleList.get(i).getName() +"\t" +Main.drawCanvas.sampleList.get(i).getName() +".tsv\t" +i +"\n");
			}
		//	namewriter.close();
			for(int i = 0; i<simplearray.length; i++) {
				simplearray[i] = 0;
			}
			double logr = 0;
			for(int i = 1 ; i<24; i++) {
				System.out.println(i);
				VarNode node = FileRead.head.getNextVisible(FileRead.head);
				while(node != null) {
					if(node.isRscode() == null || node.indel) {
						node = node.getNextVisible(node);
						continue;
					}
					for(int v = 0; v<node.vars.size(); v++) {
						entry = node.vars.get(v);
						if(Main.drawCanvas.hideNodeVar(node, entry)) {
							continue;
						}
						for(int m = 0; m<entry.getValue().size(); m++) {
							if(Main.drawCanvas.hideVar(entry.getValue().get(m), node.indel)) {
								continue;
							}
							if(entry.getValue().get(m).alleles != null) {							
								break;
							}
							
							if(entry.getKey().length() > 1) {
								continue;
							}
							logr = Math.log(entry.getValue().get(m).getCoverage()/(double)(covs.get(entry.getValue().get(m).getSample().getName())*chromNormalize[i-1]))/(double)Math.log(2);
							
							if(entry.getValue().get(m).isHomozygous()) {
							//	System.out.println(entry.getValue().get(m).getSample().getName() +"\t" +node.getChrom() +"\t" +(node.getPosition()+1) +"\t" +entry.getKey() +entry.getKey() +"\t" +entry.getValue().get(m).getAlleleFraction() +"\t0\t0\t0");
								samplewriter.get(entry.getValue().get(m).getSample().getName()).write(node.isRscode()+"\t" +node.getChrom() +"\t" +(node.getPosition()+1) +"\t" +entry.getKey() +entry.getKey() +"\t" +entry.getValue().get(m).getAlleleFraction() +"\t" +logr +"\t" +Main.getBase.get(node.getRefBase()) +"\t" +entry.getKey() +"\n");
							}
							else {
							//	System.out.println(entry.getValue().get(m).getSample().getName() +"\t" +node.getChrom() +"\t" +(node.getPosition()+1) +"\t" +Main.getBase.get(node.getRefBase()) +entry.getKey() +"\t" +entry.getValue().get(m).getAlleleFraction()+"\t0\t0\t0");
								samplewriter.get(entry.getValue().get(m).getSample().getName()).write(node.isRscode() +"\t" +node.getChrom() +"\t" +(node.getPosition()+1) +"\t" +Main.getBase.get(node.getRefBase()) +entry.getKey() +"\t" +entry.getValue().get(m).getAlleleFraction()+"\t" +logr +"\t" +Main.getBase.get(node.getRefBase()) +"\t" +entry.getKey() +"\n");
							}
							
						}
					}
					
					node = node.getNextVisible(node);
				}
				Main.nothread = true;
				Main.chromosomeDropdown.setSelectedIndex(i);
			}
			for(int i = 0 ; i<Main.drawCanvas.sampleList.size(); i++) {
				if(Main.drawCanvas.sampleList.get(i).multiVCF || Main.drawCanvas.sampleList.get(i).getTabixFile() == null) {
					continue;
				}
				samplewriter.get(Main.drawCanvas.sampleList.get(i).getName()).close();
			}
		
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	public void paint(Graphics g) {
		try {					
			
			if(array != null) {
				for(int i = 0; i<array.length; i++) {
					if(i%2 == 0) {
						g.setColor(lightgray);
						g.fillRect(samplewidth*i, 10, samplewidth, bottom);
					}
				}
				for(int i = 0;i<11;i++) {
					g.setColor(Color.LIGHT_GRAY);
					g.fillRect(0,bottom-i*30,width+30,1);
					g.setColor(Color.black);
					g.fillRect(width+20,bottom-i*30,10,1);					
					g.drawString(""+MethodLibrary.map(i*30, 0,300,0,highest), width,bottom-i*30+10 );
				}							
				g.setColor(Color.black);
				for(int j = 0; j<array[0].length; j++) {
					g.fillRect((j-VariantHandler.callSlider.getValue())*barwidth, bottom-MethodLibrary.map(array[Main.drawCanvas.selectedIndex][j],0,highest,0,300), barwidth, MethodLibrary.map(array[Main.drawCanvas.selectedIndex][j],0,highest,0,300));
					g.drawString(""+j,(j-VariantHandler.callSlider.getValue())*barwidth+2,bottom+10);					
				}				
			}
			else {
				VarNode node = FileRead.head.getNextVisible(FileRead.head);
				g.setColor(Color.black);
				Map.Entry<String, ArrayList<SampleNode>> entry;								
			
				while(node != null) {
					for(int v = 0; v<node.vars.size(); v++) {
						entry = node.vars.get(v);
						if(Main.drawCanvas.hideNodeVar(node, entry)) {
							continue;
						}
						for(int m = 0; m<entry.getValue().size(); m++) {
							if(Main.drawCanvas.hideVar(entry.getValue().get(m), node.indel)) {
								continue;
							}
							if(entry.getValue().get(m).alleles != null) {							
								break;
							}
							if(!entry.getValue().get(m).getSample().equals(Main.drawCanvas.selectedSample)) {
								continue;
							}
							if(entry.getKey().length() > 1) {
								continue;
							}
							
							g.fillOval((int)(node.getPosition()/(double)chromsize*width), (int)(this.getHeight()-entry.getValue().get(m).getAlleleFraction()*this.getHeight()), 3,3);
							
						}
					}
					
					node = node.getNextVisible(node);
				}
				/*
				g.setColor(Color.LIGHT_GRAY);
				for(int i = 0;i<simplearray.length;i++) {					
					g.fillRect(i*barwidth, bottom-MethodLibrary.map(simplearray[i],0,highest,0,300), barwidth-1, MethodLibrary.map(simplearray[i],0,highest,0,300));
				}*/
			}			
		}
		catch(Exception e) {
			ErrorLog.addError(e.getStackTrace());
			e.printStackTrace();
		}
	}
}
