package base.BasePlayer;

/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/**          +-
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Nov 20, 2009
 * Time: 3:43:04 PM
 * To change this template use File | Settings | File Templates.
 */
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;
import base.BBfile.BBFileHeader;
import base.BBfile.BBTotalSummaryBlock;
import base.BBfile.BBZoomLevelHeader;
import base.BBfile.BBZoomLevels;
import base.BBfile.BPTree;
import base.BBfile.BigBedIterator;
import base.BBfile.BigWigIterator;
import base.BBfile.RPChromosomeRegion;
import base.BBfile.RPTree;
import base.BBfile.ZoomLevelIterator;

/*
*   Broad Institute Interactive Genome Viewer Big Binary File (BBFile) Reader
*   -   File reader for UCSC BigWig and BigBed file types.
*
*   Notes:   Table entries refer to Jim Kent of UCSC's document description:
*           "BigWig and BigBed: Enabling Browsing of Large Distributed Data Sets",
*           November 2009.
*
*           The overall binary file layout is defined in Table B of the document.
*
*   BBFile Reader sequences through this binary file layout:
*
*   1) Reads in BBFile Header Table C and determine if file is a valid Big Bed or Big Wig
*      binary file type.
*
*   2) Reads in Zoom Header Tables D if zoom data is present, as defined by zoomLevels
*       in Table C, one for each zoom level.
*
*   3) Reads in  the AutoSQL block Table B if present, as referenced by autoSqlOffset in Table C.
*
*   4) Reads in Total Summary Block Table DD if present, as referenced by
*       TotalSummaryOffset in Table C.
*
*   5) Reads in B+ Tree Header Chromosome Index Header Table E, as referenced
*       by chromosomeTreeOffset in Table C.
*
*   6)Reads in B+ Tree Nodes indexing mChromosome ID's for mChromosome regions;
*       Table F for node type (leaf/child), Table G for leaf items,
*       Table H for child node items.
*
*   7) Reads in R+ Tree Chromosome ID Index Header Table K.
*
*   8) Reads in R+ Tree Nodes indexing of data arranged by mChromosome ID's;
*       Table L for node type (leaf or child), Table M for leaf items,
*       Table N for child node items.
*
*   9) Verifies Data Count of data records, as referenced by fullDataOffset in Table C
*
*   10) References data count records of data size defined in Table M of R+ Tree index
*       for all leaf items in the tree.
*
*   11) Reads in zoom level format Table O for each zoom level comprised of
*       zoom record count followed by that many Table P zoom statistics records,
*       followed by an R+ tree of zoom data locations indexed as in Tables L, M, and N.
*
*   12) Returns information on chromosome name keys and chromosome data regions.
*
*   13) Provides iterators using chromosome names and data regions to extract
*       zoom data, Wig data, and Bed data.
* 
* */

public class BBFileReader {

    public static final long BBFILE_HEADER_OFFSET = 0;

   
    private SeekableStream fis;      // BBFile input stream handle
    private long fileOffset;           // file offset for next item to be read

    static BBFileHeader fileHeader; // Big Binary file header
    private boolean isLowToHigh;       // BBFile binary data format: low to high or high to low
    private int uncompressBufSize;     // buffer byte size for data decompression; 0 for uncompressed

    // AutoSQL String defines custom BigBed formats

   
    // This section defines the zoom items if zoom data exists
    private int zoomLevelCount;       // number of zoom levels defined
    private long zoomLevelOffset;      // file offset to zoom level headers
    private BBZoomLevels zoomLevels;   // zoom level headers and data locations

    // Total Summary Block - file statistical info
    private BBTotalSummaryBlock totalSummaryBlock;

    // B+ tree
    private long chromIDTreeOffset; // file offset to mChromosome index B+ tree
    private BPTree chromosomeIDTree;     // Container for the mChromosome index B+ tree

    // R+ tree
    private long chromDataTreeOffset;  // file offset to mChromosome data R+ tree
    private RPTree chromosomeDataTree;     // Container for the mChromosome data R+ tree
    private String autoSql;


    public BBFileReader(String path, BedTrack track) throws IOException {
    	
        this(path, SeekableStreamFactory.getInstance().getStreamFor(path), track);

    }
    public BBFileReader(URL path, BedTrack track) throws IOException {
    	
        this(path.getPath(), SeekableStreamFactory.getInstance().getStreamFor(path), track);

    }
    public BBFileReader(String path, SeekableStream stream, BedTrack track) throws IOException {

    	
       
        fis = stream;
     
        if(track.getFileHeader() == null) {
        // read in file header
	        fileOffset = BBFILE_HEADER_OFFSET;
	        track.setFileHeader(new BBFileHeader(path, fis, fileOffset));
	       
        }
        Main.drawCanvas.loadbarAll = 30;
        Main.drawCanvas.loadBarSample = 30;
        fileHeader = track.getFileHeader();
        //fileHeader.print();
        if (!fileHeader.isHeaderOK()) {
          
            throw new RuntimeException("Error reading BBFile header for: " + path);
        }

        // get data characteristics
        isLowToHigh = fileHeader.isLowToHigh();
        uncompressBufSize = fileHeader.getUncompressBuffSize();

        // update file offset past BBFile header
        fileOffset += BBFileHeader.BBFILE_HEADER_SIZE;

        // get zoom level count from file header
        zoomLevelCount = fileHeader.getZoomLevels();
       
        // extract zoom level headers and zoom data records
        // Note: zoom headers Table D immediately follow the BBFile Header
        if (zoomLevelCount > 0) {

            zoomLevelOffset = fileOffset;
           if(track.getZoomlevels() == null) {
        	   track.setZoomlevels(new BBZoomLevels(fis, zoomLevelOffset, zoomLevelCount,
                    isLowToHigh, uncompressBufSize));
            
           }
           zoomLevels = track.getZoomlevels();
            // end of zoom level headers - compare with next BBFile item location
            fileOffset += zoomLevelCount * BBZoomLevelHeader.ZOOM_LEVEL_HEADER_SIZE;
        }
        if(!Main.drawCanvas.loading) {
        	return;
        }
        Main.drawCanvas.loadbarAll = 60;
        Main.drawCanvas.loadBarSample = 60;
        // get the AutoSQL custom BigBed fields
        
       
        long autoSqlOffset = fileHeader.getAutoSqlOffset();
        if (autoSqlOffset != 0) {
            fis.seek(autoSqlOffset);
            autoSql = readNullTerminatedString(fis);
        }
        
        // get the Total Summary Block (Table DD)
        fileOffset = fileHeader.getTotalSummaryOffset();
        if (fileHeader.getVersion() >= 2 && fileOffset > 0) {
            totalSummaryBlock = new BBTotalSummaryBlock(fis, fileOffset, isLowToHigh);
            fileOffset += BBTotalSummaryBlock.TOTAL_SUMMARY_BLOCK_SIZE;
        }

        // get Chromosome Data B+ Tree (Table E, F, G, H) : should always exist
        chromIDTreeOffset = fileHeader.getChromosomeTreeOffset();
        if (chromIDTreeOffset != 0) {
            fileOffset = chromIDTreeOffset;
            chromosomeIDTree = new BPTree(fis, fileOffset, isLowToHigh);
        }

        // get R+ chromosome data location tree (Tables K, L, M, N)
        chromDataTreeOffset = fileHeader.getFullIndexOffset();
        if (chromDataTreeOffset != 0) {
            fileOffset = chromDataTreeOffset;
            boolean forceDescend = false;
            chromosomeDataTree = new RPTree(fis, fileOffset, isLowToHigh, uncompressBufSize, forceDescend);
        }
        if(!Main.drawCanvas.loading) {
        	return;
        }
        Main.drawCanvas.loadbarAll = 90;
        Main.drawCanvas.loadBarSample = 90;
        
        // get number of data records indexed by the R+ chromosome data location tree
        fileOffset = fileHeader.getFullDataOffset();
    }

    /*
    *   Method returns the Big Binary File header which identifies
    *   the file type and content.
    *
    *   Returns:
    *       Big Binary File header (Table C)
    * */

    public BBFileHeader getBBFileHeader() {
        return fileHeader;
    }

    /*
    *   Method returns if the Big Binary File is BigBed.
    *
    *   Returns:
    *       Boolean identifies if Big Binary File is BigBed
    *       (recognized from magic number in file header Table C)
    * */

    public boolean isBigBedFile() {
        return fileHeader.isBigBed();
    }

    /*
    *   Method returns if the Big Binary File is BigWig
    *
    *   Returns:
    *       Boolean identifies if Big Binary File is BigWig
    *       (recognized from magic number in file header Table C)
    * */

    public boolean isBigWigFile() {
        return fileHeader.isBigWig();
    }


    /*
    *   Method returns if the Big Binary File is written with a low to high byte
    *   order for formatted data.
    *
    *   Returns:
    *       Boolean identifies if Big Binary File is low to high byte order
    *       (recognized from magic number in file header Table C); else is
    *       high to low byte order if false.
    * */

    public boolean isLowToHigh() {
        return isLowToHigh;
    }

    /*
    *   Method returns the total summary block for the Big Binary File.
    *
    *   Returns:
    *       Total summary block data statistics for the whole file (Table DD)
    * */

    public BBTotalSummaryBlock getTotalSummaryBlock() {
        return totalSummaryBlock;
    }


    /*
    *   Method returns the zoom levels in the Big Binary File.
    *
    *   Returns:
    *       Zoom level object containing zoom level headers and R+ zoom data locations tree
    *       (includes Tables  D, O)
    * */

    public BBZoomLevels getZoomLevels() {
        return zoomLevels;
    }


    /*
    *   Method finds chromosome names in the B+ chromosome index tree.
    *
    *   Returns:
    *       LIst of all chromosome key names in the B+ tree.
    * */

    public ArrayList<String> getChromosomeNames() {

        ArrayList<String> chromosomeList = chromosomeIDTree.getChromosomeNames();
        return chromosomeList;
    }

    public String getAutoSql() {
        return autoSql;
    }

    /**
     * Returns an iterator for BigBed features which occupy a chromosome selection region.
     * <p/>
     * Note: the BBFile type should be BigBed; else a null iterator is returned.
     * <p/>
     * Parameters:
     * startChromosome - name of start chromosome
     * startBase     - starting base position for features
     * endChromosome - name of end chromosome
     * endBase       - ending base position for feature
     * contained     - flag specifies bed features must be contained in the specified
     * base region if true; else can intersect the region if false
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for the requested chromosome region.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     * 2) A null object is returned if the file is not BigBed.(see isBigBedFile method)
     */
    synchronized public BigBedIterator getBigBedIterator(String startChromosome, int startBase,
                                                         String endChromosome, int endBase, boolean contained) {

        if (!isBigBedFile())
            return null;


        // go from chromosome names to chromosome ID region
        RPChromosomeRegion selectionRegion = getChromosomeBounds(startChromosome, startBase,
                endChromosome, endBase);

        // check for valid selection region
        if (selectionRegion == null)
            return new BigBedIterator();  // an empty iterator

        // compose an iterator
        BigBedIterator bedIterator = new BigBedIterator(fis, chromosomeIDTree, chromosomeDataTree,
                selectionRegion, contained);

        return bedIterator;
    }


    /**
     * Returns an iterator for BigWig values which occupy the specified startChromosome region.
     * <p/>
     * Note: the BBFile type should be BigWig; else a null iterator is returned.
     * <p/>
     * Parameters:
     * startChromosome  - name of start chromosome
     * startBase    - starting base position for features
     * endChromosome  - name of end chromosome
     * endBase      - ending base position for feature
     * contained    - flag specifies bed features must be contained in the specified
     * base region if true; else can intersect the region if false
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for the requested chromosome region.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     * 2) A null object is returned if the file is not BigWig.(see isBigWigFile method)
     */
    synchronized public BigWigIterator getBigWigIterator(String startChromosome, int startBase,
                                                         String endChromosome, int endBase, boolean contained) {


        if (!isBigWigFile())
            return null;

        // go from chromosome names to chromosome ID region
        RPChromosomeRegion selectionRegion = getChromosomeBounds(startChromosome, startBase,
                endChromosome, endBase);

        // check for valid selection region, return empty iterator if null
        if (selectionRegion == null)
            return new BigWigIterator();

        // compose an iterator
        BigWigIterator wigIterator = new BigWigIterator(fis, chromosomeIDTree, chromosomeDataTree,
                selectionRegion, contained);

        return wigIterator;
    }


    /**
     * Returns an iterator for zoom level records for the chromosome selection region.
     * <p/>
     * Note: the BBFile can be BigBed or BigWig.
     * <p/>
     * Parameters:
     * zoomLevel - zoom level for data extraction; levels start at 1
     * startChromosome - start chromosome name
     * startBase     - staring base position for features
     * endChromosome - end chromosome name
     * endBase       - ending base position for feature
     * contained     - flag specifies bed features must be contained in the
     * specified base region if true; else can intersect the region if false
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for the requested chromosome region.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     */
    synchronized public ZoomLevelIterator getZoomLevelIterator(int zoomLevel, String startChromosome, int startBase,
                                                               String endChromosome, int endBase, boolean contained) {
        // check for valid zoom level
        if (zoomLevel < 1 || zoomLevel > zoomLevelCount)
            throw new RuntimeException("Error: ZoomLevelIterator zoom level is out of range\n");

        // get the appropriate zoom level R+ zoom data index tree
        RPTree zoomDataTree = zoomLevels.getZoomLevelRPTree(zoomLevel);

        // go from chromosome names to chromosome ID region
        RPChromosomeRegion selectionRegion = getChromosomeBounds(startChromosome, startBase,
                endChromosome, endBase);

        // check for valid selection region  
      /*  if (selectionRegion == null) {
            return ZoomLevelIterator.EmptyIterator.theInstance;
        }
*/
        /// compose an iterator
        ZoomLevelIterator zoomIterator = new ZoomLevelIterator(fis, chromosomeIDTree,
                zoomDataTree, zoomLevel, selectionRegion, contained);

        return zoomIterator;
    }

    /**
     * Returns an iterator for zoom level records for all chromosome regions.
     * <p/>
     * Note: the BBFile can be BigBed or BigWig.
     * <p/>
     * Parameters:
     * zoomLevel - zoom level for data extraction; levels start at 1
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for the requested chromosome region.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     */
    synchronized public ZoomLevelIterator getZoomLevelIterator(int zoomLevel) {

        // check for valid zoom level
        if (zoomLevel < 1 || zoomLevel > zoomLevelCount)
            throw new RuntimeException("Error: ZoomLevelIterator zoom level is out of range\n");

        // get the appropriate zoom level R+ zoom data index tree
        RPTree zoomDataTree = zoomLevels.getZoomLevelRPTree(zoomLevel);

        // get all regions bounds
        RPChromosomeRegion selectionRegion = zoomDataTree.getChromosomeBounds();

        // compose an iterator
        boolean contained = true;   //all regions are contained
        ZoomLevelIterator zoomIterator = new ZoomLevelIterator(fis, chromosomeIDTree,
                zoomDataTree, zoomLevel, selectionRegion, contained);

        return zoomIterator;
    }

    /*
    *   Method generates a chromosome bounds region for the supplied chromosome region name.
    *
    *   Note: No attempt is made to verify the region exists in the file data, nor
    *   which data is being examined.
    *
    *   Parameters:
    *       startChromosome - name of start chromosome
    *       startBase - starting base position for region
    *       endChromosome - name of end chromosome
    *       endBase - ending base position for region
    *
    *   Returns:
    *       Chromosome bounds of a named chromosome region for data extraction;
    *       or null for regions not found in the B+ chromosome index tree.
    * */

    private RPChromosomeRegion getChromosomeBounds(String startChromosome, int startBase,
                                                   String endChromosome, int endBase) {

        // If the chromosome name length is > the key size we can't distinguish it
        if (startChromosome.length() > chromosomeIDTree.getKeySize()) {
            return null;
        }

        // find the chromosome ID's using the name to get a valid name key, then associated ID
        String startChromKey = chromosomeIDTree.getChromosomeKey(startChromosome);
        int startChromID = chromosomeIDTree.getChromosomeID(startChromKey);
        if (startChromID < 0)       // mChromosome not in data?
            return null;

        String endChromKey = chromosomeIDTree.getChromosomeKey(endChromosome);
        int endChromID = chromosomeIDTree.getChromosomeID(endChromKey);
        if (endChromID < 0)       // mChromosome not in data?
            return null;

        // create the bounding mChromosome region
        RPChromosomeRegion chromBounds = new RPChromosomeRegion(startChromID, startBase,
                endChromID, endBase);

        return chromBounds;
    }


    private static String readNullTerminatedString(InputStream fis) throws IOException {
        BufferedInputStream bis = new BufferedInputStream(fis, 1000);
        ByteArrayOutputStream bos = new ByteArrayOutputStream(100);
        byte b;
        while ((b = (byte) bis.read()) != 0) {
            bos.write(b);
        }
        return new String(bos.toByteArray());
    }



} // end of BBFileReader
