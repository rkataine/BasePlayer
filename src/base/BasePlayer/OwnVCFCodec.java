
/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package base.BasePlayer;

import htsjdk.tribble.TribbleException;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.util.ParsingUtils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.LazyGenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * A feature codec for the VCF 4 specification
 *
 * <p>
 * VCF is a text file format (most likely stored in a compressed manner). It contains meta-information lines, a
 * header line, and then data lines each containing information about a position in the genome.
 * </p>
 * <p>One of the main uses of next-generation sequencing is to discover variation amongst large populations
 * of related samples. Recently the format for storing next-generation read alignments has been
 * standardised by the SAM/BAM file format specification. This has significantly improved the
 * interoperability of next-generation tools for alignment, visualisation, and variant calling.
 * We propose the Variant Call Format (VCF) as a standarised format for storing the most prevalent
 * types of sequence variation, including SNPs, indels and larger structural variants, together
 * with rich annotations. VCF is usually stored in a compressed manner and can be indexed for
 * fast data retrieval of variants from a range of positions on the reference genome.
 * The format was developed for the 1000 Genomes Project, and has also been adopted by other projects
 * such as UK10K, dbSNP, or the NHLBI Exome Project. VCFtools is a software suite that implements
 * various utilities for processing VCF files, including validation, merging and comparing,
 * and also provides a general Perl and Python API.
 * The VCF specification and VCFtools are available from http://vcftools.sourceforge.net.</p>
 *
 * <p>
 * See also: @see <a href="http://vcftools.sourceforge.net/specs.html">VCF specification</a><br>
 * See also: @see <a href="http://www.ncbi.nlm.nih.gov/pubmed/21653522">VCF spec. publication</a>
 * </p>
 *
 * <h2>File format example</h2>
 * <pre>
 *     ##fileformat=VCFv4.0
 *     #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878
 *     chr1    109     .       A       T       0       PASS  AC=1    GT:AD:DP:GL:GQ  0/1:610,327:308:-316.30,-95.47,-803.03:99
 *     chr1    147     .       C       A       0       PASS  AC=1    GT:AD:DP:GL:GQ  0/1:294,49:118:-57.87,-34.96,-338.46:99
 * </pre>
 *
 * @author Mark DePristo
 * @since 2010
 */
public class OwnVCFCodec extends AbstractVCFCodec {
    // Our aim is to read in the records and convert to VariantContext as quickly as possible, relying on VariantContext to do the validation of any contradictory (or malformed) record parameters.
    public final static String VCF4_MAGIC_HEADER = "##fileformat=VCFv4";
	private static final VCFHeaderVersion VCFHeaderVersion = null;

    /**
     * Reads all of the header from the provided iterator, but no reads no further.
     * @param lineIterator the line reader to take header lines from
     * @return The parsed header
     */
    @SuppressWarnings("static-access")
	@Override
    public Object readActualHeader(final LineIterator lineIterator) {
        final List<String> headerStrings = new ArrayList<String>();

        String line;
        boolean foundHeaderVersion = false;
        while (lineIterator.hasNext()) {
            line = lineIterator.peek();
            lineNo++;
            if (line.startsWith(VCFHeader.METADATA_INDICATOR)) {
                final String[] lineFields = line.substring(2).split("=");
                if (lineFields.length == 2 && VCFHeaderVersion.isFormatString(lineFields[0]) ) {
                    if ( !VCFHeaderVersion.isVersionString(lineFields[1]) )
                        throw new TribbleException.InvalidHeader(lineFields[1] + " is not a supported version");
                    foundHeaderVersion = true;
                    version = VCFHeaderVersion.toHeaderVersion(lineFields[1]);
                    if ( ! version.isAtLeastAsRecentAs(VCFHeaderVersion.VCF4_0) )
                        throw new TribbleException.InvalidHeader("This codec is strictly for VCFv4; please use the VCF3 codec for " + lineFields[1]);
                    if ( version != VCFHeaderVersion.VCF4_0 && version != VCFHeaderVersion.VCF4_1 && version != VCFHeaderVersion.VCF4_2 )
                        throw new TribbleException.InvalidHeader("This codec is strictly for VCFv4 and does not support " + lineFields[1]);
                }
                headerStrings.add(lineIterator.next());
            }
            else if (line.startsWith(VCFHeader.HEADER_INDICATOR)) {
                if (!foundHeaderVersion) {
                    throw new TribbleException.InvalidHeader("We never saw a header line specifying VCF version");
                }
                headerStrings.add(lineIterator.next());
                super.parseHeaderFromLines(headerStrings, version);
                return this.header;
            }
            else {
                throw new TribbleException.InvalidHeader("We never saw the required CHROM header line (starting with one #) for the input VCF file");
            }

        }
        throw new TribbleException.InvalidHeader("We never saw the required CHROM header line (starting with one #) for the input VCF file");
    }
  
    private VariantContext decodeLine(final String line, final boolean includeGenotypes) {
        // the same line reader is not used for parsing the header and parsing lines, if we see a #, we've seen a header line
        if (line.startsWith(VCFHeader.HEADER_INDICATOR)) return null;

        // our header cannot be null, we need the genotype sample names and counts
        if (header == null) throw new TribbleException("VCF Header cannot be null when decoding a record");

        if (parts == null)
            parts = new String[Math.min(header.getColumnCount(), NUM_STANDARD_FIELDS+1)];

 //       final int nParts = ParsingUtils.split(line, parts, VCFConstants.FIELD_SEPARATOR_CHAR, true);

        // if we have don't have a header, or we have a header with no genotyping data check that we
        // have eight columns.  Otherwise check that we have nine (normal columns + genotyping data)
       /* if (( (header == null || !header.hasGenotypingData()) && nParts != NUM_STANDARD_FIELDS) ||
                (header != null && header.hasGenotypingData() && nParts != (NUM_STANDARD_FIELDS + 1)) )
            throw new TribbleException("Line " + lineNo + ": there aren't enough columns for line " + line + " (we expected " + (header == null ? NUM_STANDARD_FIELDS : NUM_STANDARD_FIELDS + 1) +
                    " tokens, and saw " + nParts + " )");
*/
        return parseVCFLine(parts, includeGenotypes);
    }

    /**
     * parse out the VCF line
     *
     * @param parts the parts split up
     * @return a variant context object
     */
    private VariantContext parseVCFLine(final String[] parts, final boolean includeGenotypes) {
        VariantContextBuilder builder = new VariantContextBuilder();
        builder.source(getName());

        lineNo++;

       
        final String chr = getCachedString(parts[0]);
        builder.chr(chr);
        int pos = -1;
        try {
            pos = Integer.valueOf(parts[1]);
        } catch (NumberFormatException e) {
            generateException(parts[1] + " is not a valid start position in the VCF format");
        }
        builder.start(pos);

        if ( parts[2].isEmpty() )
            generateException("The VCF specification requires a valid ID field");
        else if ( parts[2].equals(VCFConstants.EMPTY_ID_FIELD) )
            builder.noID();
        else
            builder.id(parts[2]);

        final String ref = getCachedString(parts[3].toUpperCase());
        final String alts = getCachedString(parts[4]);
        if(parts.length>5 && parts[5] != null) {
	        builder.log10PError(parseQual(parts[5]));
	
	        final List<String> filters = parseFilters(getCachedString(parts[6]));
	        if ( filters != null ) builder.filters(new HashSet<String>(filters));
	        final Map<String, Object> attrs = parseInfo(parts[7]);
	        builder.attributes(attrs);
	
	        if ( attrs.containsKey(VCFConstants.END_KEY) ) {
	            // update stop with the end key if provided
	            try {
	                builder.stop(Integer.valueOf(attrs.get(VCFConstants.END_KEY).toString()));
	            } catch (Exception e) {
	                generateException("the END value in the INFO field is not valid");
	            }
	        } else {
	            builder.stop(pos + ref.length() - 1);
	        }
	
	        // get our alleles, filters, and setup an attribute map
	        final List<Allele> alleles = parseAlleles(ref, alts, lineNo);
	        builder.alleles(alleles);
	
	        // do we have genotyping data
	        if (parts.length > NUM_STANDARD_FIELDS && includeGenotypes) {
	            final LazyGenotypesContext.LazyParser lazyParser = new LazyVCFGenotypesParser(alleles, chr, pos);
	            final int nGenotypes = header.getNGenotypeSamples();
	            LazyGenotypesContext lazy = new LazyGenotypesContext(lazyParser, parts[8], nGenotypes);
	
	            // did we resort the sample names?  If so, we need to load the genotype data
	            if ( !header.samplesWereAlreadySorted() )
	                lazy.decode();
	
	            builder.genotypesNoValidation(lazy);
	        }
        }
        else {
        	final List<Allele> alleles = parseAlleles(ref, alts, lineNo);
        	builder.alleles(alleles);
        	builder.log10PError(parseQual("20"));
        	builder.stop(pos + ref.length() - 1);
        }
        VariantContext vc = null;
        try {
            vc = builder.make();
        } catch (Exception e) {
            generateException(e.getMessage());
        }

        return vc;
    }
    /**
     * parse the filter string, first checking to see if we already have parsed it in a previous attempt
     *
     * @param filterString the string to parse
     * @return a set of the filters applied or null if filters were not applied to the record (e.g. as per the missing value in a VCF)
     */
    protected List<String> parseFilters(final String filterString) {
        // null for unfiltered
        if ( filterString.equals(VCFConstants.UNFILTERED) )
            return null;

        if ( filterString.equals(VCFConstants.PASSES_FILTERS_v4) )
            return Collections.emptyList();
        if ( filterString.equals(VCFConstants.PASSES_FILTERS_v3) )
            generateException(VCFConstants.PASSES_FILTERS_v3 + " is an invalid filter name in vcf4", lineNo);
        if (filterString.isEmpty())
            generateException("The VCF specification requires a valid filter status: filter was " + filterString, lineNo);

        // do we have the filter string cached?
        if ( filterHash.containsKey(filterString) )
            return filterHash.get(filterString);

        // empty set for passes filters
        final List<String> fFields = new LinkedList<String>();
        // otherwise we have to parse and cache the value
        if ( !filterString.contains(VCFConstants.FILTER_CODE_SEPARATOR) )
            fFields.add(filterString);
        else
            fFields.addAll(Arrays.asList(filterString.split(VCFConstants.FILTER_CODE_SEPARATOR)));

        filterHash.put(filterString, Collections.unmodifiableList(fFields));

        return fFields;
    }
    private Map<String, Object> parseInfo(String infoField) {
        Map<String, Object> attributes = new HashMap<String, Object>();

        if ( infoField.isEmpty() )
            generateException("The VCF specification requires a valid (non-zero length) info field");

        if ( !infoField.equals(VCFConstants.EMPTY_INFO_FIELD) ) {
        	
            if ( infoField.indexOf('\t') != -1 || infoField.indexOf(' ') != -1 )
                infoField = infoField.split("\\s+")[0];
            	//generateException("The VCF specification does not allow for whitespace in the INFO field. Offending field value was \"" + infoField + "\"");

            List<String> infoFields = ParsingUtils.split(infoField, VCFConstants.INFO_FIELD_SEPARATOR_CHAR);
            for (int i = 0; i < infoFields.size(); i++) {
                String key;
                Object value;

                int eqI = infoFields.get(i).indexOf("=");
                if ( eqI != -1 ) {
                    key = infoFields.get(i).substring(0, eqI);
                    String valueString = infoFields.get(i).substring(eqI + 1);

                    // split on the INFO field separator
                    List<String> infoValueSplit = ParsingUtils.split(valueString, VCFConstants.INFO_FIELD_ARRAY_SEPARATOR_CHAR);
                    if ( infoValueSplit.size() == 1 ) {
                        value = infoValueSplit.get(0);
                        final VCFInfoHeaderLine headerLine = header.getInfoHeaderLine(key);
                        if ( headerLine != null && headerLine.getType() == VCFHeaderLineType.Flag && value.equals("0") ) {
                            // deal with the case where a flag field has =0, such as DB=0, by skipping the add
                            continue;
                        }
                    } else {
                        value = infoValueSplit;
                    }
                } else {
                    key = infoFields.get(i);
                    final VCFInfoHeaderLine headerLine = header.getInfoHeaderLine(key);
                    if ( headerLine != null && headerLine.getType() != VCFHeaderLineType.Flag ) {
                      /*  if ( GeneralUtils.DEBUG_MODE_ENABLED && ! warnedAboutNoEqualsForNonFlag ) {
                            System.err.println("Found info key " + key + " without a = value, but the header says the field is of type "
                                               + headerLine.getType() + " but this construct is only value for FLAG type fields");
                            warnedAboutNoEqualsForNonFlag = true;
                        }
*/
                        value = VCFConstants.MISSING_VALUE_v4;
                    } else {
                        value = true;
                    }
                }

                // this line ensures that key/value pairs that look like key=; are parsed correctly as MISSING
                if ( "".equals(value) ) value = VCFConstants.MISSING_VALUE_v4;

                attributes.put(key, value);
            }
        }

        return attributes;
    }
    class LazyVCFGenotypesParser implements LazyGenotypesContext.LazyParser {
        final List<Allele> alleles;
        final String contig;
        final int start;

        LazyVCFGenotypesParser(final List<Allele> alleles, final String contig, final int start) {
            this.alleles = alleles;
            this.contig = contig;
            this.start = start;
        }

        @Override
        public LazyGenotypesContext.LazyData parse(final Object data) {
            //System.out.printf("Loading genotypes... %s:%d%n", contig, start);
            return createGenotypeMap((String) data, alleles, contig, start);
        }
    }
    public VariantContext decode(String line) {
        return decodeLine(line, true);
    }
   

    @Override
    public boolean canDecode(final String potentialInput) {
        return canDecodeFile(potentialInput, VCF4_MAGIC_HEADER);
    }
}

   
