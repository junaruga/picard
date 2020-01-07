package picard.arrays;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.arrays.illumina.InfiniumVcfFields;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * A simple program to combine multiple genotyping array VCF files into one VCF
 *
 */
@CommandLineProgramProperties(
        summary = CombineGenotypingArrayVcfs.USAGE_DETAILS,
        oneLineSummary = "Program to combine multiple genotyping array VCF files into one VCF.",
        programGroup = picard.cmdline.programgroups.GenotypingArraysProgramGroup.class
)
@DocumentedFeature
public class CombineGenotypingArrayVcfs extends CommandLineProgram {
    static final String USAGE_DETAILS =
            "CombineGenotypingArrayVcfs takes one or more VCF files, as generated by GtcToVcf " +
                    "and combines them into a single VCF. " +
                    "The input VCFs must have the same sequence dictionary and same list of variant loci. " +
                    "The input VCFs must not share sample Ids. " +
                    "<h4>Usage example:</h4>" +
                    "<pre>" +
                    "java -jar picard.jar VcfToAdpc \\<br />" +
                    "      INPUT=input1.vcf \\<br />" +
                    "      INPUT=input2.vcf \\<br />" +
                    "      OUTPUT=output.vcf" +
                    "</pre>";

    @Argument(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME,  doc="Input VCF file(s).")
    public List<File> INPUT;

    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output VCF file.")
    public File OUTPUT;

    private final Log log = Log.getInstance(CombineGenotypingArrayVcfs.class);

    final private ProgressLogger progressLogger = new ProgressLogger(log, 10000);

    public CombineGenotypingArrayVcfs() {
        CREATE_INDEX = true;
    }

    // These items will be removed from the merged header
    private final Set<String> sampleSpecificHeaders = new HashSet<String>(Arrays.asList(
            InfiniumVcfFields.ANALYSIS_VERSION_NUMBER,
            InfiniumVcfFields.AUTOCALL_DATE,
            InfiniumVcfFields.AUTOCALL_GENDER,
            InfiniumVcfFields.CHIP_WELL_BARCODE,
            InfiniumVcfFields.EXPECTED_GENDER,
            InfiniumVcfFields.FINGERPRINT_GENDER,
            InfiniumVcfFields.IMAGING_DATE,
            InfiniumVcfFields.P_95_GREEN,
            InfiniumVcfFields.P_95_RED,
            InfiniumVcfFields.SAMPLE_ALIAS,
            InfiniumVcfFields.SCANNER_NAME,
            "Biotin(Bgnd)", "Biotin(High)",
            "DNP(Bgnd)", "DNP(High)", "Extension(A)", "Extension(C)", "Extension(G)", "Extension(T)",
            "Hyb(High)", "Hyb(Low)", "Hyb(Medium)", "NP(A)", "NP(C)", "NP(G)", "NP(T)",
            "NSB(Bgnd)Blue", "NSB(Bgnd)Green", "NSB(Bgnd)Purple", "NSB(Bgnd)Red", "Restore",
            "String(MM)", "String(PM)", "TargetRemoval",
            "fileDate"));


    @Override
    public int doWork() {
        log.info("Checking inputs.");
        final List<File> UNROLLED_INPUT = IOUtil.unrollFiles(INPUT, IOUtil.VCF_EXTENSIONS);
        for (final File f: UNROLLED_INPUT) IOUtil.assertFileIsReadable(f);
        IOUtil.assertFileIsWritable(OUTPUT);

        final SAMSequenceDictionary sequenceDictionary = VCFFileReader.getSequenceDictionary(UNROLLED_INPUT.get(0));

        final List<String> sampleList = new ArrayList<String>();
        final Collection<CloseableIterator<VariantContext>> iteratorCollection = new ArrayList<>(UNROLLED_INPUT.size());
        final Collection<VCFHeader> headers = new HashSet<VCFHeader>(UNROLLED_INPUT.size());

        Set<String> sampleNames = new HashSet<>();

        for (final File file : UNROLLED_INPUT) {
            final VCFFileReader fileReader = new VCFFileReader(file, false);
            final VCFHeader fileHeader = fileReader.getFileHeader();

            for (final String sampleName : fileHeader.getSampleNamesInOrder()) {
                if (!sampleNames.add(sampleName)) {
                    throw new IllegalArgumentException("Input file " + file.getAbsolutePath() + " contains a sample entry (" + sampleName + ") that appears in another input file.");
                }
                sampleList.add(sampleName);
            }

            headers.add(fileHeader);
            iteratorCollection.add(fileReader.iterator());
        }

        if (CREATE_INDEX && sequenceDictionary == null) {
            throw new PicardException("A sequence dictionary must be available (either through the input file or by setting it explicitly) when creating indexed output.");
        }

        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setOutputFile(OUTPUT)
                .setReferenceDictionary(sequenceDictionary);
        if (CREATE_INDEX) {
            builder.setOption(Options.INDEX_ON_THE_FLY);
        }
        final VariantContextWriter writer = builder.build();

        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(headers, false);
        headerLines.removeIf(line -> sampleSpecificHeaders.contains(line.getKey()));
        writer.writeHeader(new VCFHeader(headerLines, sampleList));

        int closedIteratorCount = 0;
        while (closedIteratorCount == 0) {
            List<VariantContext> variantContexts = new ArrayList<>();
            for (final CloseableIterator<VariantContext> iterator: iteratorCollection) {
                if (iterator.hasNext()) {
                    variantContexts.add(iterator.next());
                } else {
                    closedIteratorCount++;
                }
            }
            if (closedIteratorCount == 0) {
                progressLogger.record(variantContexts.get(0).getContig(), variantContexts.get(0).getStart());
                writer.add(merge(variantContexts));
            }
        }
        if (closedIteratorCount != iteratorCollection.size()) {
            throw new PicardException("Mismatch in number of variants among input VCFs");
        }
        writer.close();
        return 0;
    }

    /**
     * Merges multiple VariantContexts all for the same locus into a single hybrid.
     *
     * @param variantContexts           list of VCs
     * @return new VariantContext       representing the merge of variantContexts
     */
    public static VariantContext merge(final List<VariantContext> variantContexts) {
        if ( variantContexts == null || variantContexts.isEmpty() )
            return null;

        // establish the baseline info from the first VC
        final VariantContext first = variantContexts.get(0);
        final String name = first.getSource();

        final Set<String> filters = new HashSet<>();

        int depth = 0;
        double log10PError = CommonInfo.NO_LOG10_PERROR;
        boolean anyVCHadFiltersApplied = false;
        GenotypesContext genotypes = GenotypesContext.create();

        // Go through all the VCs, verify that the loci and ID and other attributes agree.
        final Map<String, Object> firstAttributes = first.getAttributes();
        for (final VariantContext vc : variantContexts ) {
            if ((vc.getStart() != first.getStart()) || (!vc.getContig().equals(first.getContig()))) {
                throw new PicardException("Mismatch in loci among input VCFs");
            }
            if (!vc.getID().equals(first.getID())) {
                throw new PicardException("Mismatch in ID field among input VCFs");
            }
            if (!vc.getReference().equals(first.getReference())) {
                throw new PicardException("Mismatch in REF allele among input VCFs");
            }
            checkThatAllelesMatch(vc, first);

            for (final Genotype g : vc.getGenotypes()) {
                genotypes.add(g);
            }

            // We always take the QUAL of the first VC with a non-MISSING qual for the combined value
            if ( log10PError == CommonInfo.NO_LOG10_PERROR )
                log10PError =  vc.getLog10PError();

            filters.addAll(vc.getFilters());
            anyVCHadFiltersApplied |= vc.filtersWereApplied();

            // add attributes
            // special case DP (add it up)
            if (vc.hasAttribute(VCFConstants.DEPTH_KEY))
                depth += vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);

            // Go through all attributes - if any differ (other than AC, AF, or AN) that's an error.  (recalc AC, AF, and AN later)
            for (final Map.Entry<String, Object> p : vc.getAttributes().entrySet()) {
                final String key = p.getKey();
                if ((!key.equals("AC")) && (!key.equals("AF")) && (!key.equals("AN")) && (!key.equals("devX_AB")) && (!key.equals("devY_AB"))) {
                    final Object value = p.getValue();
                    final Object extantValue = firstAttributes.get(key);
                    if (extantValue == null) {
                        // attribute in one VCF but not another.  Die!
                        throw new PicardException("Attribute '" + key + "' not found in all VCFs");
                    }
                    else if (!extantValue.equals(value)) {
                        // Attribute disagrees in value between one VCF Die! (if not AC, AF, nor AN, skipped above)
                        throw new PicardException("Values for attribute '" + key + "' disagrees among input VCFs");
                    }
                }
            }
        }

        if (depth > 0)
            firstAttributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));

        final VariantContextBuilder builder = new VariantContextBuilder().source(name).id(first.getID());
        builder.loc(first.getContig(), first.getStart(), first.getEnd());
        builder.alleles(first.getAlleles());
        builder.genotypes(genotypes);

        builder.attributes(new TreeMap<>(firstAttributes));
        // AC, AF, and AN are recalculated here
        VariantContextUtils.calculateChromosomeCounts(builder, false);

        builder.log10PError(log10PError);
        if ( anyVCHadFiltersApplied ) {
            builder.filters(filters.isEmpty() ? filters : new TreeSet<>(filters));
        }

        return builder.make();
    }

    private static void checkThatAllelesMatch(final VariantContext vc1, final VariantContext vc2) {
        if (!vc1.getReference().equals(vc2.getReference())) {
            throw new PicardException("Mismatch in REF allele among input VCFs");
        }
        if (vc1.getAlternateAlleles().size() != vc2.getAlternateAlleles().size()) {
            throw new PicardException("Mismatch in ALT allele count among input VCFs");
        }
        for (int i = 0; i < vc1.getAlternateAlleles().size(); i++) {
            if (!vc1.getAlternateAllele(i).equals(vc2.getAlternateAllele(i))) {
                throw new PicardException("Mismatch in ALT allele among input VCFs for " + vc1.getContig() + "." + vc1.getStart());
            }
        }
    }
}
