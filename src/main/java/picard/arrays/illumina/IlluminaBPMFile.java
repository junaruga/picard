package picard.arrays.illumina;

import picard.PicardException;

import java.io.DataInputStream;
import java.io.IOException;

/**
 * A class to parse the contents of an Illumina Bead Pool Manifest (BPM) file
 *
 * A BPM file contains metadata about the content on an Illumina Genotyping Array
 * Each type of genotyping array has a specific BPM for that array.
 *
 * The BPM contains information including alleles, mapping and normalization
 *
 * This class will parse the binary BPM file format
 */
public class IlluminaBPMFile extends InfiniumDataFile {
    private static final String BPM_IDENTIFIER = "BPM";
    private static final int BPM_IDENTIFIER_LENGTH = 3;

    private String controlConfig = null;
    private IlluminaBPMLocusEntry[] locusEntries = null;
    // TODO - Finish reading the rest of the bpm.  Make other fields accessible.
    // TODO - With this, you don't need to read the bpm.csv into GtcToVcf anymore
    // TODO - only Gtc uses the whole table of contents thing.  Maybe break that out into a base class?

    IlluminaBPMFile(DataInputStream stream) throws IOException {
        super(stream, true);
        parse();
    }

    /**
     * Main parsing method.
     *
     * @throws IOException thrown when there is a problem reading the stream.
     */
    private void parse() throws IOException {

        stream.mark(0);

        try {
            final byte[] curIdentifier = new byte[BPM_IDENTIFIER_LENGTH];
            for (int i = 0; i < curIdentifier.length; i++) {
                curIdentifier[i] = parseByte();
            }

            final String identifier = new String(curIdentifier);
            setIdentifier(identifier);
            if (!identifier.equals(BPM_IDENTIFIER)) {
                throw new PicardException("Invalid identifier '" + identifier + "' for BPM file");
            }
            setFileVersion(parseByte());
            if (getFileVersion() != 1) {
                throw new PicardException("Unknown BPM version (" + getFileVersion() + ")");
            }
            int version = parseInt();
            final int versionFlag = 0x1000;
            if ((version & versionFlag) == versionFlag) {
                version ^= versionFlag;
            }
            if (version > 5 || version < 3) {
                throw new PicardException("Unsupported BPM version (" + version + ")");
            }
            final String manifestName = parseString();
            if (version > 1) {
                controlConfig = parseString();
            }
            readData();
        } finally {
            stream.close();
        }
    }

    private void readData() throws IOException {
        int numLoci = parseInt();

        int[] indices = new int[numLoci];

        // Read the index block
        stream.skipBytes(4 * numLoci);
//        for (int i = 0; i < numLoci; i++) {
//            indices[i] = parseInt();
//        }

        // Read the names
        String[] names = new String[numLoci];
        for (int i = 0; i < numLoci; i++) {
            names[i] = parseString();
        }

        // Read the normalization ids.
        int[] normIds = parseByteArrayAsInts(numLoci);

        // Initialize the locus entries
        locusEntries = new IlluminaBPMLocusEntry[numLoci];
        for (int i = 0; i < numLoci; i++) {
            locusEntries[i] = null;
        }

        // Read the locus entries.
        for (int i = 0; i < numLoci; i++) {
            IlluminaBPMLocusEntry locusEntry = new IlluminaBPMLocusEntry();
            final int version = parseInt();

            locusEntry.ilmnId = parseString();
            locusEntry.name = parseString();
            parseString();
            parseString();
            parseString();
            locusEntry.index = parseInt() - 1;
            parseString();
            locusEntry.ilmnStrand = parseString();
            locusEntry.snp = parseString();
            locusEntry.chrom = parseString();
            locusEntry.ploidy = parseString();
            locusEntry.species = parseString();
            locusEntry.mapInfo = Integer.parseInt(parseString());
            parseString();
            locusEntry.customerStrand = parseString();
            locusEntry.addressA = parseInt();
            locusEntry.addressB = parseInt();
            parseString();
            parseString();
            locusEntry.genomeBuild = parseString();
            locusEntry.source = parseString();
            locusEntry.sourceVersion = parseString();
            locusEntry.sourceStrand = parseString();
            parseString();
            parseByte();
            locusEntry.expClusters = parseByte();
            locusEntry.intensityOnly = parseByte();
            locusEntry.assayType = parseByte();

            if (version >= 7) {
                locusEntry.fracA = parseFloat();
                locusEntry.fracC = parseFloat();
                locusEntry.fracT = parseFloat();
                locusEntry.fracG = parseFloat();
            }
            if (version == 8) {
                locusEntry.refStrand = parseString();
            }

            if (locusEntry.assayType < 0 || locusEntry.assayType > 2) {
                throw new PicardException("Invalid assay_type '" + locusEntry.assayType + "' in BPM file");
            }
            if ((locusEntry.assayType != 0 && locusEntry.addressB == 0) || (locusEntry.assayType == 0 && locusEntry.addressB != 0)) {
                throw new PicardException("Invalid assay_type '" + locusEntry.assayType + "' for address B '" + locusEntry.addressB + "' in BPM file");
            }
            int normId = normIds[locusEntry.index];
            if (normId > 100) {
                throw new PicardException("Invalid normalization ID: " + normId + " for name: " + locusEntry.name);
            }
            locusEntry.normalizationId = normId + 100 * locusEntry.assayType;

            if (locusEntries[locusEntry.index] != null) {
                throw new PicardException("Duplicate locus entry for index: " + locusEntry.index + " '" + locusEntry.name);
            }
            locusEntries[locusEntry.index] = locusEntry;
            if (!names[locusEntry.index].equals(locusEntry.name)) {
                throw new PicardException("Mismatch in names at index: " + locusEntry.index);
            }
        }
    }

    public IlluminaBPMLocusEntry[] getLocusEntries() {
        return locusEntries;
    }
}
