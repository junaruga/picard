package picard.fingerprint;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.util.MathUtil;
import picard.vcf.VcfTestUtils;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import static picard.util.TestNGUtil.compareDoubleWithAccuracy;

/**
 * Created by farjoun on 8/27/15.
 */
public class FingerprintCheckerTest {

    private final double maf = 0.4;
    private final Snp snp = new Snp("test", "chr1", 1, (byte) 'A', (byte) 'C', maf, Collections.singletonList("dummy"));
    private final HaplotypeBlock hb = new HaplotypeBlock(maf);

    private static final double DELTA = 1e-6;

    private static final File TEST_DATA_DIR = new File("testdata/picard/fingerprint/");
    private static final File SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING = new File(TEST_DATA_DIR, "Homo_sapiens_assembly19.haplotype_database.subset.txt");

    @BeforeClass
    public void setup() {
        hb.addSnp(snp);
    }



    @DataProvider(name = "pLoH")
    public Iterator<Object[]> pLohData() {
        final List<Object[]> listOfDoubles = new ArrayList<>();

        for (int i = 1; i < 20; i++) {
            listOfDoubles.add(new Object[]{i / 40D});
        }
        return listOfDoubles.iterator();
    }

    @Test(dataProvider = "pLoH")
    public void testMatchResults(final double pLoH) {

        final Fingerprint fpObserved = new Fingerprint("test", null, "noop");
        final Fingerprint fpExpected = new Fingerprint("test", null, "noop");

        final HaplotypeProbabilities hpHet = new HaplotypeProbabilitiesFromGenotype(snp, hb, 0.0001, 1.0, 0.0001);
        final HaplotypeProbabilities hpHomRef = new HaplotypeProbabilitiesFromGenotype(snp, hb, 1.0, 0.00001, 0.000000001);

        // Expected is a het
        fpExpected.add(hpHet);

        // Observed is a hom, so possible scenario is that observed is tumor, and expected is normal
        fpObserved.add(hpHomRef);

        // get match results using pLOD
        final MatchResults mr = FingerprintChecker.calculateMatchResults(fpObserved, fpExpected, 0.01, pLoH);

        // make sure that it's more likely to be the same sample, if the observed is "tumor" and the expected is "normal"
        Assert.assertTrue(mr.getLodTN() > mr.getLOD());

        // make sure that the regular LOD is negative (we're comparing a HET to a HOM)
        Assert.assertTrue(mr.getLOD() < 0);

        // make sure that it's more likely to be tumor/normal rather than normal/tumor
        // (a hom normal isn't expected to be measured as a het in the tumor)
        Assert.assertTrue(mr.getLodTN() > mr.getLodNT());
    }

    @DataProvider(name = "checkFingerprintsVcfDataProvider")
    public Object[][] testCheckFingerprintsVcfDataProvider() {
        return new Object[][]{
                {new File(TEST_DATA_DIR, "NA12891.vcf"), new File(TEST_DATA_DIR, "NA12891.fp.vcf"), "NA12891", "NA12891", -1.048021, -2.053484,  1.005462},
                {new File(TEST_DATA_DIR, "NA12891.vcf"), new File(TEST_DATA_DIR, "NA12891.g.vcf"),  "NA12891", "NA12891", -1.037564, -2.049586,  1.012022},
                {new File(TEST_DATA_DIR, "NA12892.vcf"), new File(TEST_DATA_DIR, "NA12892.fp.vcf"), "NA12892", "NA12892", -1.105025, -2.166160,  1.061135},
                {new File(TEST_DATA_DIR, "NA12892.vcf"), new File(TEST_DATA_DIR, "NA12892.g.vcf"),  "NA12892", "NA12892", -1.094570, -2.162798,  1.068227},
                {new File(TEST_DATA_DIR, "NA12891.vcf"), new File(TEST_DATA_DIR, "NA12892.fp.vcf"), "NA12891", "NA12892", -7.024770, -2.109822, -4.914948},
                {new File(TEST_DATA_DIR, "NA12891.vcf"), new File(TEST_DATA_DIR, "NA12892.g.vcf"),  "NA12891", "NA12892", -7.718515, -2.106459, -5.612055},
                {new File(TEST_DATA_DIR, "NA12892.vcf"), new File(TEST_DATA_DIR, "NA12891.fp.vcf"), "NA12892", "NA12891", -7.024770, -2.109822, -4.914948},
                {new File(TEST_DATA_DIR, "NA12892.vcf"), new File(TEST_DATA_DIR, "NA12891.g.vcf"),  "NA12892", "NA12891", -7.679670, -2.105924, -5.573746},
                {new File(TEST_DATA_DIR, "emptyNA12892.vcf"), new File(TEST_DATA_DIR, "NA12891.g.vcf"),  "NA12892", "NA12891", 0, 0, 0},
        };
    }

    @Test(dataProvider = "checkFingerprintsVcfDataProvider")
    public void testCheckFingerprintsVcf(final File vcfFile, final File genotypesFile, final String observedSampleAlias, final String expectedSampleAlias,
                                      final double llExpectedSample, final double llRandomSample, final double lodExpectedSample) throws IOException {
        final Path indexedInputVcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(vcfFile, "fingerprintcheckertest.tmp.").toPath();
        final Path indexedGenotypesVcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(genotypesFile, "fingerprintcheckertest.tmp.").toPath();

        final FingerprintChecker fpChecker = new FingerprintChecker(SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING);
        final List<FingerprintResults> results = fpChecker.checkFingerprintsFromPaths(Collections.singletonList(indexedInputVcf),
                Collections.singletonList(indexedGenotypesVcf),
                observedSampleAlias,
                expectedSampleAlias);
        Assert.assertEquals(results.size(), 1);
        final FingerprintResults fpr = results.get(0);
        Assert.assertNull(fpr.getReadGroup());
        Assert.assertEquals(fpr.getSampleAlias(), observedSampleAlias);
        final MatchResults mr = fpr.getMatchResults().first();
        Assert.assertEquals(mr.getSample(), expectedSampleAlias);
        Assert.assertEquals(mr.getSampleLikelihood(), llExpectedSample, DELTA);
        Assert.assertEquals(mr.getPopulationLikelihood(), llRandomSample, DELTA);
        Assert.assertEquals(mr.getLOD(), lodExpectedSample, DELTA);
    }

    @Test(dataProvider = "checkFingerprintsVcfDataProvider")
    public void testFingerprintVcf(final File vcfFile, final File genotypesFile, final String observedSampleAlias, final String expectedSampleAlias,
                                   final double llExpectedSample, final double llRandomSample, final double lodExpectedSample) {
        final FingerprintChecker fpChecker = new FingerprintChecker(SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING);
        final Map<FingerprintIdDetails, Fingerprint> fp1 = fpChecker.fingerprintVcf(vcfFile.toPath());

        Assert.assertFalse(fp1.isEmpty());
    }


    @Test(dataProvider = "checkFingerprintsVcfDataProvider")
    public void testFingerprintSwapEqual(final File vcfFile, final File genotypesFile, final String observedSampleAlias, final String expectedSampleAlias,
                                   final double llExpectedSample, final double llRandomSample, final double lodExpectedSample) {
        final FingerprintChecker fpChecker = new FingerprintChecker(SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING);
        final Map<FingerprintIdDetails, Fingerprint> fp1Map = fpChecker.fingerprintVcf(vcfFile.toPath());
        final Map<FingerprintIdDetails, Fingerprint> fp2Map = fpChecker.fingerprintVcf(genotypesFile.toPath());

        for(Fingerprint fp1:fp1Map.values()){
            for(Fingerprint fp2:fp2Map.values()) {
                final MatchResults matchResults12 = FingerprintChecker.calculateMatchResults(fp1, fp2);
                final MatchResults matchResults21 = FingerprintChecker.calculateMatchResults(fp2, fp1);
                compareDoubleWithAccuracy(matchResults12.getLOD(),matchResults21.getLOD(),1e-10);
            }
        }
    }

    @Test(expectedExceptions = PicardException.class)
    public void testTerminateOnBadFile() {
        final FingerprintChecker fpChecker = new FingerprintChecker(SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING);
        final File badSam = new File(TEST_DATA_DIR, "aligned_queryname_sorted.sam");
        fpChecker.fingerprintFiles(Collections.singletonList(badSam.toPath()), 1, 1, TimeUnit.DAYS);
    }

    @DataProvider(name = "checkFingerprintsSamDataProvider")
    public Object[][] testCheckFingerprintsSamDataProvider() {
        final File na12891_r1 = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r1.sam");
        final File na12891_r2 = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r2.sam");
        final File na12892_r1 = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r1.sam");
        final File na12892_r2 = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r2.sam");

        final File na12891_noRg = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.noRgTag.sam");

        return new Object[][]{
                {na12891_r1, na12891_r2, true, true},
                {na12892_r1, na12892_r2, true, true},
                {na12892_r1, na12891_r2, false, true},
                {na12891_r1, na12891_noRg, true, true},
                {na12892_r1, na12891_noRg, false, true},

                {na12891_r2, na12891_r2, true, false},
                {na12892_r2, na12892_r2, true, false},
                {na12892_r2, na12891_r2, false, false},
                {na12891_r2, na12891_noRg, true, false},
                {na12892_r2, na12891_noRg, false, false},
        };
    }

    @Test(dataProvider = "checkFingerprintsSamDataProvider")
    public void testCheckFingerprintsSam(final File samFile1, final File samFile2, final boolean expectedMatch, final boolean silent) throws IOException {

        final File metricsFile = File.createTempFile("crosscheck",".crosscheck_metrics");
        metricsFile.deleteOnExit();

        final String[] args = {
                "EXPECT_ALL_GROUPS_TO_MATCH=true",
                "LOD_THRESHOLD=-1",
                "H=" + SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING.getAbsolutePath(),
                "I=" + samFile1.getAbsolutePath(),
                "I=" + samFile2.getAbsolutePath(),
                "VALIDATION_STRINGENCY=" + (silent ? "SILENT" : "LENIENT"),
                "CROSSCHECK_BY=FILE",
                "OUTPUT="+metricsFile.getAbsolutePath()
        };

        Assert.assertEquals(new CrosscheckFingerprints().instanceMain(args), expectedMatch ? 0 : 1);

        // this part checks that the results are symmetric (i.e LOD x,y == LOD y,x)
        final MetricsFile<CrosscheckMetric, Double> metricsFileReader = new MetricsFile<>();
        metricsFileReader.read(new FileReader(metricsFile));
        final List<CrosscheckMetric> metrics = metricsFileReader.getMetrics();

        final Map<Set<String>, Set<String>> collected = metrics.stream()
                .collect(Collectors.groupingBy(s -> CollectionUtil.makeSet(s.LEFT_GROUP_VALUE, s.RIGHT_GROUP_VALUE), Collectors.mapping(s -> s.LOD_SCORE.toString(), Collectors.toSet())));

        for (Map.Entry<Set<String>, Set<String>> entry : collected.entrySet()) {
            if (entry.getValue().size() > 1) {

                final List<CrosscheckMetric> mismatchingMetrics = metrics.stream()
                        .filter(s -> CollectionUtil.makeSet(s.LEFT_GROUP_VALUE,s.RIGHT_GROUP_VALUE).equals(entry.getKey())).collect(Collectors.toList());

                Assert.fail("Metrics disagree: LOD scores are: \n[" + String.join(", ", entry.getValue()) +
                        "],\n from the following metrics: \n" + mismatchingMetrics.get(0) + mismatchingMetrics.get(1));
            }
        }
    }

    @DataProvider(name = "checkFingerprintsSamDataProviderFail")
    public Object[][] testCheckFingerprintsSamDataProviderFail() {
        final File na12891_r1 = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r1.sam");
        final File na12892_r1 = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r1.sam");
        final File na12891_noRg = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.noRgTag.sam");

        return new Object[][]{
                {na12892_r1, na12891_noRg, false},
                {na12891_r1, na12891_noRg, true},
        };
    }

    @Test(dataProvider = "checkFingerprintsSamDataProviderFail", expectedExceptions = PicardException.class)
    public void testCheckFingerprintsFail(final File samFile1, final File samFile2, final boolean expectedMatch) {
        final String[] args = {
                "EXPECT_ALL_GROUPS_TO_MATCH=true",
                "LOD_THRESHOLD=-1",
                "H=" + SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING.getAbsolutePath(),
                "I=" + samFile1.getAbsolutePath(),
                "I=" + samFile2.getAbsolutePath(),
                "VALIDATION_STRINGENCY=STRICT",
                "CROSSCHECK_BY=FILE",
        };

        new CrosscheckFingerprints().instanceMain(args);
    }

    @DataProvider(name = "queryableData")
    public Iterator<Object[]> queryableData() throws IOException {
        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{new File(TEST_DATA_DIR, "NA12891.fp.vcf"), false});
        tests.add(new Object[]{new File(TEST_DATA_DIR, "NA12891.vcf"), false});
        tests.add(new Object[]{VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891.vcf"), "fingerprintcheckertest.tmp."), true});
        tests.add(new Object[]{VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891.vcf.gz"), "fingerprintcheckertest.tmp."), true});

        return tests.iterator();
    }

    @Test(dataProvider = "queryableData")
    public void testQueryable(final File vcf, boolean expectedQueryable) {
        try(VCFFileReader reader = new VCFFileReader(vcf, false)) {
            Assert.assertEquals(reader.isQueryable(), expectedQueryable);
        }
    }

    @Test
    public void testWriteFingerprint() throws IOException {
        final File haplotype_db = new File(TEST_DATA_DIR, "haplotypeMap_small.vcf");
        final File vcfInput = new File(TEST_DATA_DIR, "testSample_small.vcf");
        final File fasta = new File(TEST_DATA_DIR, "reference.fasta");
        final File vcfExpected = new File(TEST_DATA_DIR, "expectedFingerprint_small.vcf");
        final FingerprintChecker fpchecker = new FingerprintChecker(haplotype_db);
        final Fingerprint fp = fpchecker.fingerprintVcf(vcfInput.toPath()).values().iterator().next();

        final File vcfOutput = File.createTempFile("fingerprint", ".vcf");
        FingerprintUtils.writeFingerPrint(fp, vcfOutput, fasta, "Dummy", "Testing");

        VcfTestUtils.assertVcfFilesAreEqual(vcfOutput, vcfExpected);
    }
}