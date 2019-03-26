/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.arrays;

import picard.arrays.illumina.ArraysControlInfo;
import picard.arrays.illumina.ExtendedIlluminaManifest;
import picard.arrays.illumina.ExtendedIlluminaManifestRecord;
import picard.arrays.illumina.IlluminaManifestRecord;
import picard.arrays.illumina.InfiniumEGTFile;
import picard.arrays.illumina.InfiniumGTCFile;
import picard.arrays.illumina.InfiniumNormalizationManifest;
import picard.arrays.illumina.InfiniumVcfFields;
import picard.pedigree.Sex;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFRecordCodec;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

/**
 * Class to convert a GTC file and a BPM file to a VCF file.
 */
@CommandLineProgramProperties(
        summary = GtcToVcf.USAGE_DETAILS,
        oneLineSummary = "Program to convert a GTC file to a VCF",
        programGroup = picard.cmdline.programgroups.GenotypingArraysProgramGroup.class
)
public class GtcToVcf extends CommandLineProgram {


    static final String USAGE_DETAILS =
            "<p>GtcToVcf takes an Illumina GTC file and converts it to a VCF file using several supporting files.</p>" +
                    "<p>A GTC file is an abbreviation for Genotype Call file, a tab-delimited text file containing data about gene expression.</p>" +
                    "<p>A VCF, aka Variant Calling Format, is a text file for storing how a sequenced sample differs from the reference genome. " +
                    "<a href='http://software.broadinstitute.org/software/igv/book/export/html/184'></a></p>" +
                    "<p>The name of the INPUT must match the name of the EXTENDED_ILLUMINA_MANIFEST, sans file extensions.</p>" +
                    "<h3>Usage example:</h3>" +
                    "<pre>" +
                    "GtcToVcf\n" +
                    "\t--INPUT input.gtc\n" +
                    "\t--REFERENCE_SEQUENCE reference.fasta\n" +
                    "\t--OUTPUT output.vcf\n" +
                    "\t--EXTENDED_ILLUMINA_MANIFEST input.csv\n" +
                    "\t--CLUSTER_FILE cluster.egt\n" +
                    "\t--ILLUMINA_NORMALIZATION_MANIFEST normalization_manifest.bpm.csv\n" +
                    "\t--SAMPLE_ALIAS my_sample_alias\n" +
                    "\t--ANALYSIS_VERSION_NUMBER 0\n" +
                    "\t--EXPECTED_GENDER UNKNOWN\n" +
                    "\t--GENDER_GTC gender.gtc\n" +
                    "\t--ZCALL_THRESHOLDS_FILE thresholds.\n" +
                    "\t--FINGERPRINT_GENOTYPES_VCF_FILE fingerprint_genotypes.vcf" +
                    "</pre><hr />";

    private final Log log = Log.getInstance(GtcToVcf.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "GTC file to be converted")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output VCF file to write.")
    public File OUTPUT;

    @Argument(shortName = "MANIFEST", doc = "An Extended Illumina Manifest file (csv)")
    public File EXTENDED_ILLUMINA_MANIFEST;

    @Argument(shortName = "CF", doc = "An Infinium cluster file (egt)")
    public File CLUSTER_FILE;

    @Argument(shortName = "NORM_MANIFEST", doc = "An Infinium manifest containing the normalization ids (bpm.csv)")
    public File ILLUMINA_NORMALIZATION_MANIFEST;

    @Argument(shortName = "E_GENDER", doc = "The expected gender for this sample.")
    public String EXPECTED_GENDER;

    @Argument(doc = "The sample alias")
    public String SAMPLE_ALIAS;

    @Argument(doc = "The analysis version of the data used to generate this VCF")
    public Integer ANALYSIS_VERSION_NUMBER;

    @Argument(shortName = "G_GTC", doc = "GTC file that was called using the gender cluster file.", optional = true)
    public File GENDER_GTC;

    @Argument(shortName = "ZCALL_T_FILE", doc = "The zcall thresholds file.", optional = true)
    public File ZCALL_THRESHOLDS_FILE = null;

    @Argument(shortName = "FP_VCF", doc = "The fingerprint VCF for this sample", optional = true)
    public File FINGERPRINT_GENOTYPES_VCF_FILE;

    private static final List<Allele> NO_CALL_ALLELES = Collections.unmodifiableList(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL));

    // This file gets initialized during customCommandLineValidation.
    // It is a static member so we don't have to parse the file twice.
    private static InfiniumGTCFile infiniumGTCFile;

    private static InfiniumEGTFile infiniumEGTFile;

    private static String gtcGender = null;

    private static HashMap<String, String[]> zCallThresholds = new HashMap<>();

    private static Sex fingerprintGender;

    private static final DecimalFormat df = new DecimalFormat();

    private static final String dot = ".";

    static {
        df.setMaximumFractionDigits(3);
    }


    @Override
    protected boolean requiresReference() {
        return true;
    }

    @Override
    protected int doWork() {
        final ExtendedIlluminaManifest manifest = setupAndGetManifest();

        try (final ReferenceSequenceFile refSeq = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE)) {

            final VCFHeader vcfHeader = createVCFHeader(manifest, infiniumGTCFile, gtcGender, CLUSTER_FILE,
                    REFERENCE_SEQUENCE, refSeq.getSequenceDictionary());

            // Setup a collection that will sort contexts properly
            // Necessary because input GTC file is not sorted
            final SortingCollection<VariantContext> contexts =
                    SortingCollection.newInstance(
                            VariantContext.class,
                            new VCFRecordCodec(vcfHeader),
                            new VariantContextComparator(refSeq.getSequenceDictionary()),
                            MAX_RECORDS_IN_RAM,
                            TMP_DIR.stream().map(File::toPath).toArray(Path[]::new));

            // fill the sorting collection
            fillContexts(contexts, infiniumGTCFile, manifest, infiniumEGTFile);

            writeVcf(contexts, OUTPUT, refSeq.getSequenceDictionary(), vcfHeader);

        } catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        }

        return 0;
    }

    @Override
    protected String[] customCommandLineValidation() {

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(EXTENDED_ILLUMINA_MANIFEST);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (FINGERPRINT_GENOTYPES_VCF_FILE != null) {
            IOUtil.assertFileIsReadable(FINGERPRINT_GENOTYPES_VCF_FILE);
        }
        if (GENDER_GTC != null) {
            IOUtil.assertFileIsReadable(GENDER_GTC);
        }
        if (ZCALL_THRESHOLDS_FILE != null) {
            IOUtil.assertFileIsReadable(ZCALL_THRESHOLDS_FILE);
        }

        return super.customCommandLineValidation();
    }

    private ExtendedIlluminaManifest setupAndGetManifest() {
        fingerprintGender = getFingerprintSex(FINGERPRINT_GENOTYPES_VCF_FILE);
        final InfiniumNormalizationManifest infiniumNormalizationManifest = new InfiniumNormalizationManifest(ILLUMINA_NORMALIZATION_MANIFEST);
        try (final DataInputStream gtcInputStream = new DataInputStream(new FileInputStream(INPUT))) {

            infiniumEGTFile = new InfiniumEGTFile(CLUSTER_FILE);
            infiniumGTCFile = new InfiniumGTCFile(gtcInputStream, infiniumNormalizationManifest);
            final ExtendedIlluminaManifest manifest = new ExtendedIlluminaManifest(EXTENDED_ILLUMINA_MANIFEST);

            if (GENDER_GTC != null) {
                try (DataInputStream genderGtcStream = new DataInputStream(new FileInputStream(GENDER_GTC))) {
                    gtcGender = new InfiniumGTCFile(genderGtcStream, infiniumNormalizationManifest).getGender();
                }
            }

            if (ZCALL_THRESHOLDS_FILE != null) {
                parseZCallThresholds();
            }

            final String gtcManifestName = FilenameUtils.removeExtension(infiniumGTCFile.getSnpManifest());
            final String illuminaManifestName = FilenameUtils.removeExtension(manifest.getDescriptorFileName());

            if (!gtcManifestName.equalsIgnoreCase(illuminaManifestName)) {
                throw new PicardException("The GTC's manifest name " + gtcManifestName +
                        " does not match the Illumina manifest name " + illuminaManifestName);
            }

            if (infiniumGTCFile.getNumberOfSnps() != manifest.getNumAssays()) {
                log.warn("The number of SNPs in the GTC file: " + infiniumGTCFile.getNumberOfSnps() +
                        " does not equal the number of SNPs in the Illumina manifest file: " + manifest.getNumAssays());
            }
            return manifest;
        } catch (IOException e) {
            throw new PicardException("Failed to setup");
        }
    }

    private void parseZCallThresholds() {
        final String notApplicable = "NA";
        try (Stream<String> stream = Files.lines(ZCALL_THRESHOLDS_FILE.toPath())) {

            stream.forEach(line -> {
                String[] tokens = line.split("\t");
                if ((!tokens[1].equals(notApplicable)) && (!tokens[2].equals(notApplicable))) {
                    zCallThresholds.put(tokens[0], new String[]{tokens[1], tokens[2]});
                } else {
                    zCallThresholds.put(tokens[0], new String[]{dot, dot});
                }
            });

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    Sex getFingerprintSex(final File file) {
        if (file == null) {
            return Sex.Unknown;
        } else {
            try (VCFFileReader reader = new VCFFileReader(file, false)) {
                final VCFHeader header = reader.getFileHeader();
                final VCFHeaderLine gender = header.getMetaDataLine("gender");
                if (gender != null) {
                    return Sex.valueOf(gender.getValue());
                } else {
                    return Sex.Unknown;
                }
            }
        }
    }

    private void fillContexts(final SortingCollection<VariantContext> contexts, final InfiniumGTCFile gtcFile,
                              final ExtendedIlluminaManifest manifest, final InfiniumEGTFile egtFile) {
        final ProgressLogger progress = new ProgressLogger(log, 100000, "sorted");

        final Iterator<ExtendedIlluminaManifestRecord> iterator = manifest.extendedIterator();
        int gtcIndex = 0;

        int numVariantsWritten = 0;

        while (iterator.hasNext()) {
            final ExtendedIlluminaManifestRecord record = iterator.next();

            if (record.isBad()) {
                continue;
            }
            // If the record is not flagged as errant in the manifest we include it in the VCF
            Allele A = record.getAlleleA();
            Allele B = record.getAlleleB();
            Allele ref = record.getRefAllele();

            //if A, B and ref are all . then parse the alleles from the probeSeq
            if (A.isNoCall() && B.isNoCall() && ref.isNoCall()) {
                final String sourceSeq = record.getSourceSeq();
                final String alleles = sourceSeq.substring(sourceSeq.indexOf('['), sourceSeq.lastIndexOf(']') + 1);
                if (alleles.charAt(1) == '-') {
                    A = Allele.SPAN_DEL;
                    B = ref = Allele.create(alleles.substring(alleles.indexOf("/") + 1, alleles.lastIndexOf(']')), true);
                } else {
                    A = ref = Allele.create(alleles.substring(alleles.indexOf("[") + 1, alleles.lastIndexOf('/')), true);
                    B = Allele.create(alleles.substring(alleles.indexOf("/") + 1, alleles.lastIndexOf(']')));
                }
            }

            final String chr = record.getB37Chr();
            final Integer position = record.getB37Pos();
            final Integer endPosition = position + ref.length() - 1;

            progress.record(chr, position);

            // Create list of unique alleles
            final List<Allele> assayAlleles = new ArrayList<>();
            assayAlleles.add(ref);

            if (!ref.equals(A, true)) {
                assayAlleles.add(A);
            }

            if (!ref.equals(B, true)) {
                assayAlleles.add(B);
            }

            final Genotype genotype = getGenotype(gtcFile, record, gtcIndex, A, B);

            final VariantContextBuilder builder = new VariantContextBuilder();

            builder.source(record.getName());
            builder.chr(chr);
            builder.start(position);
            builder.stop(endPosition);
            builder.alleles(assayAlleles);
            builder.log10PError(VariantContext.NO_LOG10_PERROR);
            builder.id(record.getName());
            builder.genotypes(genotype);

            VariantContextUtils.calculateChromosomeCounts(builder, false);

            //custom info fields
            builder.attribute(InfiniumVcfFields.ALLELE_A, record.getAlleleA());
            builder.attribute(InfiniumVcfFields.ALLELE_B, record.getAlleleB());
            builder.attribute(InfiniumVcfFields.ILLUMINA_STRAND, record.getIlmnStrand());
            builder.attribute(InfiniumVcfFields.PROBE_A, record.getAlleleAProbeSeq());
            builder.attribute(InfiniumVcfFields.PROBE_B, record.getAlleleBProbeSeq());
            builder.attribute(InfiniumVcfFields.BEADSET_ID, record.getBeadSetId());
            builder.attribute(InfiniumVcfFields.ILLUMINA_CHR, record.getChr());
            builder.attribute(InfiniumVcfFields.ILLUMINA_POS, record.getPosition());
            builder.attribute(InfiniumVcfFields.ILLUMINA_BUILD, record.getGenomeBuild());
            builder.attribute(InfiniumVcfFields.SOURCE, record.getSource().replace(' ', '_'));
            builder.attribute(InfiniumVcfFields.GC_SCORE, formatFloatForVcf(egtFile.totalScore[gtcIndex]));
            builder.attribute(InfiniumVcfFields.N_AA, egtFile.nAA[gtcIndex]);
            builder.attribute(InfiniumVcfFields.N_AB, egtFile.nAB[gtcIndex]);
            builder.attribute(InfiniumVcfFields.N_BB, egtFile.nBB[gtcIndex]);
            builder.attribute(InfiniumVcfFields.DEV_R_AA, formatFloatForVcf(egtFile.devRAA[gtcIndex]));
            builder.attribute(InfiniumVcfFields.DEV_R_AB, formatFloatForVcf(egtFile.devRAB[gtcIndex]));
            builder.attribute(InfiniumVcfFields.DEV_R_BB, formatFloatForVcf(egtFile.devRBB[gtcIndex]));
            builder.attribute(InfiniumVcfFields.MEAN_R_AA, formatFloatForVcf(egtFile.meanRAA[gtcIndex]));
            builder.attribute(InfiniumVcfFields.MEAN_R_AB, formatFloatForVcf(egtFile.meanRAB[gtcIndex]));
            builder.attribute(InfiniumVcfFields.MEAN_R_BB, formatFloatForVcf(egtFile.meanRBB[gtcIndex]));
            builder.attribute(InfiniumVcfFields.DEV_THETA_AA, formatFloatForVcf(egtFile.devThetaAA[gtcIndex]));
            builder.attribute(InfiniumVcfFields.DEV_THETA_AB, formatFloatForVcf(egtFile.devThetaAB[gtcIndex]));
            builder.attribute(InfiniumVcfFields.DEV_THETA_BB, formatFloatForVcf(egtFile.devThetaBB[gtcIndex]));
            builder.attribute(InfiniumVcfFields.MEAN_THETA_AA, formatFloatForVcf(egtFile.meanThetaAA[gtcIndex]));
            builder.attribute(InfiniumVcfFields.MEAN_THETA_AB, formatFloatForVcf(egtFile.meanThetaAB[gtcIndex]));
            builder.attribute(InfiniumVcfFields.MEAN_THETA_BB, formatFloatForVcf(egtFile.meanThetaBB[gtcIndex]));

            EuclideanValues aaVals = polarToEuclidean(egtFile.meanRAA[gtcIndex], egtFile.devRAA[gtcIndex],
                    egtFile.meanThetaAA[gtcIndex], egtFile.devThetaAA[gtcIndex]);
            EuclideanValues abVals = polarToEuclidean(egtFile.meanRAB[gtcIndex], egtFile.devRAB[gtcIndex],
                    egtFile.meanThetaAB[gtcIndex], egtFile.devThetaAB[gtcIndex]);
            EuclideanValues bbVals = polarToEuclidean(egtFile.meanRBB[gtcIndex], egtFile.devRBB[gtcIndex],
                    egtFile.meanThetaBB[gtcIndex], egtFile.devThetaBB[gtcIndex]);

            builder.attribute(InfiniumVcfFields.DEV_X_AA, formatFloatForVcf(aaVals.devX));
            builder.attribute(InfiniumVcfFields.DEV_X_AB, formatFloatForVcf(abVals.devX));
            builder.attribute(InfiniumVcfFields.DEV_X_BB, formatFloatForVcf(bbVals.devX));
            builder.attribute(InfiniumVcfFields.MEAN_X_AA, formatFloatForVcf(aaVals.meanX));
            builder.attribute(InfiniumVcfFields.MEAN_X_AB, formatFloatForVcf(abVals.meanX));
            builder.attribute(InfiniumVcfFields.MEAN_X_BB, formatFloatForVcf(bbVals.meanX));
            builder.attribute(InfiniumVcfFields.DEV_Y_AA, formatFloatForVcf(aaVals.devY));
            builder.attribute(InfiniumVcfFields.DEV_Y_AB, formatFloatForVcf(abVals.devY));
            builder.attribute(InfiniumVcfFields.DEV_Y_BB, formatFloatForVcf(bbVals.devY));
            builder.attribute(InfiniumVcfFields.MEAN_Y_AA, formatFloatForVcf(aaVals.meanY));
            builder.attribute(InfiniumVcfFields.MEAN_Y_AB, formatFloatForVcf(abVals.meanY));
            builder.attribute(InfiniumVcfFields.MEAN_Y_BB, formatFloatForVcf(bbVals.meanY));
            if (zCallThresholds.containsKey(egtFile.rsNames[gtcIndex])) {
                String[] zThresh = zCallThresholds.get(egtFile.rsNames[gtcIndex]);
                builder.attribute(InfiniumVcfFields.ZTHRESH_X, zThresh[0]);
                builder.attribute(InfiniumVcfFields.ZTHRESH_Y, zThresh[1]);
            }
            final String rsid = record.getRsId();
            if (StringUtils.isNotEmpty(rsid)) {
                builder.attribute(InfiniumVcfFields.RS_ID, rsid);
            }

            if (record.isDupe()) {
                builder.filter(InfiniumVcfFields.DUPE);
            }

            numVariantsWritten++;
            contexts.add(builder.make());


            gtcIndex++;
        }

        log.info(numVariantsWritten + " Variants were written to file");
        log.info(gtcFile.getNumberOfSnps() + " SNPs in the GTC file");
        log.info(manifest.getNumAssays() + " Variants on the " + manifest.getDescriptorFileName() + " genotyping array manifest file");
    }

    //Uses Manhattan distance conversion
    private EuclideanValues polarToEuclidean(float r, float rDeviation, float theta, float thetaDeviation) {
        //calculate variance (deviation^2)
        final double thetaVariance = Math.pow(thetaDeviation, 2.0);
        final double rVariance = Math.pow(rDeviation, 2.0);

        final double halfPi = Math.PI / 2;
        final double normalizedTheta = halfPi * theta;
        final double rOverX = (1 + Math.tan(normalizedTheta));

        //calculate X and Y variances from R and Theta variances
        final double thetaVarianceFactorX = -1 * (halfPi * r) * Math.pow((1 + Math.tan(normalizedTheta)), -2) * (1 / Math.pow(Math.cos(normalizedTheta), 2));
        final double rVarianceFactorX = 1 / (1 + Math.tan(normalizedTheta));
        final double varianceX = (Math.pow(thetaVarianceFactorX, 2) * thetaVariance) + (Math.pow(rVarianceFactorX, 2) * rVariance);
        final double thetaVarianceFactorY = -1 * thetaVarianceFactorX;
        final double rVarianceFactorY = 1 - rVarianceFactorX;
        final double varianceY = (Math.pow(thetaVarianceFactorY, 2) * thetaVariance) + (Math.pow(rVarianceFactorY, 2) * rVariance);

        /*
            Theta quantifies the relative amount of signal measured by the A and B intensities, defined by the equation:
            1/(2pi) * arctan(1/XY). R is a measurement of the total intensity observed from the A and B signals, defined as: R = A+B
            Illumina uses Manhattan distance https://en.wikipedia.org/wiki/Taxicab_geometry which is why R is A+B and not sqrt(A^2 + B^2)
            So Theta = 1/(2pi) * arctan(X/Y) and R = X + Y
         */

        final double meanX = r / rOverX;
        final double meanY = r - meanX;
        final double devX = Math.pow(varianceX, 0.5);
        final double devY = Math.pow(varianceY, 0.5);

        return new EuclideanValues((float) meanX, (float) meanY, (float) devX, (float) devY);
    }

    private class EuclideanValues {
        private final float meanX, meanY, devX, devY;

        EuclideanValues(float meanX, float meanY, float devX, float devY) {
            this.meanX = meanX;
            this.meanY = meanY;
            this.devX = devX;
            this.devY = devY;
        }

    }

    private Genotype getGenotype(final InfiniumGTCFile infiniumGTCFile,
                                 final IlluminaManifestRecord record,
                                 final int index,
                                 final Allele A,
                                 final Allele B) {
        final byte[] genotypeByteArray = infiniumGTCFile.getGenotype(index);

        // The Sample Alleles
        final List<Allele> alleles;

        if (genotypeByteArray[0] == 'N' && genotypeByteArray[1] == 'C') alleles = NO_CALL_ALLELES;
        else if (genotypeByteArray[0] == 'A' && genotypeByteArray[1] == 'A') alleles = Arrays.asList(A, A);
        else if (genotypeByteArray[0] == 'A' && genotypeByteArray[1] == 'B') alleles = Arrays.asList(A, B);
        else if (genotypeByteArray[0] == 'B' && genotypeByteArray[1] == 'B') alleles = Arrays.asList(B, B);
        else {
            throw new PicardException("Unexpected genotype call [" + genotypeByteArray[0] + "]" +
                    "[" + genotypeByteArray[1] + "]" + " for SNP: " + record.getName());
        }

        final Map<String, Object> attributes = new HashMap<>();
        attributes.put(InfiniumVcfFields.IGC, formatFloatForVcf(infiniumGTCFile.getGenotypeScore(index)));
        attributes.put(InfiniumVcfFields.X, infiniumGTCFile.getRawXIntensity(index));
        attributes.put(InfiniumVcfFields.Y, infiniumGTCFile.getRawYIntensity(index));
        attributes.put(InfiniumVcfFields.NORMX, formatFloatForVcf(infiniumGTCFile.getNormalizedXIntensity(index)));
        attributes.put(InfiniumVcfFields.NORMY, formatFloatForVcf(infiniumGTCFile.getNormalizedYIntensity(index)));
        attributes.put(InfiniumVcfFields.R, formatFloatForVcf(infiniumGTCFile.getRIlmn(index)));
        attributes.put(InfiniumVcfFields.THETA, formatFloatForVcf(infiniumGTCFile.getThetaIlmn(index)));
        attributes.put(InfiniumVcfFields.BAF, formatFloatForVcf(infiniumGTCFile.getBAlleleFreq(index)));
        attributes.put(InfiniumVcfFields.LRR, formatFloatForVcf(infiniumGTCFile.getLogRratiosIlmn(index)));

        final String sampleName = FilenameUtils.removeExtension(INPUT.getName());

        return GenotypeBuilder.create(sampleName, alleles, attributes);
    }

    private String formatFloatForVcf(final float value) {
        if (Float.isNaN(value)) {
            return dot;
        }
        return df.format(value);
    }

    /**
     * Writes out the VariantContext objects in the order presented to the supplied output file
     * in VCF format.
     */
    private void writeVcf(final SortingCollection<VariantContext> variants,
                          final File output,
                          final SAMSequenceDictionary dict,
                          final VCFHeader vcfHeader) {

        try (final VariantContextWriter writer = new VariantContextWriterBuilder()
                .setOutputFile(output)
                .setReferenceDictionary(dict)
                .setOptions(VariantContextWriterBuilder.DEFAULT_OPTIONS)
                .build()) {

            writer.writeHeader(vcfHeader);

            for (final VariantContext variant : variants) {
                if (variant.getAlternateAlleles().size() > 1) {
                    variant.getCommonInfo().addFilter(InfiniumVcfFields.TRIALLELIC);
                }

                writer.add(variant);
            }
        }
    }

    private VCFHeader createVCFHeader(final ExtendedIlluminaManifest manifest,
                                      final InfiniumGTCFile gtcFile,
                                      final String gtcGender,
                                      final File clusterFile,
                                      final File reference,
                                      final SAMSequenceDictionary dict) {
        final String inputName = INPUT.getName();
        final String chipWellBarcode = inputName.substring(0, inputName.lastIndexOf('.'));

        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("fileDate", new Date().toString()));
        lines.add(new VCFHeaderLine("source", "BPM file"));
        final String descriptorFileName = manifest.getDescriptorFileName();
        lines.add(new VCFHeaderLine(InfiniumVcfFields.ARRAY_TYPE, descriptorFileName.substring(0, descriptorFileName.lastIndexOf(dot))));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.EXTENDED_ILLUMINA_MANIFEST_FILE, EXTENDED_ILLUMINA_MANIFEST.getName()));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.EXTENDED_ILLUMINA_MANIFEST_VERSION, manifest.getExtendedManifestVersion()));

        lines.add(new VCFHeaderLine(InfiniumVcfFields.CHIP_WELL_BARCODE, chipWellBarcode));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.ANALYSIS_VERSION_NUMBER, ANALYSIS_VERSION_NUMBER.toString()));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.SAMPLE_ALIAS, SAMPLE_ALIAS));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.EXPECTED_GENDER, EXPECTED_GENDER));
        //add control codes
        final int measurementCount = gtcFile.getRawControlXIntensities().length / ArraysControlInfo.CONTROL_INFO.length;
        for (int i = 0; i < ArraysControlInfo.CONTROL_INFO.length; i++) {
            final int offset = i * measurementCount;
            ArraysControlInfo controlInfo = ArraysControlInfo.CONTROL_INFO[i];
            final int redIntensity = gtcFile.getRawControlXIntensity(offset);
            final int greenIntensity = gtcFile.getRawControlYIntensity(offset);
            lines.add(new VCFHeaderLine(controlInfo.getControl(), controlInfo.toString() + "|" + redIntensity + "|" + greenIntensity));
        }
        lines.add(new VCFHeaderLine(InfiniumVcfFields.FINGERPRINT_GENDER, fingerprintGender.name()));
        if (gtcGender != null) {
            lines.add(new VCFHeaderLine(InfiniumVcfFields.AUTOCALL_GENDER, gtcGender));
        } else {
            lines.add(new VCFHeaderLine(InfiniumVcfFields.AUTOCALL_GENDER, gtcFile.getGender()));
        }
        lines.add(new VCFHeaderLine(InfiniumVcfFields.AUTOCALL_DATE, gtcFile.getAutoCallDate()));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.IMAGING_DATE, gtcFile.getImagingDate()));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.CLUSTER_FILE, clusterFile.getName()));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.MANIFEST_FILE, descriptorFileName));
        lines.add(new VCFHeaderLine("content", manifest.getManifestFile().getName()));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.AUTOCALL_VERSION, gtcFile.getAutoCallVersion()));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.ZCALL_VERSION, "1.0.0.0"));
        if (ZCALL_THRESHOLDS_FILE != null)
            lines.add(new VCFHeaderLine(InfiniumVcfFields.ZCALL_THRESHOLDS, ZCALL_THRESHOLDS_FILE.getName()));
        lines.add(new VCFHeaderLine("reference", reference.getAbsolutePath()));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.GENOME_BUILD, "HG19"));
        lines.add(new VCFHeaderLine("picardVersion", this.getVersion()));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.P_95_RED, String.valueOf(gtcFile.getP95Red())));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.P_95_GREEN, String.valueOf(gtcFile.getP95Green())));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.SCANNER_NAME, gtcFile.getScannerName()));

        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
        lines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
        lines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
        lines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));

        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.IGC, 1, VCFHeaderLineType.Float, "Illumina GenCall Confidence Score"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.X, 1, VCFHeaderLineType.Integer, "Raw X intensity"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.Y, 1, VCFHeaderLineType.Integer, "Raw Y intensity"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.NORMX, 1, VCFHeaderLineType.Float, "Normalized X intensity"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.NORMY, 1, VCFHeaderLineType.Float, "Normalized Y intensity"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.R, 1, VCFHeaderLineType.Float, "Normalized R value"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.THETA, 1, VCFHeaderLineType.Float, "Normalized Theta value"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.BAF, 1, VCFHeaderLineType.Float, "B Allele Frequency"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.LRR, 1, VCFHeaderLineType.Float, "Log R Ratio"));

        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.ALLELE_A, 1, VCFHeaderLineType.String, "A allele"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.ALLELE_B, 1, VCFHeaderLineType.String, "B allele"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.ILLUMINA_STRAND, 1, VCFHeaderLineType.String, "Probe strand"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.PROBE_A, 1, VCFHeaderLineType.String, "Probe base pair sequence"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.PROBE_B, 1, VCFHeaderLineType.String, "Probe base pair sequence; not missing for strand-ambiguous SNPs"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.BEADSET_ID, 1, VCFHeaderLineType.Integer, "Bead set ID for normalization"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.ILLUMINA_CHR, 1, VCFHeaderLineType.String, "Chromosome in Illumina manifest"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.ILLUMINA_POS, 1, VCFHeaderLineType.Integer, "Position in Illumina manifest"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.ILLUMINA_BUILD, 1, VCFHeaderLineType.String, "Genome Build in Illumina manifest"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.SOURCE, 1, VCFHeaderLineType.String, "Probe source"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.GC_SCORE, 1, VCFHeaderLineType.Float, "Gentrain Score"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.N_AA, 1, VCFHeaderLineType.Integer, "Number of AA calls in training set"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.N_AB, 1, VCFHeaderLineType.Integer, "Number of AB calls in training set"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.N_BB, 1, VCFHeaderLineType.Integer, "Number of BB calls in training set"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_R_AA, 1, VCFHeaderLineType.Float, "Standard deviation of normalized R for AA cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_R_AB, 1, VCFHeaderLineType.Float, "Standard deviation of normalized R for AB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_R_BB, 1, VCFHeaderLineType.Float, "Standard deviation of normalized R for BB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_THETA_AA, 1, VCFHeaderLineType.Float, "Standard deviation of normalized THETA for AA cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_THETA_AB, 1, VCFHeaderLineType.Float, "Standard deviation of normalized THETA for AB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_THETA_BB, 1, VCFHeaderLineType.Float, "Standard deviation of normalized THETA for BB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_X_AA, 1, VCFHeaderLineType.Float, "Standard deviation of normalized X for AA cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_X_AB, 1, VCFHeaderLineType.Float, "Standard deviation of normalized X for AB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_X_BB, 1, VCFHeaderLineType.Float, "Standard deviation of normalized X for BB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_Y_AA, 1, VCFHeaderLineType.Float, "Standard deviation of normalized Y for AA cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_Y_AB, 1, VCFHeaderLineType.Float, "Standard deviation of normalized Y for AB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_Y_BB, 1, VCFHeaderLineType.Float, "Standard deviation of normalized Y for BB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_R_AA, 1, VCFHeaderLineType.Float, "Mean of normalized R for AA cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_R_AB, 1, VCFHeaderLineType.Float, "Mean of normalized R for AB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_R_BB, 1, VCFHeaderLineType.Float, "Mean of normalized R for BB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_THETA_AA, 1, VCFHeaderLineType.Float, "Mean of normalized THETA for AA cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_THETA_AB, 1, VCFHeaderLineType.Float, "Mean of normalized THETA for AB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_THETA_BB, 1, VCFHeaderLineType.Float, "Mean of normalized THETA for BB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_X_AA, 1, VCFHeaderLineType.Float, "Mean of normalized X for AA cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_X_AB, 1, VCFHeaderLineType.Float, "Mean of normalized X for AB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_X_BB, 1, VCFHeaderLineType.Float, "Mean of normalized X for BB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_Y_AA, 1, VCFHeaderLineType.Float, "Mean of normalized Y for AA cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_Y_AB, 1, VCFHeaderLineType.Float, "Mean of normalized Y for AB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_Y_BB, 1, VCFHeaderLineType.Float, "Mean of normalized Y for BB cluster"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.ZTHRESH_X, 1, VCFHeaderLineType.Float, "zCall X threshold"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.ZTHRESH_Y, 1, VCFHeaderLineType.Float, "zCall Y threshold"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.RS_ID, 1, VCFHeaderLineType.String, "dbSNP rs ID"));

        lines.add(new VCFFilterHeaderLine(InfiniumVcfFields.DUPE, "Duplicate assays position."));
        lines.add(new VCFFilterHeaderLine(InfiniumVcfFields.TRIALLELIC, "Tri-allelic assay."));
        lines.add(new VCFFilterHeaderLine(InfiniumVcfFields.FAIL_REF, "Assay failed to map to reference."));

        final VCFHeader header = new VCFHeader(lines, Collections.singletonList(chipWellBarcode));
        header.setSequenceDictionary(dict);
        return header;
    }
}