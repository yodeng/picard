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

package picard.arrays.illumina;


import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang.StringUtils;
import picard.PicardException;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A class to represent a record (line) from a Extended Illumina Manifest block
 */
public class ExtendedIlluminaManifestRecord extends IlluminaManifestRecord {
    protected enum Flag {
        ILLUMINA_FLAGGED,   // Illumina flagged
        LIFTOVER_FAILED,
        UNSUPPORTED_GENOME_BUILD,
        SEQUENCE_MISMATCH,      // mismatch between what the manifest claims is the reference vs. the actual reference.
        INDEL_SEQ_MISMATCH,     // Unable to reconcile indel situation
        INDEL_EXTENSION_ERROR,  // extension base conflict.
        DUPE,
        PASS,
    }

    private String b37Chr;
    private Integer b37Pos;
    private String snpRefAllele;
    private String snpAlleleA;
    private String snpAlleleB;
    private String rsId;
    private Flag flag = Flag.PASS;

    private Allele A;
    private Allele B;
    private Allele ref;

    private Strand calculatedStrand;

    private final Log log = Log.getInstance(ExtendedIlluminaManifestRecord.class);

    private static final String SRC_SEQ_REGEX = "([AGTCagtc]*)\\[([AGTCagtc-])\\/([AGTC]*)\\]([agtcAGTC]*)";
    // Symbolics for the regex groups...
    private static final int SRC_SEQ_SEQUENCE_BEFORE_VARIANT = 1;
    private static final int SRC_SEQ_FIRST_VARIANT_SEQUENCE  = 2;
    private static final int SRC_SEQ_SECOND_VARIANT_SEQUENCE = 3;
    private static final int SRC_SEQ_SEQUENCE_AFTER_VARIANT  = 4;


    private static String build36 = "36";
    private static String build37 = "37";
    public static final Pattern pattern = Pattern.compile(SRC_SEQ_REGEX);

    /**
     * This constructor is used to read records from an already created ExtendedIlluminaManifestRecord file.
     * It does not work to set the Extended-specific fields
     */
    ExtendedIlluminaManifestRecord(final Map<String, Integer> columnNameToIndex, final String[] line, final int index) {
        super(columnNameToIndex, line, index);

        final int end = line.length;
        flag = Flag.valueOf(line[end - 1]);

        if (!isBad()) {
            b37Chr = line[end - 7];
            b37Pos = parseIntOrNull(line[end - 6]);
            snpRefAllele = line[end - 5];
            snpAlleleA = line[end - 4];
            snpAlleleB = line[end - 3];
            rsId = line[end - 2];

            A = Allele.create(snpAlleleA, snpAlleleA.equals(snpRefAllele));
            B = Allele.create(snpAlleleB, snpAlleleB.equals(snpRefAllele));
            ref = Allele.create(snpRefAllele, true);
        } else {
            b37Chr = "0";
            b37Pos = 0;
            snpRefAllele = "";
            snpAlleleA = "";
            snpAlleleB = "";
            rsId = "";

            A = Allele.NO_CALL;
            B = Allele.NO_CALL;
            ref = Allele.NO_CALL;
        }
    }

    /**
     * This constructor is used to take a record from an Illumina Manifest and sets the Extended-specific fields
     * in preparation for writing out the ExtendedIlluminaManifestRecord to file (or otherwise using it)
     */
    ExtendedIlluminaManifestRecord(final IlluminaManifestRecord record,
                                   final Map<String, ReferenceSequenceFile> referenceFilesMap,
                                   final Map<String, File> chainFilesMap,
                                   final boolean dupe,
                                   final String passedRsId) {
        super(record);

        validate(record, dupe);

        if (!isBad()) {
            liftOverToBuild37(record, chainFilesMap);
        }

        final ReferenceSequenceFile refFile = referenceFilesMap.get(record.getMajorGenomeBuild());

        Strand strand = Strand.INVALID;
        if (!isBad()) {     // getStrand may flag a record as bad.
            strand = getStrand(record, refFile);
        }
        calculatedStrand = strand;

        if (!isBad()) {         // getStrand may flag a record as bad
            rsId = passedRsId == null ? "" : passedRsId;
            if (record.isSnp()) {
                populateSnpAlleles(record, refFile, strand);
            } else {
                populateIndelAlleles(record, refFile, strand);
            }
        }

        if (!isBad()) {     // Note that populateSnp/IndelAlleles may flag a record as bad
            A = Allele.create(snpAlleleA, snpAlleleA.equals(snpRefAllele));
            B = Allele.create(snpAlleleB, snpAlleleB.equals(snpRefAllele));
            ref = Allele.create(snpRefAllele, true);
        } else {
            A = Allele.NO_CALL;
            B = Allele.NO_CALL;
            ref = Allele.NO_CALL;
        }
    }

    public Allele getAlleleA() {
        return A;
    }

    public Allele getAlleleB() {
        return B;
    }

    public Allele getRefAllele() {
        return ref;
    }

    public Strand getCalculatedStrand() {
        return calculatedStrand;
    }

    public String getB37Chr() {
        return b37Chr;
    }

    public Integer getB37Pos() {
        return b37Pos;
    }

    public String getRsId() { return rsId; }

    public Boolean isBad() {
        return flag != Flag.DUPE && flag != Flag.PASS;
    }

    public Boolean isDupe() {
        return flag == Flag.DUPE;
    }

    public Flag getFlag() {
        return flag;
    }

    private void validate(final IlluminaManifestRecord r, boolean dupe) {
        //set dupe first so it can be overridden by fail flags
        if (dupe) flag = Flag.DUPE;

        // short circuit on non-existant values
        if (r.getChr().trim().equals("0")) {
            flag = Flag.ILLUMINA_FLAGGED;
        }

        if (!r.getMajorGenomeBuild().trim().equals(build36) && !r.getMajorGenomeBuild().trim().equals(build37)) {
            flag = Flag.UNSUPPORTED_GENOME_BUILD;
        }
    }

    /**
     * Determines the chromosome and position of the record on Build 37.
     */
    private void liftOverToBuild37(final IlluminaManifestRecord r, final Map<String, File> chainFilesMap) {

        // no liftover needed
        if (r.getMajorGenomeBuild().trim().equals(build37)) {
            b37Chr = r.getChr();
            b37Pos = r.getPosition();

        } else {
            final File chainFileToBuild37 = chainFilesMap.get(r.getMajorGenomeBuild());
            final LiftOver liftOver = new LiftOver(chainFileToBuild37);
            final Interval interval = new Interval(r.getChr(), r.getPosition(), r.getPosition());
            final Interval b37Interval = liftOver.liftOver(interval);

            log.info("Performing a liftover:");
            log.info("Old chr, pos: " + r.getChr() + ":" + r.getPosition());

            if (b37Interval != null) {
                b37Chr = b37Interval.getContig();
                b37Pos = b37Interval.getStart();

                log.info("New chr, pos: " + b37Chr + ":" + b37Pos);
            } else {
                flag = Flag.LIFTOVER_FAILED;
                log.warn("Could not perform liftover from build " + r.getGenomeBuild() + " " + r.getChr() + ":" + r.getPosition() + " to build37.");
            }
        }
    }

    /**
     * Uses source sequence to determine snpAlleleA, snpAlleleB, and snpRefAllele.
     * <p>
     * In this context, variant is either a SNP or an indel.
     */
    private void populateSnpAlleles(final IlluminaManifestRecord record, final ReferenceSequenceFile refFile, final Strand strand) {

        snpAlleleA = record.getSnp().substring(1, 2);
        snpAlleleB = record.getSnp().substring(3, 4);

        if (strand == Strand.NEGATIVE) {
            snpAlleleA = SequenceUtil.reverseComplement(snpAlleleA);
            snpAlleleB = SequenceUtil.reverseComplement(snpAlleleB);
        }

        //extra validation for ambiguous snps
        if (isAmbiguous()) {
            if (record.getAlleleBProbeSeq() != null) {
                String probeAAllele = record.getAlleleAProbeSeq().substring(record.getAlleleAProbeSeq().length() - 1);
                String probeBAllele = record.getAlleleBProbeSeq().substring(record.getAlleleBProbeSeq().length() - 1);
                if (!probeAAllele.equals(snpAlleleA) && !probeBAllele.equals(snpAlleleB) && (strand == Strand.POSITIVE)) {
                    snpAlleleA = probeAAllele;
                    snpAlleleB = probeBAllele;
                }
            } else {
                // This manifest contains no Allele B Probe Sequence.  We (currently) need this for validating/trusting
                // these ambiguous SNPs, so we are flagging it.
                log.warn("Ambiguous probe without alleleBProbeSeq!!!  Record at " + record.getChr() + ":" + record.getPosition());
            }
        }

        snpRefAllele = getSequenceAt(refFile, record.getChr(), record.getPosition(), record.getPosition());
    }

    /**
     * Uses source sequence to determine snpAlleleA, snpAlleleB, and snpRefAllele.
     * <p>
     * In this context, variant is either a SNP or an indel.
     */
    private void populateIndelAlleles(final IlluminaManifestRecord record, final ReferenceSequenceFile refFile, final Strand strand) {
        // Verify that alleleAProbeSeq is at the location it's supposed to be, and that it is found in the source sequence
        final Strand alleleAProbeSeqStrand = getAlleleAProbeSeqStrand(record, refFile);
        if (isBad()) {      // Couldn't find alleleA probe sequence in reference.
            return;
        }
        final Matcher matcher = parseSourceSeq(record.getSourceSeq());

        String sequenceBeforeInsertion;
        String insertion = matcher.group(SRC_SEQ_SECOND_VARIANT_SEQUENCE).toUpperCase();
        String sequenceAfterInsertion;

        if (strand == Strand.NEGATIVE) {
            sequenceBeforeInsertion = SequenceUtil.reverseComplement(matcher.group(SRC_SEQ_SEQUENCE_AFTER_VARIANT).toUpperCase());
            insertion = SequenceUtil.reverseComplement(insertion);
            sequenceAfterInsertion = SequenceUtil.reverseComplement(matcher.group(SRC_SEQ_SEQUENCE_BEFORE_VARIANT).toUpperCase());
        } else {
            sequenceBeforeInsertion = matcher.group(SRC_SEQ_SEQUENCE_BEFORE_VARIANT).toUpperCase();
            sequenceAfterInsertion = matcher.group(SRC_SEQ_SEQUENCE_AFTER_VARIANT).toUpperCase();
        }
        final String baseBeforeInsertion = sequenceBeforeInsertion.substring(sequenceBeforeInsertion.length() - 1);
        final String sequenceWithoutInsertion = sequenceBeforeInsertion + sequenceAfterInsertion;
        final String sequenceWithInsertion = sequenceBeforeInsertion + insertion + sequenceAfterInsertion;

        // A note on determining which is alleleA and which is alleleB.
        // align the alleleA probe sequence against the source sequence with and without the insertion.
        // Then look at the base immediately after these two alignments.
        // This is the single base extension probe.
        // If it is an A or T, then that one is alleleA
        // If it is a  C or G, then that one is alleleB
        // They must differ or that's an error.
        // from Illumina: "That extension base is the determining factor in which allele gets genotype A and which gets genotype B (A and T are genotype A, C and G are genotype B)."
        String alleleAProbeSeq = getAlleleAProbeSeq();
        if (alleleAProbeSeqStrand == Strand.NEGATIVE) {
            alleleAProbeSeq = SequenceUtil.reverseComplement(alleleAProbeSeq);
        }

        // Note.  Source Sequence can be long and probe sequence is short.  Have found several cases where probe seq is found multiple times in source sequence,
        //        We are interested in where probe sequence overlaps the middle of source sequence (where the insertion is).  So that's why we use the 'fromIndex' parameter of indexOf
        int fromIndex = sequenceBeforeInsertion.length() - alleleAProbeSeq.length() - 5;      // fudge factor of 5
        int locnWith = sequenceWithInsertion.indexOf(alleleAProbeSeq, fromIndex);
        int locnWithout = sequenceWithoutInsertion.indexOf(alleleAProbeSeq, fromIndex);
        if ((locnWith < 0) || (locnWithout < 0)) {
            flag = Flag.INDEL_EXTENSION_ERROR;
            log.warn("Error in indel processing.  Record at " + record.getChr() + ":" + record.getPosition() + " " + record.getSnp() + " src strand:" + strand + " probe strand:" + alleleAProbeSeqStrand + " " + flag.toString());
            log.warn("Problem find alleleAProbeSeq in source sequence with and without insertion");
            return;
        }
        String extensionBaseWithInsertion;
        String extensionBaseWithoutInsertion;
        if (alleleAProbeSeqStrand == Strand.POSITIVE) {
            extensionBaseWithInsertion = sequenceWithInsertion.substring(locnWith + alleleAProbeSeq.length(), locnWith + alleleAProbeSeq.length() + 1);
            extensionBaseWithoutInsertion = sequenceWithoutInsertion.substring(locnWithout + alleleAProbeSeq.length(), locnWithout + alleleAProbeSeq.length() + 1);
        } else {        // alleleA probe is on the negative strand
            extensionBaseWithInsertion = sequenceWithInsertion.substring(locnWith - 1, locnWith);
            extensionBaseWithoutInsertion = sequenceWithoutInsertion.substring(locnWithout - 1, locnWithout);
        }

        if (!validateAlleleABBases(extensionBaseWithInsertion, extensionBaseWithoutInsertion)) {
            flag = Flag.INDEL_EXTENSION_ERROR;
            log.warn("Error in indel processing.  Record at " + record.getChr() + ":" + record.getPosition() + " " + record.getSnp() + " src strand:" + strand + " probe strand:" + alleleAProbeSeqStrand + " " + flag.toString());
            log.warn("Extension base conflict: " + extensionBaseWithInsertion + ", " + extensionBaseWithoutInsertion);
            return;
        }
        // Now determine whether the reference with insertion matches sequence with insertion or reference without matches sequence without.
        // Use that to determine whether the insertion in the record exists in reference (in which case we are testing for a deletion)
        //                 of if the insertion in the record does not exist in reference (in which case we are testing for the insertion)
        final String referenceWithoutInsertion = getSequenceAt(refFile, record.getChr(),
                record.getPosition() - sequenceBeforeInsertion.length(),
                record.getPosition() + sequenceAfterInsertion.length() - 1);
        final String referenceWithInsertion = getSequenceAt(refFile, record.getChr(),
                record.getPosition() - sequenceBeforeInsertion.length(),
                record.getPosition() + insertion.length() + sequenceAfterInsertion.length() - 1);
        boolean insFound    = referenceWithInsertion.equals(sequenceWithInsertion);
        boolean insNotFound = referenceWithoutInsertion.equals(sequenceWithoutInsertion);
        if (insFound != insNotFound) {
            // One or the other condition
            // Whether the condition is an insertion or deletion is only reflected in the refAllele.
            snpRefAllele = insFound ? baseBeforeInsertion + insertion : baseBeforeInsertion;
            if (isAorT(extensionBaseWithInsertion)) {
                snpAlleleA = baseBeforeInsertion + insertion;
                snpAlleleB = baseBeforeInsertion;
            } else {        // extensionBaseWithInsertion is C or G
                snpAlleleA = baseBeforeInsertion;
                snpAlleleB = baseBeforeInsertion + insertion;
            }
        } else {
            flag = Flag.INDEL_SEQ_MISMATCH;
            log.warn("Error in indel processing.  Record at " + record.getChr() + ":" + record.getPosition() + " " + record.getSnp() + " src strand:" + strand + " probe strand:" + alleleAProbeSeqStrand + " " + flag.toString());
            log.warn("Cannot match reference and sequence with or without insertion (or matches both conditions...).");
            log.warn("Reference with Insertion:    " + referenceWithInsertion);
            log.warn("Sequence with Insertion:     " + sequenceWithInsertion);
            log.warn("Sequence without Insertion:  " + sequenceWithoutInsertion);
            log.warn("Reference without Insertion: " + referenceWithoutInsertion);
            return;
        }

        // correcting indel coordinates for VCF format.
        b37Pos = b37Pos - 1;
    }

    private boolean validateAlleleABBases(String allele1, String allele2) {
        return ((isAorT(allele1) && isCorG(allele2)) || (isCorG(allele1) && isAorT(allele2)));
    }

    private boolean isAorT(String allele) {
        final String localAllele = allele.toUpperCase();
        return (localAllele.equals("A") || localAllele.equals("T"));
    }

    private boolean isCorG(String allele) {
        final String localAllele = allele.toUpperCase();
        return (localAllele.equals("C") || localAllele.equals("G"));
    }


    private Strand getAlleleAProbeSeqStrand(final IlluminaManifestRecord record, final ReferenceSequenceFile refFile) {
        return getStrand(record, refFile, getAlleleAProbeSeq());
    }

    /**
     * Compares probeASeq to reference to determine strand.
     */
    private Strand getStrand(final IlluminaManifestRecord record, final ReferenceSequenceFile refFile) {
        String probeSeq;

        if (isIndel()) {
            Matcher matcher = parseSourceSeq(record.getSourceSeq());
            probeSeq = matcher.group(SRC_SEQ_SEQUENCE_BEFORE_VARIANT).toUpperCase();
        } else if (isAmbiguous()) {
            //ambiguous snps contain the probed base so we need to truncate the string
            probeSeq = record.getAlleleAProbeSeq().substring(0, record.getAlleleAProbeSeq().length() - 1);
        } else {
            probeSeq = record.getAlleleAProbeSeq();
        }
        return getStrand(record, refFile, probeSeq);
    }

    private Strand getStrand(final IlluminaManifestRecord record, final ReferenceSequenceFile refFile, final String probeSeq) {
        Strand strand = Strand.INVALID;

        final String reference = getSequenceAt(refFile, record.getChr(), record.getPosition() - probeSeq.length(), record.getPosition() - 1);
        final String reverseReference = SequenceUtil.reverseComplement(getSequenceAt(refFile, record.getChr(), record.getPosition() + 1, record.getPosition() + probeSeq.length()));

        if (reference.equals(probeSeq)) {
            strand = Strand.POSITIVE;
        } else if (reverseReference.equals(probeSeq)) {
            strand = Strand.NEGATIVE;
        } else {
            flag = Flag.SEQUENCE_MISMATCH;
            log.warn("Error in getStrand.  Record at " + record.getChr() + ":" + record.getPosition() + record.getSnp() + " strand:" + strand + " " + flag.toString());
            log.warn("Reference:        " + reference);
            log.warn("probeSeq:         " + probeSeq);
            log.warn("ReverseReference: " + reverseReference);
        }

        return strand;
    }

    /**
     * Find the sequence at a given position in a reference sequence file.
     */
    private static String getSequenceAt(final ReferenceSequenceFile refFile, final String chr, final int startPos, final int endPos) {
        final int contigLength = refFile.getSequenceDictionary().getSequence(chr).getSequenceLength();
        int usedEndPos = Math.min(endPos, contigLength);
        return new String(refFile.getSubsequenceAt(chr, startPos, usedEndPos).getBases()).toUpperCase();
    }

    /**
     * Use regex to capture the insertion sequence and the sequence after the indel.
     * <p>
     * The source sequence is of the from V[W/X]Y, where V,W,X,Y are sequences.
     * <p>
     * - A SNP example looks like:    AGGGAGTC[A/G]GGTTGCGA
     * V     W X    Y
     * <p>
     * - A InDel example looks like:  AGCCTCGA[-/CGAA]TCACC
     * V     W   X   Y
     */
    private static Matcher parseSourceSeq(final String sourceSeq) {
        final Matcher matcher = pattern.matcher(sourceSeq);
        if (matcher.find()) {
            return matcher;
        } else {
            throw new PicardException("Could not find the pattern V[W/X]Y in the SourceSeq: " + sourceSeq);
        }
    }

    @Override
    public String getLine() {
        final String originalLine = super.getLine();

        final List<String> extensions = new ArrayList<>();
        extensions.add(b37Chr);
        extensions.add(b37Pos != null ? b37Pos.toString() : null);
        extensions.add(snpRefAllele);
        extensions.add(snpAlleleA);
        extensions.add(snpAlleleB);
        extensions.add(rsId);
        extensions.add(flag.name());

        return originalLine + "," + StringUtils.join(extensions, ",");
    }

}