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

import htsjdk.samtools.util.Iso8601Date;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.text.ParseException;
import java.text.SimpleDateFormat;

public class InfiniumVcfFields {
    // Header Fields
    public static final String ARRAY_TYPE = "arrayType";

    public static final String EXTENDED_ILLUMINA_MANIFEST_VERSION = "extendedIlluminaManifestVersion";
    public static final String CHIP_WELL_BARCODE = "chipWellBarcode";
    public static final String ANALYSIS_VERSION_NUMBER = "analysisVersionNumber";
    public static final String SAMPLE_ALIAS = "sampleAlias";

    public static final String GENOME_BUILD = "genomeBuild";
    public static final String EXPECTED_GENDER = "expectedGender";
    public static final String FINGERPRINT_GENDER = "fingerprintGender";
    public static final String AUTOCALL_GENDER = "autocallGender";
    public static final String AUTOCALL_DATE = "autocallDate";
    public static final String IMAGING_DATE = "imagingDate";
    public static final String CLUSTER_FILE = "clusterFile";
    public static final String MANIFEST_FILE = "manifestFile";
    public static final String EXTENDED_ILLUMINA_MANIFEST_FILE = "extendedManifestFile";
    public static final String AUTOCALL_VERSION = "autocallVersion";
    public static final String ZCALL_VERSION = "zcallVersion";
    public static final String ZCALL_THRESHOLDS = "zcallThresholds";
    public static final String P_95_RED = "p95Red";
    public static final String P_95_GREEN = "p95Green";
    public static final String SCANNER_NAME = "scannerName";

    //FORMAT Fields
    public static final String X = "X";
    public static final String Y = "Y";
    public static final String NORMX = "NORMX";
    public static final String NORMY = "NORMY";
    public static final String R = "R";
    public static final String THETA = "THETA";
    public static final String BAF = "BAF";
    public static final String LRR = "LRR";
    public static final String IGC = "IGC";
    public static final String GTA = "GTA";
    public static final String GTZ = "GTZ";

    //INFO Fields
    public static final String ALLELE_A = "ALLELE_A";
    public static final String ALLELE_B = "ALLELE_B";
    public static final String ILLUMINA_STRAND = "ILLUMINA_STRAND";
    public static final String PROBE_A = "PROBE_A";
    public static final String PROBE_B = "PROBE_B";
    public static final String BEADSET_ID = "BEADSET_ID";
    public static final String ILLUMINA_CHR = "ILLUMINA_CHR";
    public static final String ILLUMINA_POS = "ILLUMINA_POS";
    public static final String ILLUMINA_BUILD = "ILLUMINA_BUILD";
    public static final String SOURCE = "SOURCE";
    public static final String GC_SCORE = "GC_SCORE";
    public static final String N_AA = "N_AA";
    public static final String N_AB = "N_AB";
    public static final String N_BB = "N_BB";
    public static final String DEV_R_AA = "devR_AA";
    public static final String DEV_R_AB = "devR_AB";
    public static final String DEV_R_BB = "devR_BB";
    public static final String MEAN_R_AA = "meanR_AA";
    public static final String MEAN_R_AB = "meanR_AB";
    public static final String MEAN_R_BB = "meanR_BB";
    public static final String DEV_THETA_AA = "devTHETA_AA";
    public static final String DEV_THETA_AB = "devTHETA_AB";
    public static final String DEV_THETA_BB = "devTHETA_BB";
    public static final String MEAN_THETA_AA = "meanTHETA_AA";
    public static final String MEAN_THETA_AB = "meanTHETA_AB";
    public static final String MEAN_THETA_BB = "meanTHETA_BB";
    public static final String DEV_X_AA = "devX_AA";
    public static final String DEV_X_AB = "devX_AB";
    public static final String DEV_X_BB = "devX_BB";
    public static final String MEAN_X_AA = "meanX_AA";
    public static final String MEAN_X_AB = "meanX_AB";
    public static final String MEAN_X_BB = "meanX_BB";
    public static final String DEV_Y_AA = "devY_AA";
    public static final String DEV_Y_AB = "devY_AB";
    public static final String DEV_Y_BB = "devY_BB";
    public static final String MEAN_Y_AA = "meanY_AA";
    public static final String MEAN_Y_AB = "meanY_AB";
    public static final String MEAN_Y_BB = "meanY_BB";
    public static final String ZTHRESH_X = "zthresh_X";
    public static final String ZTHRESH_Y = "zthresh_Y";
    public static final String RS_ID = "refSNP";

    //FILTER Fields
    public static final String TRIALLELIC = "TRIALLELIC";
    public static final String DUPE = "DUPE";
    public static final String FAIL_REF = "FAIL_REF";
    public static final String ZCALL_DIFF = "ZCALL_DIFF";

    public static String  getValueFromVcfOtherHeaderLine(final VCFHeader vcfHeader, final String keyName) {
        VCFHeaderLine otherHeaderLine = vcfHeader.getOtherHeaderLine(keyName);
        if (otherHeaderLine != null) {
            return otherHeaderLine.getValue();
        } else {
            throw new IllegalArgumentException("Input VCF file is missing header line of type '" + keyName + "'");
        }
    }

    public static Iso8601Date getDateFromVcfOtherHeaderLine(final VCFHeader vcfHeader, final String keyName, final SimpleDateFormat dateformat) {
        String dateString = InfiniumVcfFields.getValueFromVcfOtherHeaderLine(vcfHeader, keyName);
        try {
            return new Iso8601Date(dateformat.parse(dateString));
        } catch (ParseException pe) {
            throw new IllegalArgumentException("Unrecognized date for '" + keyName + "' in VCF header (" + dateString + ")");
        }
    }

    public static Integer getIntegerFromVcfOtherHeaderLine(final VCFHeader vcfHeader, final String keyName) {
        VCFHeaderLine otherHeaderLine = vcfHeader.getOtherHeaderLine(keyName);
        if (otherHeaderLine != null) {
            return Integer.valueOf(otherHeaderLine.getValue());
        } else {
            throw new IllegalArgumentException("Input VCF file is missing header line of type '" + keyName + "'");
        }
    }

    public static String getOptionalValueFromVcfOtherHeaderLine(final VCFHeader vcfHeader, final String keyName) {
        VCFHeaderLine otherHeaderLine = vcfHeader.getOtherHeaderLine(keyName);
        if (otherHeaderLine != null) {
            return otherHeaderLine.getValue();
        } else {
            return null;
        }
    }
}
