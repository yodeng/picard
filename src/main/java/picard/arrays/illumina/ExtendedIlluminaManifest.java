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

import org.apache.commons.lang.ArrayUtils;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

/**
 * The Extended Illumina Manifest has additional BUILD37 columns.
 *
 * Reads the header, stores the contents, and then provides an iterator to allow
 * access to the ManifestRecords (the assay records - that's all for now).
 */
public class ExtendedIlluminaManifest extends IlluminaManifest {

    private static final String BUILD37_CHR_HEADER_NAME = "build37Chr";
    private static final String BUILD37_POS_HEADER_NAME = "build37Pos";
    private static final String BUILD37_REF_ALLELE_HEADER_NAME = "build37RefAllele";
    private static final String BUILD37_ALLELE_A_HEADER_NAME = "build37AlleleA";
    private static final String BUILD37_ALLELE_B_HEADER_NAME = "build37AlleleB";
    private static final String BUILD37_RSID_NAME = "build37Rsid";
    private static final String BUILD37_FLAG_HEADER_NAME = "build37Flag";

    static final String[] EXTENDED_MANIFEST_HEADERS = {
            BUILD37_CHR_HEADER_NAME,
            BUILD37_POS_HEADER_NAME,
            BUILD37_REF_ALLELE_HEADER_NAME,
            BUILD37_ALLELE_A_HEADER_NAME,
            BUILD37_ALLELE_B_HEADER_NAME,
            BUILD37_RSID_NAME,
            BUILD37_FLAG_HEADER_NAME
    };

    static final String EXTENDED_MANIFEST_VERSION_HEADER_NAME = "CreateExtendedIlluminaManifest.version";
    static final String EXTENDED_MANIFEST_TARGET_BUILD_HEADER_NAME = "Target Build";
    static final String EXTENDED_MANIFEST_TARGET_REFERENCE_HEADER_NAME = "Target Reference File";
    static final String EXTENDED_MANIFEST_CLUSTER_FILE_HEADER_NAME = "Cluster File";
    static final String EXTENDED_MANIFEST_DBSNP_FILE_HEADER_NAME = "dbSNP File";
    static final String EXTENDED_MANIFEST_SUPPORTED_BUILD_HEADER_NAME = "Supported Build";
    static final String EXTENDED_MANIFEST_SUPPORTED_REFERENCE_HEADER_NAME = "Supported Reference File";
    static final String EXTENDED_MANIFEST_SUPPORTED_CHAIN_FILE_HEADER_NAME = "Supported Chain File";

    @Override
    public String[] getAllPossibleHeaderNames() {
        return (String[])ArrayUtils.addAll(super.getAllPossibleHeaderNames(), EXTENDED_MANIFEST_HEADERS);
    }

    public ExtendedIlluminaManifest(final File manifestFile) throws IOException {
        super(manifestFile);
    }

    public Iterator<ExtendedIlluminaManifestRecord> extendedIterator() {

        return new Iterator<ExtendedIlluminaManifestRecord>() {
            private int assayCount = 0;
            private int numAssays = getNumAssays();

            public boolean hasNext() {
                return (assayCount < numAssays) && (manifestFileParser.hasNext()) ;
            }

            public ExtendedIlluminaManifestRecord next() {
                assayCount++;
                return new ExtendedIlluminaManifestRecord(getAssayHeaderNameToIndex(), manifestFileParser.next(), assayCount - 1);
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    public String getExtendedManifestVersion() {
        String version = "?";
        for (String[] headerLine: getHeaderContents()) {
            if (headerLine[0].equals(EXTENDED_MANIFEST_VERSION_HEADER_NAME)) {
                version = headerLine[1];
            }
        }
        return version;
    }
}
