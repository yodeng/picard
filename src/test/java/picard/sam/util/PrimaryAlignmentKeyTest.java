package picard.sam.util;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Tuple;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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

public class PrimaryAlignmentKeyTest {

    final SAMFileHeader samHeader = new SAMFileHeader();

    @DataProvider(name="positiveTestData")
    public Object[][] positiveTestData() {
        return new Object[][]{
                // unpaired < first
                { makeRecordPair("r1", 0,
                        "r1", SAMFlag.READ_PAIRED.intValue() | SAMFlag.FIRST_OF_PAIR.intValue()), -1 },
                // unpaired < second
                { makeRecordPair("r1", 0,
                        "r1", SAMFlag.READ_PAIRED.intValue() | SAMFlag.SECOND_OF_PAIR.intValue()), -2 },
                // first < second
                { makeRecordPair("r1", SAMFlag.READ_PAIRED.intValue() | SAMFlag.FIRST_OF_PAIR.intValue(),
                        "r1", SAMFlag.READ_PAIRED.intValue() | SAMFlag.SECOND_OF_PAIR.intValue()), -1 },

                // first > unpaired
                { makeRecordPair("r1", SAMFlag.READ_PAIRED.intValue() | SAMFlag.FIRST_OF_PAIR.intValue(),
                        "r1", 0), 1 },
                // second > unpaired
                { makeRecordPair("r1", SAMFlag.READ_PAIRED.intValue() | SAMFlag.SECOND_OF_PAIR.intValue(),
                        "r1", 0), 2 },
                // second > first
                { makeRecordPair("r1", SAMFlag.READ_PAIRED.intValue() | SAMFlag.SECOND_OF_PAIR.intValue(),
                        "r1", SAMFlag.READ_PAIRED.intValue() | SAMFlag.FIRST_OF_PAIR.intValue()), 1 },

                // unpaired == unpaired
                { makeRecordPair("r1", 0, "r1", 0), 0 },
                // first == first
                { makeRecordPair("r1", SAMFlag.READ_PAIRED.intValue() | SAMFlag.FIRST_OF_PAIR.intValue(),
                        "r1", SAMFlag.READ_PAIRED.intValue() | SAMFlag.FIRST_OF_PAIR.intValue()), 0 },
                // second == second
                { makeRecordPair("r1", SAMFlag.READ_PAIRED.intValue() | SAMFlag.SECOND_OF_PAIR.intValue(),
                        "r1", SAMFlag.READ_PAIRED.intValue() | SAMFlag.SECOND_OF_PAIR.intValue()), 0 },

                // different read names, "r1" < "r2"
                { makeRecordPair("r1", 0, "r2", 0), -1 }
        };
    }

    @Test(dataProvider="positiveTestData")
    public void testPrimaryAlignmentKeyCompare(
            Tuple<SAMRecord, SAMRecord> recPair,
            final int expectedCompare) {
        Assert.assertEquals(new PrimaryAlignmentKey(recPair.a).compareTo(new PrimaryAlignmentKey(recPair.b)), expectedCompare);
    }

    @DataProvider(name="negativeTestData")
    public Object[][] negativeTestData() {
        return new Object[][]{
                { makeSAMRecord("r1", SAMFlag.NOT_PRIMARY_ALIGNMENT.intValue()) },  // secondary!
                { makeSAMRecord("r1", SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue()) }, //supplementary
        };
    }

    // reject secondary and supplementary
    @Test(dataProvider="negativeTestData", expectedExceptions = IllegalArgumentException.class)
    public void testRejectNonPrimary(final SAMRecord rec) {
        new PrimaryAlignmentKey(rec);
    }

    private Tuple<SAMRecord, SAMRecord> makeRecordPair(
            final String firstReadName,
            final int firstReadFlags,
            final String secondReadName,
            final int secondReadFlags) {
        return new Tuple(
                makeSAMRecord(firstReadName, firstReadFlags),
                makeSAMRecord(secondReadName, secondReadFlags)
        );
    }

    private SAMRecord makeSAMRecord(final String readName, final int readFlags)
    {
        final SAMRecord rec = new SAMRecord(samHeader);
        rec.setReadName(readName);
        rec.setFlags(readFlags);
        return rec;
    }

}
