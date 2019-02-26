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

import picard.PicardException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.Set;
import java.util.TreeSet;

public class InfiniumNormalizationManifest {
    private long[] positions;
    private byte[] chromosomes;
    private int[] normIds;
    private Integer[] allNormIds;

    public InfiniumNormalizationManifest(File illuminaNormalizationManifest) {
        parse(illuminaNormalizationManifest);
    }

    private void parse(File illuminaNormalizationManifest) {
        try (LineNumberReader lnr = new LineNumberReader(new FileReader(illuminaNormalizationManifest))) {
            Set<Integer> allNormIdSet = new TreeSet<>();
            lnr.skip(Long.MAX_VALUE);
            int numberOfSnps = lnr.getLineNumber() - 1;

            normIds = new int[numberOfSnps];
            positions = new long[numberOfSnps];
            chromosomes = new byte[numberOfSnps];
            try (BufferedReader reader = new BufferedReader(new FileReader(illuminaNormalizationManifest))) {
                boolean headerRead = false;
                int count = 0;
                while (reader.ready()) {
                    if (!headerRead) {
                        reader.readLine();
                        headerRead = true;
                    }

                    String line = reader.readLine();
                    String[] tokens = line.split(",");

                    allNormIdSet.add(new Integer(tokens[8].trim()));
                    normIds[count] = new Integer(tokens[8].trim());
                    String chrom = tokens[2].trim();
                    Byte chromByte = convertChromosomeToByte(chrom);
                    chromosomes[count] = chromByte;
                    positions[count] = new Long(tokens[3].trim());
                    count++;
                }
                allNormIds = allNormIdSet.toArray(new Integer[allNormIdSet.size()]);
                reader.close();
            }
        } catch (IOException e) {
            throw new PicardException("Error parsing Infinium normalization manifest", e);
        }
    }

    Integer[] getAllNormIds() {
        return allNormIds;
    }

    int[] getNormIds() {
        return normIds;
    }

    public byte[] getChromosomes() {
        return chromosomes;
    }

    public void setChromosomes(byte[] chromosomes) {
        this.chromosomes = chromosomes;
    }

    public long[] getPositions() {
        return positions;
    }

    public void setPositions(long[] positions) {
        this.positions = positions;
    }

    private Byte convertChromosomeToByte(String chrom) {
        switch (chrom) {
            case "X":
                chrom = "23";
                break;
            case "Y":
                chrom = "24";
                break;
            case "XY":
                chrom = "23";
                break;
            case "MT":
                chrom = "25";
                break;
        }
        return new Byte(chrom);
    }
}

