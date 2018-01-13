package picard.util;

import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import picard.PicardException;
import picard.sam.CreateSequenceDictionary;

import java.io.*;
import java.util.Iterator;

/**
 * Created by farjoun on 1/12/18.
 */
public class SequenceDictionaryUtils {


    /** Encodes a sequence dictionary
     *
     *
     * @param writer a Buffered writer into which the dictinoary will be written
     * @param samSequenceRecordIterator an iterator that produces SAMSequenceRecords
     *
     * @throws IllegalArgumentException if the iterator produces two SAMSequenceRecord with the same name
     */
    public static void encodeDictionary(final BufferedWriter writer, Iterator<SAMSequenceRecord> samSequenceRecordIterator) {
        final SAMSequenceDictionaryCodec samDictCodec = new SAMSequenceDictionaryCodec(writer);

        samDictCodec.encodeHeaderLine(false);
        // SortingCollection is used to check uniqueness of sequence names
        final SortingCollection<String> sequenceNames = makeSortingCollection();

        // read reference sequence one by one and write its metadata
        while (samSequenceRecordIterator.hasNext()) {
            final SAMSequenceRecord samSequenceRecord = samSequenceRecordIterator.next();
            samDictCodec.encodeSequenceRecord(samSequenceRecord);
            sequenceNames.add(samSequenceRecord.getSequenceName());
        }

        // check uniqueness of sequences names
        final CloseableIterator<String> iterator = sequenceNames.iterator();

        if(!iterator.hasNext()) return;

        String current = iterator.next();
        while (iterator.hasNext()) {
            final String next = iterator.next();
            if (current.equals(next)) {
                throw new PicardException("Sequence name " + current +
                        " appears more than once in reference file");
            }
            current = next;
        }
    }

    public static SortingCollection<String> makeSortingCollection() {
        final File tmpDir = IOUtil.createTempDir("SamDictionaryNames", null);
        tmpDir.deleteOnExit();
        // 256 byte for one name, and 1/10 part of all memory for this, rough estimate
        long maxNamesInRam = Runtime.getRuntime().maxMemory() / 256 / 10;
        return SortingCollection.newInstance(
                String.class,
                new StringCodec(),
                String::compareTo,
                (int) Math.min(maxNamesInRam, Integer.MAX_VALUE),
                tmpDir
        );
    }

    private static class StringCodec implements SortingCollection.Codec<String> {
        private DataInputStream dis;
        private DataOutputStream dos;

        public StringCodec clone() {
            return new StringCodec();
        }

        public void setOutputStream(final OutputStream os) {
            dos = new DataOutputStream(os);
        }

        public void setInputStream(final InputStream is) {
            dis = new DataInputStream(is);
        }

        public void encode(final String str) {
            try {
                dos.writeUTF(str);
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            }
        }

        public String decode() {
            try {
                return dis.readUTF();
            } catch (EOFException e) {
                return null;
            } catch (IOException e) {
                throw new PicardException("Exception reading sequence name from temporary file.", e);
            }
        }
    }
}
