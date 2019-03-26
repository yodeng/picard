package picard.pedigree;

import picard.PicardException;

import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Represents the sex of an individual.
 */
public enum Sex {
    Male(1, "M"), Female(2, "F"), Unknown(-9, "U"), Not_Reported(-10, "N");

    /** The integer code used when reading/writing ped files. */
    private final int code;

    /** The String symbol used when reading/writing vcf/gtc files */
    private final String symbol;

    /** Private constructor that takes the pedigree code for sex. */
    private Sex(final int code, final String symbol) {
        this.code = code;
        this.symbol = symbol;
    }

    /** Returns the code used to encode this sex in a ped/fam file. */
    public int toCode() { return this.code;}

    /** Decodes the Sex from a numeric code. Note that any value other than 1 or 2 will return Unknown. */
    public static Sex fromCode(final int code) {
        if (code == Male.code) return Male;
        else if (code == Female.code) return Female;
        else return Unknown;
    }

    public String toSymbol() {
        return this.symbol;
    }

    /**
     * Throw or return the Sex matching sexString.
     *
     * @param sexString  the string to match
     * @return the matched Sex
     */
    public static Sex fromString(String sexString) {
        final Predicate<Sex> match =
                s -> sexString.equalsIgnoreCase(s.symbol)
                        ||   sexString.equalsIgnoreCase(s.name());
        final List<Sex> genders = Stream.of(Sex.values()).filter(match).collect(Collectors.toList());
        if (genders.size() == 1) return genders.get(0);
        throw new PicardException("Unrecognized sex string: " + sexString);
    }
}
