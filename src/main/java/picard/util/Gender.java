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

package picard.util;

import picard.PicardException;

import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public enum Gender {
    MALE("M", "Male"),
    FEMALE("F", "Female"),
    UNKNOWN("U", "Unknown"),
    NOT_REPORTED("N", "Not Reported");

    private final String symbol;
    private final String name;

    Gender(String symbol, String name) {
        this.symbol = symbol;
        this.name = name;
    }

    /**
     * Throw or return the Gender matching genderString.
     *
     * @param genderString  the string to match
     * @return the matched Gender
     */
    public static Gender fromString(String genderString) {
        final Predicate<Gender> match =
                g -> genderString.equalsIgnoreCase(g.symbol)
                        ||   genderString.equalsIgnoreCase(g.name);
        final List<Gender> genders = Stream.of(Gender.values()).filter(match).collect(Collectors.toList());
        if (genders.size() == 1) return genders.get(0);
        throw new PicardException("Unrecognized gender string: " + genderString);
    }

    public String getSymbol() {
        return symbol;
    }

    public String getName() {
        return name;
    }
}

