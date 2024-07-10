#!/usr/bin/awk -f

# this script merges male and female samples from the sync file

{
    # Print the first three columns (chromosome, position, and reference allele)
    # with a tab between each
    printf "%s\t%s\t%s\t", $1, $2, $3;

    # Iterate over pairs of columns starting from the 4th column up to the 34th
    for (i = 4; i <= 34; i += 2) {
        # Split the allele counts of the first column in the pair into an array
        split($i, countsA, ":");

        # Split the allele counts of the second column in the pair into another array
        split($(i+1), countsB, ":");

        # Initialize a variable to store the sum for each allele type
        sumCounts = "";

        # Sum the allele counts for the pair
        for (j = 1; j <= length(countsA); j++) {
            countsA[j] += countsB[j];
            sumCounts = sumCounts countsA[j]; 

            # Add a colon after each count, except for the last one
            if (j < length(countsA)) sumCounts = sumCounts ":";
        }

        # Print the summed counts for this pair
        printf "%s", sumCounts;

        # Add a tab separator between pairs, except after the last pair
        if (i < 34) printf "\t";
    }
    printf "\n";
}

