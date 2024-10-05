import pandas as pd
import sys

genome_version = sys.argv[1]

output = sys.argv[2]
input = sys.argv[3]

# Import transcript dataframe
dfTranscript = pd.read_csv(input, sep="\t")

# Depending on genome_version, use different rules to generate new 'Gene_name' column, based on existing 'Transcript_name' column:
if genome_version == "Zmays_Ensembl_18":
    dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 3) == "GRM", "Gene_name"
    ] = dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 3) == "GRM", "Transcript_name"
    ].str.slice(
        0, 13
    )
    dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 1) == "A", "Gene_name"
    ] = (
        dfTranscript.loc[
            dfTranscript.Transcript_name.str.slice(0, 1) == "A", "Transcript_name"
        ].str.slice(0, 13)
        + "P"
        + dfTranscript.loc[
            dfTranscript.Transcript_name.str.slice(0, 1) == "A", "Transcript_name"
        ].str.slice(
            13,
        )
    )
    dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 1) == "E", "Gene_name"
    ] = (
        dfTranscript.loc[
            dfTranscript.Transcript_name.str.slice(0, 1) == "E", "Transcript_name"
        ].str.slice(0, 13)
        + "P"
        + dfTranscript.loc[
            dfTranscript.Transcript_name.str.slice(0, 1) == "E", "Transcript_name"
        ].str.slice(
            13,
        )
    )

elif genome_version == "Osativa_323_v7.0":
    dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 3) == "LOC", "Gene_name"
    ] = dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 3) == "LOC", "Transcript_name"
    ].str.slice(
        0, 14
    )
    dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 3) == "Chr", "Gene_name"
    ] = dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 3) == "Chr", "Transcript_name"
    ]

elif genome_version == "Pvirgatum_516_v5.1":
    dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 7) != "Pavir.J", "Gene_name"
    ] = dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 7) != "Pavir.J", "Transcript_name"
    ].str.slice(
        0, 15
    )
    dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 7) == "Pavir.J", "Gene_name"
    ] = dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 7) == "Pavir.J", "Transcript_name"
    ].str.slice(
        0, 13
    )

elif genome_version == "Sbicolor_454_v3.1.1":
    dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 7) == "Sobic.0", "Gene_name"
    ] = dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 7) == "Sobic.0", "Transcript_name"
    ].str.slice(
        0, 16
    )
    dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 7) == "Sobic.K", "Gene_name"
    ] = dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 7) == "Sobic.0", "Transcript_name"
    ].str.slice(
        0, 13
    )

elif genome_version == "Sviridis_500_v2.1":
    dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 6) == "Sevir.", "Gene_name"
    ] = dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 6) == "Sevir.", "Transcript_name"
    ].str.slice(
        0, 14
    )

elif genome_version == "Athaliana_TAIR10":
    dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 2) == "AT", "Gene_name"
    ] = dfTranscript.loc[
        dfTranscript.Transcript_name.str.slice(0, 2) == "AT", "Transcript_name"
    ].str.slice(
        0, 9
    )

else:
    print(
        "Did not recognize genome version. Make sure spelling is correct. Currently supported genome versions:"
    )
    print(" Zmays_Ensembl_18")
    print(" Osativa_323_v7.0")
    print(" Pvirgatum_516_v5.1")
    print(" Sbicolor_454_v3.1.1")
    print(" Sviridis_500_v2.1")
    print(" Athaliana_TAIR10")

# Generate Gene dataframe by summing Transcript dataframe based on Gene column
dfGene = dfTranscript.groupby(["Gene_name"], as_index=False).sum()
dfGene.columns = dfGene.columns.str.replace("Transcript", "Gene")

# Write to CSV
dfGene.to_csv(output, sep="\t", mode="w", index=False)
