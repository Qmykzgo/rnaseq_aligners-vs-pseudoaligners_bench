# extract.R
# Load and inspect the qPCR dataset from RNAontheBENCH and export it for downstream analysis.

load("/home/yer_kanat/RNAontheBENCH/data/qpcr.rda")
print("Objects loaded from qpcr.rda:")
print(ls())

# Print structure and head of the object to see its format
if (exists("qpcr")) {
    print("Structure of qpcr object:")
    print(str(qpcr))
    print("Head of qpcr:")
    print(head(qpcr))
    
    # Save to the expected location
    write.table(qpcr, file="/home/yer_kanat/Downloads/rnqbench/data/truth_tables/GSE83402_qPCR_normalized.txt", sep="\t", quote=FALSE, row.names=TRUE)
    print("Successfully exported qPCR data to GSE83402_qPCR_normalized.txt")
} else {
    print("Error: qpcr object not found in the loaded workspace.")
}
