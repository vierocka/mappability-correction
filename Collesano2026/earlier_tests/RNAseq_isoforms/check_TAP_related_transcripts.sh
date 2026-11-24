python3 dedup_kallisto_accessions.py HCC1937/kalisto/abundance.tsv > HCC1937.accession_dedup.tsv
python3 dedup_kallisto_accessions.py JY/abundance.tsv > JY.accession_dedup.tsv
python3 dedup_kallisto_accessions.py HCT116/abundance.tsv > HCT116.accession_dedup.tsv

# check TAP-related transcripts:

grep -H -E "NP_000584\.3|NP_001278951\.1|NP_000535\.3|NP_001276972\.1|NP_061313\.2|NP_001397804\.1|NP_003181\.3|NP_757345\.2|NP_757346\.2|NP_004039\.1|NP_004334\.1|NP_005304\.3|NP_001336173\.1|NP_057526\.3|NP_001185470\.1|NP_001035548\.1|NP_071745\.1|NP_001123612\.1|NP_001424731\.1|NP_001425687\.1|NP_001316158\.1|NP_001316162\.1|NP_004150\.1|NP_683720\.2|NP_002791\.1|NP_002784\.1|NP_001338284\.1" \
JY.accession_dedup.tsv \
HCC1937.accession_dedup.tsv  \
HCT116.accession_dedup.tsv
