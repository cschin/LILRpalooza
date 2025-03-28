# LILRpalooza

A imple API server that a user can
submit a LILR/LAIR sequence to get GFE annotation form 
https://feature.b12x.org

## build
Build the docker image locally:
```
cd docker/
bash build.sh
```
This will create an image tagged with `lilr_palooza`. 

## pre build image
You can pull down a pre-built image using `docker pull`
```
docker pull cschin/lilr_palooza-arm64 # for Apple silicon
docker pull cschin/lilr_palooza-amd64 # for amd64 architechure
```

## run the server
This start a API server running on port 8080, please ensure you have
the port free for the service. If not, you need to map it another
free port. 
```
docker run -p 8080:8080 lilr_palooza
```
or
```
docker run -p 8080:8080 cschin/lilr_palooza-arm64 # for Apple silicon
docker run -p 8080:8080 cschin/lilr_palooza-amd64 # for amd64 architechure
```

## example of submit a sequence for annotation
An example using `curl`
```
cd server/tests
bash test1.sh
```

This is the content of the test1.sh. You need construct 
a json object of three fields `gene_name`, `seq_name` and
`sequence`. The sequence should be from the start codon to
the stop codon including all introns of a gene's DNA sequence.
The current supported genes are "CDC42EP5", "LILRA1", "LILRB1", "LILRB4", "LAIR1", "LENG8", "LILRA2", "LILRA5", "LILRB2",
"LILRB5", "LAIR2", "LENG9", "LILRA3", "LILRA6", "LILRB3", "TTYH1".

```
curl -X POST http://127.0.0.1:8080/api/dna2gfe \
     -H "Content-Type: application/json" \
     -d '{"gene_name": "LILRB5", "seq_name": "HG03486.paternal.f1_assembly_v2_genbank::HG03486#1#JAHEOQ010000020.1_7673_303909_1:139353-145759 ID=MP000009;Rank=1;Identity=0.9983;Positive=0.9983;Target=LILRB5", "sequence": "ATGACCCTCACCCTCTCAGTCCTGATTTGCCTCGGTGAGGTTTGAAGAGGGGGAAGGAAGGTCCCCGTCTTGGAGGGAGCTCACTCTAAAGCGAGGCTCTGGTCTATCAGAGAATCTGGTCTATCAGAGGCTCCGAGGGAGGAGAGGAACTGCTGGGGCTTCCAGGGGCAAATCCCTCACAGGGAACTCTCTTCCAGGGCTGAGTGTGGGCCCCAGGACCTGCGTGCAGGCAGGTGAGTCTGTCCCCAGCTGTCCCAGGTCCCTTCTTCTCACTGGGGACAAGGGCCCAACCCCGGGCAGCTGGGGGTGGAGATAGCTGTTCTGGGCTGACTGATGGGGACGTCTGGAGGGTCCTGGGGCTGAGAACTGGAATCTGAGGGATGGGGATGTCTTGGGATCCAGCCTCTGATTCCATTCTAGGCACCCTCCCCAAACCCACCCTCTGGGCTGAGCCAGCCTCTGTGATAGCTCGGGGGAAGCCCGTGACCCTCTGGTGTCAGGGGCCCCTGGAGACTGAGGAGTACCGTCTGGATAAGGAGGGACTCCCATGGGCCCGGAAGAGACAGAACCCACTGGAGCCTGGAGCCAAGGCCAAGTTCCACATTCCATCCACGGTGTATGACAGTGCAGGGCGATACCGCTGCTACTATGAGACCCCTGCAGGCTGGTCAGAGCCCAGTGACCCCCTGGAGCTGGTGGCGACAGGTGAGAGGACACTCAGGGGTCCCAGCCCCAGGCTCTGCCCTCAGGAAGAGGGTCGGCTCTTAGGGACGTCTACCTCTCACAGCCCAGCCCTGGGGATGATGTGGGAGGTCGGAGCCCCACTTAAGACGTGCCTCCTTCTCTGCTAGGATTCTATGCAGAACCCACTCTTTTAGCCCTGCCGAGTCCTGTGGTGGCCTCAGGAGGAAATGTGACCCTCCAGTGTGATACACTGGACGGACTTCTCACGTTTGTTCTTGTTGAGGAAGAACAGAAGCTCCCCAGGACCCTGTACTCACAGAAGCTCCCCAAAGGGCCATCCCAGGCCCTGTTCCCTGTGGGTCCCGTGACCCCCAGCTGCAGGTGGAGGTTCAGATGCTATTACTATTACAGGAAAAACCCTCAGGTGTGGTCGAACCCCAGTGACCTCCTGGAGATTCTGGTCCCAGGTGAAAAAGCCACCACACTTCTTTATATAATTTTGGGGAACCAGATAGGTTGTTGGGAGTTTGGTTGATGACTGATCATGGCAAGGACCCCAGAAGGATGTGTTGATGGATGGGCTGAAGGCGTGAGGAAGACCCCACGGGGAGGCTCAGATGGGGAAACAGGAGCCTGAGTCACCCTCACCTGGAAGGGGTCGACTCAGGAAGGCAATGGGTGTATTTGCTGCAATTTCCTGTCCCTCGATGAGGAGAGGACAGACCAGACAGACAGTGGCCAGGAGTCAGAGAGACACTATCGGTCTGGAACTACTCCAAGACAGACCCAGGTGAGAAGGAGGCCCCGGGATCCGAGACACAGAGCGTGAGAGACAGTGAGACCTGCAGGGCCAGGACGCCAGGACGGGAGAAGGAAGGGGCGGGGGAGGAACCAGCCTTCCAAGTCCCAATTCCTCTTTCCCTCCAGGCGTGTCTAGGAAGCCCTCCCTCCTGATCCCGCAGGGCTCTGTCGTGGCCCGCGGAGGCAGCCTGACCCTGCAGTGTCGCTCTGATGTCGGCTATGACATATTCGTTCTGTACAAGGAGGGGGAACATGACCTCGTCCAGGGCTCTGGCCAGCAGCCCCAGGCTGGGCTCTCCCAGGCCAACTTCACCCTGGGCCCTGTGAGCCGCTCCCACGGGGGCCAGTACAGATGCTACGGTGCACACAACCTCTCCCCTAGGTGGTCGGCCCCCAGCGACCCCCTGGACATCCTGATCGCAGGTGAGGAGCCCAGCGGGTTCAGTCAGGGACCCAGGCTCTGCACAGGTCCTGCCGGGGGAATCCAATTAGTGATGGCCGGGATGAGGCGGGGGGGTGGTCCCAAGGGAGGGAGAGACAGACAGAGACAGGGGATGGGTGGGGAGGGGAAGACTCAGAGAAAACAGAGACAGAGGCTCCTAGAGAGGCCTGGGGAGGTCTCAGCTCAGAGCAAGGTGGGGCAGCCCCTCACCCATCCTTCTTCTCTCCAGGACTGATCCCTGACATACCCGCCCTCTCGGTGCAGCCGGGCCCCAAGGTGGCCTCAGGAGAGAACGTGACCCTGCTGTGTCAGTCATGGCATCAGATAGACACTTTCTTTTTGACCAAGGAGGGGGCAGCCCATCCCCCGCTGTGTCTAAAGTCAAAGTACCAGTCTTATAGACACCAGGCTGAATTCTCCATGAGTCCTGTGACCTCAGCCCAGGGTGGAACCTACCGATGCTACAGCGCAATCAGGTCCTACCCCTACCTGCTGTCCAGCCCTAGTTACCCCCAGGAGCTCGTGGTCTCAGGTGAGAGCCCTGACCCTGTCCTGTCCAAGCTCAAAGGCTCAGCTCAGGCCCTGCCCCCAGGAGAGCTCTGGGCTGGGATGGAGTCGCGGTGCGGGGGGGAGGGTTTGAGGGGGGCTCAGCCAGAGGGAGACTCACCCCTCAGAGGGGAGGAGGACAACGGGGGCTCCCCAGGCATGCCCACACTTGGCCCCATCTCCTGGGATGCAAATGGTGAAAGGTGAGCAGAAGAAAGTTTCCAGAGAAGCCACGGGCAGGTGGAGGGACGGGTTTCCTCACTCAGCACCAAAGCGCCTCGCTCCCTTTCTGTGCTTATTCCCAGGACCCTCTGGGGATCCCAGCCTCTCACCTACAGGCTCCACCCCCACACCTGGTGAGTCACTGAGGCCTCTGGGCTCGGAGGGAGCGTGGTCTCCCCCCAGGCAGCCCTGAGTCTCCCCGAGGATCCTATTCCCCTCAAAGACTCAAGCGGGAGCTTCCCTCCAGGGAGCTGGGCAGAGCCAGAGGAGGGGCCACAGGCTCCCCGGGGCTCTGAGGCTGGGCCGGTGAGGGGGCGGGCGTCGAGGCAGAGAGAGATGTTGGGTGTTGGGGCCCAGCCTGGGGGAGGAGCAGCCGGGCTGATGTGGGGAGCAGGGCAGCCCCAGCCCTCACCTCCCCGTCCTGACCCAGCAGGCCCTGAGGACCAGCCCCTCACCCCCACGGGGTTGGATCCCCAGAGTGGTGAGTGAGGGGCTCTGAGTGGGAGGTGGGCGGAGACCAGGGGAGGCAGGGGTGGGTTCTGTCGTAGGTTCAGGCTCCTCTGGAGATGGTGAAGTGCACAAGCCCTTCCCCTGCCTGGGCCTCAGTTTCTCCAAGTGTAAAGGAGAGAGGCCTGCATTGATGGGATTCTTCAGGGGACTGTCCTGTCCCACCGGCAGCAGTGACGGTGACCTGGGGCAGGGGAGGGGAGCAGGGCCGTGGTTTGGGGCATTCAGGCTCTTTCCCTGCAGCTCCGGGGCTCCGCTCAGGTGCAGAGAACAAGGGCTGCGGGTCAGACTCCTGGGTTCACTTCCCAGCTCTGCCGCATCCCACCGTGGGCCCAGGCAGGTCAACTTTCTACTCTGACTCAGTTTCAGCAGCTGTAAACTGGCTCAGTCCCATCCAGCTCACAGAACTGCTGCGAGGCGTAAGCAAAATCATGGGACCTGGCCCTGTACACAGCTCGGCAGGGGCACCGTCCTCCTGCTACCCTCAGCCCTTCCCAGATACACACAGAGCCCCTATCCAGACAGGTTCTGCATGGGAGTATGGGAACTTGGCAGAGTGGGAAACGGACCTGGCTGAGCTGGGAGTGAGAGCAATGCAGGGTCCGTCCTGCACAACCCACTCCCTCTCCCAGGCCCTGCTGTGCTGGGGAAGGGAGGATCCTAAGAAGGACACCAGCCCCAGATGGAGACACTAGGACAGGCCCCTCCTGTCAAATAGGAAACAGGTGTGCACCTGGTGGGGCAGCAGGAGGACAGCTGGGGAAAACACAAAGTCCCTGGTTCCCTTCCCAGACCTGCTTCTTCCAGGCTGAGGAGCCTGGGGCAGGCGATTCCCCTCTCTGAGCCTCAGTTTGCTCCTCTGTGAATTGGGGGGTTGGCAATCCCATGTTGCACAACTGCTGTGAGGGTTGGAGCTCATGAAGGAAAGACCTAGCTCGCGCCTGCACACAGAAGGTGCTCACATCAATGACGTCATCCCCATTCCCAACGTCATCACGCTCAAGGTCTGGGAAGGCACCTGGGGGTTGTGACTGGGGTCTCAGTGGCCTTCGTCCTGCTGCTGTTCCTCCTCCTCTTCCTCCTCCTCCGACATCGGCATCAGAGCAAACACAGGACATCGGGTGAGTAGGGAATGGGGGAACCCGTGGGCCGACCGAGGGTGGGCTCGGGGCACCAGCCAGAGGGAAACCAAACACAAAGGAAAGTCAGCTTAGAAAAACTGCTCCAGAAATTCCCAGGTGAAAAATCGATCGAGAAAGAAGAGAATAAATGTGAGCATGTGTGGAAGTGCTTGATTCTTCTGATTTTACTTTAAACTTACGACGTATTTAAAGCCTCAGTGCCAGTGGGCCTCCAGGTTTCCTTCTTTCCGCTCGAGTTGTGTGTGCAGGGCAGCTGGTTCGAATTCTCCCAGGCCTGACCCTCTGTCCATCTCTGTCCAGCCCATTTCTACCGTCCTGCAGGGGCTGCGGGGCCAGAGCCCAAGGACCAGGGCCTGCAGAAGAGGTAATTCTGCATGAAGACCCAAGACTCCCATCCACCCGCACAGCCCCCTCACTGCCCCTCACACTCCCGTGTCCTCCCCCAGGGCCAGCCCAGTTGCTGACATCCAGGAGGAAATTCTCAGTGAGTGACTAGAAGCGGAGGGCACCTGGGGTGGGCAAGGGAGCACCAAAGTTTCTGTAGCAATGGGGGCAGGAGCACAGGCTGGGAAGGGTCTGGGGCCAAGGGGGAGGTGGTCTGAACCCACACTGTGGGACCTCAGGGACATCACAGTCCCTCCCTGGATCTCAGCCACCCTAGTGGGAACAGGGCAAGGGCTGGCAGGACTGAGAAGTCTCAGAGAACCTTCCCAGGAGATGAACCCCTTGCTCTGGCCCAGCAGATGCTGCTGTGAAGGACACACAGCCCAAGGACGGGGTGGAGATGGATGCTCGGGTGAGGCCCCGCCCCTGTCCCGGGCACCAAAGGCCTCCTGGTGCCAGATCTAATCCAGCAGGACTTCTCTGTCCTCCTTCCCCCGGCTCTCAGCATCGTCACGGTGGACCCCTCCTTGTCCAGCACGCTGCCTCCCGCCTGCTGCGACCTCACTCTCTTCTGCTGTCCTGGGACCTCATGGGCCTCCTCCCGGGTCCCCTTCCTGCTCCTCATCCTCTGTTTGGCCGTCTGGTTGTTAGAGCTCTCCCCAGGCCTCAGGAGGATGAGGAATAAATGAACCAACCCGGTCCCCCAGGCTCCCCTTCATTCATTCAACCAGCGAGTGTTCCCAGGGAGCTCACTGTGGATGGGGCTCCCCATGGGAGCTGCAGACACAGCAGGGAGCAAAGCCGCCCCCGCCTCCTGAGCTCACCTCATGGTGGGAGACAAAATGCAAATAAATGCATCGTGTCCAGGAGTGCAACGTGCTGTAAGGAACATAAACCAGGGAAAGGGCAGAGAGTGTGGGGCAGTGGGGCCAGTCTGAATGGAAGGGGAGGGCTGTCTGCTCAGCTGTCATCTGAGAAGCCTGGACAGAGTGGGGCACACGATCCTCTGATGGACGAGCCCCTGCAGGCAGAGGAAACAGCCATGCAAAGGCCCCGAGGCAGCAGCGAGCTCTTGCGGGAAGGCCCGTGAGGCTGCAGCCAAATGGGCAAGGTCAGAGTGAGGAGCAGAGGCCAGAACCACAGGGAGGGAGCGGCCAGACCCTCCACGGCCTTAGGGCGTCCCTGAGATTCCATCGGGAAAGGGATGTAATCGGATCACCCCGGGAACAGTGAGGAAAATTGACTCCAGGAGGTCAGGGGGACTCAAGGACACCCCCCACCACTGTCTCTCTCCAGCAGAGCCCACACGATGAAGACCCCCAGGCAGTGACGTATGCCGAGGTGAAACACTCCAGACCTAGGAGAGAAATGGCCTCTCCTCCCTCCCCACTGTCTGGGGAATTCCTGGACACAAAGGACAGACAGGCAGAAGAGGACAGACAGATGGACACTGAGAGAGTCCTTTCCTCTCCAGGCCCCCAGGCCTCCCCACCCCCACCACGTTCCTTCCCTCTCACTCTCCCCCGCTGCAGGCTGCTGCATCTGAAGCCCCCCAGGATGTGACCTACGCCCAGCTGCACAGCTTGACCCTTAGACGGGAGGCAACTGAGCCTCCTCCATCCCAGGAAAGGGAACCTCCAGCTGAACCCAGCATCTACGCCCCCCTGGCCATCCACTAG"}'
```

## Method 
The submitted gene sequnece are aligned to a reference set of protein
sequences (in `server/ref_data/`) to identify the exons and introns with `miniprot` with the following parameter:
```
./miniprot dna_sequence.fa protein_ref.fa -j 2 --trans --aln --max-intron-out 20000 -G 20000 --outs=0.975 --outc=0.8 --gff
```

In the output from `miniprot`, we use the `##ATN` line to identify
the exons and introns (annotated with upper case and lower case letter respectively). The exon and intron sequences are sent to `https://feature.b12x.org:443/features` to get annotated GFE features. The output
json records of the server are concatenated and returne by the LILRpalooza API server. Here is an output example:

```
{
  "valid_input": true,
  "seq_name": "HG03486.paternal.f1_assembly_v2_genbank::HG03486#1#JAHEOQ010000020.1_7673_303909_1:139353-145759 ID=MP000009;Rank=1;Identity=0.9983;Positive=0.9983;Target=LILRB5",
  "data": [
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "exon",
      "rank": 1,
      "sequence": "ATGACCCTCACCCTCTCAGTCCTGATTTGCCTCG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "intron",
      "rank": 1,
      "sequence": "GTGAGGTTTGAAGAGGGGGAAGGAAGGTCCCCGTCTTGGAGGGAGCTCACTCTAAAGCGAGGCTCTGGTCTATCAGAGAATCTGGTCTATCAGAGGCTCCGAGGGAGGAGAGGAACTGCTGGGGCTTCCAGGGGCAAATCCCTCACAGGGAACTCTCTTCCAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "exon",
      "rank": 2,
      "sequence": "GGCTGAGTGTGGGCCCCAGGACCTGCGTGCAGGCAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "intron",
      "rank": 2,
      "sequence": "GTGAGTCTGTCCCCAGCTGTCCCAGGTCCCTTCTTCTCACTGGGGACAAGGGCCCAACCCCGGGCAGCTGGGGGTGGAGATAGCTGTTCTGGGCTGACTGATGGGGACGTCTGGAGGGTCCTGGGGCTGAGAACTGGAATCTGAGGGATGGGGATGTCTTGGGATCCAGCCTCTGATTCCATTCTAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "exon",
      "rank": 3,
      "sequence": "GCACCCTCCCCAAACCCACCCTCTGGGCTGAGCCAGCCTCTGTGATAGCTCGGGGGAAGCCCGTGACCCTCTGGTGTCAGGGGCCCCTGGAGACTGAGGAGTACCGTCTGGATAAGGAGGGACTCCCATGGGCCCGGAAGAGACAGAACCCACTGGAGCCTGGAGCCAAGGCCAAGTTCCACATTCCATCCACGGTGTATGACAGTGCAGGGCGATACCGCTGCTACTATGAGACCCCTGCAGGCTGGTCAGAGCCCAGTGACCCCCTGGAGCTGGTGGCGACAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "intron",
      "rank": 3,
      "sequence": "GTGAGAGGACACTCAGGGGTCCCAGCCCCAGGCTCTGCCCTCAGGAAGAGGGTCGGCTCTTAGGGACGTCTACCTCTCACAGCCCAGCCCTGGGGATGATGTGGGAGGTCGGAGCCCCACTTAAGACGTGCCTCCTTCTCTGCTAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "exon",
      "rank": 4,
      "sequence": "GATTCTATGCAGAACCCACTCTTTTAGCCCTGCCGAGTCCTGTGGTGGCCTCAGGAGGAAATGTGACCCTCCAGTGTGATACACTGGACGGACTTCTCACGTTTGTTCTTGTTGAGGAAGAACAGAAGCTCCCCAGGACCCTGTACTCACAGAAGCTCCCCAAAGGGCCATCCCAGGCCCTGTTCCCTGTGGGTCCCGTGACCCCCAGCTGCAGGTGGAGGTTCAGATGCTATTACTATTACAGGAAAAACCCTCAGGTGTGGTCGAACCCCAGTGACCTCCTGGAGATTCTGGTCCCAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "intron",
      "rank": 4,
      "sequence": "GTGAAAAAGCCACCACACTTCTTTATATAATTTTGGGGAACCAGATAGGTTGTTGGGAGTTTGGTTGATGACTGATCATGGCAAGGACCCCAGAAGGATGTGTTGATGGATGGGCTGAAGGCGTGAGGAAGACCCCACGGGGAGGCTCAGATGGGGAAACAGGAGCCTGAGTCACCCTCACCTGGAAGGGGTCGACTCAGGAAGGCAATGGGTGTATTTGCTGCAATTTCCTGTCCCTCGATGAGGAGAGGACAGACCAGACAGACAGTGGCCAGGAGTCAGAGAGACACTATCGGTCTGGAACTACTCCAAGACAGACCCAGGTGAGAAGGAGGCCCCGGGATCCGAGACACAGAGCGTGAGAGACAGTGAGACCTGCAGGGCCAGGACGCCAGGACGGGAGAAGGAAGGGGCGGGGGAGGAACCAGCCTTCCAAGTCCCAATTCCTCTTTCCCTCCAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "exon",
      "rank": 5,
      "sequence": "GCGTGTCTAGGAAGCCCTCCCTCCTGATCCCGCAGGGCTCTGTCGTGGCCCGCGGAGGCAGCCTGACCCTGCAGTGTCGCTCTGATGTCGGCTATGACATATTCGTTCTGTACAAGGAGGGGGAACATGACCTCGTCCAGGGCTCTGGCCAGCAGCCCCAGGCTGGGCTCTCCCAGGCCAACTTCACCCTGGGCCCTGTGAGCCGCTCCCACGGGGGCCAGTACAGATGCTACGGTGCACACAACCTCTCCCCTAGGTGGTCGGCCCCCAGCGACCCCCTGGACATCCTGATCGCAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "intron",
      "rank": 5,
      "sequence": "GTGAGGAGCCCAGCGGGTTCAGTCAGGGACCCAGGCTCTGCACAGGTCCTGCCGGGGGAATCCAATTAGTGATGGCCGGGATGAGGCGGGGGGGTGGTCCCAAGGGAGGGAGAGACAGACAGAGACAGGGGATGGGTGGGGAGGGGAAGACTCAGAGAAAACAGAGACAGAGGCTCCTAGAGAGGCCTGGGGAGGTCTCAGCTCAGAGCAAGGTGGGGCAGCCCCTCACCCATCCTTCTTCTCTCCAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "exon",
      "rank": 6,
      "sequence": "GACTGATCCCTGACATACCCGCCCTCTCGGTGCAGCCGGGCCCCAAGGTGGCCTCAGGAGAGAACGTGACCCTGCTGTGTCAGTCATGGCATCAGATAGACACTTTCTTTTTGACCAAGGAGGGGGCAGCCCATCCCCCGCTGTGTCTAAAGTCAAAGTACCAGTCTTATAGACACCAGGCTGAATTCTCCATGAGTCCTGTGACCTCAGCCCAGGGTGGAACCTACCGATGCTACAGCGCAATCAGGTCCTACCCCTACCTGCTGTCCAGCCCTAGTTACCCCCAGGAGCTCGTGGTCTCAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "intron",
      "rank": 6,
      "sequence": "GTGAGAGCCCTGACCCTGTCCTGTCCAAGCTCAAAGGCTCAGCTCAGGCCCTGCCCCCAGGAGAGCTCTGGGCTGGGATGGAGTCGCGGTGCGGGGGGGAGGGTTTGAGGGGGGCTCAGCCAGAGGGAGACTCACCCCTCAGAGGGGAGGAGGACAACGGGGGCTCCCCAGGCATGCCCACACTTGGCCCCATCTCCTGGGATGCAAATGGTGAAAGGTGAGCAGAAGAAAGTTTCCAGAGAAGCCACGGGCAGGTGGAGGGACGGGTTTCCTCACTCAGCACCAAAGCGCCTCGCTCCCTTTCTGTGCTTATTCCCAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "exon",
      "rank": 7,
      "sequence": "GACCCTCTGGGGATCCCAGCCTCTCACCTACAGGCTCCACCCCCACACCTG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "intron",
      "rank": 7,
      "sequence": "GTGAGTCACTGAGGCCTCTGGGCTCGGAGGGAGCGTGGTCTCCCCCCAGGCAGCCCTGAGTCTCCCCGAGGATCCTATTCCCCTCAAAGACTCAAGCGGGAGCTTCCCTCCAGGGAGCTGGGCAGAGCCAGAGGAGGGGCCACAGGCTCCCCGGGGCTCTGAGGCTGGGCCGGTGAGGGGGCGGGCGTCGAGGCAGAGAGAGATGTTGGGTGTTGGGGCCCAGCCTGGGGGAGGAGCAGCCGGGCTGATGTGGGGAGCAGGGCAGCCCCAGCCCTCACCTCCCCGTCCTGACCCAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "exon",
      "rank": 8,
      "sequence": "CAGGCCCTGAGGACCAGCCCCTCACCCCCACGGGGTTGGATCCCCAGAGTG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "intron",
      "rank": 8,
      "sequence": "GTGAGTGAGGGGCTCTGAGTGGGAGGTGGGCGGAGACCAGGGGAGGCAGGGGTGGGTTCTGTCGTAGGTTCAGGCTCCTCTGGAGATGGTGAAGTGCACAAGCCCTTCCCCTGCCTGGGCCTCAGTTTCTCCAAGTGTAAAGGAGAGAGGCCTGCATTGATGGGATTCTTCAGGGGACTGTCCTGTCCCACCGGCAGCAGTGACGGTGACCTGGGGCAGGGGAGGGGAGCAGGGCCGTGGTTTGGGGCATTCAGGCTCTTTCCCTGCAGCTCCGGGGCTCCGCTCAGGTGCAGAGAACAAGGGCTGCGGGTCAGACTCCTGGGTTCACTTCCCAGCTCTGCCGCATCCCACCGTGGGCCCAGGCAGGTCAACTTTCTACTCTGACTCAGTTTCAGCAGCTGTAAACTGGCTCAGTCCCATCCAGCTCACAGAACTGCTGCGAGGCGTAAGCAAAATCATGGGACCTGGCCCTGTACACAGCTCGGCAGGGGCACCGTCCTCCTGCTACCCTCAGCCCTTCCCAGATACACACAGAGCCCCTATCCAGACAGGTTCTGCATGGGAGTATGGGAACTTGGCAGAGTGGGAAACGGACCTGGCTGAGCTGGGAGTGAGAGCAATGCAGGGTCCGTCCTGCACAACCCACTCCCTCTCCCAGGCCCTGCTGTGCTGGGGAAGGGAGGATCCTAAGAAGGACACCAGCCCCAGATGGAGACACTAGGACAGGCCCCTCCTGTCAAATAGGAAACAGGTGTGCACCTGGTGGGGCAGCAGGAGGACAGCTGGGGAAAACACAAAGTCCCTGGTTCCCTTCCCAGACCTGCTTCTTCCAGGCTGAGGAGCCTGGGGCAGGCGATTCCCCTCTCTGAGCCTCAGTTTGCTCCTCTGTGAATTGGGGGGTTGGCAATCCCATGTTGCACAACTGCTGTGAGGGTTGGAGCTCATGAAGGAAAGACCTAGCTCGCGCCTGCACACAGAAGGTGCTCACATCAATGACGTCATCCCCATTCCCAACGTCATCACGCTCAAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "exon",
      "rank": 9,
      "sequence": "GTCTGGGAAGGCACCTGGGGGTTGTGACTGGGGTCTCAGTGGCCTTCGTCCTGCTGCTGTTCCTCCTCCTCTTCCTCCTCCTCCGACATCGGCATCAGAGCAAACACAGGACATCGG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "intron",
      "rank": 9,
      "sequence": "GTGAGTAGGGAATGGGGGAACCCGTGGGCCGACCGAGGGTGGGCTCGGGGCACCAGCCAGAGGGAAACCAAACACAAAGGAAAGTCAGCTTAGAAAAACTGCTCCAGAAATTCCCAGGTGAAAAATCGATCGAGAAAGAAGAGAATAAATGTGAGCATGTGTGGAAGTGCTTGATTCTTCTGATTTTACTTTAAACTTACGACGTATTTAAAGCCTCAGTGCCAGTGGGCCTCCAGGTTTCCTTCTTTCCGCTCGAGTTGTGTGTGCAGGGCAGCTGGTTCGAATTCTCCCAGGCCTGACCCTCTGTCCATCTCTGTCCAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "exon",
      "rank": 10,
      "sequence": "CCCATTTCTACCGTCCTGCAGGGGCTGCGGGGCCAGAGCCCAAGGACCAGGGCCTGCAGAAGAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "intron",
      "rank": 10,
      "sequence": "GTAATTCTGCATGAAGACCCAAGACTCCCATCCACCCGCACAGCCCCCTCACTGCCCCTCACACTCCCGTGTCCTCCCCCAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "exon",
      "rank": 11,
      "sequence": "GGCCAGCCCAGTTGCTGACATCCAGGAGGAAATTCTCA"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "intron",
      "rank": 11,
      "sequence": "GTGAGTGACTAGAAGCGGAGGGCACCTGGGGTGGGCAAGGGAGCACCAAAGTTTCTGTAGCAATGGGGGCAGGAGCACAGGCTGGGAAGGGTCTGGGGCCAAGGGGGAGGTGGTCTGAACCCACACTGTGGGACCTCAGGGACATCACAGTCCCTCCCTGGATCTCAGCCACCCTAGTGGGAACAGGGCAAGGGCTGGCAGGACTGAGAAGTCTCAGAGAACCTTCCCAGGAGATGAACCCCTTGCTCTGGCCCAGCAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "exon",
      "rank": 12,
      "sequence": "ATGCTGCTGTGAAGGACACACAGCCCAAGGACGGGGTGGAGATGGATGCTCGG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "intron",
      "rank": 12,
      "sequence": "GTGAGGCCCCGCCCCTGTCCCGGGCACCAAAGGCCTCCTGGTGCCAGATCTAATCCAGCAGGACTTCTCTGTCCTCCTTCCCCCGGCTCTCAGCATCGTCACGGTGGACCCCTCCTTGTCCAGCACGCTGCCTCCCGCCTGCTGCGACCTCACTCTCTTCTGCTGTCCTGGGACCTCATGGGCCTCCTCCCGGGTCCCCTTCCTGCTCCTCATCCTCTGTTTGGCCGTCTGGTTGTTAGAGCTCTCCCCAGGCCTCAGGAGGATGAGGAATAAATGAACCAACCCGGTCCCCCAGGCTCCCCTTCATTCATTCAACCAGCGAGTGTTCCCAGGGAGCTCACTGTGGATGGGGCTCCCCATGGGAGCTGCAGACACAGCAGGGAGCAAAGCCGCCCCCGCCTCCTGAGCTCACCTCATGGTGGGAGACAAAATGCAAATAAATGCATCGTGTCCAGGAGTGCAACGTGCTGTAAGGAACATAAACCAGGGAAAGGGCAGAGAGTGTGGGGCAGTGGGGCCAGTCTGAATGGAAGGGGAGGGCTGTCTGCTCAGCTGTCATCTGAGAAGCCTGGACAGAGTGGGGCACACGATCCTCTGATGGACGAGCCCCTGCAGGCAGAGGAAACAGCCATGCAAAGGCCCCGAGGCAGCAGCGAGCTCTTGCGGGAAGGCCCGTGAGGCTGCAGCCAAATGGGCAAGGTCAGAGTGAGGAGCAGAGGCCAGAACCACAGGGAGGGAGCGGCCAGACCCTCCACGGCCTTAGGGCGTCCCTGAGATTCCATCGGGAAAGGGATGTAATCGGATCACCCCGGGAACAGTGAGGAAAATTGACTCCAGGAGGTCAGGGGGACTCAAGGACACCCCCCACCACTGTCTCTCTCCAGCAGAGCCCACACGATGAAGACCCCCAGGCAGTGACGTATGCCGAGGTGAAACACTCCAGACCTAGGAGAGAAATGGCCTCTCCTCCCTCCCCACTGTCTGGGGAATTCCTGGACACAAAGGACAGACAGGCAGAAGAGGACAGACAGATGGACACTGAGAGAGTCCTTTCCTCTCCAGGCCCCCAGGCCTCCCCACCCCCACCACGTTCCTTCCCTCTCACTCTCCCCCGCTGCAG"
    },
    {
      "locus": "LILRB5",
      "accession": 1,
      "term": "exon",
      "rank": 13,
      "sequence": "GCTGCTGCATCTGAAGCCCCCCAGGATGTGACCTACGCCCAGCTGCACAGCTTGACCCTTAGACGGGAGGCAACTGAGCCTCCTCCATCCCAGGAAAGGGAACCTCCAGCTGAACCCAGCATCTACGCCCCCCTGGCCATCCACTAG"
    }
  ]
}
```

## Issues
This is just from a couple of hours coding for the hackathon #dash16.
It has very limited testing and minimum error handling.


