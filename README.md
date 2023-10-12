# PrDb
A protein sequence fetching and database-building tool

### Introduction

The PRDB (Protein Retrieval Database) program is designed to fetch protein sequences of Influenza A Virus (IAV) from the NCBI database using taxonomic IDs. This manual provides instructions on how to use the program effectively.

### Prerequisites

- Python 3.x
- Biopython library (install using `pip install biopython`)

### Usage

```
python prdb.py [-h] (--email EMAIL{STR}) (--txid TXID [TXID ...][INT]) [--db DB [DB ...]{STR}] 
    [--sub SUB [SUB ...]{STR}] [--gene GENE [GENE ...]{STR}] [--min MIN][INT] [--max MAX][INT] 
    [--meta <filename>]{STR} [--F <filename>]{STR}
```

#### Options:

- `-h`: Show the help message and exit.
- `--email`: Your email ID to access the NCBI database. Required.
- `--txid`: Taxonomic ID of the organism(s) of interest. Can provide multiple IDs. Required.
- `--db`: Search specific database(s). Can provide multiple database names.
- `--sub`: Search for specific subtypes. Can provide multiple subtype names.
- `--gene`: Search for specific gene(s) encoding the desired protein(s). Can provide multiple gene names.
- `--min`: Threshold for the minimum length of the sequence.
- `--max`: Threshold for the maximum length of the sequence.
- `--full`: Only search for complete sequences.
- `--meta`: Output filename to store the metadata for the sequences found.
- `--F`: Output filename to store the sequences found.

### Steps to Operate

1. Provide Your Email: Start by entering your email ID using the `--email` option. This is necessary to access the NCBI database.

2. Taxonomic IDs (txid): Use the `--txid` option to specify the taxonomic IDs of the organisms you want to search for. You can provide multiple IDs separated by spaces. ** Use 11320 for IAV viral proteomes. It is required to input the taxonomical ID to initiate the search.

3. Search Criteria:
   - Database (db): If you want to search in specific databases, use the `--db` option and provide one or more database names. e.g. genbank, refseq or swissprot.
   - Subtype (sub): Use the `--sub` option to search for specific subtypes. Provide subtype names.
   - Gene Name (gene): To search for sequences associated with specific gene names, use the `--gene` option and provide gene names.

4. Sequence Length Criteria:
   - Minimum Length (min): Use the `--min` option to specify the minimum sequence length.
   - Maximum Length (max): Use the `--max` option to specify the maximum sequence length.

5. Complete Sequences Only: If you're interested in complete sequences only, use the `--full` flag.

6. Output Files:
   - Sequences File (--F): To save the searched protein sequences, use the `--F` option followed by the desired output filename. It outputs a fasta file.
   - Metadata File (--meta): To save associated metadata, use the `--meta` option followed by the desired output filename. It outputs a CSV file.

### Example Usage

1. Fetch protein sequences of Influenza A Virus subtype H1N1 from taxonomic ID 11320:
   ```
   python prdb.py --email your@email.com --txid 11320 --sub H1N1 --F sequences_output
   ```

2. Fetch sequences from taxonomic IDs 11320, search for subtype H1N1 and H3N2, and save metadata:
   ```
   python prdb.py --email your@email.com --txid 11320 --sub H1N1 H3N2 --meta metadata_output
   ```
