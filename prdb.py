# PRDB
# NOTE - some strains don't have host, but it is not human, it may be written in the strain itself
# Import modules
import urllib.error
from Bio import Entrez
import time
import pandas as pd
import argparse
import re
import sys
import logging

# Implement argparse
custom_usage = """
python prdb.py [-h] (--email EMAIL{STR}) (--txid TXID [TXID ...][INT]) [--db DB [DB ...]{STR}] 
    [--sub SUB [SUB ...]{STR}] [--gene GENE [GENE ...]{STR}] [--min MIN][INT] [--max MAX][INT] 
    [--meta <filename>]{STR} [--F <filename>]{STR}
        
Description:
    Fetch protein sequences of IAV from the NCBI database using taxonomic ID (Version 1.6)

Notes:
    1. Options inside () are required, while [] are optional.
    2. Options can be called in any order.
    3. Options with [--Option [OPTION]] indicates it can take in a list of variables separated by 
        whitespaces, e.g. --sub H1N1 H3N2 H10N7 ... HnNn.
    4. <filename> has to be replaced with the name of the file for the specific Option.
    5. {STR} indicates the Option takes in only strings as input.    
    6. {INT} indicates the Option takes in only integers as input.    
"""
parser = argparse.ArgumentParser(usage=custom_usage)
# Adding arguments for the search parameters as well as the output files
parser.add_argument("--email", default=None, type=str, required=True, help="Email ID to access NCBI database")
parser.add_argument("--txid", default=None, type=int,
                    nargs="+", action="append", help="Taxonomic ID")
parser.add_argument("--db", default=None, type=str,
                    nargs="+", action="append", help="Search specific Database/s")
parser.add_argument("--sub", default=None, type=str,
                    nargs="+", action="append", help="Search for specific subtypes")
parser.add_argument("--gene", default=None, type=str,
                    nargs="+", action="append", help="Search for specific gene/s encoding the desired protein/s")
parser.add_argument("--min", default=None, type=int, help="Threshold for minimum length of sequence")
parser.add_argument("--max", default=None, type=int, help="Threshold for maximum length of sequence")
parser.add_argument("--full", default=False, action="store_true", help="Only search for the complete sequences")
parser.add_argument("--meta", default=None, type=str,
                    help="Output filename to store the metadata for the sequences found")
parser.add_argument("--F", default=None, type=str, help="Output filename to store the sequences found")
args = parser.parse_args()

# Implement logger
# Create logging level of error
logger = logging.getLogger()
logger.setLevel(logging.ERROR)
# Add logging handler to print messages onto the console
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(levelname)s - %(message)s'))
logger.addHandler(ch)


# Functions

# Validating email ID:

def is_valid_email(email):
    pattern = r'^[\w\.-]+@[\w\.-]+\.\w+$'  # Pattern for verifying email address
    if re.match(pattern, email):
        return True
    else:
        return False


# Search field formatting

def Filter(term):  # Adds the filter criteria for NCBI searching. Essential for specifying database
    input = term.copy()
    keyword = ""
    if len(input) >= 1:
        if len(input) == 1:
            input = input[0] + "[filter]"
            keyword += input
        else:
            for i in range(len(input)):
                input[i] = input[i] + "[filter]"
            keyword += " OR ".join(input)  # Joins the search terms with OR functionality to indicate optionality
        keyword = "(" + keyword + ")"
    else:
        keyword = None
    return keyword


# Organism criteria function

def Organism(term):  # Adds the Organism criteria for specifying the term is associated with organism of interest
    input = term.copy()
    keyword = ""
    if len(input) >= 1:
        if len(input) == 1:
            input = "txid" + str(input[0]) + "[Organism]"
            keyword += input
        else:
            for i in range(len(input)):
                input[i] = "txid" + str(input[i]) + "[Organism]"
            keyword += " OR ".join(input)  # Joins the search terms with OR functionality to indicate optionality
        keyword = "(" + keyword + ")"
    else:
        keyword = None
    return keyword


# Subtype criteria function

def Subtype(term):  # This field is unique to viruses and may not work for other organisms
    input = term.copy()
    keyword = ""
    if len(input) >= 1:
        if len(input) == 1:
            input = input[0] + " subtype[Organism]"  # Specifies the subtype to look for in the database
            keyword += input
        else:
            for i in range(len(input)):
                input[i] = input[i] + " subtype[Organism]"
            keyword += " OR ".join(input)  # Joins the search terms with OR functionality to indicate optionality
        keyword = "(" + keyword + ")"
    else:
        keyword = None
    return keyword


# Gene Name criteria function

def Gene_name(term):  # Specifies the gene encoding the protein of interest
    input = term.copy()
    keyword = ""
    if len(input) >= 1:
        if len(input) == 1:
            input = input[0] + "[Gene Name]"
            keyword += input
        else:
            for i in range(len(input)):
                input[i] = input[i] + "[Gene Name]"
            keyword += " OR ".join(input)  # Joins the search terms with OR functionality to indicate optionality
        keyword = "(" + keyword + ")"
    else:
        keyword = None
    return keyword


def minlength(term):  # Minimum length filtering criteria for fetching specific lengths of sequences
    if term is not None:
        input = int(term) - 1
        keyword = f"(0:{input}[Sequence Length])"
    else:
        keyword = None
    return keyword


def maxlength(term, min):  # Maximum length filtering criteria for fetching specific lengths of sequences
    keyword = ""
    if term is not None:
        if min is not None:
            keyword += f"({min}:{term}[Sequence Length])"
        if min is None:
            keyword += f"(0:{term}[Sequence Length])"
    else:
        keyword = None
    return keyword


def search_term(org, db, subtype, gene, min, max):
    # Gets all the search inputs and concatenates them with AND or NOT functionality to search the NCBI database
    keyword = ""
    search_list = []
    if org is not None:
        search_list.append(org)
    if db is not None:
        search_list.append(db)
    if subtype is not None:
        search_list.append(subtype)
    if gene is not None:
        search_list.append(gene)
    if max is not None:
        search_list.append(max)
    keyword += " AND ".join(search_list)
    if min is not None:
        keyword += f" NOT {min}"  # The minimum sequence length requires a NOT criteria due to its unique usage
    return keyword


# Spelling check

# E search

def Esearch(term):  # Connects to the NCBI Entrez Direct API to search for the sequences
    try:
        handle = Entrez.esearch(db="protein", term=f"{term}",
                                usehistory="y", idtype="acc")
        record = Entrez.read(handle)
        handle.close()
        return record
    except urllib.error.HTTPError as e:
        # Handle HTTPError
        logger.error("An HTTP error occurred:", e)
    except urllib.error.URLError as e:
        # Handle URLError
        logger.error("A URL error occurred:", e)
    except TimeoutError:
        # Handle TimeoutError
        logger.error("The request timed out.")
    except Exception as e:
        # Handle other exceptions
        logger.error("An error occurred:", e)


# E fetch

def Fetch_sequences_in_batches(record, file):  # Connects to the NCBI Entrez Direct API to fetch the sequences
    batch_size = 9999  # Using batches to compensate for the limitation of retmax.
    # Fetch sequence in numbers higher than 10,000
    count = int(record["Count"])
    try:
        out_handle = open(fr"{file}.fasta", "a")
        for start in range(0, count, batch_size):
            end = min(start + batch_size - 1, count - 1)
            print(f"Downloading sequence {start + 1} to {end + 1}")
            try:
                handle_fetch = Entrez.efetch(
                    db="protein",
                    rettype="fasta",
                    retmode="text",
                    retstart=start,
                    retmax=int(end - start + 1),
                    webenv=record["WebEnv"],
                    query_key=record["QueryKey"],
                    idtype="acc")
                data = handle_fetch.read()
                handle_fetch.close()
                out_handle.write(data)
            except urllib.error.HTTPError as e:
                # Handle HTTPError
                logger.error("An HTTP error occurred:", e)
                pass
            except urllib.error.URLError as e:
                # Handle URLError
                logger.error("A URL error occurred:", e)
                pass
            except TimeoutError:
                # Handle TimeoutError
                logger.error("The request timed out.")
                pass
            except Exception as e:
                # Handle other exceptions
                logger.error("An error occurred:", e)
                pass
        out_handle.close()
    except FileNotFoundError as e:
        logger.error(f"File not found\n{e}")
    except PermissionError as e:
        logger.error(f"Insufficient permissions to write to file\n{e}")
    except IOError as e:
        logger.error(f"An error occurred while writing to file: {e}")
    return


# Fetch Metadata
def Split_data_table(stype, sname):
    data_list = []  # List to store dictionaries
    type_list = stype.rsplit("|")
    name_list = sname.rsplit("|")
    dict_meta = {}
    for y in range(len(type_list)):
        dict_meta[f"{type_list[y]}"] = f"{name_list[y]}"
    data_list.append(dict_meta)  # Append the dictionary to the list
    df = pd.DataFrame(data_list)  # Create DataFrame using the list of dictionaries
    column_order = ['serotype', 'host', 'country', 'collection_date']  # Define the column order
    new_dict = {}
    for column_name in column_order:
        if column_name in df.columns:
            new_dict[column_name] = df[column_name]
        else:
            new_dict[column_name] = pd.Series(dtype="object")  # Add an empty series if column is missing
    new_df = pd.DataFrame(new_dict)
    return new_df


def output_metadata(df1, df2, file):  # Write out the metadata dataframe into a csv file
    try:
        df = pd.concat([df1, df2], axis=1)
        df.to_csv(fr"{file}_metadata.csv", index=False, mode="a", header=False)
    except FileNotFoundError as e:
        logger.error(f"File not found\n{e}")
    except PermissionError as e:
        logger.error(f"Insufficient permissions to write to file\n{e}")
    except IOError as e:
        logger.error(f"An error occurred while writing to file: {e}")
    return


def add_headers(file):  # Add the appropriate headers to the csv file generated above
    try:
        df = pd.read_csv(fr"{file}_metadata.csv", header=None)
        df_c = pd.DataFrame(
            columns=["AccessionVersion", "Gi", "Length", "Moltype", "Strain", "Serotype", "Host", "Country",
                     "Collection_date"])
        df.columns = df_c.columns
        df.to_csv(fr"{file}_metadata.csv", index=False)
    except FileNotFoundError as e:
        logger.error(f"File not found\n{e}")
    except PermissionError as e:
        logger.error(f"Insufficient permissions to write to file\n{e}")
    except IOError as e:
        logger.error(f"An error occurred while writing to file: {e}")
    return


# E summary

def Esummary(record, file):  # Connects to the NCBI Entrez Direct API to fetch the sequence metadata
    batch_size = 9999
    count = int(record["Count"])
    query_key = record["QueryKey"]
    webenv = record["WebEnv"]
    for start in range(0, count, batch_size):
        end = min(start + batch_size - 1, count - 1)
        print(f"Downloading sequence metadata {start + 1} to {end + 1}")
        try:
            handle = Entrez.esummary(db="protein", version="2.0", query_key=query_key, webenv=webenv, retstart=start,
                                     retmax=int(end - start + 1))
            record_summary = Entrez.read(handle)
            handle.close()
            proteins = record_summary['DocumentSummarySet']
            # The record contains a dictionary and list elements with our data of interest
            # Parse through the dictionary to extract what we need
            for protein in proteins:
                acc_list = protein[0]
                gi_list = protein[2]
                length_list = protein[8]
                mol_list = protein[10]
                strain = protein[25]
                subtype = protein[16]
                subname = protein[17]
                df2 = Split_data_table(subtype, subname)
                df1 = pd.DataFrame({
                    'AccessionVersion': [acc_list],
                    'Gi': [gi_list],
                    'Length': [length_list],
                    'Moltype': [mol_list],
                    'Strain': [strain]
                })
                output_metadata(df1, df2, file)
        except urllib.error.HTTPError as e:
            # Handle HTTPError
            logger.error("An HTTP error occurred:", e)
        except urllib.error.URLError as e:
            # Handle URLError
            logger.error("A URL error occurred:", e)
        except TimeoutError:
            # Handle TimeoutError
            logger.error("The request timed out.")
        except Exception as e:
            # Handle other exceptions
            logger.error("An error occurred:", e)
    add_headers(file)
    return


# Start the time for running code
start_time = time.time()

# Input lists for search parameters, contains all the terms used by the user

taxid = []

db_fil = []

sub = []

gname = []

# Always tell NCBI who you are, as the NCBI API will not accept a connection request
if args.email is not None:
    if is_valid_email(args.email) is True:
        Entrez.email = args.email
    else:
        print("Please use a valid email or check the email input.")
        sys.exit()
else:
    print("Please input email id!")
    sys.exit()

# Iterating through input to check whether terms are present or not
# If terms are present then iterate through them as a list
if args.txid is not None:
    for Input in args.txid:
        for parameters in Input:
            taxid.append(parameters)

if args.db is not None:
    for Input in args.db:
        for parameters in Input:
            db_fil.append(parameters)

if args.sub is not None:
    for Input in args.sub:
        for parameters in Input:
            sub.append(parameters)

if args.gene is not None:
    for Input in args.gene:
        for parameters in Input:
            gname.append(parameters)

if args.min is not None:
    minlen = args.min
else:
    minlen = None

if args.max is not None:
    maxlen = args.max
else:
    maxlen = None

# Output file for writing sequences
output_file_seq = f"{args.F}" + ".aa"
output_file_meta = f"{args.meta}"

# Check for tax ID input which is the main search criteria for the program
if args.txid is not None:
    # Search Parameters
    organism = Organism(taxid)
    database_filter = Filter(db_fil)
    serotype = Subtype(sub)
    gene_name = Gene_name(gname)
    min_length = minlength(minlen)
    max_length = maxlength(maxlen, minlen)
    if args.full is False:
        terms = search_term(organism, database_filter, serotype, gene_name, min_length, max_length)
    else:
        terms = search_term(organism, database_filter, serotype, gene_name, min_length,
                            max_length) + " NOT (partial[Properties])"
    records = Esearch(terms)
    if int(records["Count"]) > 0:
        # Usage of program based on optionality of the arguments used by the user
        if args.F is None and args.meta is None:
            print(f"\n{records['Count']} protein sequences found in NCBI database.\n")
            print(f"Search query: {terms}\n")
            print(
                "To fetch the sequences found above please use the --F argument with the output filename along with the"
                " same search parameters.\n")
            print(
                "To get the associated metadata for the sequences found above please use the --meta argument with the "
                "output filename along with the"
                "same search parameters.\n")
        elif args.F is not None and args.meta is None:
            print("Fetching Sequences\n")
            Fetch_sequences_in_batches(records, output_file_seq)
        elif args.F is None and args.meta is not None:
            print("\nDownloading sequence metadata\n")
            Esummary(records, output_file_meta)
        elif args.F is not None and args.meta is not None:
            print("Fetching Sequences\n")
            Fetch_sequences_in_batches(records, output_file_seq)
            print("\nDownloading sequence metadata\n")
            Esummary(records, output_file_meta)
    else:
        print(f"{records['Count']} sequences found!\nPlease recheck your input")
else:
    print("Please input taxonomic id to search for sequences!")

# Get ending time
end_time = time.time()  # Get the current time after the code finishes
elapsed_time = round(end_time - start_time, ndigits=2)  # Calculate the elapsed time
print(f"Elapsed Time: {elapsed_time} seconds")
