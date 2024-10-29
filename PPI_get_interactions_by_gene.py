#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fetch interactions for a specific gene or gene list
"""

import requests
import json
import argparse
from core import config as cfg

parser = argparse.ArgumentParser(description="Fetch interactions for a gene list from a file.")
parser.add_argument("gene_file", help="Path to the file containing the gene list (one gene per line)")
parser.add_argument("output_file", help="Path to the file where results will be saved")
args = parser.parse_args()

# Read gene list from the provided file
with open(args.gene_file, "r") as file:
    content = file.read().strip() 
    geneList = "[" + content + "]"
print("Gene list loaded from file:", geneList)
request_url = cfg.BASE_URL + "/interactions"
#evidenceList=[]

# These parameters can be modified to match any search criteria following
# the rules outlined in the Wiki: https://wiki.thebiogrid.org/doku.php/biogridrest
params = {
    "accesskey": "129f271d7737a5a522d4e7785a333034",
    "format": "tab2",  # Return results in TAB2 format
    "geneList": "|".join(geneList),  # Must be | separated
    "searchNames": "true",  # Search against official names
    "includeInteractors": "true",  # Set to true to get any interaction involving EITHER gene, set to false to get interactions between genes
    "taxId": 9606,  # Limit to Homo sapiens
    #"evidenceList": "|".join(evidenceList),  # Exclude these two evidence types
    #"includeEvidence": "false",  # If false "evidenceList" is evidence to exclude, if true "evidenceList" is evidence to show
    "includeHeader": "true",
}

# Additional options to try, you can uncomment them as necessary
# params["start"] = 5 # Specify where to start fetching results from if > 10,000 results being returned
# params["max"] = 10 # Specify the number of results to return, max is 10,000
params["interSpeciesExcluded"] = "true" # true or false, If ‘true’, interactions with interactors from different species will be excluded (ex. no Human -> Mouse interactions)
params["selfInteractionsExcluded"] = "true" # true or false, If ‘true’, interactions with one interactor will be excluded. (ex. no STE11 -> STE11 interactions)
#params["searchIds"] = "true" # true or false, If ‘true’, ENTREZ_GENE, ORDERED LOCUS and SYSTEMATIC_NAME (orf) will be examined for a match with the geneList
# params["searchSynonyms"] = "false" # true or false, If ‘true’, SYNONYMS will be examined for a match with the geneList
# params["searchBiogridIds"] = "false" # true or false, If ‘true’, BIOGRID INTERNAL IDS will be examined for a match with the geneList
# params["excludeGenes"] = "false" # true or false, If 'true' the geneList becomes a list of genes to EXCLUDE rather than to INCLUDE
#params["includeInteractorInteractions"] = "true" # true or false, If ‘true’ interactions between the geneList’s first order interactors will be included. Ignored if includeInteractors is ‘false’ or if excludeGenes is set to ‘true’.
# params["htpThreshold"] = 50 # Any publication with more than this many interactions will be excluded
# params["throughputTag"] = "any" # any, low, high. If set to low, only `low throughput` interactions will be returned, if set to high, only `high throughput` interactions will be returned
# params["additionalIdentifierTypes"] = "SGD|FLYBASE|REFSEQ" # You can specify a | separated list of additional identifier types to search against (see get_identifier_types.py)

r = requests.get(request_url, params=params)
interactions = r.text

# Save interactions to the specified output file
with open(args.output_file, "w") as output_file:
    output_file.write(interactions)

# Optionally, print a message indicating success
print(f"Interactions saved to {args.output_file}")
