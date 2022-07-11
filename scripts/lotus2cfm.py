#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This script prepares LOTUS structures for CFM-predict"""

import pandas
from requests_html import HTMLSession

# Get last version from Zenodo
s = HTMLSession()
r = s.get("https://doi.org/10.5281/zenodo.6378223")
urls = [l for l in r.html.absolute_links if "structure_metadata.tsv.gz" in l]
if not urls:
    raise ValueError(f"can't get the link for structure_metadata.tsv.gz")
print(urls[0])

last_lotus_structures = pandas.read_csv(
    filepath_or_buffer=urls[0],
    sep="\t",
    compression="gzip"
)

# Minimal cleaning and formatting
lotus_4cfm = last_lotus_structures[~last_lotus_structures['structureCleaned_smiles2D'].str.contains(
    "He")].query(
    "structureCleaned_exactMass > 50 & structureCleaned_exactMass < 2000").rename(
    columns={
        "structureCleaned_inchikey2D": "structure_inchikey_2D",
        "structureCleaned_smiles2D": "structure_smiles_2D"
    }
)[[
    "structure_inchikey_2D",
    "structure_smiles_2D"
]].drop_duplicates()

# Export
lotus_4cfm.to_csv(
    path_or_buf="smiles4cfm.txt",
    sep= " ",
    header=False,
    index=False
)
