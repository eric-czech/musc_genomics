"""
Wrapper for HUGO gene nomenclature centralization service

See http://www.genenames.org/ for more details
"""

import pandas as pd


HUGO_URL = 'https://s3-us-west-2.amazonaws.com/eric.a.czech/Public/musc_genomics/hugo_genes.csv'


def get_huge_genes():
    """ Return HUGO Gene List (and related information)
    :return: DataFrame containing gene ID's as well as past names and related aliases
    """
    return pd.read_csv(HUGO_URL)
