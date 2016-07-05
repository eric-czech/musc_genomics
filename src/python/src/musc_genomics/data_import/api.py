
from musc_genomics.data_import import hugo
from musc_genomics import common

class DataAPI(object):

    def __init__(self, cache_dir=common.get_data_dir('')):
    def get_gene_ids():
        d = hugo.get_huge_genes()
        d = d[~d['Approved Name'].str.lower().str.contains('symbol withdrawn|entry withdrawn')]
        return d['Approved Symbol'].unique()

    def get_genetic_data(genetic_profile)

