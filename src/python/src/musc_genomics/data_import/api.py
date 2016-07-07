
from musc_genomics.data_import import hugo
from musc_genomics import data
from musc_genomics.data_import import cgds


def get_gene_ids():
    """ Get Gene Names (HUGO) """
    key = ('materialized', 'hugo')
    if data.exists(*key):
        return data.load(*key)
    d = hugo.get_huge_genes()
    d = d[~d['Approved Name'].str.lower().str.contains('symbol withdrawn|entry withdrawn')]
    d = d['Approved Symbol'].unique()
    data.save(*key, d)
    return d


def get_gp_data(gene_ids, genetic_profile_id, key, batch_size=50):
    """ Get Genetic Profile Data (CGDS)
    :param gene_ids: List of Gene Ids to collect data for
    :param genetic_profile_id: Genetic profile name (eg "cellline_ccle_broad_mutations")
    :param key: Two-item key specifying data type and name for local storage
    :param batch_size: Number of gene ids to collect data for at once
    :return: DataFrame
    """
    d = data.load(*key, raise_on_absent=False)
    if d is not None:
        return d
    d = cgds.get_genetic_profile_data(
        cgds.CCLE_CASE_LIST_ID, genetic_profile_id,
        gene_ids, gene_id_batch_size=batch_size
    )
    data.save(*key, d)
    return d