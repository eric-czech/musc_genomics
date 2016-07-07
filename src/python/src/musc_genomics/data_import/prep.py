
import pandas as pd
import numpy as np


def prep_raw_cgds_genetic_data(d):

    # Drop completely duplicated records
    d = d[~d.duplicated()]

    # Ensure that only one GENE_ID (an int at this point) is present for every "COMMON"
    # gene symbol (a string identifier)
    assert d.groupby('COMMON')['GENE_ID'].apply(lambda x: len(x.unique())).max() == 1

    # Drop the numeric GENE_ID and redefine as the string identifier
    d = d.drop('GENE_ID', axis=1).rename(columns={'COMMON': 'GENE_ID'}).copy()

    # Ensure that there is now only one record per GENE_ID
    assert d.groupby('GENE_ID').size().max() == 1

    # Melt to long format, with all tumor + gene combinations in rows
    d = pd.melt(d, id_vars='GENE_ID', var_name='TUMOR_ID', value_name='VALUE')

    # Remove any rows with NA values to get to a more efficient, sparse representation
    d = d[d['VALUE'].notnull()]

    # Ensure no tumor + gene combinations are duplicated
    assert np.all(~d[['TUMOR_ID', 'GENE_ID']].duplicated())

    return d.reset_index(drop=True)


def split_mutation_value(v):
    assert len(v) == 1
    mu = v.iloc[0].split(',')
    res = pd.Series(np.repeat(1, len(mu)), index=mu)
    res.index.name = 'MUTATION'
    return res


def get_impact_score(v):
    if pd.isnull(v):
        return None
    v = v.upper().strip()
    if v == 'L':
        return 'low'
    elif v == 'M':
        return 'med'
    elif v == 'N':
        return 'neutral'
    elif v == 'H':
        return 'high'
    return None


def add_mutation_metadata(d_mu, d_mm):
    # Project out relevant clinical data and rename for consistency with other data sets
    d_mm_ft = d_mm[['case_id', 'amino_acid_change', 'mutation_type', 'functional_impact_score']].copy()
    d_mm_ft.columns = ['TUMOR_ID', 'MUTATION', 'MUTATION_TYPE', 'IMPACT_SCORE']

    # Map the impact score value to a more informative one
    d_mm_ft['IMPACT_SCORE'] = d_mm_ft['IMPACT_SCORE'].apply(get_impact_score)

    # Fill missing values and convert all clinical feature values to lower case
    d_mm_ft[['MUTATION_TYPE', 'IMPACT_SCORE']] = d_mm_ft[['MUTATION_TYPE', 'IMPACT_SCORE']]\
        .fillna('unknown').applymap(lambda x: x.lower())

    # Merge mutation data on clinical data
    d_mu_all = pd.merge(d_mu, d_mm_ft, how='left', on=['TUMOR_ID', 'MUTATION'])

    # Ensure that all mutations were able to joined on metadata (this should contain all
    # mutations present and at TOW, no records were ever not joinable)
    assert np.all(d_mu_all.notnull()), 'Found null values in joined mutation features'

    return d_mu_all


def get_gene_mutation_features(g):
    gene = g['GENE_ID'].iloc[0]

    # Create either a single item list containing the one mutation associated with this
    # gene, or a list containing all such mutations as well as one extra item equal
    # to the combination of each mutation
    mu = list(g['MUTATION'].unique())
    if len(mu) > 1:
        mu = mu + [','.join(sorted(mu))]

    # Return series with gene + mutation combo in index and 1 as value
    res = pd.Series(np.repeat(1, len(mu)), index=['AAC:' + gene + ':' + m for m in mu])

    # Add features indicating frequency of each mutation type
    # (note that the gene id itself is not associated with these counts so
    # that they may become specific to the tumor id / case, but not the individual genes)
    ct = g['MUTATION_TYPE'].value_counts()
    ct.index = ['TYP:' + x for x in ct.index]
    res = res.append(ct)

    # Similarly, add the frequency of each impact score value
    ct = g['IMPACT_SCORE'].value_counts()
    ct.index = ['IMP:' + x for x in ct.index]
    res = res.append(ct)

    # Rename the resulting axes
    res.index.name = 'FEATURE'
    res.name = 'VALUE'
    return res
