
import pandas as pd
import numpy as np
from sklearn.preprocessing import Imputer
import logging
logger = logging.getLogger(__name__)


def get_modeling_data(d_cosmic, d_cgds, drug_names, drug_unit='IC_50'):
    """
    Returns merged view of COSMIC and CGDS data for the given drug names

    :param d_cosmic: COSMIC drug sensitivity data
    :param d_cgds: CGDS genetic data
    :param drug_names: Names of drugs to subset cosmic data to
    :param drug_unit: Unit of measurement for which sensitivity was assessed (e.g. 'IC_50')
    :return: Data frame with drug sensitivities prefixed by 'RES:VAL:<drug name>', genetic features
        with original names, and tumor id in index (data is aligned based on this)
    """

    # Subset COSMIC sensitivity data to these drugs and unit
    d_drug = d_cosmic[(d_cosmic['DrugName'].isin(drug_names)) & (d_cosmic['ValueUnit'] == drug_unit)].copy()

    # Ensure that at least one record was found
    assert not d_drug.empty, 'Failed to find COSMIC data for drug names = {}, value unit = {}"'\
        .format(drug_names, drug_unit)

    # Ensure that no cell line sensitivity numbers are duplicated
    assert d_drug.groupby(['CellLine', 'DrugName']).size().max() == 1

    # Conform COSMIC cell line IDs to those in CGDS data
    d_drug['CellLine'] = d_drug['CellLine'].str.replace('-', '')

    # Conform COSMIC metadata values to lower case representation
    d_drug['Tissue'] = d_drug['Tissue'].str.lower()
    d_drug['CancerType'] = d_drug['CancerType'].str.lower()

    # Pivot response values for each drug into columns with cell line / tumor id in index
    d_drug = d_drug.set_index(['CellLine', 'CancerType', 'Tissue', 'DrugName'])['Value']\
        .unstack().rename(columns=lambda c: 'VAL:' + c).reset_index()\
        .set_index('CellLine').rename(columns=lambda c: 'RES:' + c.upper())
    d_drug.columns.name = 'FEATURE'
    d_drug.index.name = 'TUMOR_ID'

    # Merge sensitivity data and genomic features
    d = pd.merge(d_cgds.reset_index(), d_drug.reset_index(), on='TUMOR_ID', how='inner').set_index('TUMOR_ID')
    return d


def prep_modeling_data(d, response=None, impute=True):

    if response is not None:
        mask = d[response].isnull()
        if np.any(mask):
            logger.warning(
                'Removing {} records due to missing response value (response = "{}")'\
                .format(mask.sum(), response)
            )
            d = d[~mask.values]

    # Isolate drug concentration fields not related to target response
    c_res = d.filter(regex='^RES:').columns.tolist()

    if response is not None:
        c_res = [c for c in c_res if c != response]

    # Ensure that all non-response fields are non-null
    assert np.all(d[[c for c in d if c not in c_res]].notnull()), \
        'Found NA values in non-response features (which is not allowed at this point)'

    # Apply mean value imputation to all non-target response fields
    d_res = d.copy()

    imp_summary = {}
    imputers = {}
    for c in c_res:
        imp_summary[c] = d_res[c].isnull().value_counts()
        imputers[c] = Imputer().fit(d_res[[c]])
        if impute:
            d_res[c] = imputers[c].transform(d_res[[c]])[:, 0]
    imp_summary = pd.DataFrame(imp_summary)
    d_res.response_feature_imputers = imputers

    # Ensure no values are NA at this point
    if impute:
        assert np.all(d_res.notnull())
    else:
        assert np.all(d_res[[c for c in d_res if c not in c_res]].notnull())

    return d_res, imp_summary
