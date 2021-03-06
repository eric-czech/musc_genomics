

import pandas as pd
import numpy as np


def get_feature_group(c):
    """
    Returns the group prefix assigned to a feature (e.g. CN:A1BG-AS1 -> CN)
    :param c: Feature to extract group from
    :return: Group name as a string
    """
    parts = c.split(':')
    if parts[0] in ['GE', 'CN']:
        return parts[0]
    if c.startswith('RES:VAL'):
        return ':'.join(parts[:3])
    return ':'.join(parts[:2])


def get_feature_type_counts(features):
    """
    Returns frequency of different feature types within a list of feature names

    Features are grouped into types based on their prefix (e.g. 'GE:', 'CN:', 'CL:GENDER', etc) and
    the frequency of those types are returned.  This mostly a convenience on analyzing metadata for
    very wide frames.

    :param features: List of features
    :return: Frequency of feature types as series
    """
    return pd.Series([get_feature_group(c) for c in features]).value_counts().sort_index()


def filter_na_values(d, ge_thresh=0., cn_thresh=0., mu_thresh=1., col_thresh=1.):
    """
    Filters a dataset containing genomic and drug response data based on various row and
    columnwise NA frequency thresholds.

    :param d: Data frame to be filtered
    :param ge_thresh: Row-wise NA threshold for gene expression features (in [0, 1])
    :param cn_thresh: Row-wise NA threshold for copy number features (in [0, 1])
    :param mu_thresh: Row-wise NA threshold for mutation features (in [0, 1])
    :param col_thresh: Col-wise NA threshold for all of the above (in [0, 1])
    :return: (d_filtered, na_summary) -- Filtered data frame and summary of NA mask values
    """
    # Create complete NA mask for entire data frame
    d_na = d.isnull()

    # ##### Mask Generation ##### #

    # Filter above mask to various sets of features and compute rowwise
    # masks of records over NA frequency thresholds
    na_ge = d_na.filter(regex='^GE:').mean(axis=1) > ge_thresh
    na_ge.name = 'GE'
    na_cn = d_na.filter(regex='^CN:').mean(axis=1) > cn_thresh
    na_cn.name = 'CN'
    na_mu = d_na.filter(regex='^MU:AAC:').mean(axis=1) >= mu_thresh
    na_mu.name = 'MU:AAC'

    # Filter global mask to similar set of features but instead compute column-wise
    # masks of features over NA frequency threshold
    na_col = d_na.filter(regex='^GE:|^CN:|^MU:AAC:|^CL:').mean(axis=0) >= col_thresh

    # Create row and column masks that will ultimately be used to filter the given dataset
    # *Note: these must by numpy boolean arrays, not indexed series, or alignment issues may occur
    row_mask = (na_ge | na_cn | na_mu).values
    col_mask = d.columns.isin(na_col[na_col].index.values)
    assert not hasattr(row_mask, 'index'), 'Row mask cannot be indexed (must be numpy array)'
    assert not hasattr(col_mask, 'index'), 'Col mask cannot be indexed (must be numpy array)'

    # ##### Mask Summarization ##### #

    # Compute mask tabulation for row-wise filters
    d_row = pd.concat([na_ge.value_counts(), na_cn.value_counts(), na_mu.value_counts()], axis=1)\
        .rename_axis('FEATURE_GROUP', axis=1)\
        .rename_axis('OVER_THRESHOLD', axis=0)\
        .fillna(0).astype(np.int64)

    # Compute mask tabulation for col-wise filters
    d_col = na_col\
        .rename('OVER_THRESHOLD').reset_index()\
        .assign(FEATURE_GROUP=lambda x: x['FEATURE'].apply(get_feature_group))\
        .groupby(['OVER_THRESHOLD', 'FEATURE_GROUP']).size().unstack()\
        .fillna(0).astype(np.int64)

    # Combine tabulations of number of rows/cols over NA threshold
    d_na = d_row.T.assign(Axis='rowwise')\
        .append(d_col.T.assign(Axis='colwise'))

    # Append tabulations of NA mask values over all rows and cols
    d_na = d_na.append(
        pd.Series(col_mask).value_counts().to_frame().T\
            .assign(Axis='colwise', FEATURE_GROUP='Overall')\
            .set_index('FEATURE_GROUP')
    )
    d_na = d_na.append(
        pd.Series(row_mask).value_counts().to_frame().T\
            .assign(Axis='rowwise', FEATURE_GROUP='Overall')\
            .set_index('FEATURE_GROUP')
    )

    # Add percentage calculations
    d_na = d_na\
        .rename_axis('OVER_THRESHOLD', axis=1)\
        .assign(PctTrue=lambda x: 100. * x[True] / (x[True] + x[False]))

    # ##### Mask Application ##### #

    return d.loc[~row_mask, ~col_mask], d_na


def fill_na_values(d):
    """
    Fills NA values in genetic feature data with appropriate values.

    Currently, this operation will only fill in missing values in mutation features and an error
    will be thrown if any other features outside of these contain NA values (these must be handled by the caller)

    :param d: Data frame to fill missing values within
    :return: (d_fill, fill_summary) -- Data frame containing filled values, number of cells filled by feature type
    """

    # Create complete NA mask for entire data frame
    d_na = d.isnull()

    # Get count by feature group of all fields containing NA values
    ct = d_na.sum(axis=0)
    na_cts = get_feature_type_counts(ct[ct > 0].index)

    # Get all mutation features containing NA values
    ct = d_na.filter(regex='^MU:').sum(axis=0)
    ct = ct[ct > 0]
    fill_summary = get_feature_type_counts(ct.index)
    fill_summary = pd.concat([get_feature_type_counts(d.columns), fill_summary], axis=1)
    fill_summary.columns = ['NumFeatures', 'NumFilled']
    c_fill = ct.index.values

    # Replace missing mutation feature values with 0
    d_fill = d.copy()
    d_fill[c_fill] = d_fill[c_fill].fillna(0)

    # Ensure there are no more NA values at this point.  If there are, then they
    # will be for non-mutation features and these types of imputations are not yet supported
    # (so thrown an error)
    if not np.all(d_fill.notnull()):
        msg = 'NA values were found in non-mutation features (and imputation for this is not yet supported).'  \
            'Summary of features with NA values grouped by type:\n{}'\
            .format(na_cts)
        raise AssertionError(msg)

    return d_fill, fill_summary


def remove_sparse_responses(d, threshold=.25):
    """
    Removes response values that don't meet a certain threshold for NA values
    :param d: Modeling data frame
    :param threshold: Maximum percentage of cases to allow NA values for
    :return:
    """
    d_res = d.filter(regex='^RES:')
    cts = d_res.isnull().mean()

    ax = cts.hist()
    ax.set_title('NA Count Distribution Across Responses')

    # Set fields to remove as those with at least {threshold}% null values
    c_rm = cts[cts >= threshold].index.values

    print('Removing {} response fields of {} due to sparsity'.format(len(c_rm), len(cts)))
    return d.drop(c_rm, axis=1)


def filter_responses(d, drug_names):
    """
    Removes response/drug features not matching the given drug names (any non-response features are left intact)
    :param d: Data frame containing response data (may contain more fields than just responses)
    :param drug_names: List of drugs to filter to
    :return: DataFrame like d with unnecessary response fields removed
    """
    pattern = '|'.join(drug_names)
    c_other = [c for c in d if not c.startswith('RES:')]
    c_drug = d.filter(regex='^RES:.*({}).*$'.format(pattern)).columns.tolist()
    return d[c_other + c_drug]
