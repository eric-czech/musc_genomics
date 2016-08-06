
import pandas as pd
from musc_genomics import data
from musc_genomics.data_modeling import modeling, features


def get_modeling_data_01(encoding_set):
    """
    Returns modeling data with preset operations for preparation

    :param encoding_set: Name of feature encoding set to load
    :return: Modeling data (no NA values)
    """
    d = data.load('features', 'encode_{}'.format(encoding_set))

    d, na_summary = features.filter_na_values(d)
    d_res = d.filter(regex='^RES:')
    d_feat = d[[c for c in d if c not in d_res]]
    d_feat, fill_summary = features.fill_na_values(d_feat)

    d_fill = pd.concat([d_feat, d_res], axis=1)
    assert d_fill.shape == d.shape

    d_fill, imp_summary = modeling.prep_modeling_data(d_fill)

    return d_fill, na_summary, fill_summary, imp_summary


