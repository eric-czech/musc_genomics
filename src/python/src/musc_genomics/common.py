
from research.project import manager
import pandas as pd
import os
import re


ALPHANUM = re.compile('[^a-zA-Z0-9]')

DATA_DIR = os.path.expanduser('~/data/research/musc_genomics')
PROJ_DIR = os.path.expanduser('~/repos/musc_genomics/src/python')
PROJ_MANAGER = manager.ProjectManager(data_dir=DATA_DIR, proj_dir=PROJ_DIR)


SRC_COSMIC5 = 'cosmic5'
SRC_COSMIC6 = 'cosmic6'
SRC_CGDS = 'cgds'

TUMOR_ID_TT_OESOPH = 'TTOESOPH'
TUMOR_ID_TT_THYR = 'TTTHYR'
TUMOR_ID_KMH_BLOOD = 'KMH2BLOOD'
TUMOR_ID_KMH_THYR = 'KMH2THYR'


def _normalize_tumor_id(id, source):
    assert not pd.isnull(id), 'Tumor ID cannot be null'

    # Deal with specific cases for tumor id "TT"
    if source == SRC_COSMIC6 and id.upper().strip() == 'T-T':
        id = TUMOR_ID_TT_OESOPH
    if source == SRC_COSMIC6 and id.upper().strip() == 'TT':
        id = TUMOR_ID_TT_THYR

    assert id.upper().strip() != 'TT', \
        'Resolution of annoying tumor ID "TT" not supported yet for source "{}"'.format(source)

    # Deal with specific cases for tumor id "KMH2"
    if source == SRC_COSMIC6 and id.upper().strip() == 'KM-H2':
        id = TUMOR_ID_KMH_BLOOD
    if source == SRC_COSMIC6 and id.upper().strip() == 'KMH-2':
        id = TUMOR_ID_KMH_THYR
    if source == SRC_CGDS and id.upper().strip() == 'KMH2':
        id = TUMOR_ID_KMH_BLOOD

    # Remove all non-alphanumeric characters from ID
    return ALPHANUM.sub('', id).upper()


def normalize_tumor_ids(ids, source, verify_integrity=True):
    """
    Normalizes Tumor ID values by replacing non-alphanumeric characters within them.

    :param ids: List of ids to normalize
    :param source: Source from which IDs came
    :param verify_integrity: Flag indicating whether or not an error should be thrown if the normalize list
        of IDs contains duplicates
    :return:
    """
    # Make sure the source of these identifiers is one that has been seen before
    assert source in [SRC_COSMIC6, SRC_CGDS], 'Tumor ID source "{}" is not yet supported'.format(source)

    ids = pd.Series(list(ids))
    res = pd.Series(ids.apply(_normalize_tumor_id, source=source).values, index=ids.values)
    cts = res.value_counts()
    assert not verify_integrity or cts.max() == 1, \
        'The following tumor ids resulted in duplicates after normalization:\n{}'\
        .format(res[res.isin(cts[cts > 1].index.values)])
    return res.values


# def get_data_dir(data_type):
#     path = os.path.join(DATA_DIR, data_type)
#     if not os.path.exists(path):
#         os.mkdir(path)
#     return path

