"""
Wrapper API for CGDS (cbioportal) Web Services

See http://www.cbioportal.org/web_api.jsp for more details.
"""
import pandas as pd
from py_utils.collection_utils import to_batches
import functools
import logging
from urllib import parse
from IPython.core.debugger import Tracer
logger = logging.getLogger(__name__)

BASE_URL = 'http://www.cbioportal.org/public-portal/webservice.do?'


def _to_url(cmd, data=None):
    url = '{}cmd={}'.format(BASE_URL, cmd)
    if data:
        url += '&' + '&'.join(['{}={}'.format(k, v) for k, v in data.items()])
    return url


def _get(cmd, data=None):
    url = _to_url(cmd, data)
    return pd.read_csv(url, sep='\t', comment='#')


def _is_id(idv):
    return isinstance(idv, str) and len(idv) > 0


def _is_iterable(seq, check_empty=True):
    try:
        _ = (e for e in seq)
        if check_empty and len(seq) == 0:
            return False
        return True
    except TypeError:
        return False


def get_cancer_studies():
    return _get('getCancerStudies')


def get_genetic_profiles(cancer_study_id):
    assert _is_id(cancer_study_id), \
        'Cancer Study ID must be a non-empty string (e.g. "cellline_ccle_broad_all")'
    return _get('getGeneticProfiles', {'cancer_study_id': cancer_study_id})


def get_case_lists(cancer_study_id):
    assert _is_id(cancer_study_id), \
        'Cancer Study ID must be a non-empty string (e.g. "cellline_ccle_broad_all")'
    return _get('getCaseLists', {'cancer_study_id': cancer_study_id})


def get_genetic_profile_data(case_list_id, genetic_profile_id, gene_ids,
                             gene_id_batch_size=50, print_progress=True):
    assert _is_id(case_list_id), \
        'Case list ID must be a non-empty string (e.g. "cellline_ccle_broad_all")'
    assert _is_id(genetic_profile_id), \
        'Genetic Profile ID must be a non-empty string (e.g. "cellline_ccle_broad_log2CNA")'
    assert _is_iterable(gene_ids), 'Gene IDs must be iterable and non-empty'

    # Split given gene list into batches and combine results from each batch
    res = None
    gene_id_batches = list(to_batches(gene_ids, gene_id_batch_size))
    for i, gene_ids in enumerate(gene_id_batches):
        if print_progress:
            logger.info('Processing batch {} of {}'.format(i + 1, len(gene_id_batches)))
        data = {
            'id_type': 'gene_symbol',
            'case_set_id': case_list_id,
            'genetic_profile_id': genetic_profile_id,
            'gene_list': ','.join(gene_ids)
        }
        # Tracer()()
        part = _get('getProfileData', data)
        if res is None:
            res = part
        else:
            res = res.append(part)

    # Return results as single, combined data frame
    return res
