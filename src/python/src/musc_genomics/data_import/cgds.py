"""
Wrapper API for CGDS (cbioportal) Web Services

See http://www.cbioportal.org/web_api.jsp for more details.
"""
import pandas as pd
from py_utils.collection_utils import to_batches
from urllib.error import HTTPError
import logging
import time
from IPython.core.debugger import Tracer
logger = logging.getLogger(__name__)

CCLE_CANCER_STUDY = 'cellline_ccle_broad'
CCLE_CASE_LIST_ID = 'cellline_ccle_broad_all'
GEN_PROF_COPY_NUMBER = 'cellline_ccle_broad_log2CNA'
GEN_PROF_GENE_EXPRESSION = 'cellline_ccle_broad_mrna_median_Zscores'
GEN_PROF_MUTATION = 'cellline_ccle_broad_mutations'

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


def get_cancer_types():
    return _get('getTypesOfCancer')


def get_genetic_profiles(cancer_study_id):
    assert _is_id(cancer_study_id), \
        'Cancer Study ID must be a non-empty string (e.g. "cellline_ccle_broad")'
    return _get('getGeneticProfiles', {'cancer_study_id': cancer_study_id})


def get_case_lists(cancer_study_id):
    assert _is_id(cancer_study_id), \
        'Cancer Study ID must be a non-empty string (e.g. "cellline_ccle_broad")'
    return _get('getCaseLists', {'cancer_study_id': cancer_study_id})


def get_clinical_data(case_list_id):
    assert _is_id(case_list_id), \
        'Case list ID must be a non-empty string (e.g. "cellline_ccle_broad_all")'
    return _get('getClinicalData', {'case_set_id': case_list_id})


def _get_batch_batch_result(gene_ids, batch_size, cmd, args, print_progress=True):
    res = None
    gene_id_batches = list(to_batches(gene_ids, batch_size))
    n = len(gene_id_batches)
    m = max(int(n / 10.), 1)
    for i, gene_ids in enumerate(gene_id_batches):
        if print_progress and i % m == 0:
            logger.info('Processing batch {} of {}'.format(i + 1, n))
        data = dict(args)
        data['gene_list'] = ','.join(gene_ids)
        # Tracer()()

        attempts = 0
        while True:
            attempts += 1
            try:
                part = _get(cmd, data)
                break
            except:
                if attempts > 5:
                    raise
                logger.warn('An http error occurred.  Will try again in 30 seconds ...')
                time.sleep(30)

        if res is None:
            res = part
        else:
            res = res.append(part)
    return res


def get_genetic_profile_data(case_list_id, genetic_profile_id, gene_ids,
                             gene_id_batch_size=50, print_progress=True):
    assert _is_id(case_list_id), \
        'Case list ID must be a non-empty string (e.g. "cellline_ccle_broad_all")'
    assert _is_id(genetic_profile_id), \
        'Genetic Profile ID must be a non-empty string (e.g. "cellline_ccle_broad_log2CNA")'
    assert _is_iterable(gene_ids), 'Gene IDs must be iterable and non-empty'

    # Split given gene list into batches and combine results from each batch
    args = {
        'id_type': 'gene_symbol',
        'case_set_id': case_list_id,
        'genetic_profile_id': genetic_profile_id
    }
    return _get_batch_batch_result(gene_ids, gene_id_batch_size, 'getProfileData', args, print_progress=print_progress)


def get_mutation_data(case_list_id, genetic_profile_id, gene_ids,
                      gene_id_batch_size=50, print_progress=True):
    assert _is_id(case_list_id), \
        'Case list ID must be a non-empty string (e.g. "cellline_ccle_broad_all")'
    assert _is_id(genetic_profile_id), \
        'Genetic Profile ID must be a non-empty string (e.g. "cellline_ccle_broad_log2CNA")'
    assert _is_iterable(gene_ids), 'Gene IDs must be iterable and non-empty'

    args = {
        'id_type': 'gene_symbol',
        'case_set_id': case_list_id,
        'genetic_profile_id': genetic_profile_id
    }
    return _get_batch_batch_result(gene_ids, gene_id_batch_size, 'getMutationData', args, print_progress=print_progress)


def get_protein_array_info(cancer_study_id, array_type, gene_ids, gene_id_batch_size=500, print_progress=True):
    assert _is_id(cancer_study_id), \
        'Cancer Study ID must be a non-empty string (e.g. "cellline_ccle_broad")'
    array_types = ['protein_level', 'phosphorylation']
    assert array_type in array_types, \
        'Protein array type must be one of {}'.format(array_types)
    assert _is_iterable(gene_ids), 'Gene IDs must be iterable and non-empty'

    args = {
        'id_type': 'gene_symbol',
        'cancer_study_id': cancer_study_id,
        'protein_array_type': array_type
    }
    return _get_batch_batch_result(gene_ids, gene_id_batch_size, 'getProteinArrayInfo',
                                   args, print_progress=print_progress)


def get_protein_array_data(case_list_id, array_info):
    assert _is_id(case_list_id), \
        'Case list ID must be a non-empty string (e.g. "cellline_ccle_broad_all")'
    array_types = ['0', '1']
    assert array_info in array_types, \
        'Protein array info flag must be one of {}'.format(array_types)

    data = {
        'case_set_id': case_list_id,
        'array_info': array_info
    }
    return _get('getProteinArrayData', data)
