

from musc_genomics import common
import logging
logger = logging.getLogger(__name__)


def save(datatype, dataset, d):
    return common.PROJ_MANAGER.save(datatype, dataset, d)


def exists(datatype, dataset):
    return common.PROJ_MANAGER.exists(datatype, dataset)


def load(datatype, dataset):
    return common.PROJ_MANAGER.load(datatype, dataset)
