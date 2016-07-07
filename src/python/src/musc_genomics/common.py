
from research.project import manager
import os

DATA_DIR = os.path.expanduser('~/data/research/musc_genomics')
PROJ_DIR = os.path.expanduser('~/repos/musc_genomics/src/python')
PROJ_MANAGER = manager.ProjectManager(data_dir=DATA_DIR, proj_dir=PROJ_DIR)


# def get_data_dir(data_type):
#     path = os.path.join(DATA_DIR, data_type)
#     if not os.path.exists(path):
#         os.mkdir(path)
#     return path
