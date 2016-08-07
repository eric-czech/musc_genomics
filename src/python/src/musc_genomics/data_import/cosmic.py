
import pandas as pd
import os
from urllib.request import urlretrieve
from musc_genomics import common

FTP_URL = 'ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_manova_input_w5.csv'


def get_raw_cosmic_data(local_file):
    if not os.path.exists(local_file):
        urlretrieve(FTP_URL, local_file)
    return pd.read_csv(local_file, header=[0, 1, 2, 3, 4, 5])


def get_cosmic_rel6_data(filename):
    path = os.path.join(common.DATA_DIR, 'raw', 'cosmic_6.0', filename)
    return pd.read_excel(path, sheetname=None)
