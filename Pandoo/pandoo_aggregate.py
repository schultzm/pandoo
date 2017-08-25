'''
    Uses python3.

    This script:
        1. aggregates all genes found into a column
        2. aggregates all 'maybe found' genes into a dataframe column
        3. rules based on CARalert guidelines

    Copyright (C) 2017 Mark B Schultz
    https://github.com/schultzm/pandoo
    email: dr.mark.schultz@gmail.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.
    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import pandas as pd
import re

# GENE_CLASSES = {'CPE': ['kpc', 'ndm', 'vim', 'imi', 'imp', 'oxa', 'ges', 'sme'],
#                 '16S_methylase': ['armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'rmtG', 'rmtH', 'npmA'],
#                 }

def create_summary_df(dframe, db_prefix):
    '''
    Returns a dataframe, with one row and a column each for all genes found
    and genes to be confirmed.
    '''
    pregx = re.compile('_[^_]+$')
    for idx in dframe.index.values:
        genes = [item for item in
                 list(zip(dframe.loc[idx,dframe.columns.to_series().str \
                                 .contains(db_prefix.lower())].index.values,
                 dframe.loc[idx,dframe.columns.to_series().str. \
                        contains(db_prefix.lower())].values))
                 if isinstance(item[1], str)]
        genes_2 = {db_prefix+'all_genes':
                   ':'.join([pregx.sub('', gene[0].replace(db_prefix.lower(), ''))
                                 for gene in genes]),
                   db_prefix+'genes_to_be_confirmed':
                   ':'.join([pregx.sub('', gene[0].replace(db_prefix.lower(), ''))
                             for gene in genes if gene[1]=='maybe'])}
        df2 = pd.DataFrame([genes_2], index=[idx])
        return df2
