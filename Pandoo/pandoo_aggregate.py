#!/usr/bin/env python3

'''
    Uses python3.

    This script:
        1. summarises a pandoo output table for LIMS
        2. rules based on CARalert guidelines

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
from pathlib import Path
from io import StringIO
import sys
# import numpy as np
import re

INFILE = ('/home/schultzm/jobs/mdu/AMR_ongoing/rerun_20170808_results/fromstdin_metadataAll_simplified.csv')
PANDAS_INDEX_LABEL = 'Isolate'

def read_pandas_df(infile):
    '''
    Will read in the csv and return a Pandas dataframe.
    '''
    contents = Path(infile).read_text()
    df1 = pd.read_csv(StringIO(contents), header=0,
                      converters={PANDAS_INDEX_LABEL: str}) \
                                 .set_index(PANDAS_INDEX_LABEL)
    return df1

def create_summary_df(df):
    db_prefix = 'abricate_resfinder_'
    p = re.compile('_[^_]+$')
    for idx in df.index.values:
        genes = [item for item in
                 list(zip(df.loc[idx,df.columns.to_series().str \
                                 .contains(db_prefix)].index.values,
                 df.loc[idx,df.columns.to_series().str. \
                        contains(db_prefix)].values))
                 if isinstance(item[1], str)]
        
        genes_2 = {'All genes': [p.sub('', gene[0].replace(db_prefix, ''))
                                 for gene in genes],
                   'Genes to be confirmed': [p.sub('', gene[0] \
                                                   .replace(db_prefix, ''))
                                             for gene in genes if
                                             gene[1]=='maybe']}
        df2 = pd.DataFrame([genes_2], index=[idx])
        return df2
    
    
if __name__ == '__main__':
    out_df = create_summary_df(read_pandas_df(INFILE))
    print(out_df)