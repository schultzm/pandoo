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

INFILE = ('/home/schultzm/jobs/mdu/AMR_ongoing/rerun_results/fromstdin_metadataAll_simplified.csv')
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


# def aggregate_yes(colname, rowname, string_target, df):
#     subset = df[df[rowname].str.contains(string_target)]
#     print(subset)
#     gene_yes = df.loc[colname, rowname]
    
if __name__ == '__main__':
    df = read_pandas_df(INFILE)
    extract_from = df[df['2017-01708'].str.contains('abricate_resfinder')]
    print(extract_from)
#     print(extract_from)
#     genes = df.loc['2017-01708', 'Sp_krkn_FinalCall']
#     print(genes)