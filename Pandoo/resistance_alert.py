import pandas as pd
from collections import defaultdict

FTP_URL = 'ftp://ftp.ncbi.nlm.nih.gov/pathogen/betalactamases/Allele.tab'
DFRAME = pd.read_table(FTP_URL, header=0, index_col='#allele name')

def gene_list(dframe, k_word, columnname):
    '''
    Parse DFRAME to retrieve gene list matching k_word and columnname criteria.
    '''
    return [dframe.loc[value, 'gene name'] for value in dframe.index.values
            if k_word in dframe.loc[value, columnname]]


METALLO_BETALACTAMASES = gene_list(DFRAME, 'metallo-beta-lactamase',
                                   'curated gene product name')
CARBAPENEMASES = gene_list(DFRAME, 'carbapenem-hydrolyzing',
                           'curated gene product name')
FULL_CPE_GENES = METALLO_BETALACTAMASES + CARBAPENEMASES
ESBL_GENES = gene_list(DFRAME, 'extended-spectrum beta-lactamase',
                       'curated gene product name')
RMTASE_GENES = ['armA',
                'npmA',
                'rmtA',
                'rmtA_1',
                'rmtB',
                'rmtB1',
                'rmtB2',
                'rmtB2_1',
                'rmtB3',
                'rmtB4',
                'rmtB_1',
                'rmtC',
                'rmtC_1',
                'rmtD',
                'rmtD1',
                'rmtD2',
                'rmtD2_1',
                'rmtD_1',
                'rmtE',
                'rmtE1',
                'rmtE2',
                'rmtE_1',
                'rmtF',
                'rmtF1',
                'rmtF2',
                'rmtG',
                'rmtG_1',
                'rmtH',
                'rmtH_1',
                'rmtf_1']

COLISTIN_GENES = ['mcr-1',
                  'mcr-1.1',
                  'mcr-1.2',
                  'mcr-1.2_1',
                  'mcr-1.3',
                  'mcr-1.3_1',
                  'mcr-1.4',
                  'mcr-1.4_1',
                  'mcr-1.5',
                  'mcr-1.5_1',
                  'mcr-1.6',
                  'mcr-1.6_1',
                  'mcr-1.7',
                  'mcr-1.7_1',
                  'mcr-1.8',
                  'mcr-1.9_1',
                  'mcr-1_1',
                  'mcr-2',
                  'mcr-2_1',
                  'mcr-3_1',
                  'mcr_1.8_1']

ENTEROBACTERIACEAE = ['Cedecea',
                      'Citrobacter',
                      'Morganella',
                      'Pantoea',
                      'Cronobacter',
                      'Edwardsiella',
                      'Enterobacter',
                      'Escherichia',
                      'Ewingella',
                      'Hafnia alvei',
                      'Hafnia paralvei',
                      'Klebsiella',
                      'Kluyvera',
                      'Leclercia',
                      'Plesiomonas',
                      'Proteus',
                      'Providencia',
                      'Raoultella',
                      'Salmonella',
                      'Serratia',
                      'Shigella',
                      'Yersinia']
# RULES:
# CARALERT NOTIFICATION
# Enterobacteriaceae with 16S_Methyltransferase_gene.
# CPE alert
# ENTEROBACTERIACEAE with CPE_gene

GENES_DICT = {'GENES': {'CPE_genes': FULL_CPE_GENES,
              'ESBL_genes': ESBL_GENES,
              '16S_Methyltransferase_genes': RMTASE_GENES,
              'Colistin_genes': COLISTIN_GENES},
              'Enterobacteriaceae': ENTEROBACTERIACEAE}
