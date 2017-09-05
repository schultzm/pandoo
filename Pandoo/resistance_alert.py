#this code is development
import pandas as pd
from collections import defaultdict


DFRAME = pd.read_table('ftp://ftp.ncbi.nlm.nih.gov/pathogen/betalactamases/Allele.tab', header=0, index_col='#allele name')
# print(DFRAME)

def gene_list(dframe, k_word, columnname):
    '''
    Parse DFRAME to retrieve gene list matching k_word and columnname criteria.
    '''
    return [dframe.loc[value, 'gene name'] for value in dframe.index.values if k_word in dframe.loc[value, columnname]]


METALLO_BETALACTAMASES = gene_list(DFRAME, 'metallo-beta-lactamase', 'curated gene product name')
CARBAPENEMASES = gene_list(DFRAME, 'carbapenem-hydrolyzing', 'curated gene product name')
FULL_CPE_GENES = METALLO_BETALACTAMASES + CARBAPENEMASES
ESBL_GENES = gene_list(DFRAME, 'extended-spectrum beta-lactamase', 'curated gene product name')
RMTASE_GENES = ['armA',
                'rmtA',
                'rmtB',
                'rmtC',
                'rmtD',
                'rmtE',
                'rmtF',
                'rmtG',
                'rmtH',
                'npmA']

COLISTIN_GENES = ['mcr-1']#.1',
#                   'mcr-1.2',
#                   'mcr-1.3',
#                   'mcr-1.4',
#                   'mcr-1.5',
#                   'mcr-1.6',
#                   'mcr-1.7',
#                   'mcr-1.8',
#                   'mcr-2.1',
#                   'mcr-3.1',
#                   'mcr-4.1']
# print(METALLO_BETALACTAMASES)
# print(CARBAPENEMASES)
# print('carbapenemases', FULL_CPE_GENES)
# print(ESBL_GENES)


#CARALERT NOTIFICATION
#Enterobacteriaceae with 16S rRNA methylase gene detected.

#CPE alert
#ENTEROBACTERIACEAE with full_cpe_gene
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

GENES_DICT = {'GENES': {'CPE_genes': FULL_CPE_GENES,
              'ESBL_genes': ESBL_GENES,
              '16S_Methyltransferase_genes': RMTASE_GENES,
              'Colistin_genes': COLISTIN_GENES},
              'Enterobacteriaceae': ENTEROBACTERIACEAE}
# print(GENES_DICT)

# print(genes_dict)
# class Isolate:
#     def __init__(name, genus_epithet, specific_epithet, amr_profile):
#         '''
#         Isolate name, genus, species, and resistance genes.
#         '''
#         self.name = name
#         self.genus_epithet = genus_epithet
#         self.specific_epithet = specific_epithet
#         self.amr_profile = amr_profile
#     def x():
#         pass
# 
    