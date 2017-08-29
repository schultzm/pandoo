#this code is development
import pandas as pd
# Read the ncbibetalactamase table from ncbi and store as a pandas dframe
dframe = pd.read_table('ftp://ftp.ncbi.nlm.nih.gov/pathogen/betalactamases/Allele.tab', header=0, index_col='#allele name')
print(dframe)

def gene_list(dframe, k_word, columnname):
    '''
    Parse dframe to retrieve gene list matching k_word and columnname criteria.
    '''
    return [dframe.loc[value, 'gene name'] for value in dframe.index.values if k_word in dframe.loc[value, columnname]]

metallo_betalactamases = gene_list(dframe, 'metallo-beta-lactamase', 'curated gene product name')
carbapenemases = gene_list(dframe, 'carbapenem-hydrolyzing', 'curated gene product name')
full_CPE_genes = metallo_betalactamases + carbapenemases
esbl_genes = gene_list(dframe, 'extended-spectrum beta-lactamase', 'curated gene product name')
rmtase_genes = ['armA',
                'rmtA',
                'rmtB',
                'rmtC',
                'rmtD',
                'rmtE',
                'rmtF',
                'rmtG',
                'rmtH',
                'npmA']

colistin_genes = ['mcr-1.1',
                  'mcr-1.2',
                  'mcr-1.3',
                  'mcr-1.4',
                  'mcr-1.5',
                  'mcr-1.6',
                  'mcr-1.7',
                  'mcr-1.8',
                  'mcr-2.1',
                  'mcr-3.1',
                  'mcr-4.1']
# print(metallo_betalactamases)
# print(carbapenemases)
print(full_CPE_genes)
print(esbl_genes)


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