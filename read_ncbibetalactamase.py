#this code is development
import pandas as pd
# Read the ncbibetalactamase table from ncbi and store as a pandas dframe
dframe = pd.read_table('ftp://ftp.ncbi.nlm.nih.gov/pathogen/betalactamases/Allele.tab', header=0, index_col='#allele name')


def gene_list(dframe, k_word, columnname):
    '''
    Parse dframe to retrieve gene list matching k_word and columnname criteria.
    '''
    return [value for value in dframe.index.values if k_word in dframe.loc[value, columnname]]

metallo_betalactamases = gene_list(dframe, 'metallo-beta-lactamase', 'curated gene product name')
carbapenemases = gene_list(dframe, 'carbapenem-hydrolyzing', 'curated gene product name')
full_CPE_genes = metallo_betalactamases + carbapenemases
# print(metallo_betalactamases)
# print(carbapenemases)
print(full_CPE_genes)


ENTEROBACTERIACEAE = ["Cedecea",
                      "Citrobacter",
                      "Morganella",
                      "Pantoea",
                      "Cronobacter",
                      "Edwardsiella",
                      "Enterobacter",
                      "Escherichia",
                      "Ewingella",
                      "Hafnia alvei",
                      "Klebsiella",
                      "Kluyvera",
                      "Leclercia",
                      "Plesiomonas",
                      "Proteus",
                      "Providencia",
                      "Raoultella",
                      "Salmonella",
                      "Serratia",
                      "Shigella",
                      "Yersinia"]