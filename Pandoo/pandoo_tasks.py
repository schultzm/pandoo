#!/usr/bin/env python3

'''
    Uses python3.

    This script contains the task functions for pando

    Copyright (C) 2017 Mark B Schultz
    https://github.com/schultzm/
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


import os
import math
from multiprocessing import cpu_count
import shlex
from subprocess import Popen, PIPE
import sys
from ete3 import Tree
import numpy as np
import pandas as pd


def calc_threads(n_isos, ncores):
    '''
    Based on the number of isolates, calculate the number of threads for
    multithread jobs.
    '''
    header = ('isolates	threads	sec_per_job	ncores	concurrent_jobs	' +
              'blocks	Total_sec').split('\t')
    if cpu_count() >= 4:
        data_rows = [('\t'+str(i)+'\t'*5).split('\t')
                     for i in range(4, cpu_count()+1)]

        def block_calc(iso, concjobs):
            '''
            Calculate the number of analysis 'lots'.
            '''
            if iso < concjobs:
                return int(1)
            else:
                return int(math.ceil(1/(concjobs/iso)))

        df1 = pd.DataFrame(data_rows, columns=header)
        df1 = df1[df1.columns[:]].apply(pd.to_numeric)
        df1['isolates'] = n_isos
        df1['ncores'] = ncores
        df1['concurrent_jobs'] = df1['ncores']/df1['threads']
        # The following formula was calculated using a power curve.
        df1['sec_per_job'] = 441.61*df1['threads']**-0.591
        df1['blocks'] = np.vectorize(block_calc)(df1['isolates'],
                                                 df1['concurrent_jobs'])
        df1['Total_sec'] = df1['blocks']*df1['sec_per_job']
        df1 = df1[df1.threads <= ncores]
        return int(df1.loc[df1['Total_sec'].idxmin()].threads)
    else:
        return cpu_count()


# This value sets the name of first column in all output tables.
PANDAS_INDEX_LABEL = 'Isolate'


def get_paths(infile):
    '''
    Reads in the tab-delimited file containing:
    Column 1 - isolate IDs
    Column 2 - paths to contigs.fa
    Column 3 - paths to read1.fq.gz
    Column 4 - paths to read2.fq.gz
    '''
    while True:
        try:
            df1 = pd.read_csv(infile, header=None, sep='\t', index_col=0)
            df1.columns = ['pathContigs', 'pathReads1', 'pathReads2']
            df1.index.rename(PANDAS_INDEX_LABEL, inplace=True)
            # Don't convert the index to str.  Notify user though.
            # The user then will have to
            # go into their .tab file and put an X in front of the offending
            # index (e.g., 14E7 should be X14E7 or '14E7', NOT 14E7.).
            # Don't do  df1.index = df1.index.map(str)
            return df1
        except ValueError:
            print(infile+' contains no data')
            sys.exit()


def write_pandas_df(outfile, dframe):
    '''
    Write the Pandas dataframe to file.
    '''
    with open(outfile, 'w') as output:
        dframe.to_csv(output, sep=',', mode='w', index=True,
                      index_label=PANDAS_INDEX_LABEL)


def read_pandas_df(infile):
    '''
    Will read in the csv and return a Pandas dataframe.
    '''
    # Old code for converting to string:
    # converters={'Isolate': lambda x: str(x)}
    df1 = pd.read_csv(infile, index_col=0, header=0, dtype=str)
    return df1


def create_pandas_df(dictionary, isolate):
    '''
    Creates a Pandas dataframe from a dictionary of results.
    '''
    return pd.DataFrame([dictionary], index=[isolate])


def run_abricate(infile, outfile, outfile_simple, isolate, dbase, cutoff,
                 identity):
    '''
    Runs abricate on the infile.
    '''
    if len(infile) == 0:
        ab_results = {}
        ab_results_simplified = {}
        abricate_result = create_pandas_df(ab_results, isolate)
        abricate_result = create_pandas_df(ab_results_simplified, isolate)
    if len(infile) == 1:
        infile = infile[0]
        # The full path to the database comes in under 'dbase' variable, which
        # allows the user to point abricate to their own database.  When the
        # default dbs are being used, the end of the dbase path will be 
        # 'default'.  Here we extract 'default' as a keyword which allows
        # correct building of the abricate run command.
        dbases_path = os.path.split(dbase[1])
        if os.path.split(dbase[1])[1] == 'default':
            print('Using pre-packaged abricate database '+dbase[0])
            cmd = 'abricate --db '+dbase[0] +\
                  ' --minid '+str(identity)+' '+infile+' > '+outfile
        else:
            if not os.path.exists(dbase[1]):
                sys.exit('Print unable to find '+dbase[1])
            else:
                print('Using custom db '+ dbase[0]+' at '+dbase[1])
                cmd = 'abricate --db '+dbase[0]+' --datadir '+dbase[1] +\
                      ' --minid '+str(identity)+' '+infile+' > '+outfile
        os.system(cmd)
        ab_data = pd.read_table(outfile, sep='\t', header=0)
        ab_results_df_list = []
        ab_results_simplified = {}
        for i in ab_data.index.values:
            ab_results = {}
            # Generate the simplified dict
            simplifiedtable_key = 'abricate_'+dbase[0] +\
                                  '_'+ab_data.loc[i, 'GENE']
            if ab_data.loc[i, '%COVERAGE'] >= \
            cutoff and ab_data.loc[i, '%IDENTITY'] == 100:
                if simplifiedtable_key not in ab_results_simplified.items():
                    ab_results_simplified[simplifiedtable_key] = 'yes'
                else:
                    ab_results_simplified[simplifiedtable_key] = 'maybe'
            if ab_data.loc[i, '%COVERAGE'] < cutoff or \
            ab_data.loc[i, '%IDENTITY'] < 100:
                ab_results_simplified[simplifiedtable_key] = 'maybe'

            # Generate the complex (close to raw data) dict
            # Bind the coverage, ID and gaps into a single cell value.
            # Bind the dfs so that if there is a duplicate gene, there will be
            # Genecopy.1 and genecopy.2 in the final table.
            ab_results['abricate_' +
                       dbase[0]+'_' +
                       ab_data.loc[i, 'GENE']] = 'COV_' +\
                      str(ab_data.loc[i, '%COVERAGE']) +\
                      '_ID_'+str(ab_data.loc[i, '%IDENTITY'])
            ab_results_df_list.append(create_pandas_df(ab_results, isolate))
        if len(ab_results_df_list) > 0:
            abricate_result = pd.concat(ab_results_df_list, axis=1)
        else:
            ab_results = {}
            abricate_result = create_pandas_df(ab_results, isolate)

    def get_abricate_version():
        '''
        Get the Abricate software version.
        '''
        args = shlex.split('abricate -v')
        proc = Popen(args, stdout=PIPE)
        version = proc.stdout.read().decode('UTF-8').rstrip().split('\n')[0]
        return {'softwareAbricateVersion_'+dbase[0]: version,
                'softwareAbricateDB_'+dbase[0]: dbase[1]+'/'+dbase[0],
                'softwareAbricateSettings_'+dbase[0]: 'COV'+str(cutoff) +\
                '_ID'+str(identity)}

    # Bind the version and the results dataframes.
    # NB: writing of the outfile below will overwrite the full 
    # abricate results table as it was written when redirected from stdout.
    # write both simple ('yes', 'maybe') and complex ('near raw') dataframes.
    version_df = create_pandas_df(get_abricate_version(), isolate)
    abricate_df = pd.concat([abricate_result, version_df], axis=1)
    write_pandas_df(outfile, abricate_df)
    abricate_df_simple = pd.concat([create_pandas_df(ab_results_simplified,
                                                     isolate), version_df],
                                   axis=1)
    write_pandas_df(outfile_simple, abricate_df_simple)



def seqtk_version():
    '''
    Get the version of SeqTK,
    '''
    proc = Popen(shlex.split('seqtk'), stderr=PIPE)
    return proc.stderr.read().decode('UTF-8').rstrip().split('\n')[2]


def run_seqtk_comp(infile, outfile, isolate):
    '''
    Runs SeqTK on the contigs.fa (assemblies).  Gathers the fasta metrics.
    '''
    if len(infile) == 0:
        metrics = {}
    else:
        os.system('seqtk comp '+''.join(infile)+' > '+outfile)
        df1 = pd.read_csv(outfile, header=None, index_col=0, sep='\t',
                          names=['chr', 'length', '#A', '#C', '#G', '#T', '#2',
                                 '#3', '#4', '#CpG', '#tv', '#ts', '#CpG-ts'])
        contig_lengths = df1['length'].tolist()
        contig_lengths.sort(), contig_lengths.reverse()
        contig_lengths_prime = [[i]*i for i in contig_lengths]
        metrics = {}
        pfx = 'metricsContigs_'
        metrics[pfx+'N50'] = int(np.median([i for j in contig_lengths_prime
                                            for i in j]))
        bps = sum(contig_lengths)
        nns = sum(df1['#4'].tolist())
        metrics[pfx+'Ns'] = nns
        metrics[pfx+'no'] = len(df1.index.values)
        metrics[pfx+'bp'] = bps
        metrics[pfx+'avg'] = bps/len(contig_lengths)
        metrics[pfx+'Max'] = max(contig_lengths)
        metrics[pfx+'Min'] = min(contig_lengths)
        metrics[pfx+'ok'] = bps - nns
    metrics['sotwareSeqTKversion_comp'] = seqtk_version()
    seqtk_comp_df = create_pandas_df(metrics, isolate)
    # NB: this will overwrite the outfile that was read in at the first step.
    write_pandas_df(outfile, seqtk_comp_df)


def run_seqtk_fqchk(infiles, outfile, isolate):
    '''
    Runs SeqTK fqchk on the reads.  Gathers fastq metrics.
    '''
    pfx = 'metricsReads_'
    if len(infiles) == 0:
        metrics = {}
    else:
        os.system('seqtk fqchk '+' '.join(infiles)+' > '+outfile)
        df_all = pd.read_csv(outfile, index_col=0, sep='\t',
                             skiprows=1)
        with open(outfile, 'r') as outf:
            metrics = outf.readline()
            metrics = dict([j.split(':') for j in [i.replace(' ', '') for i in
                                                   metrics.split(';')[0:3]]])
            for key, value in list(metrics.items()):
                metrics[pfx+key] = metrics.pop(key)
        metrics[pfx+'AvgQual'] = round(df_all.loc['ALL', 'avgQ'], 2)
        metrics[pfx+'GeeCee'] = round(df_all.loc['ALL', '%C'] +
                                df_all.loc['ALL', '%G'], 2)
        metrics[pfx+'Yield'] = int(df_all.loc['ALL', '#bases']*2)

        # Count the number of reads in the infiles.
        cmd = 'zcat '+' '.join(infiles)
        cmd2 = 'wc -l'
        args = shlex.split(cmd)
        args2 = shlex.split(cmd2)
        proc = Popen(args, stdout=PIPE, stderr=PIPE)
        proc2 = Popen(args2, stdin=proc.stdout, stdout=PIPE)
        output = proc2.stdout.read().decode('UTF-8')
        metrics[pfx+'nReads'] = int(int(output)/4)

        # Get the mode read length.
        n_bases = df_all['#bases'].values[1:]
        lengths = []
        pos = 0
        while pos < len(n_bases)-1:
            lengths.append(n_bases[pos]-n_bases[pos+1])
            pos += 1
        else:
            lengths.append(n_bases[pos])
        metrics[pfx+'ModeLen'] = lengths.index(max(lengths))+1
    # Add the version, create the pandas object, write it to file.
    metrics['softwareSeqTKversion_fqchk'] = seqtk_version()
    metrics_df = create_pandas_df(metrics, isolate)
    # NB: this will overwrite the outfile that was read in at the first step.
    write_pandas_df(outfile, metrics_df)


def run_kraken(infile, outfile, fmt, isolate, dbase, threads):
    '''
    Run Kraken on the infile.
    '''
    def do_kraken(cmd_kraken):
        '''
        This function exists so do_kraken() can be bypassed if there are no
        infiles.
        '''
        cmd_krk_r = 'kraken-report'
        cmd_grep = "grep -P '\tS\t'"
        cmd_sort = 'sort -k 1 -r'
        cmd_head = 'head -3'
        # Split the cmds using shlex, store in args.
        args_kraken = shlex.split(cmd_kraken)
        args_krk_report = shlex.split(cmd_krk_r)
        args_grep = shlex.split(cmd_grep)
        args_sort = shlex.split(cmd_sort)
        args_head = shlex.split(cmd_head)

        # Pipe the output of one args to another.
        proc1 = Popen(args_kraken, stdout=PIPE)
        proc2 = Popen(args_krk_report, stdin=proc1.stdout, stdout=PIPE,
                      stderr=PIPE)
        proc3 = Popen(args_grep, stdin=proc2.stdout, stdout=PIPE, stderr=PIPE)
        proc4 = Popen(args_sort, stdin=proc3.stdout, stdout=PIPE, stderr=PIPE)
        proc5 = Popen(args_head, stdin=proc4.stdout, stdout=PIPE, stderr=PIPE)

        output = proc5.stdout.read().decode('UTF-8')
        kraken = output.rstrip().split('\n')
        kraken = [line.strip().split('\t') for line
                  in [_f for _f in kraken if _f]]
        return kraken

    if fmt == 'reads':
        if len(infile) == 2:
            infiles = ' '.join(infile)
            compression = ''
            for i in infile:
                # Compression test based on file extension using linux 'file'.
                args = shlex.split('file '+i)
                proc = Popen(args, stdout=PIPE)
                f_fmt = proc.stdout.read().decode('UTF-8').rstrip().split()
                if 'gzip' in f_fmt:
                    compression = '--gzip-compressed '
                    break
                if 'bzip2' in f_fmt:
                    compression = '--bzip2-compressed '
                    break

            cmd_kraken = 'kraken --threads '+str(threads)+' --db '+dbase +\
                         ' --fastq-input '+compression +\
                         '--paired --check-names '+infiles
            kraken = do_kraken(cmd_kraken)
        else:
            # If no read pairs in list, kraken is an empty list.
            kraken = []

    if fmt == 'contigs':
        if len(infile) == 1:
            infile = ''.join(infile)
            cmd_kraken = 'kraken --threads '+str(threads)+' --db '+dbase +\
                         ' --fasta-input '+infile
            kraken = do_kraken(cmd_kraken)
        else:
            kraken = []

    def kraken_results_df_creator(kraken_hits, fmt):
        '''
        Take the 2D array from a kraken search and return result as a
        dictionary.
        '''
        dict_hits = {}
        k = 1
        for i in range(0, len(kraken_hits)):
            dict_hits['Sp_krkn_'+fmt+'_'+str(k)] =\
                      kraken_hits[i][5].lstrip()
            dict_hits['Sp_krkn_'+fmt+'_'+str(k)+'_pc'] = kraken_hits[i][0]
            k += 1
        return dict_hits
    if len(kraken) > 0:
        krk_df = kraken_results_df_creator(kraken, fmt)
    else:
        krk_df = {}

    def get_krkn_version():
        '''
        Get the Kraken software version and Kraken Database path.
        '''
        args = shlex.split('kraken -v')
        proc = Popen(args, stdout=PIPE)
        version = proc.stdout.read().decode('UTF-8').rstrip().split('\n')[0]
        return {'softwareKrakenVersion_'+fmt: version,
                'softwareKrakenDB_'+fmt: dbase}

    krk_vers_no = create_pandas_df(get_krkn_version(), isolate)
    krk_result = create_pandas_df(krk_df, isolate)
    kraken_df = pd.concat([krk_result, krk_vers_no], axis=1)
    write_pandas_df(outfile, kraken_df)


def run_mlst(assembly, outfile, isolate, species):
    '''
    Run Torsten's MLST program on the assembly.
    '''
    if len(assembly) == 0:
        mlst_formatted_dict = {'MLST_Scheme': None, 'MLST_ST': None}
        k = 1
        for i in range(3, 10):
            mlst_formatted_dict['MLST_Locus'+str(k)] = None
            k += 1
    if len(assembly) == 1:
        assembly = assembly[0]
        sp_scheme = None
        if species in FORCE_MLST_SCHEME:
            sp_scheme = FORCE_MLST_SCHEME[species]
        elif species.split(' ')[0] in FORCE_MLST_SCHEME:
            sp_scheme = FORCE_MLST_SCHEME[species.split(' ')[0]]
        if sp_scheme is not None:
            print("species scheme is:", sp_scheme)
            cmd = 'mlst --scheme '+sp_scheme+' --quiet ' +\
                   assembly
            args_mlst = shlex.split(cmd)
            proc = Popen(args_mlst, stdout=PIPE)
            output = proc.stdout.read().decode('UTF-8')
            mlst = output.rstrip().split('\n')
            mlst = [line.strip().split('\t') for line in [_f for _f in
                                                          mlst if _f]]
            header = mlst[0]
            data = mlst[1]
            ncol = len(header)
            mlst_formatted_dict = {'MLST_Scheme': data[1], 'MLST_ST': data[2]}
            k = 1
            for i in range(3, ncol):
                locus = header[i]+'('+data[i]+')'
                mlst_formatted_dict['MLST_Locus'+str(k)] = locus
                k += 1
        else:
            cmd = 'mlst --quiet '+assembly
            args_mlst = shlex.split(cmd)
            proc = Popen(args_mlst, stdout=PIPE)
            output = proc.stdout.read().decode('UTF-8')
            out = output.rstrip().split('\t')[1:]
            ncol = len(out)
            mlst_formatted_dict = {'MLST_Scheme': out[0],
                                   'MLST_ST': out[1]}
            k = 1
            for i in range(3, ncol):
                mlst_formatted_dict['MLST_Locus'+str(k)] = out[i]
                k += 1

    def get_mlst_version():
        '''
        Get the MLST software version.
        '''
        args = shlex.split('mlst -v')
        proc = Popen(args, stdout=PIPE)
        version = proc.stdout.read().decode('UTF-8').rstrip().split('\n')[0]
        return {'softwareMLSTversion': version}

    mlst_version = create_pandas_df(get_mlst_version(), isolate)
    mlst_result = create_pandas_df(mlst_formatted_dict, isolate)
    mlst_df = pd.concat([mlst_result, mlst_version], axis=1)
    write_pandas_df(outfile, mlst_df)


def run_ariba(infiles, outfile, isolate, dbase, result_basedir):
    '''
    Run Ariba on the reads.
    '''
    def ariba_version():
        '''
        Get the Ariba software version.
        '''
        args = shlex.split('ariba version')
        proc = Popen(args, stdout=PIPE)
        version = proc.stdout.read().decode('UTF-8').rstrip().split('\n')[0]
        return {'softwareAribaVersion_'+dbase[0]: version,
                'softwareAribaDB_'+dbase[0]: dbase[1]}

    if len(infiles) < 2:
        ariba_summary_dict = {}
    else:
        cmd = 'ariba run --force '+dbase[1]+' '+' '.join(infiles)+' ' +\
              outfile  # --threads generates an error
        print(cmd)
        os.system(cmd)
        if os.path.exists(os.path.join(outfile, 'report.tsv')):
            cmd_sum = 'ariba summary --cluster_cols assembled,known_var,' +\
                      'match,ref_seq,novel_var,pct_id ' +\
                      os.path.join(result_basedir, 'ariba_summary') +\
                      ' '+os.path.join(outfile, 'report.tsv')
            print(cmd_sum)
            os.system(cmd_sum)
            ariba_summary_dict = {}
            ariba_data = pd.read_table(os.path.join(result_basedir,
                                                    'ariba_summary.csv'),
                                       sep=',', header=0)
            col_chunks = list(set([i.split('.')[0]
                                   for i in list(ariba_data.columns.values)
                                   if i != 'name']))
            # Uncomment next line to get possible 'extensions' for ariba
            # columns
            # col_chunk_extensions = list(set([i.split('.')[-1]
            #                                  for i in
            #                                  list(ariba_data.columns.values)
            #                                  if i != 'name']))
            for i in col_chunks:
                if i+'.match' in ariba_data.columns.values:
                    ariba_summary_dict['ariba_'+dbase[0] +
                                       '_'+ariba_data.loc[0, i+'.ref_seq']] = \
                     'ASSMBLD_'+ariba_data.loc[0, i+'.assembled'] +\
                     '_PCTID_'+str(ariba_data.loc[0, i+'.pct_id']) +\
                     '_MTCH_'+ariba_data.loc[0, i+'.match']
                else:
                    ariba_summary_dict['ariba_'+dbase[0] +
                                       '_'+ariba_data.loc[0, i+'.ref_seq']] = \
                     'ASSMBLD_'+ariba_data.loc[0, i+'.assembled'] +\
                     '_PCTID_'+str(ariba_data.loc[0, i+'.pct_id'])

        else:
            ariba_summary_dict = {'ariba_error_'+dbase[0]: 'returned_error'}
    ariba_df = pd.concat([create_pandas_df(ariba_version(), isolate),
                          create_pandas_df(ariba_summary_dict, isolate)],
                         axis=1)
    write_pandas_df(os.path.join(result_basedir, 'ariba_'+dbase[0]+\
                                 '_summary_melted.txt'),
                    ariba_df)


def run_quicktree(infile, outfile):
    '''
    Run quicktree, replace the quotes in the stdout, write the tree to file
    '''
    cmd = 'quicktree -in m -out t '+infile
    args = shlex.split(cmd)
    proc = Popen(args, stdout=PIPE)
    treestring = proc.stdout.read().decode('UTF-8') \
                 .rstrip().replace("'", "").replace("\n", "")
    tre = Tree(treestring, format=1)
    # This redundant print statement is to make use of (to trick pylint) the
    # dummy outfile variable (required as input by @files decorator)
    print('Got tree for writing to '+outfile)
    return tre


def relabel_tree_tips(tree, dframe, out, matrix):
    '''
    Take a tree, traverse it, swapping out the tip labels with new ones in
    df['tempID'].
    '''
    for leaf in tree.traverse():
        iso = dframe.loc[dframe['tempID'] == leaf.name].index.values
        if len(iso) > 0:
            leaf.name = iso[0]
    tree.set_outgroup(tree.get_midpoint_outgroup())
    print(tree)
    print('If tip labels are not as expected, delete ' +
          matrix+' then rerun the analysis')
    tree.write(outfile=out, format=1)
    print('Tree file at '+out)


def symlink_contigs(infile, outfile):
    '''
    Create symlinks
    '''
    cmd = 'ln -s '+infile+' '+outfile
    os.system(cmd)


def run_andi(infiles, outfile, model, cpus):
    '''
    Run Andi phylo.
    Todo: if a distance has already been calculated, don't recalculate it on
    a re-run.
    '''
    andi_c = 'andi -j -m '+model+' -t ' +\
             str(cpus)+' '+' '.join(infiles)+' > '+outfile
    print(andi_c)
    os.system(andi_c)

# From https://github.com/tseemann/mlst/blob/master/db/species_scheme_map.tab
FORCE_MLST_SCHEME = {"Acinetobacter baumannii": "abaumannii_2",
                     "Achromobacter": "achromobacter",
                     "Aeromonas": "aeromonas",
                     "Aspergillus afumigatus": "afumigatus",
                     "Anaplasma aphagocytophilum": "aphagocytophilum",
                     "Arcobacter": "arcobacter",
                     "Borrelia burgdorferi": "bburgdorferi",
                     "Burkholderia cepacia": "bcc",
                     "Bacillus cereus": "bcereus",
                     "Brachyspira hyodysenteriae": "bhyodysenteriae",
                     "Bifidobacterium bifidobacterium": "bifidobacterium",
                     "Brachyspira intermedia": "bintermedia",
                     "Bacillus licheniformis": "blicheniformis",
                     "Bordetella pertussis": "bordetella",
                     "Brachyspira pilosicoli": "bpilosicoli",
                     "Burkholderia pseudomallei": "bpseudomallei",
                     "Brachyspira": "brachyspira",
                     "Candida albicans": "calbicans",
                     "Campylobacter jejuni": "campylobacter",
                     "Campylobacter coli": "campylobacter",
                     "Clostridium botulinum": "cbotulinum",
                     "Campylobacter concisus": "cconcisus",
                     "Peptoclostridium difficile": "cdifficile",
                     "Clostridium difficile": "cdifficile",
                     "Corynebacterium diphtheriae": "cdiphtheriae",
                     "Campylobacter fetus": "cfetus",
                     "Candida glabrata": "cglabrata",
                     "Campylobacter helveticus": "chelveticus",
                     "Chlamydia": "chlamydiales",
                     "Campylobacter hyointestinalis": "chyointestinalis",
                     "Campylobacter insulaenigrae": "cinsulaenigrae",
                     "Candida krusei": "ckrusei",
                     "Campylobacter lanienae": "clanienae",
                     "Campylobacter lari": "clari",
                     "Cryptococcus neoformans": "cneoformans",
                     "Cronobacter": "cronobacter",
                     "Clostridium septicum": "csepticum",
                     "Clonorchis sinensis": "csinensis",
                     "Campylobacter sputorum": "csputorum",
                     "Candida tropicalis": "ctropicalis",
                     "Campylobacter upsaliensis": "cupsaliensis",
                     "Enterobacter cloacae": "ecloacae",
                     "Escherichia": "ecoli",
                     "Shigella": "ecoli",
                     "Enterococcus faecalis": "efaecalis",
                     "Enterococcus faecium": "efaecium",
                     "Flavobacterium psychrophilum": "fpsychrophilum",
                     "Haemophilus": "haemophilus",
                     "Helicobacter cinaedi": "hcinaedi",
                     "Haemophilus parasuis": "hparasuis",
                     "Helicobacter pylori": "hpylori",
                     "Haematopinus suis": "hsuis",
                     "Klebsiella oxytoca": "koxytoca",
                     "Klebsiella pneumoniae": "kpneumoniae",
                     "Lactobacillus casei": "lcasei",
                   # "Legionella": "legionella", #STOP. Omp locus problem
                     "Leptospira": "leptospira",
                     "Listeria monocytogenes": "lmonocytogenes",
                     "Lactobacillus salivarius": "lsalivarius",
                     "Mycobacterium abscessus": "mabscessus",
                     "Mycoplasma agalactiae": "magalactiae",
                     "Moraxells catarrhalis": "mcatarrhalis",
                     "Mannheimia haemolytica": "mhaemolytica",
                     "Mycoplasma hyorhinis": "mhyorhinis",
                     "Mycobacterium massiliense": "mmassiliense",
                     "Melissococcus plutonius": "mplutonius",
                     "Neisseria": "neisseria",
                     "Propionibacterium acnes": "pacnes",
                     "Pseudomonas aeruginosa": "paeruginosa",
                     "Pantoea agglomerans": "pagglomerans",
                     "Pseudomonas fluorescens": "pfluorescens",
                     "Propionibacterium freudenreichii": "pfreudenreichii",
                     "Porphyromonas gingivalis": "pgingivalis",
                     "Pasteurella multocida": "pmultocida_multihost",
                     "Pediococcus pentosaceus": "ppentosaceus",
                     "Plesiomonas shigelloides": "pshigelloides",
                     "Streptococcus agalactiae": "sagalactiae",
                     "Staphylococcus aureus": "saureus",
                     "Streptococcus canis": "scanis",
                     "Streptococcus dysgalactiae": "sdysgalactiae",
                     "Salmonella enterica": "senterica",
                     "Staphylococcus epidermidis": "sepidermidis",
                     "Streptococcus gallolyticus": "sgallolyticus",
                     "Sinorhizobium": "sinorhizobium",
                     "Stenotrophomonas maltophilia": "smaltophilia",
                     "Streptococcus oralis": "soralis",
                     "Streptococcus pneumoniae": "spneumoniae",
                     "Staphylococcus pseudintermedius": "spseudintermedius",
                     "Streptococcus pyogenes": "spyogenes",
                     "Streptococcus suis": "ssuis",
                     "Streptococcus thermophilus": "sthermophilus",
                     "Streptomyces": "streptomyces",
                     "Streptococcus uberis": "suberis",
                     "Streptococcus equi": "szooepidemicus",
                     "Taylorella": "taylorella",
                     "Vibrio cholerae": "vcholerae",
                     "Vibrio": "vibrio",
                     "Vibrio parahaemolyticus": "vparahaemolyticus",
                     "Vibrio tapetis": "vtapetis",
                     "Vibrio vulnificus": "vvulnificus",
                     "Wolbachia": "wolbachia",
                     "Xylella fastidiosa": "xfastidiosa",
                     "Yersinia": "yersinia",
                     "Yersinia pseudotuberculosis": "ypseudotuberculosis",
                     "Yersinia ruckeri": "yruckeri"}