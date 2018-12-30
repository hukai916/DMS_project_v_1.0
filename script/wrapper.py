"""
A Python wrapper of DMS pipeline.
"""
from pathlib import Path
from toolbox_kai import fastq_splitter_wrapper, quality_control_library, parse_bamfile, wt_codon
from toolbox_enrich2 import enrich2_json_encoder, enrich2_hdf5_extractor, enrich2_hdf5_extractor2, tsv_plot_output, tsv_plot_output_aa, enrich2_count_wrapper
from toolbox_sam import md2list
from ast import literal_eval as make_tuple
import shutil # complement pathlib with copy file utilities
import subprocess
import pandas as pd
import numpy as np
import pysam
import sys
import os
import glob # Unix style pathname pattern expansion, this seems to work better then Path().glob()
from Bio import SeqIO
from constants import CodonList, AA2CODON, AAList
from Bio import SeqIO			

class ConfigParam(): # store the configuration paramters
    def __init__(self, configFile):
        self.parseConfig(configFile)
    def parseConfig(self, filename):
        for line in open(filename):
            if line.startswith("NGS data storage(if NCBI)"):
                self.ngs_data_ncbi  = line.split(":")[1].split()
            elif line.startswith("NGS data storage(if local)"):
                self.ngs_data_local = line.split(":")[1].strip().split()
            elif line.startswith("WT fasta"):
                self.wtfile = workdir.joinpath(line.split(":")[1].strip())
            elif line.startswith("Condition"):
                self.condition = line.split("Condition:")[1].strip().split()
            elif line.startswith("Amplicon locations:"):
                self.amplicon = list(map(make_tuple, line.split("Amplicon locations:")[1].strip().split()))
            elif line.startswith("Mutation coordination in WT:"):
                self.mut_pos  = list(map(make_tuple, line.split("Mutation coordination in WT:")[1].strip().split()))
        self.mut_list = [item for subitem in self.mut_pos for item in subitem]

def get_ngs_data_ncbi(param, folder_data_sra): # retrieve NGS data from NCBI, save them to data_sra folder
    error = 0
    if not len(glob.glob(folder_data_sra.as_posix() + "/*")) == 0: # check to see if file already downloaded
        print("Step0: sra already downloaded!")
    else:
        Path(folder_data_sra).mkdir(parents=True, exist_ok=True) # create folder to store sra
        folder_sra = Path(folder_data_sra).as_posix()
        print("Step0: downloading sra from NCBI ...")
        for sra in param.ngs_data_ncbi:
            command = ' '.join(['fastq-dump', '--split-files', sra, '-O', folder_sra])
            try:
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'))
            except:
                print("ERROR in Step0! Fail to download data from NCBI!!")
                error = 1
                break
        if not error: print("Step0: download succeeded!")

    # Rename sra file to match temperature: NrdF_28_R1.fastq etc
    if not error:
        print("Step0: renaming ...")
        file_pair = list(zip(param.ngs_data_ncbi, param.condition))
        try:
            for pair in file_pair:
                p = Path(folder_data_sra.joinpath(Path(pair[0]+"_1.fastq"))).as_posix()
                target = Path(folder_data_sra.joinpath(Path(pair[1]+'_R1.fastq'))).as_posix()
                shutil.copy(p, target)
                p = Path(folder_data_sra.joinpath(Path(pair[0]+"_2.fastq"))).as_posix()
                target = Path(folder_data_sra.joinpath(Path(pair[1]+'_R2.fastq'))).as_posix()
                shutil.copy(p, target)
        except:
            print("Step0: rename failed")
            error = 1
        if not error: print('Step0: rename succeeded!')

def get_ngs_data_local(param, folder_data_sra): # retrieve NGS data from local storage
    error = 0
    if not len(glob.glob(folder_data_sra.as_posix() + "/*")) == 0: # check to see if file already downloaded
        print("Step0: sra already downloaded!")
    else:
        Path(folder_data_sra).mkdir(parents=True, exist_ok=True) # create folder to store sra
        folder_sra = Path(folder_data_sra).as_posix()
        print("Step0: copying raw ngs files from local storage ...")
        for file in param.ngs_data_local:
            command = ' '.join(['cp', file, folder_sra])
            try:
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'))
            except:
                print("ERROR in Step0! Fail to copy data from local storage!!")
                error = 1
                break
        if not error: print("Step0: data copied!")

    # Rename sra file to match temperature: NrdF_28_R1.fastq etc
    if not error:
        print("Step0: renaming ...")
        name_list = [Path(file) for file in param.ngs_data_local]
        #name_list = [item[1] for item in enumerate(name_list) if item[0]%2 == 0]
        file_pair = list(zip(name_list, param.condition))
        try:
            for i in range(len(param.condition)):
                prior_file  = name_list[i*2]
                target_file = folder_data_sra.joinpath(param.condition[i] + '_R1.fastq')
                shutil.copy(prior_file, target_file)

                prior_file  = name_list[i*2 +1]
                target_file = folder_data_sra.joinpath(param.condition[i] + '_R2.fastq')
                shutil.copy(prior_file, target_file)

        except:
            print("Step0: rename failed")
            error = 1
        if not error: print('Step0: rename succeeded!')
    pass

def prep_ref(param, folder_ref): # prepare ref seq by copying it into ref folder and renaming seq name to match file name
    filename = param.wtfile.name
    Path(folder_ref).mkdir(parents=True, exist_ok=True)
    wt_filename = folder_ref.joinpath(filename)

    Path(wt_filename).touch(exist_ok=True)
    file   = open(wt_filename, 'w+')
    _wtfile = SeqIO.parse(param.wtfile.as_posix(), 'fasta')
    for record in _wtfile:
        file.write('>' + filename + '\n')
        file.write(str(record.seq))
    file.close()
    param.wtfile = wt_filename

def fastq_merger(param, folder_data_sra, folder_merge):
    error = 0
    if not len(glob.glob(folder_merge.as_posix() + "/*.fastq")) == 0:
        print("Step1: merged file already in merge folder!")
        error = 1
    else:
        Path(folder_merge).mkdir(parents=True, exist_ok=True) # create folder to store merged data
        print("Step1: merging ...")
        try:
            for cond in param.condition:
                infile1 = folder_data_sra.joinpath(Path(cond+'_R1.fastq'))
                infile2 = folder_data_sra.joinpath(Path(cond+'_R2.fastq'))
                outfile = folder_merge.joinpath(Path(cond+'.fastq'))
                command = ' '.join(['bbmerge.sh', 'in1='+infile1.as_posix(),
                'in2='+infile2.as_posix(), 'out='+outfile.as_posix(), '-Xmx250m'])
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb')) # stderr is to silence bbmerge.sh
        except:
            print("Step1: merging failed!")
            error = 1
            # Note that if bbmerge.sh is not installed, the error will not be set to 1 because all error infor are ouptut to os.devnull.
    if not error:
        print("Step1: merge succeeded!")

def first_mapper(param, folder_merge, folder_first_mapper): # map merged reads onto reference
    error = 0
    try: # indexing reference
        command = ' '.join(['bwa index', param.wtfile.as_posix()])
        subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
    except:
        print("Step2: indexing ref failed!")
        error = 1

    if not len(glob.glob(folder_first_mapper.as_posix() + "/*.csi")) == 0:
        print("Step2: first-round mapped files already in first_mapper folder!")
        error = 1
    else:
        Path(folder_first_mapper).mkdir(parents=True, exist_ok=True) # create folder to store first-round mapped files
        print("Step2: first-round mapping ...")
        try: # bwa mem
            for cond in param.condition:
                # mapping
                infile  = Path(folder_merge).joinpath(cond+'.fastq')
                outfile = Path(folder_first_mapper).joinpath(cond+'_bwa.sam')
                command = ' '.join(['bwa mem', param.wtfile.as_posix(), infile.as_posix(),
                 '>', outfile.as_posix()])
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

                # sorting
                infile  = Path(folder_first_mapper).joinpath(cond+'_bwa.sam')
                outfile = Path(folder_first_mapper).joinpath(cond+'_bwa.bam')
                command = ' '.join(['samtools sort -m 250M', infile.as_posix(),
                '>', outfile.as_posix()])
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

                # indexing
                infile  = Path(folder_first_mapper).joinpath(cond+'_bwa.bam')
                command = ' '.join(['samtools index -m 250M', infile.as_posix()])
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        except:
            print("Step2: first-round mapping failed!")
            error = 1

    if not error:
            print("Step2: first-round mapping succeeded!")

def quality_control_INDEL(param, workdir, folder_qc_INDEL): # quality control step, filter out unqualified reads with INDELs
    error = 0
    if not len(glob.glob(folder_qc_INDEL.as_posix() + "/*.fastq")) == 0:
        print("Step3: qc INDEL already performed!")
        error = 1
    else:
        try:
            Path(folder_qc_INDEL).mkdir(parents=True, exist_ok=True) # create folder to store qc INDEL files.
            for cond in param.condition:
                #print(cond, test, "!!!!")
                samfile = pysam.AlignmentFile(workdir.joinpath('first_mapper/' + cond + "_bwa.bam"), 'rb')
                wt_name = param.wtfile.name
                itr = list(samfile.fetch(wt_name)) # Note the ref seq name must match its file name.
                                                   # itr can only be referenced once, so use list to unpack it first.

                target = [[item.cigartuples, item.query_name, item.query_sequence] for item in itr]
                df = pd.DataFrame(target, columns=['Cigar', 'Query_name', 'Query_seq'])
                df['Cigar_len'] = df.apply(lambda row:len(row['Cigar']), axis=1)
                df['Cigar_first'] = df.apply(lambda row:row['Cigar'][0][0], axis=1)
                df['Query_name']  = [item.query_name for item in itr]

                read_count_without_INDEL = df[(df['Cigar_len'] == 1)
                                    & (df['Cigar_first'] ==0)]
                read_count_without_INDEL_df     = df[(df['Cigar_len'] == 1)
                                                & (df['Cigar_first'] ==0)]

                read_count_without_INDEL_count  = read_count_without_INDEL_df.shape[0]
                read_count_with_INDEL_count = df[~((df['Cigar_len'] == 1)
                                        & (df['Cigar_first'] == 0))].count()[0]

                print("Step3: sequences with INDELs removed: " + str(read_count_without_INDEL_count)
                      + ' qualified while ' + str(read_count_with_INDEL_count) + " with INDELs"
                      + ''.join(["(", cond, ")."]))

                # Below is to split valid fastq files out and save them into folder_qc_INDEL folder.
                read_count_without_INDEL_query_name_file = workdir.joinpath('TemFolder/tem.txt')
                np.savetxt(read_count_without_INDEL_query_name_file,
                            read_count_without_INDEL_df['Query_name'].values, fmt='%s')
                raw_fastq_file = workdir.joinpath('merge/' + cond + ".fastq")
                out_fastq_file = folder_qc_INDEL.joinpath('qc_INDEL_' + cond + '.fastq')
                fastq_splitter_wrapper(raw_fastq_file, read_count_without_INDEL_query_name_file, out_fastq_file)
        except:
            error = 1
            print("Step3: qc step: INDEL removal failed!")
    if not error:
        print("Step3: qc step: INDEL removed!")

def quality_control_library_wrapper(param, folder_qc_INDEL, folder_qc_library):
    error = 0
    if not len(glob.glob(folder_qc_library.as_posix() + "/*.fastq")) == 0:
        print("Step4: library seqs already in qc_library folder!")
        error = 1
    else:
        folder_qc_library.mkdir(parents=True, exist_ok=True)
        assert (len(param.amplicon) == len(param.mut_pos)), "Amplicon # doesn't match mut_pos #!"
        library_seq = quality_control_library(param.wtfile, param.amplicon, param.mut_pos)
        for cond in param.condition:
            INDEL_list = open(Path(folder_qc_INDEL.joinpath('qc_INDEL_' + cond + '.fastq'))).readlines()
            library_list = INDEL_list + library_seq
            library_list = [x.strip() for x in library_list]
            library_seq_df = pd.DataFrame({'fastq': library_list})
            library_file_name = Path(folder_qc_library.joinpath(Path('qc_library_' + cond + '.fastq')))
            Path(library_file_name).touch(exist_ok=True)
            np.savetxt(library_file_name, library_seq_df['fastq'].values, fmt='%s')
    if not error:
        print("Step4: qc step: library seqs added!")

def second_mapper(param, folder_qc_library, folder_second_mapper):
    """
    Performs essentially the same operation on quality-controlled
    sequences (removal of INDELs and addition of library queries).
    """
    error = 0
    if not len(glob.glob(folder_second_mapper.as_posix() + "/*.csi")) == 0:
        print("Step5: second-round mapped files already in second_mapper folder!")
        error = 1
    else:
        Path(folder_second_mapper).mkdir(parents=True, exist_ok=True) # create folder to store second-round mapped files
        print("Step5: second-round mapping ...")
        try: # bwa mem
            for cond in param.condition:
                # mapping
                infile  = Path(folder_qc_library).joinpath('qc_library_' + cond+'.fastq')
                outfile = Path(folder_second_mapper).joinpath(cond+'_bwa.sam')
                command = ' '.join(['bwa mem', param.wtfile.as_posix(), infile.as_posix(),
                 '>', outfile.as_posix()])
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

                # sorting
                infile  = Path(folder_second_mapper).joinpath(cond+'_bwa.sam')
                outfile = Path(folder_second_mapper).joinpath(cond+'_bwa.bam')
                command = ' '.join(['samtools sort -m 250M', infile.as_posix(),
                '>', outfile.as_posix()])
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

                # indexing
                infile  = Path(folder_second_mapper).joinpath(cond+'_bwa.bam')
                command = ' '.join(['samtools index -m 250M', infile.as_posix()])
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        except:
            print("Step5: second-round mapping failed!")
            error = 1

    if not error:
            print("Step5: second-round mapping succeeded!")

def bam2enrich_wrapper(param, folder_second_mapper, folder_enrich2_input):
    """
    This function converts bam files into Enrich2 required input format
    and stores the output into folder enrich2_input folder.
    """
    error = 0
    if not len(glob.glob(folder_enrich2_input.as_posix() + "/*.fastq")) == 0:
        print("Step6: formatted Enrich2_input files already in enrich2_input folder!")
        error = 1
    else:
        try:
            print("Step6: converting bam files into Enrich2 input format ...")
            folder_enrich2_input.mkdir(parents=True, exist_ok=True)
            for cond in param.condition:
                bamfile = folder_second_mapper.joinpath(cond + '_bwa.bam')
                df = parse_bamfile(bamfile, param.wtfile)

                """
                df columns:
                           'Cigar'
                           'Query_name'
                           'Query_seq'
                           'Query_quality'
                           'Ref_pos' # list object, mapped query coordinates onto reference.
                           'Query_md'
                           'Cigar_len'
                           'Cigar_first'
                           'Query_name'
                           'Mut_pos' # list object, mapped query mutation site onto reference.
                           ‘Which_amplicon’ # stores the amplicon id, starting from 1
                """

                count_INDEL = df[~((df['Cigar_len'] == 1) & (df['Cigar_first'] == 0))].count()

                def which_amplicon(row, amplicon_range_list): # Assign the amplicon id for each of the query
                    """
                    The amplicon which overlaps the most to the given sequence is
                    the correct amplicon that the sequence comes from.
                    """
                    set1 = set(row['Ref_pos'])
                    max_value = 0
                    which_amplicon = 0
                    for i in range(len(amplicon_range_list)):
                        tem  = len(set1.intersection(
                                set(np.arange(amplicon_range_list[i][0],amplicon_range_list[i][1]))
                                ))/len(set1)
                        if tem > max_value:
                            max_value = tem
                            which_amplicon = i + 1
                    return(which_amplicon)

                def which_mut_site(row, mut_pos):
                    which_mut_site = '*' # default is *
                    for i in range(len(mut_pos)): # control for amplicon
                        for pos in mut_pos[i]:
                            start, end = pos, pos + 2
                            set_real    = set(row['Mut_pos'])
                            set_design  = set(np.arange(start, end+1))
                            if set_real.union(set_design) == set_design and row['Which_amplicon'] == i + 1:
                                if not len(set_real) == 0:
                                    which_mut_site = start
                                else: which_mut_site = 'wt'
                    return(which_mut_site)

                df['Which_amplicon'] = df.apply(which_amplicon, axis=1, amplicon_range_list=param.amplicon)
                df['Which_mut_site'] = df.apply(which_mut_site, axis=1, mut_pos=param.mut_pos)

                wtseq = SeqIO.parse(param.wtfile.as_posix(), 'fasta')
                wtseq = next(wtseq).seq
                for i in range(len(param.mut_pos)):
                    for site in param.mut_pos[i]:
                        df_subset = df[(df["Which_mut_site"]==site) |
                                       (df["Which_mut_site"]=='wt') &
                                       (df["Which_amplicon"]==i+1)].copy()
                       # df_subset keeps only the sequences that have a single mutated aa at the specified site or wt.
                       # This, however, will also keep wt fragment that don't cover target mutation site. This can cause trouble.
                       # Use _fragment_filter to get rid of those fragments.

                        def _fragment_filter(row): # to remove fragments that don't cover mutation site
                            for pos in param.mut_pos[i]:
                                if not pos in row['Ref_pos']:
                                    return("Yes")
                            return("No")

                        df_subset['Fragment'] = df_subset.apply(lambda row: _fragment_filter(row), axis=1)
                        df_subset = df_subset[df_subset['Fragment'] == "No"]


                        df_subset["Which_mut_site2"] = site

                        df_subset['Ref_index'] = df_subset.apply(lambda row: row['Ref_pos'].index(row['Which_mut_site2']), axis=1)

                        df_subset['Mut_seq']   = df_subset.apply(lambda row: row['Query_seq'][row['Ref_index']-1:row['Ref_index']+2], axis=1)

                        df_subset['Mut_seq_quality'] = df_subset.apply(lambda row: row['Query_quality'][row['Ref_index']-1:row['Ref_index']+2], axis=1)

                        df_subset['Mut_seq_fastq']   = df_subset.apply(lambda row: '\n'.join(['@'+row['Query_name'], row['Mut_seq'], '+', ''.join([chr(x+33) for x in row['Mut_seq_quality']])]), axis=1)

                        outfilename = '_'.join([cond, str(site), str(site+2), str(wtseq[site-1:site+2])])
                        enrich_infile_name = Path(folder_enrich2_input.joinpath(outfilename + '.fastq'))
                        Path(enrich_infile_name).touch(exist_ok=True)
                        np.savetxt(enrich_infile_name, df_subset['Mut_seq_fastq'].values, fmt='%s')
                        # Should add one info to log.txt stating the number of reads that fall into each category.
        except:
            error = 1
        if not error:
            print("Step6: Enrich2 input preparation succeeded!")

def enrich2_json_encoder_wrapper(param, folder_enrich2_json, folder_enrich2_input, folder_enrich2_output):
    error = 0
    cond_input = param.condition[-1]
    condition_list = param.condition[:-1]

    if not len(glob.glob(folder_enrich2_json.as_posix() + "/*.json")) == 0:
        print("Step7: configuration files already in enrich2_json folder!")
        error = 1
    else:
        print("Step7: building configuration json files ...")
        folder_enrich2_json.mkdir(parents=True, exist_ok=True)
        json_commandfile_name = folder_enrich2_json.joinpath("json.sh")
        json_commandfile_name.touch(exist_ok=True)
        json_commandfile = open(json_commandfile_name, 'w+')
        json_commandfile.write("source activate py2-test\n")
        try:
            for cond in condition_list:
                for mut in param.mut_list:
                    #file_prefix = '_'.join([cond, str(mut), str(mut+2)])
                    wt_code   = wt_codon(param.wtfile, mut)
                    #print("out: ", cond, cond_input, mut, wt_code)
                    res_json = enrich2_json_encoder(cond, cond_input, folder_enrich2_input, mut, wt_code)
                    jsonfilename = Path(folder_enrich2_json.joinpath("_".join([cond, str(mut), str(mut+2), str(wt_code)+".json"])))
                    jsonfilename.touch(exist_ok=True)
                    jsonfile = open(jsonfilename, 'w+')
                    for x in res_json: jsonfile.write(x)
                    jsonfile.close()

                    json_command = ' '.join(['enrich_cmd', jsonfilename.as_posix(),
                                        "ratios complete --no-plots --output-dir ",
                                        folder_enrich2_output.as_posix() + "/" + cond + "_" + "_".join([str(mut), str(mut+2), str(wt_code)])])
                    json_commandfile.write(json_command + '\n')
            json_commandfile.close()
        except: error = 1
    if not error:
        print("Step7: configuration jsons created in enrich2_json folder!")

def enrich2_wrapper(json_bashfile, folder_enrich2_output):
    error = 0
    if not len(glob.glob(folder_enrich2_output.as_posix() + "/*")) == 0:
        print("Step8: Enrich2 already performed, results saved in enrich2_output folder!")
        error = 1
    else:
        print("Step8: performing Enrich2 ...")
        if 1:
            folder_enrich2_output.mkdir(parents=True, exist_ok=True)
            command = "bash " + json_bashfile.as_posix()
            subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        else:
            error = 1
    if not error:
        print("Step8: Enrich2 performed, results saved in enrich2_output folder!")

def enrich2_hdf5_extractor_wrapper(param, folder_enrich2_output, _type='codon'):
    """
    For each target mutation site, extract and combine information from
    raw_count file and main_score file. Detailed rational, see notes in
    enrich2_hdf5_extractor function. The returned DataFrame contains the
    following columns:
    count0, codon, count1, score, SE, logratio, variance, mut_aa, wt_aa
    The type parameter specifies the extraction type: either codon or aa.
    codon will use default enrich2_hdf5_extractor function whereas aa uses
    a slightly modified enrich2_hdf5_extractor2 to retrieve aa data.
    use _type instead of type, because type is built-in keyword.
    """

    cond_input = param.condition[-1]
    condition_list = param.condition[:-1]

    error = 0

    # Check type, if codon, save results into tsv folder, if aa, save to aa_tsv folder.
    if not _type == 'codon':
        folder_enrich2_output_tsv = folder_enrich2_output.joinpath('aa_tsv')
    else: folder_enrich2_output_tsv = folder_enrich2_output.joinpath('tsv')

    if not len(glob.glob(folder_enrich2_output_tsv.as_posix() + "/*")) == 0:
        print("Step9: extracted tsv/aa_tsv files already in enrich2_output/" +
              folder_enrich2_output_tsv.name + " folder!")
        error = 1
    else:
        if _type == 'codon':
            print("Step9: extracting codon data from Enrich2 output hdf5 files ...")
        else: print("Step9: extracting aa data by re-running Enrich2 via count mode ...")
        if 1:
            folder_enrich2_output_tsv.mkdir(parents=True, exist_ok=True)
            for cond in condition_list:
                for mut in param.mut_list:
                    wt_code = str(wt_codon(param.wtfile, mut))
                    file_folder = folder_enrich2_output.joinpath(
                    "_".join([cond, str(mut), str(mut+2), wt_code]))

                    raw_count_file1 = file_folder.joinpath(cond_input + "_lib.h5")
                    raw_count_file2 = file_folder.joinpath(cond + "_lib.h5")
                    score_file      = file_folder.joinpath(cond + "_sel.h5") # In fact, score_file.h5 also contain raw counts data

                    if _type == 'codon':
                        df = enrich2_hdf5_extractor(raw_count_file1, raw_count_file2, wt_code, score_file)
                    else:
                        df = enrich2_hdf5_extractor2(wt_code, score_file)

                    outfilename = folder_enrich2_output_tsv.joinpath("_".join([cond, str(mut), str(mut+2), wt_code + '.tsv']))

                    df_original = np.concatenate((df.columns.values.reshape((1,-1)), df.values), axis = 0)
                    np.savetxt(outfilename, df_original, fmt="%s", delimiter='\t')

                    """
                    Above saves the original df into tsv folder.
                    """

                    """
                    If type == aa, combine the SY and WT into one single row using the formula:

                    """
                    if _type == 'aa': # combine SY and WT into one row (if specified analysis type is 'aa') by using count file to perform Enrich2 again
                        df_filtered = df.dropna() # Enrich2 will keep rows without NaNs, and use that to peform calculation.
                        c0 = df_filtered['count0'].sum()
                        c1 = df_filtered['count1'].sum()
                        c0_wt = df_filtered[df_filtered['wt_aa'] == 'WT']['count0'].sum()
                        c1_wt = df_filtered[df_filtered['wt_aa'] == 'WT']['count1'].sum()
                        res = enrich2_count_wrapper(c0_wt, c1_wt, c0, c1)
                        df = df_filtered[~(df_filtered['wt_aa']=='WT')]
                        res['count0'] = c0_wt
                        res['count1'] = c1_wt
                        res['wt_aa']  = df_filtered.iloc[-1, -2]
                        res['mut_aa'] = res['wt_aa']
                        df = df.append(res, ignore_index=True)

                    df2 = df[(df['count0']/df['count0'].sum() > 0.005) | (df['count1']/df['count1'].sum() > 0.005)].copy()
                    # Keep the rows where at least one count is no less than 0.005 of the total counts

                    df2['score'] = df2.apply(lambda row: row['score'] if row['count1']/df2['count1'].sum() > 0.005 else "deplete", axis=1)
                    # Assign "deplete" to rows with count1 below 0.005

                    outfilename = folder_enrich2_output_tsv.joinpath("_".join([cond, str(mut), str(mut+2), wt_code + '.preprocessed.tsv']))

                    df_preprocessed = np.concatenate((df2.columns.values.reshape((1,-1)), df2.values), axis = 0)
                    np.savetxt(outfilename, df_preprocessed, fmt="%s", delimiter='\t')
                    """
                    Above save preprocessed df into tsv folder.
                    """
        else:
            error = 1
            print("Step9: information extraction failed!")

    if not error:
        print("Step9: information extracted from Enrich2 output!")

def tsv_plot_input(param, folder_plot_input, folder_enrich2_output):
    condition_list = param.condition[:-1]

    error = 0
    if not len(glob.glob(folder_plot_input.as_posix() + "/*")) == 0:
        print("Step10: plot input files already in plot_input folder!")
        error = 1
    else:
        print("Step10: preparing plot input from extracted tsv files ...")
        folder_plot_input.mkdir(parents=True, exist_ok=True)
        CodonList.append('pos')
        #AAList = [key for key in AA2CODON]
        AAList.append('pos')

        for _type in ('tsv', 'aa_tsv'):
            folder_enrich2_output_tsv = folder_enrich2_output.joinpath(_type)

            if _type == 'tsv':
                index = 'codon'
                ColumnList = CodonList
            else:
                index = 'mut_aa'
                ColumnList = AAList
            try: # change to try
                for cond in condition_list:
                    df = pd.DataFrame(columns=ColumnList)
                    df_se = df.copy()
                    for mut in param.mut_list:
                        row = pd.Series(index=df.columns)
                        tsvfile = "_".join([cond, str(mut), str(mut+2), str(wt_codon(param.wtfile, mut))+ ".preprocessed.tsv"])
                        df_tsv = pd.read_csv(folder_enrich2_output_tsv.joinpath(tsvfile), delimiter = '\t')

                        df_tsv.set_index([index], inplace=True)
                        #df_tsv['score_se'] = df_tsv.apply(lambda row: (row['score'], row['SE']), axis=1)
                        row_se = row.copy()

                        row.update(df_tsv['score']) # only keep Enrich2 score
                        row['pos'] = int((mut+2)/3)
                        df = df.append(row, ignore_index=True)

                        row_se.update(df_tsv['SE'])
                        row_se['pos'] = int((mut+2)/3)
                        df_se = df_se.append(row_se, ignore_index=True)

                    if _type == "aa_tsv":
                        filler = '_aa'
                        df = df.rename(columns={'Ter': 'Stop'})
                        df_se = df_se.rename(columns={'Ter': 'Stop'})
                    else: filler = ''

                    savefile = folder_plot_input.joinpath(cond+filler+".tsv")
                    savefile.touch(exist_ok=True)

                    df2array = np.concatenate((df.columns.values.reshape((1,-1)), df.values), axis = 0)
                    np.savetxt(savefile, df2array, fmt='%s', delimiter='\t')

                    savefile = folder_plot_input.joinpath(cond+filler+".se.tsv")
                    savefile.touch(exist_ok=True)

                    df2array = np.concatenate((df_se.columns.values.reshape((1,-1)), df_se.values), axis = 0)
                    np.savetxt(savefile, df2array, fmt='%s', delimiter='\t')

            except:
                error = 1
    if not error:
        print("Step10: tsv to plot_input step succeeded!")

def tsv_plot_output_wrapper(param, folder_plot_input, folder_plot_output_codon):
    """
    Generate multiple plots using tsv_plot_output function.
    Output will be saved in plot_output folder.
    Multiple plots will be created:
    1. raw: contains every information: all codons (64), SE
    2. simple1: all-NaN columns removed, SE
    3. simple2: all-NaN columns removed, SE removed
    """
    error = 0
    if not len(glob.glob(folder_plot_output_codon.as_posix() + "/*")) == 0:
        print("Step11: plot output files already in plot_output/codon folder!")
        error = 1
    else:
        print("Step11: preparing plot output files (codon version) ...")
        folder_plot_output_codon.mkdir(parents=True, exist_ok=True)
        if 1:
            # Plot1:
            for cond in param.condition[:-1]:
                df1      = pd.read_csv(folder_plot_input.joinpath(cond+'.tsv'), delimiter='\t')
                df_se    = pd.read_csv(folder_plot_input.joinpath(cond+'.se.tsv'), delimiter='\t')
                outfile1 = folder_plot_output_codon.joinpath(cond + '_raw.pdf')
                tsv_plot_output(param.wtfile, cond, df1, df_se=df_se, outfilename=outfile1)
            # Plot2:
            for cond in param.condition[:-1]:
                df2      = pd.read_csv(folder_plot_input.joinpath(cond+'.tsv'), delimiter='\t')
                df_se    = pd.read_csv(folder_plot_input.joinpath(cond+'.se.tsv'), delimiter='\t')
                df2.replace(to_replace='deplete', value=1000, inplace=True)
                df2 = df2.apply(lambda col: pd.to_numeric(col), axis=0)
                for col in df2.columns:
                    if sum(np.isnan(df2[col])) == len(df2.index):
                        df2.drop(columns=col, inplace=True)
                            # Note tat np.nan != np.nan
                outfile2 = folder_plot_output_codon.joinpath(cond + '_simple1.pdf')
                tsv_plot_output(param.wtfile, cond, df2, df_se=df_se, outfilename=outfile2, version=2)
            # Plot3:
                outfile3 = folder_plot_output_codon.joinpath(cond + '_simple2.pdf')
                tsv_plot_output(param.wtfile, cond, df2, outfilename=outfile3, version=2)
        else:
            error = 1
    if not error:
        print("Step11: plots (codon version) created!")

def tsv_plot_output_aa_wrapper(param, folder_plot_input, folder_plot_output_aa):
    """
    Generate multiple plots using tsv_plot_output function.
    Output will be saved in plot_output folder.
    Multiple plots will be created:
    1. raw: contains every information: all codons (64), SE
    2. simple1: all-NaN columns removed, SE
    3. simple2: all-NaN columns removed, SE removed
    """
    error = 0
    if not len(glob.glob(folder_plot_output_aa.as_posix() + "/*")) == 0:
        print("Step11: plot output files already in plot_output/aa folder!")
        error = 1
    try:
        print("Step11: preparing plot output files (aa version) ...")
        folder_plot_output_aa.mkdir(parents=True, exist_ok=True)
        try:
            # Plot1:
            for cond in param.condition[:-1]:
                df1      = pd.read_csv(folder_plot_input.joinpath(cond+'_aa.tsv'), delimiter='\t')
                df_se    = pd.read_csv(folder_plot_input.joinpath(cond+'_aa.se.tsv'), delimiter='\t')
                outfile1 = folder_plot_output_aa.joinpath(cond + '_raw.pdf')
                tsv_plot_output_aa(param.wtfile, cond, df1, df_se=df_se, outfilename=outfile1)
            # Plot2:
            for cond in param.condition[:-1]:
                df2      = pd.read_csv(folder_plot_input.joinpath(cond+'_aa.tsv'), delimiter='\t')
                df_se    = pd.read_csv(folder_plot_input.joinpath(cond+'_aa.se.tsv'), delimiter='\t')
                df2.replace(to_replace='deplete', value=1000, inplace=True)
                df2 = df2.apply(lambda col: pd.to_numeric(col), axis=0)
                for col in df2.columns:
                    if sum(np.isnan(df2[col])) == len(df2.index):
                        df2.drop(columns=col, inplace=True)
                            # Note tat np.nan != np.nan
                outfile2 = folder_plot_output_aa.joinpath(cond + '_simple1.pdf')
                tsv_plot_output_aa(param.wtfile, cond, df2, df_se=df_se, outfilename=outfile2)
            # Plot3:
                outfile3 = folder_plot_output_aa.joinpath(cond + '_simple2.pdf')
                tsv_plot_output_aa(param.wtfile, cond, df2, outfilename=outfile3)
        except:
            error = 1
    except: error = 1
    if not error:
        print("Step11: plots (aa version) created!")

def func_wrapper(param, workdir):
    error = 0
    if 1:
        folder_data_sra = workdir.joinpath('data_sra')
        if not param.ngs_data_local[0] == 'to': # If both NCBI sra and local storage provided, prefer to use local copy.
            get_ngs_data_local(param, folder_data_sra) # Step0: download and rename sra files
        else:
            get_ngs_data_ncbi(param, folder_data_sra)

        folder_ref      = workdir.joinpath('ref')
        prep_ref(param, folder_ref)

        folder_merge    = workdir.joinpath('merge')
        fastq_merger(param, folder_data_sra, folder_merge)

        folder_first_mapper = workdir.joinpath('first_mapper')
        first_mapper(param, folder_merge, folder_first_mapper)

        folder_qc_INDEL     = workdir.joinpath('qc_INDEL')
        quality_control_INDEL(param, workdir, folder_qc_INDEL)

        folder_qc_library   = workdir.joinpath('qc_library')
        quality_control_library_wrapper(param, folder_qc_INDEL, folder_qc_library)

        folder_second_mapper = workdir.joinpath('second_mapper')
        second_mapper(param, folder_qc_library, folder_second_mapper)

        folder_enrich2_input = workdir.joinpath('enrich2_input')
        bam2enrich_wrapper(param, folder_second_mapper, folder_enrich2_input)

        folder_enrich2_json = workdir.joinpath('enrich2_json')
        folder_enrich2_output = workdir.joinpath('enrich2_output')
        enrich2_json_encoder_wrapper(param, folder_enrich2_json, folder_enrich2_input, folder_enrich2_output)
        enrich2_wrapper(folder_enrich2_json.joinpath('json.sh'), folder_enrich2_output)

        enrich2_hdf5_extractor_wrapper(param, folder_enrich2_output, _type='codon',)
        enrich2_hdf5_extractor_wrapper(param, folder_enrich2_output, _type='aa', )

        folder_plot_input = workdir.joinpath('plot_input')
        tsv_plot_input(param, folder_plot_input, folder_enrich2_output)

        folder_plot_output_codon = workdir.joinpath('plot_output/codon')
        tsv_plot_output_wrapper(param, folder_plot_input, folder_plot_output_codon)

        folder_plot_output_aa = workdir.joinpath('plot_output/aa')
        tsv_plot_output_aa_wrapper(param, folder_plot_input, folder_plot_output_aa)
    else:
        error = 1

    if not error:
        print("Congrats! All analyisis finished! Check plot_output folder!")
    else:
        print("Failed!")


if __name__ == '__main__':
    configFile = sys.argv[1]
    workdir    = Path(Path.cwd()).parents[0]
    param = ConfigParam(configFile)

    Path(workdir.joinpath('TemFolder')).mkdir(parents=True, exist_ok=True) # create a temperate folder to contain tem files.
    func_wrapper(param, workdir)
