from pathlib import Path
import json
import pandas as pd
import numpy as np
import pylab
import subprocess
import os
import shutil # use to remove non empty folders
from toolbox_kai import wt_codon
from constants import AA2CODON, CODON_GROUPS, CODON2AA, Property2AA, AA_GROUPS
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from toolbox_matplotlib import recentered_cmap

## Note that pd.HDFStore requires: conda install -c conda-forge pytables

def enrich2_json_encoder(cond, cond_input, fastq_folder, mut, wt_code):
    """
    Customize enrich2 json configuration file by specifying:
    condition, mutation position, wt_codon and output folder name.
    """
    wt_code = str(wt_code)

    json_root = {"libraries": [],
                 "name"     : cond,
                 "output directory": ""
                 }

    json_library_lib    = {"fastq":
                                {
                                "filters":{},
                                "reads": Path(fastq_folder.joinpath('_'.join(
                                [cond_input, str(mut), str(mut+2), wt_code + '.fastq']))).as_posix(),
                                "reverse": False
                                },
                            "name": cond_input,
                            "report filtered reads": False,
                            "timepoint": 0,
                            "variants":
                                {
                                "use aligner": False,
                                "wild type":
                                    {
                                    "coding": True,
                                    "reference offset": 0,
                                    "sequence": wt_code
                                    }
                                }
                            }

    json_library_cond    = {"fastq":
                                {
                                "filters":{},
                                "reads": Path(fastq_folder.joinpath('_'.join(
                                [cond, str(mut), str(mut+2), wt_code + '.fastq']))).as_posix(),
                                "reverse": False
                                },
                            "name": cond,
                            "report filtered reads": False,
                            "timepoint": 1,
                            "variants":
                                {
                                "use aligner": False,
                                "wild type":
                                    {
                                    "coding": True,
                                    "reference offset": 0,
                                    "sequence": wt_code
                                    }
                                }
                            }

    json_root["libraries"] = [json_library_lib, json_library_cond]

    return(json.dumps(json_root, indent=4))

def enrich2_json_encoder_count(countfile_pathlist, output_folder):
    json_root = {"libraries": [],
                 "name"     : "Enrich2_count_wrapper",
                 "output directory": output_folder
                 }

    json_library_lib    = {"counts file": countfile_pathlist[0],
                           "fastq":
                                {
                                "filters":{},
                                "length": 'null',
                                "reads": 'null',
                                "reverse": False
                                },
                            "name": 'time0',
                            "report filtered reads": False,
                            "timepoint": 0,
                            "variants":
                                {
                                "use aligner": False,
                                "wild type":
                                    {
                                    "coding": False,
                                    "reference offset": 0,
                                    "sequence": 'ATG'
                                    }
                                }
                            }

    json_library_cond    = {"counts file": countfile_pathlist[1],
                           "fastq":
                                {
                                "filters":{},
                                "length": 'null',
                                "reads": 'null',
                                "reverse": False
                                },
                            "name": 'time1',
                            "report filtered reads": False,
                            "timepoint": 1,
                            "variants":
                                {
                                "use aligner": False,
                                "wild type":
                                    {
                                    "coding": False,
                                    "reference offset": 0,
                                    "sequence": 'ATG'
                                    }
                                }
                            }

    json_root["libraries"] = [json_library_lib, json_library_cond]

    return(json.dumps(json_root, indent=4))

def enrich2_count_wrapper(c0_t, c1_t, c0, c1):
    """
    input 4 count integers and retrieve the Enrich2 output.
    c0_t: target count at time 0
    c1_t: target count at time 1
    c0: all count at time 0 including c0_t
    c1: all count at time 1 including c1_t
    """
    #print("inside count_wrapper: \n", c0_t, c1_t, c0, c1)

    workdir    = Path(Path.cwd()).parents[0]
    folder_json_tem   = workdir.joinpath("json_tem")
    folder_json_tem.mkdir(parents=True, exist_ok=True)

    file_count0_tem   = folder_json_tem.joinpath('c0.tsv')
    file_count0_tem.touch(exist_ok=True)
    count0 = open(file_count0_tem, 'w+')
    count0.write("\tcount\n")
    count0.write('c_t\t' +  str(c0_t) + "\n")
    count0.write('c_other\t' + str(c0-c0_t))
    count0.close()

    file_count1_tem   = folder_json_tem.joinpath('c1.tsv')
    file_count1_tem.touch(exist_ok=True)
    count1 = open(file_count1_tem, 'w+')
    count1.write("\tcount\n")
    count1.write('c_t\t' + str(c1_t) + '\n')
    count1.write('c_other\t' +  str(c1-c1_t))
    count1.close()

    res_json = enrich2_json_encoder_count([file_count0_tem.as_posix(),
                                           file_count1_tem.as_posix()],
                                           folder_json_tem.as_posix())

    file_json_tem     = folder_json_tem.joinpath('tem.json')

    file_json_tem.touch(exist_ok=True)
    jsonfile = open(file_json_tem, 'w+')
    for x in res_json: jsonfile.write(x)
    jsonfile.close()

    json_command = ' '.join(['enrich_cmd', file_json_tem.as_posix(),
                        "ratios complete --no-plots --output-dir ",
                        folder_json_tem.as_posix()])

    json_commandfile = open(folder_json_tem.joinpath('json.sh'), 'w+')
    json_commandfile.write("source activate py2\n")
    json_commandfile.write(json_command)
    json_commandfile.close()

    command = "bash " + folder_json_tem.joinpath('json.sh').as_posix()
    subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

    my_store = pd.HDFStore(folder_json_tem.joinpath('Enrich2_count_wrapper_sel.h5'))
    df = my_store.select("/main/variants/scores")
    my_store.close()

    shutil.rmtree(folder_json_tem) # clean the tem folder
    return(df.iloc[0,:])

def extract_codon(sentence, wt_codon):
    wt_codon = list(wt_codon)
    tem = sentence.split()
    tem = ''.join(tem[::2])
    for i in range(len(tem)):
        if tem[i] == str(1):
            wt_codon[0] = tem[i+3]
        elif tem[i] == str(2):
            wt_codon[1] = tem[i+3]
        elif tem[i] == str(3):
            wt_codon[2] = tem[i+3]
    return("".join(wt_codon))

def enrich2_hdf5_extractor(raw_count1, raw_count2, wt_codon, score_file):
    """
    This function generates a big DataFrame that combines data from two raw_counts file and
    one score_file. We need to incorporate "raw count files" because "main" files of Enrich2 output
    only keep items that show in both raw files.
    """

    my_store = pd.HDFStore(raw_count1)
    df_raw_count1 = my_store.select("/raw/variants/counts")
    df_raw_count1.reset_index(inplace=True)
    my_store.close()

    my_store = pd.HDFStore(raw_count2)
    df_raw_count2 = my_store.select("/raw/variants/counts")
    df_raw_count2.reset_index(inplace=True)
    my_store.close()

    my_store = pd.HDFStore(score_file)
    df_score_file = my_store.select("/main/variants/scores")
    df_score_file.reset_index(inplace=True)
    my_store.close()

    df_raw_count1['codon'] = df_raw_count1.apply(lambda row: extract_codon(row['index'], wt_codon), axis=1)
    df_raw_count1.drop(['index'], axis=1, inplace=True)
    df_raw_count1.rename(columns={'count': 'count0'}, inplace=True)

    df_raw_count2['codon'] = df_raw_count2.apply(lambda row: extract_codon(row['index'], wt_codon), axis=1)
    df_raw_count2.drop(['index'], axis=1, inplace=True)
    df_raw_count2.rename(columns={'count': 'count1'}, inplace=True)

    df_score_file['codon'] = df_score_file.apply(lambda row: extract_codon(row['index'], wt_codon), axis=1)
    df_score_file.drop(['index'], axis=1, inplace=True)

    df_raw = pd.merge(df_raw_count1, df_raw_count2, how='outer', left_on='codon', right_on='codon')
    df_all = pd.merge(df_raw, df_score_file, how='outer', left_on='codon', right_on='codon')

    df_all['mut_aa'] = df_all.apply(lambda row: CODON2AA[row['codon']][0] if row['codon'] in CODON2AA else np.nan, axis=1)
    df_all['wt_aa']  = df_all.apply(lambda row: CODON2AA[wt_codon][0], axis=1)

    return(df_all)

def extract_aa(sentence):
    """
    wt_codon = list(wt_codon)
    tem = sentence.split()
    tem = ''.join(tem[::2])
    for i in range(len(tem)):
        if tem[i] == str(1):
            wt_codon[0] = tem[i+3]
        elif tem[i] == str(2):
            wt_codon[1] = tem[i+3]
        elif tem[i] == str(3):
            wt_codon[2] = tem[i+3]
    """
    wt_aa = sentence
    mut_aa = sentence
    if sentence == '_sy':
        wt_aa = 'WT'
        mut_aa = 'SY'
    elif sentence == '_wt':
        wt_aa = 'WT'
        mut_aa = 'WT'
    else:
        wt_aa  = sentence[2:5]
        mut_aa = sentence[-3:]



    return([wt_aa, mut_aa])

def enrich2_hdf5_extractor2(wt_codon, score_file):
    """
    Perform a similar operation and generate amino acid based data extraction.
    We can actually obtain all required info from _sel.h5 file instead of separately
    retrieving them from _lib.h5 files. _sel.h5 might be actually based on the two _lib.h5
    files.
    Score_file is the sole _sel.h5 file.
    """
    my_store = pd.HDFStore(score_file)
    df_raw_count = my_store.select('/main/synonymous/counts_unfiltered')
    my_store.close()

    my_store = pd.HDFStore(score_file)
    df_raw_score =  my_store.select('/main/synonymous/scores')
    my_store.close()

    df_raw = pd.merge(df_raw_count, df_raw_score, how='outer', left_index=True, right_index=True)
    df_raw.rename(columns={'c_0':'count0', 'c_1':'count1'}, inplace=True)

    df_raw.reset_index(inplace=True)

    df_raw['aa_list'] = df_raw.apply(lambda row: extract_aa(row['index']), axis=1)
    df_raw['wt_aa']   = df_raw.apply(lambda row: row['aa_list'][0], axis=1)
    df_raw['mut_aa']   = df_raw.apply(lambda row: row['aa_list'][1], axis=1)
    df_raw.drop(['aa_list', 'index'], axis=1, inplace=True)
    return(df_raw)

def tsv_plot_output(wtfile, cond, df, df_se=pd.DataFrame(), outfilename='test.pdf'):
    """
    Read in score and se dataframes, output plots accordingly.
    wtfile: WT sequence file, used for y-labeling
    cond: The experiment condition, used for title info
    df: enrich score dataframe, from plot_input folder, can also used customized df
    df_se: enrich score dataframe, from plot_input folder, will match df (if customized) automatically
    outfilename: outputfile path + filename, where to store the output
    """

    row_number, col_number = df.shape

    # Create Figure object that consists of two axes:
    grid = GridSpec(2,1, height_ratios=[row_number,1], hspace=1/(row_number+1)*2)
    fig  = plt.figure()
    fig.set_size_inches((col_number+8)*0.2, (10+row_number)*0.2)
    # Later, will use subplots_adjust to make sure each square is 1*0.2^2 inches^2


    ### Up to now, there is no axes object created, therefore, no plot will be shown.

    # Create two axes object subplot.
    mesh_ax = plt.subplot(grid[0])
    cbar_ax = plt.subplot(grid[1])
    # Adjust subplots to make room for add-on information:
        ## top: (title + amino acid grouping): 5 unit (each unit is 1*0.2 inch)
        ## padding between main figure and bar: 1 unit
        ## bottom: 3 units
        ## left 2.5 units
        ## right: 3.5 units
    plt.gcf().subplots_adjust(top=(row_number+5)/(row_number+10),
                              bottom=3/(row_number+10),
                              left=3/(col_number+8),
                              right=(col_number+3)/(col_number+8))

    # Replace 'deplete' with an arbituraliy large value, say, 1000
    df_raw = df.set_index(['pos'])
    for col in df_raw.columns:
        df_raw[col] = df_raw.apply(lambda row: float(row[col])
                                   if not row[col] == 'deplete'
                                   else 1000, axis=1)
    # Find score range:
    ls = df_raw.values.reshape(1,-1).tolist()[0]
    ls = [x for x in ls if not np.isnan(x)]
    vmin = sorted(set(ls))[0]
    vmax = sorted(set(ls))[-2]
    #print(vmin, vmax)

    # Prepare a numpy array and mask NaNs for later plotting:
    arr_masked = np.ma.array(df_raw, mask=np.isnan(df_raw))

    # Get color map and set colors:
    cmap = plt.get_cmap("RdBu_r")
    colorlim = max(-vmin, vmax)
    # Recenter color map:
    cmap = recentered_cmap(cmap, -colorlim, colorlim)
    # Set grey to NaN values, and black to totally depleted ones:
    cmap.set_bad("#808080",1)
    cmap.set_over("black")

    # Plot the heatmap: return the value as mesn_pcolor as a mappable object in order to add color_bar:
    mesh_pcolor = mesh_ax.pcolormesh(arr_masked, cmap=cmap, vmin=-colorlim, vmax=colorlim)


    ## Below are modifications to plotting:

    # Add in the color bar
    cbar = fig.colorbar(mesh_pcolor, cax=cbar_ax, orientation='horizontal')
    cbar.set_label("Enrich2 Score")

    # Set mesh_ax y labels:
    mesh_wt_ax = mesh_ax.twinx()
    mesh_ax.set_ylabel("Position in WT")
    mesh_wt_ax.set_ylabel("WT sequence")

    # Set mesh_ax title:
    mesh_ax.set_title("Codon Distribution Map for Experiment: " + cond, pad=50)
    # pad uses points, not sure the relationship between point and inch.
    # 50 works nicely it seems.

    # Add in column information:
    for i, x in enumerate(list(df_raw.columns)):
        mesh_ax.text(i + 0.5, len(df_raw.index) + 1, x,
        horizontalalignment="center",
        verticalalignment="center",
        rotation = 90)

    # Add in amino acid grouping information:
    new_CODON_GROUPS = []
    cstart = 0
    cend   = -1
    for i in range(0, len(CODON_GROUPS)):
        aaName  = CODON_GROUPS[i][0]
        count = 0
        for item in df.columns.values:
            if item in AA2CODON[aaName]:
                count = count + 1
        if count == 1: cstart = cend = cend + 1
        else:
            gap = count - 1
            cstart = cend + 1
            cend = cstart + gap
        if cend >= cstart:
            newTuple = (aaName, cstart, cend)
            new_CODON_GROUPS.append(newTuple)

    for codon, start, end in new_CODON_GROUPS:
        mesh_ax.text((end - start + 1) / 2 + start,
            row_number + 2.5, codon,
            horizontalalignment="center",
            verticalalignment="center")

        bar = Line2D([start + 0.125, end + 1 - 0.125],
                [row_number + 2, row_number + 2], color="black")
        bar.set_clip_on(False)
        mesh_ax.add_line(bar)

        # Add in the deliminator for amino acid groups:
        delimBar = Line2D([end + 1, end + 1],
                [0, len(df_raw.index) + 1], color="white")
        delimBar.set_clip_on(False)
        mesh_ax.add_line(delimBar)

    WT_codon  = [str(wt_codon(wtfile, int(x*3-2))) for x in df['pos'].values]
    wt_aa     = [str(CODON2AA[codon][0]) for codon in WT_codon]

    # Add in mesh_ax y-label: aa coordinates in WT sequence:
    ypos = np.arange(len(df['pos'])) + 0.5
    mesh_ax.set_yticks(ypos)
    mesh_ax.set_yticklabels(list(map(int, df['pos'])), ha='right')

    # Add in wt_ax y-label: WT aa and codon:
    mesh_wt_ax.set_ylim(0, len(df['pos']))
    mesh_wt_ax.set_yticks(ypos)
    label = ['-'.join(element) for element in zip(WT_codon, wt_aa)]
    mesh_wt_ax.set_yticklabels(label, ha='left')

    # Add in WT label onto corresponding cell:
    x_coordinate = [list(df.columns).index(x) if x in df.columns else np.nan for x in WT_codon]
    y_coordinate = list(range(df.shape[1]))
    wt_xy        = zip(x_coordinate, y_coordinate)
    for x, y in wt_xy:
        if not x == np.nan:
            mesh_ax.add_patch(Circle((x + 0.5, y + 0.5), .1666,
            fill=True, facecolor="black",
            edgecolor="none", alpha=0.8))

    # Make the figure cleaner by removing ticks:
    mesh_ax.tick_params(bottom=False, left=False)
    mesh_wt_ax.tick_params(right=False)
    mesh_ax.get_xaxis().set_visible(False)
    cbar_ax.tick_params(bottom=False)

    # Add in SE if specified:
    if not df_se.shape == pd.DataFrame().shape: # Add in SE
        # Subset df_se to match df that might be truncated:
        for col in df_se.columns:
            if not col in df.columns: df_se.drop(columns=col, inplace=True)
        df_se = df_se[df_se['pos'].isin(df['pos'])]

        for row in range(len(df_se['pos'])):
            for col in range(len(df_se.columns)-1):
                se_value = df_se.iloc[row, col]
                corner_dist = (1-se_value)/2
                diag = Line2D([col + corner_dist, col + 1 - corner_dist],
                        [row + corner_dist, row + 1 - corner_dist], color="black")

                if se_value > 0.02 and df_raw.iloc[row, col] != 1000: # se_value below 0.02 will not be displayed so as totally depleted ones
                    mesh_ax.add_line(diag)

    pylab.savefig(outfilename)

def tsv_plot_output_aa(wtfile, cond, df, df_se=pd.DataFrame(), outfilename='test.pdf'):
    """
    Read in score and se dataframes, output plots accordingly.
    wtfile: WT sequence file, used for y-labeling
    cond: The experiment condition, used for title info
    df: enrich score dataframe, from plot_input folder, can also used customized df
    df_se: enrich score dataframe, from plot_input folder, will match df (if customized) automatically
    outfilename: outputfile path + filename, where to store the output
    This function is customized for aa.
    """

    row_number, col_number = df.shape

    # Create Figure object that consists of two axes:
    grid = GridSpec(2,1, height_ratios=[row_number,1], hspace=1/(row_number+1)*2)
    fig  = plt.figure()
    fig.set_size_inches((col_number+8)*0.2, (10+row_number)*0.2)
    # Later, will use subplots_adjust to make sure each square is 1*0.2^2 inches^2


    ### Up to now, there is no axes object created, therefore, no plot will be shown.

    # Create two axes object subplot.
    mesh_ax = plt.subplot(grid[0])
    cbar_ax = plt.subplot(grid[1])
    # Adjust subplots to make room for add-on information:
        ## top: (title + amino acid grouping): 5 unit (each unit is 1*0.2 inch)
        ## padding between main figure and bar: 1 unit
        ## bottom: 3 units
        ## left 2.5 units
        ## right: 3.5 units
    plt.gcf().subplots_adjust(top=(row_number+5)/(row_number+10),
                              bottom=3/(row_number+10),
                              left=5/(col_number+8),
                              right=(col_number+5)/(col_number+8))

    # Replace 'deplete' with an arbituraliy large value, say, 1000
    df_raw = df.set_index(['pos'])
    for col in df_raw.columns:
        df_raw[col] = df_raw.apply(lambda row: float(row[col])
                                   if not row[col] == 'deplete'
                                   else 1000, axis=1)
    # Find score range:
    ls = df_raw.values.reshape(1,-1).tolist()[0]
    ls = [x for x in ls if not np.isnan(x)]
    vmin = sorted(set(ls))[0]
    vmax = sorted(set(ls))[-2]

    # Prepare a numpy array and mask NaNs for later plotting:
    arr_masked = np.ma.array(df_raw, mask=np.isnan(df_raw))

    # Get color map and set colors:
    cmap = plt.get_cmap("RdBu_r")
    colorlim = max(-vmin, vmax)
    # Recenter color map:
    cmap = recentered_cmap(cmap, -colorlim, colorlim)

    colors = [cmap(i) for i in range(65,200)]  # R -> G -> B
    cmap_new = LinearSegmentedColormap.from_list('place_holder', colors)
    # The new cmap_new will take advantage of old cmap, "shrink" the spectrum to make it lighter for better printing.

    # Set grey to NaN values, and black to totally depleted ones:
    cmap_new.set_bad("#808080",1)
    rgba = cmap_new(0) # The darkest in the new cmap_new, equal to cmap(65)
    cmap_new.set_over(rgba)

    # Plot the heatmap: return the value as mesn_pcolor as a mappable object in order to add color_bar:
    mesh_pcolor = mesh_ax.pcolormesh(arr_masked, cmap=cmap_new, vmin=-colorlim, vmax=colorlim)

    ## Below are modifications to plotting:

    # Add in the color bar
    cbar = fig.colorbar(mesh_pcolor, cax=cbar_ax, orientation='horizontal')
    cbar.set_label("Enrich2 Score")

    # Set mesh_ax y labels:
    mesh_ax.set_ylabel("Position in WT")

    # Set mesh_ax title:
    mesh_ax.set_title("Amino Acid Distribution Map for Experiment: " + cond, pad=50)
    # pad uses points, not sure the relationship between point and inch.
    # 50 works nicely it seems.

    # Add in column information:
    for i, x in enumerate(list(df_raw.columns)):
        mesh_ax.text(i + 0.5, len(df_raw.index) + 1, x,
        horizontalalignment="center",
        verticalalignment="center",
        rotation = 90)

    # Add in amino acid grouping information:
    new_AA_GROUPS = []
    cstart = 0
    cend   = -1
    for i in range(0, len(AA_GROUPS)):
        aaProperty  = AA_GROUPS[i][0]
        count = 0
        for item in df.columns.values:
            if item in Property2AA[aaProperty]:
                count = count + 1
        if count == 1: cstart = cend = cend + 1
        else:
            gap = count - 1
            cstart = cend + 1
            cend = cstart + gap
        if cend >= cstart:
            newTuple = (aaProperty, cstart, cend)
            new_AA_GROUPS.append(newTuple)

    for codon, start, end in new_AA_GROUPS:
        mesh_ax.text((end - start + 1) / 2 + start,
            row_number + 2.5, codon,
            horizontalalignment="center",
            verticalalignment="center")

        bar = Line2D([start + 0.125, end + 1 - 0.125],
                [row_number + 2, row_number + 2], color="black")
        bar.set_clip_on(False)
        mesh_ax.add_line(bar)


    WT_codon  = [str(wt_codon(wtfile, int(x*3-2))) for x in df['pos'].values]
    wt_aa     = [str(CODON2AA[codon][0]) for codon in WT_codon]

    # Add in mesh_ax y-label: aa coordinates in WT sequence:
    ypos = np.arange(len(df['pos'])) + 0.5
    mesh_ax.set_yticks(ypos)
    labelPos = list(map(str,list(map(int, df['pos']))))
    labelAA  = wt_aa
    labelCombine = ['-'.join(item) for item in zip(labelPos, labelAA)]
    mesh_ax.set_yticklabels(labelCombine, ha='right')

    # Add in deliminator horizontally:
    for i in range(0, df_raw.shape[0]):
        delimBar = Line2D([0, df_raw.shape[1]],[i, i],
                transform=mesh_ax.transData, color="white")
        delimBar.set_clip_on(False)
        mesh_ax.add_line(delimBar)

    # Add in WT label onto corresponding cell:
    WT_aa = [CODON2AA[x][0] for x in WT_codon]
    x_coordinate = [list(df.columns).index(x) if x in df.columns else np.nan for x in WT_aa]
    y_coordinate = list(range(df.shape[1]))

    wt_xy        = zip(x_coordinate, y_coordinate)
    for x, y in wt_xy:
        if not x == np.nan:
            mesh_ax.add_patch(Circle((x + 0.5, y + 0.5), .1666,
            fill=True, facecolor="black",
            edgecolor="none", alpha=0.8))

    # Make the figure cleaner by removing ticks:
    mesh_ax.tick_params(bottom=False, left=False)
    mesh_ax.get_xaxis().set_visible(False)
    cbar_ax.tick_params(bottom=False)

    # Add in SE if specified:
    if not df_se.shape == pd.DataFrame().shape: # Add in SE
        # Subset df_se to match df that might be truncated:
        for col in df_se.columns:
            if not col in df.columns: df_se.drop(columns=col, inplace=True)
        df_se = df_se[df_se['pos'].isin(df['pos'])]

        for row in range(len(df_se['pos'])):
            for col in range(len(df_se.columns)-1):
                se_value = df_se.iloc[row, col]
                corner_dist = (1-se_value)/2
                diag = Line2D([col + corner_dist, col + 1 - corner_dist],
                        [row + corner_dist, row + 1 - corner_dist], color="black")

                if se_value > 0.02 and df_raw.iloc[row, col] != 1000: # se_value below 0.02 will not be displayed so as totally depleted ones
                    mesh_ax.add_line(diag)

    pylab.savefig(outfilename)
