def get_UMI_counting(in_df, options):
    import numpy as np
    import math
    import pandas

    log10_res = math.ceil(-math.log10(options.smooth_res))

    print("Generating UMI counting table ...", flush = True)
    df = in_df.loc[in_df['BC'].str.len() == options.BC_len, ['BC', 'UMI']].groupby('BC').nunique().sort_values(by = 'UMI', ascending = False).reset_index()
    df['idx'] = df.index.values + 1

    print("Calculating log10(slope) ...", flush = True)
    df['log10_idx'] = np.log10(df['idx'])
    df['log10_UMI'] = np.log10(df['UMI'])

    df['r_log10_idx'] = np.log10(df['idx']).round(decimals = log10_res)
    df['r_log10_UMI'] = np.log10(df['UMI']).round(decimals = (log10_res + 1))
    df['r_log10_UMI_max'] = df.groupby('r_log10_idx')['r_log10_UMI'].transform('max')
    df['r_log10_UMI_min'] = df.groupby('r_log10_idx')['r_log10_UMI'].transform('min')

    log10_idx_array     = np.asarray(df.loc[:, ["r_log10_idx"]])
    log10_UMI_array_max = np.asarray(df.loc[:, ["r_log10_UMI_max"]])
    log10_UMI_array_min = np.asarray(df.loc[:, ["r_log10_UMI_min"]])

    tmp_array = log10_idx_array[1:log10_idx_array.shape[0], 0] - log10_idx_array[0:log10_idx_array.shape[0] - 1, 0]
    tmp_array[tmp_array < options.smooth_res] = options.smooth_res

    df['log10_slope'] = np.append(np.nan, \
                        (log10_UMI_array_min[1:log10_UMI_array_min.shape[0],     0] -  \
                         log10_UMI_array_max[0:log10_UMI_array_max.shape[0] - 1, 0]) / \
                        tmp_array * -1)
    df.loc[df['UMI'] < options.min_read_no, 'log10_slope'] = 0

    return df

def estimate_cell_no(df, options):
    import numpy as np

    df = df.iloc[(options.min_cellno - 1):]

    log10_idx_ori = np.max(df.loc[df['log10_slope'] == np.nanmax(df.loc[df['UMI'] > options.min_read_no, 'log10_slope']), "log10_idx"])
    idx_ori       = np.max(df.loc[df['log10_idx'] == log10_idx_ori, "idx"])

    #===including 10% more CB===
    log10_idx     = np.max(df.loc[df["idx"] == round(idx_ori * (1 + options.CB_no_ext)), "log10_idx"])
    idx           = np.max(df.loc[df["idx"] == round(idx_ori * (1 + options.CB_no_ext)), "idx"])

    return idx_ori, idx

def batch_seq_comp(query, target, options):
    import time, distance, os

    print("Calculating Levenshtein distance on " + str(query[0]) + ": " + query[1] + " ...", flush = True)
    start_time = time.time()

# target:
#       rep_id                BC
# 0          1  CTACGAAGTGATGAGG

#   if not options.whitelist:
#       target = target.loc[target["rep_id"] > query[0], :].copy()

    target.loc[:, "id"]       = query[0]
    target.loc[:, "query"]    = query[1]
    target.loc[:, "distance"] = target.apply(lambda x: distance.levenshtein(x["BC"], x["query"]), axis = 1)

    tmp_f = os.path.join(options.tmp_dir, "assigner_tmp_") + str(query[0]) + ".tsv"

    res = target.loc[target["distance"] <= options.CB_mrg_thr, ["rep_id", "id", "distance"]].copy()
    if res.shape[0] == 1:
        res.to_csv(tmp_f, header = None, index = None, sep = "\t")

    time_elapse = time.time() - start_time
    hours   = time_elapse // 3600
    rest_t  = time_elapse % 3600
    minutes = rest_t // 60
    seconds = rest_t % 60
    print("Calc. LD on " + query[1] + " spent %d : %d : %.2f" % (hours, minutes, seconds), flush = True)

    return 1

def merge_cb(options):
    import glob
    import os
    import pandas as pd
    import subprocess

    dist_tmp = options.CB_mrg_dist + ".tmp"

    # Always create the output table with its header, even when no
    # unambiguous barcode matches were found.
    with open(dist_tmp, "wt") as fh:
        fh.write("rep_id\tNanopore_id\tdistance\n")

    tmp_files = sorted(
        glob.glob(os.path.join(options.tmp_dir, "assigner_tmp_*.tsv"))
    )

    if tmp_files:
        with open(dist_tmp, "at") as out:
            subprocess.run(["cat", *tmp_files], stdout=out, check=True)

        for tmp_file in tmp_files:
            os.remove(tmp_file)
    else:
        print(
            "Warning: no unambiguous barcode matches were found; "
            "no assigner_tmp_*.tsv files were generated.",
            flush=True,
        )

    if options.CB_mrg_dist_compression:
        subprocess.run(["gzip", "-f", dist_tmp], check=True)
        os.replace(dist_tmp + ".gz", options.CB_mrg_dist)
    else:
        os.replace(dist_tmp, options.CB_mrg_dist)

    dist_df = pd.read_csv(
        options.CB_mrg_dist,
        sep="\t",
        header=0,
        compression=options.CB_mrg_dist_compression,
    )

    return dist_df
