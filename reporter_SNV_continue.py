#! /usr/bin/env python3

import time, os, sys, glob, gzip, subprocess
import multiprocessing as mp
import pandas as pd
from Bio import bgzf
from functools import partial 
from contextlib import contextmanager
from optparse import OptionParser
from scanner_core import misc
from curator_core import curator_io

@contextmanager
def poolcontext(*args, **kwargs):
	pool = mp.Pool(*args, **kwargs)
	yield pool
	pool.terminate()

def proc_longshot(CB, options): ## fixed. Index file was not being generated for one, but it was also expecting the output in a .tbi format, which is not the default output of `bcftools index`.
    bam_pref = os.path.join(options.tmp_dir, CB)
    vcf_pref = os.path.join(options.tmp_dir, options.longshot_o + '.' + CB)
    vcf_file = vcf_pref + '.vcf'
    vcf_gz_file = vcf_pref + '.vcf.gz'
    vcf_tbi_file = vcf_gz_file + '.tbi'

    # === skip if result already exists ===
    if os.path.isfile(vcf_gz_file):
        print(f"Skip {CB}: {vcf_gz_file} exists")
        return

    longshot_time = time.time()
    print(f'Longshot: Calling SNV in {CB} ...', flush=True)
    
    # 1. Run Longshot
    cmd = (f"{options.longshot} -F "
           f"--bam {bam_pref}.curated.minimap2.bam "
           f"--ref {options.ref_genome} "
           f"--min_cov {options.longshot_min_cov} "
           f"--min_alt_count {options.longshot_min_alt_count} "
           f"--out {vcf_file} "
           f"--sample_id {CB}")

    code_msg, out_msg, err_msg = curator_io.sys_run(cmd)
    
    # Log output
    with open(vcf_pref + '.log', 'wt') as fh:
        fh.write(out_msg.decode('utf-8', errors='replace'))
        fh.write(err_msg.decode('utf-8', errors='replace'))

    # CRITICAL: Check if longshot actually succeeded
    if code_msg != 0:
        raise RuntimeError(f"Longshot failed for {CB} with exit code {code_msg}")
    
    if not os.path.isfile(vcf_file) or os.path.getsize(vcf_file) == 0:
        raise RuntimeError(f"Longshot produced empty or missing VCF: {vcf_file}")

    # 2. Compression
    compress_cmd = [
        options.bcftools, 'view', '-Oz',
        '-o', vcf_gz_file, vcf_file
    ]

    with open(vcf_pref + '.compress.log', 'w') as log:
        result = subprocess.run(compress_cmd, stdout=log, stderr=subprocess.STDOUT)

    if result.returncode != 0:
        raise RuntimeError(f"VCF compression failed (code {result.returncode}): {vcf_file}")

    if not os.path.isfile(vcf_gz_file) or os.path.getsize(vcf_gz_file) == 0:
        raise RuntimeError(f"Compressed VCF is missing or empty: {vcf_gz_file}")

    # 3. Indexing
    # FIX: Use '-t' to force a .tbi index
    index_cmd = [options.bcftools, 'index', '-t', vcf_gz_file]

    with open(vcf_pref + '.index.log', 'w') as log:
        result = subprocess.run(index_cmd, stdout=log, stderr=subprocess.STDOUT)

    if result.returncode != 0:
        raise RuntimeError(f"VCF indexing failed (code {result.returncode}): {vcf_gz_file}")

    if not os.path.isfile(vcf_tbi_file):
        raise RuntimeError(f"Indexing succeeded but {vcf_tbi_file} was not created")

    # === remove temp files ===
    if not options.keep_meta:
        try:
            os.remove(vcf_file)
        except OSError:
            pass

    hours, minutes, seconds = misc.get_time_elapse(longshot_time)
    print(f"Longshot spent {hours}:{minutes}:{seconds:.2f} on {CB}", flush=True)

def filter_by_prevalence(options):
	counter_raw, counter_pass = 0, 0
	oh = bgzf.open(os.path.join(options.o_dir, options.o_pref) + '.filtered.vcf.gz', 'wt')
	fh = bgzf.open(os.path.join(options.o_dir, options.o_pref) + '.raw.vcf.gz', 'rt')
	lh = gzip.open(os.path.join(options.o_dir, options.o_snv_l), 'wt')
	while True:
		line = fh.readline()
		if not line:
			break

		if line.startswith('#'):
			oh.write(line)
		else:
			counter_raw += 1
			line_list = line.rstrip().split("\t")

#     0,   1,  2,   3,   4,    5,      6,    7,      8,                9 :: 
# CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, AAACAGCCAACTCGCG, AAACAGCCAATAATCC, AAACAGCCACGCAACT
			SNV_no, empty_no = 0, 0
			for col_idx in range(9, len(line_list), 1):
				if line_list[col_idx].split(':')[1] == '.':
					empty_no += 1
				else:
					SNV_no += 1
			if SNV_no / (SNV_no + empty_no) >= options.prevalence:
				counter_pass += 1
				oh.write(line)
				lh.write(line_list[0] + "\t" + line_list[1] + "\n")

	lh.close()
	fh.close()
	oh.close()

	options.counter_pass = counter_pass
	options.counter_raw  = counter_raw

	if not options.keep_meta:
		cmd = 'rm ' + os.path.join(options.o_dir, options.o_pref) + '.raw.vcf.gz'
		os.system(cmd)

def merge_longshot(CB_list, options):
	#===skip if file exist===
	if os.path.isfile(os.path.join(options.o_dir, options.o_pref) + '.filtered.vcf.gz'):
		print("Skip merging longshot result: " + os.path.join(options.o_dir, options.o_pref) + ".filtered.vcf.gz exists")
		return

	vcf_files = sorted(glob.glob(os.path.join(options.tmp_dir, options.longshot_o + '.*.vcf.gz')))
	if not vcf_files:
		print(f"Warning: No longshot VCF files found in {options.tmp_dir} with prefix {options.longshot_o}")
		return
	chunk_size = 500
	for i in range(0, len(vcf_files), chunk_size):
		chunk = vcf_files[i:i + chunk_size]
		with open(os.path.join(options.tmp_dir, f'splitted_vcf_list_{i // chunk_size}'), 'w') as f:
			f.write('\n'.join(chunk) + '\n')

	counter = 0
	f_list = glob.glob(os.path.join(options.tmp_dir, 'splitted_vcf_list_*'))
	for sp_f in f_list:
		counter += 1
		print('Merging part ' + str(counter) + ' of ' + str(len(f_list)) + ' ...', end = "\r")
		cmd = options.bcftools + ' merge -0 -l ' + sp_f + \
		      ' -Oz -o ' + sp_f + '.merge.vcf.gz &> /dev/null'
		os.system(cmd)

		cmd = options.tabix + ' --verbosity 0 -p vcf ' + sp_f + '.merge.vcf.gz'
		os.system(cmd)

	print("\nGenerating final matrix...")
	merge_vcfs = sorted(glob.glob(os.path.join(options.tmp_dir, 'splitted_vcf_list_*.merge.vcf.gz')))
	with open(os.path.join(options.tmp_dir, 'final_vcf_list'), 'w') as f:
		f.write('\n'.join(merge_vcfs) + '\n')

	cmd = options.bcftools + ' merge -0 -l ' + os.path.join(options.tmp_dir, 'final_vcf_list') + ' -Ov -o ' + os.path.join(options.o_dir, options.o_pref) + '.raw.vcf.gz'
	os.system(cmd)

	cmd = 'rm ' + os.path.join(options.tmp_dir, 'splitted_vcf_list_*')
	os.system(cmd)

	cmd = 'rm ' + os.path.join(options.tmp_dir, 'final_vcf_list')
	os.system(cmd)

	print('Filtering by prevalence: ' + str(options.prevalence) + ' ...')
	filter_by_prevalence(options)
	print("\nTotal SNVs no: " + str(options.counter_raw))
	print(str(options.counter_pass) + ' SNVs pass filtering...')

	cmd = options.tabix + ' --verbosity 0 -p vcf ' + os.path.join(options.o_dir, options.o_pref) + '.filtered.vcf.gz'
	os.system(cmd)
	print("Done.\n")

def proc_annovar(options):
	vcf_o = os.path.join(options.o_dir, options.o_name)

	#---unzip vcf---
	if options.o_name.endswith('.gz'):
		o_name_list = options.o_name.split('.')
		vcf_o = os.path.join(options.o_dir, '.'.join(o_name_list[0:(len(o_name_list)-1)]))
		cmd = 'gunzip -c ' + os.path.join(options.o_dir, options.o_name) + ' > ' + vcf_o
		os.system(cmd)

	#---proc annovar---
	cmd = 'perl ' + os.path.join(options.annovar, 'table_annovar.pl') + \
	      ' ' + vcf_o + ' ' + options.annovar_db + \
	      ' -buildver ' + options.annovar_gv + \
	      ' -out ' + options.o_dir + '/annovar -remove' + \
              ' -protocol ' + options.annovar_protocol + \
	      ' -operation ' + options.annovar_operation + \
	      ' -nastring . -polish -vcfinput -otherinfo'
	if options.annovar_xref:
		cmd += ' -xref ' + options.annovar_xref
	print('Execute ANNOVAR: ' + cmd)
	os.system(cmd)

	print("Done.\n")

def compute_read_no(CB, options):
	mpileup_time = time.time()
	print('Processing ' + CB + ' ...', flush = True)

	#===skip if file exist===
	if os.path.isfile(os.path.join(options.tmp_dir, CB) + '.mpileup'):
		print("Skip " + CB + ": " + os.path.join(options.tmp_dir, CB) + ".mpileup exists")
		return

	cmd = options.samtools + ' mpileup ' + os.path.join(options.tmp_dir, CB) + '.curated.minimap2.bam' + \
	      ' --reference ' + options.ref_genome + ' -l ' + os.path.join(options.o_dir, options.o_snv_l) + \
	      ' --no-output-ins --no-output-ins --no-output-del --no-output-del' + \
	      ' -Q ' + options.min_read_quality + ' -o ' + os.path.join(options.tmp_dir, CB) + '.mpileup'
	os.system(cmd)

	hours, minutes, seconds = misc.get_time_elapse(mpileup_time)
	print("Samtools mpileup spends %d:%d:%.2f on %s" % (hours, minutes, seconds, CB), flush = True)

def merge_mpileup(CB_list, options):
	#=== init df ===
	pileup_obj = dict()
	fh = bgzf.open(os.path.join(options.o_dir, options.o_pref) + '.filtered.vcf.gz', 'rt')
	while True:
		line = fh.readline()
		if not line:
			break
		if line.startswith('##'):
			continue
		elif line.startswith('#'):
			line_list = line.rstrip().split("\t")
			line_list[0] = line_list[0].split("#")[1]
			options.vcf_header = line_list
			continue
		line_list = line.split("\t")
		pileup_obj[str(line_list[0] + ':' + line_list[1])] = {'CHROM': line_list[0], 'POS': line_list[1], 'REF': line_list[3], 'ALT': line_list[4]}
	fh.close()

	pileup_df = pd.DataFrame.from_dict(pileup_obj, orient='index')
	pileup_df['POS'] = pileup_df['POS'].apply(str)

	#=== read real DP ===
	col_list = ['CHROM', 'POS', 'REF', 'RNO', 'MATCH', 'QUALITY']
	CB_no, CB_total = 0, len(CB_list)
	for CB in CB_list:
		CB_no += 1
		print("Processing " + str(CB_no) + " of " + str(CB_total) + " ... ", end = "\r", flush = True)
		dp_f  = os.path.join(options.tmp_dir, CB) + '.mpileup'
		dp_df = pd.read_csv(dp_f, sep = '\t', names = col_list, quoting = 3)
		dp_df['CHROM'] = dp_df['CHROM'].apply(str)
		dp_df['POS']   = dp_df['POS'].apply(str)

		dp_df = dp_df.merge(pileup_df.loc[:, ['CHROM', 'POS', 'ALT']], how = 'left', on = ['CHROM', 'POS'])

		dp_df[CB] = ""
		dp_df['REF_no'] = dp_df['MATCH'].str.count("[\\.\\,]").astype('Int64')
		dp_df['A']      = dp_df['MATCH'].str.count('A').astype('Int64') + dp_df['MATCH'].str.count('a').astype('Int64')
		dp_df['T']      = dp_df['MATCH'].str.count('T').astype('Int64') + dp_df['MATCH'].str.count('t').astype('Int64')
		dp_df['C']      = dp_df['MATCH'].str.count('C').astype('Int64') + dp_df['MATCH'].str.count('c').astype('Int64')
		dp_df['G']      = dp_df['MATCH'].str.count('G').astype('Int64') + dp_df['MATCH'].str.count('g').astype('Int64')

		for idx in range(0, dp_df.shape[0]):
			dp_str = str(dp_df.loc[idx, 'REF_no'])
			alt_list = dp_df.loc[idx, 'ALT'].split(',')
			for Nu in alt_list:
				dp_str += '/' + str(dp_df.loc[idx, Nu])
			dp_df.loc[idx, CB] = dp_str

		pileup_df = pileup_df.merge(dp_df.loc[:, ['CHROM', 'POS', 'REF', 'ALT', CB]], how = 'left', on = ['CHROM', 'POS', 'REF', 'ALT'], )

		if not options.keep_meta:
			cmd = 'rm ' + os.path.join(options.tmp_dir, CB) + '.mpileup'
			os.system(cmd)

	#===separate rows having multiple ALT===
	pileup_sinalt_df = pileup_df.loc[~pileup_df.ALT.str.contains(',')]
	pileup_mulalt_df = pileup_df.loc[ pileup_df.ALT.str.contains(',')]

	#=== fill na ===
		#---singular ATL---
	pileup_sinalt_df = pileup_sinalt_df.fillna('0/0')

		#---multiple ALT---
	for i in range(0, pileup_mulalt_df.shape[0]):
		na_str = "0" + "/0" * (pileup_mulalt_df.iloc[i, :].ALT.count(",") + 1)
		pileup_mulalt_df.iloc[i, :] = pileup_mulalt_df.iloc[i, :].fillna(na_str)

	#=== merge back ===
	pileup_df = pd.concat([pileup_sinalt_df, pileup_mulalt_df]).sort_values(["CHROM", "POS"]).reset_index(drop = True)

	pileup_df.to_csv(os.path.join(options.o_dir, options.o_snv_dp), sep = '\t', header = True, index = False, compression = 'gzip')

	return pileup_df, options

def correct_vcf(pileup_df, options):
	pileup_obj = dict()
	pileup_df_idx = -1
	fh = bgzf.open(os.path.join(options.o_dir, options.o_pref) + '.filtered.vcf.gz', 'rt')
	oh = bgzf.open(os.path.join(options.o_dir, options.o_name), 'wt')
	while True:
		line = fh.readline()
		if not line:
			break
		if line.startswith('##'):
			oh.write(line)
			continue
		elif line.startswith('#'):
			oh.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth of reads passing MAPQ filter">' + "\n")
			oh.write(line)
			continue

		line_list = line.rstrip().split("\t")

		sel_row = pileup_df.loc[(pileup_df['CHROM'] == line_list[0]) & (pileup_df['POS'] == line_list[1])]
		if sel_row.shape[0] > 0:
			for idx in range(0, len(options.vcf_header)):
				if idx == 8:
					line_list[idx] += ':DP'
				elif idx > 8:
					cell_list = line_list[idx].split(':')
					GT_str = mod_GT(sel_row[options.vcf_header[idx]].values[0])
					line_list[idx] = GT_str + ':' + ':'.join(cell_list[1:len(cell_list)]) + ':' + sel_row[options.vcf_header[idx]].values[0]
			pileup_df = pileup_df.drop(sel_row.index)
			oh.write("\t".join(line_list) + "\n")
		else:
			#no data in DP
			print("\n\nno data in DP file!")
			print("VCF: " + line)

	oh.close()
	fh.close()

	cmd = options.tabix + ' --verbosity 0 -p vcf ' + os.path.join(options.o_dir, options.o_name)
	os.system(cmd)

	if not options.keep_meta:
		cmd = 'rm ' + os.path.join(options.o_dir, options.o_pref) + '.filtered.vcf.gz'
		os.system(cmd)
		cmd = 'rm ' + os.path.join(options.o_dir, options.o_pref) + '.filtered.vcf.gz.tbi'
		os.system(cmd)

def mod_GT(DP_str):
# https://samtools.github.io/hts-specs/VCFv4.2.pdf
# unphased, always use '/'
# The allele values are 0 for the reference allele (what is in the REF field), 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on.
	if (DP_str == './.') | (DP_str == '<NA>/<NA>'):
		return './.'
	elif (DP_str == '././.') | (DP_str == '<NA>/<NA>/<NA>'):
		return '././.'
	elif (DP_str == './././.') | (DP_str == '<NA>/<NA>/<NA>/<NA>'):
		return './././.'
	else:
		DP_list = DP_str.split('/')
		res_str, GT_counter = '', 0
		for DP_idx in range(0, len(DP_list)):
			if int(float(DP_list[DP_idx])) > 0:
				GT_counter += 1
				res_str += '/' + str(DP_idx)
		if GT_counter == 0:
			res_str = '/./.'
		elif GT_counter == 1:
			res_str += res_str
		return res_str[1:]

if __name__ == "__main__":

	#===get params===
	parser = OptionParser()

	parser.add_option("-d",           dest = "o_dir",     default = "scNanoGPS_res",
	                  help = "Output directory name. "
	                         "Default: scNanoGPS_res")
	parser.add_option("--tmp_dir",    dest = "tmp_dir",    default = "tmp",
	                  help = "Temporary folder name. "
	                         "Default: tmp")
	parser.add_option("--CB_file",  dest = "CB_file",    default = "filtered_barcode_list.txt",
                          help = "File name for filtered barcode list. "
                                 "Default: filtered_barcode_list.txt")
	parser.add_option("--ref_genome", dest = "ref_genome", default = None,
	                   help = "* Required ! "
	                          "File for reference genome.")
	parser.add_option("--longshot_min_cov", dest = "longshot_min_cov", default = "2",
                          help = "Minimal coverage for longshot. "
                                 "Default: 2", type = "string")
	parser.add_option("--longshot_min_alt_count", dest = "longshot_min_alt_count", default = "2",
                          help = "Minimal alternative count for longshot. "
                                 "Default: 2", type = "string")
	parser.add_option("--longshot_o", dest = "longshot_o", default = "longshot.output",
                          help = "Prefix of LongShot output VCF file. "
                                 "Default: longshot.output")
	parser.add_option("-o",           dest = "o_name",     default = "SNV_matrix.vcf.gz",
	                  help = "Result SNV matrix file name. Must be ended with .vcf.gz. "
	                         "Default: SNV_matrix.vcf.gz")
	parser.add_option("--log",        dest = "log_f_name", default = "logs/reporter_SNV.log.txt",
	                  help = "Log file name. "
	                         "Default: logs/reporter_SNV.log.txt")
	parser.add_option("--o_snv_l",    dest = "o_snv_l",    default = "filtered_SNV_position_list.tsv.gz",
                          help = "Filtered SNVs position list. "
                                 "Default: filtered_SNV_position_list.tsv.gz")
	parser.add_option("--o_snv_dp",   dest = "o_snv_dp",   default = "SNV_dp_matrix.tsv.gz",
                          help = "SNVs depth matrix. "
                                 "Default: SNV_dp_matrix.tsv.gz")
	parser.add_option("-t",           dest = "ncores",     default = 1,
	                  help = "Number of cores for program running. "
	                         "Default: 1", type = "int")
	parser.add_option("--samtools",   dest = "samtools",   default = "samtools",
                          help = "Path to samtools. "
                                 "Default: samtools")
	parser.add_option("--bcftools",   dest = "bcftools",   default = "bcftools",
                          help = "Path to bcftools. "
                                 "Default: bcftools")
	parser.add_option("--tabix",      dest = "tabix",      default = "tabix",
                          help = "Path to tabix. "
                                 "Default: tabix")
	parser.add_option("--longshot",   dest = "longshot",   default = "longshot",
	                  help = "Path to longshot. "
	                         "Default: None")
	parser.add_option("--prevalence", dest = "prevalence", default = 0.01,
                          help = "SNV prevalence rate. "
                                 "Default: 0.01", type = "float")
	parser.add_option("--min_read_quality", dest = "min_read_quality", default = "0",
                          help = "Minimal read quality. "
                                 "Default: 0", type = "string")
	parser.add_option("--annovar",    dest = "annovar",  default = None,
                          help = "Path to Annovar. "
                                 "Default: None")
	parser.add_option("--annovar_db", dest = "annovar_db", default = None,
                          help = "Path to Annovar database. "
                                 "Default: None")
	parser.add_option("--annovar_gv", dest = "annovar_gv", default = "hg38",
                          help = "Annovar database genome version. "
                                 "Default: hg38")
	parser.add_option("--annovar_protocol", dest = "annovar_protocol", default = "refGene,cytoBand,gnomad30_genome,avsnp150,dbnsfp42c",
                          help = "Protocols of Annovar database. "
                                 "Default: refGene,cytoBand,gnomad30_genome,avsnp150,dbnsfp42c")
	parser.add_option("--annovar_operation", dest = "annovar_operation", default = "gx,r,f,f,f",
                          help = "Operation of Annovar database. "
                                 "Default: gx,r,f,f,f")
	parser.add_option("--annovar_xref", dest = "annovar_xref", default = None,
                          help = "Path to Omim xref. "
                                 "Default: None")
	parser.add_option("--keep_meta",    dest = "keep_meta",   default = None,
                          help = "Keep meta files. "
                                 "Default: None")
	options, arguments = parser.parse_args()

	#===pre-check===
	termination = False
	if not options.ref_genome or \
	   not os.path.isfile(options.ref_genome):
		print("\nCannot find reference genome file: " + str(options.ref_genome) + "\n")
		termination = True
	if not options.ref_genome or \
	   not os.path.isfile(options.ref_genome + ".fai"):
		print("\nCannot find reference genome index: " + str(str(options.ref_genome) + ".fai") + "\n")
		print("Please use \"samtools faidx\" to create reference genome index file\n")
		termination = True
	if not os.path.isdir(options.o_dir):
		print("\nOutput directory is not exist: "    + options.o_dir   + "\n")
		termination = True
	if not os.path.isdir(options.tmp_dir):
		print("\nTemporary directory is not exist: " + options.tmp_dir + "\n")
		print(options.tmp_dir)
		termination = True
	if not os.path.isfile(os.path.join(options.o_dir, options.CB_file)):
		options.CB_file = None
	if not options.o_name.endswith(".vcf.gz"):
		print("\nThe given output SNV matrix file name: " + options.o_name)
		print("Output SNV matrix file name must be ended with .vcf.gz\n")
		termination = True
	o_name_list = options.o_name.split('.')
	options.o_pref = '.'.join(o_name_list[0:(len(o_name_list) - 2)])

	cmd = "which " + options.samtools
	code_msg, out_msg, err_msg = curator_io.sys_run(cmd)
	options.samtools = out_msg.decode("utf-8").rstrip()
	if not os.path.isfile(options.samtools):
		print("\nCannot find samtools: " + options.samtools)
		termination = True

	cmd = "which " + options.bcftools
	code_msg, out_msg, err_msg = curator_io.sys_run(cmd)
	options.bcftools = out_msg.decode("utf-8").rstrip()
	if not os.path.isfile(options.bcftools):
		print("\nCannot find bcftools: " + options.bcftools)
		termination = True

	cmd = "which " + options.tabix
	code_msg, out_msg, err_msg = curator_io.sys_run(cmd)
	options.tabix = out_msg.decode("utf-8").rstrip()
	if not os.path.isfile(options.tabix):
		print("\nCannot find tabix: " + options.tabix)
		termination = True

	cmd = "which " + options.longshot
	code_msg, out_msg, err_msg = curator_io.sys_run(cmd)
	options.longshot = out_msg.decode("utf-8").rstrip()
	if not os.path.isfile(options.longshot):
		print("\nCannot find longshot: " + options.longshot)
		termination = True

	options.run_annovar = True
	if not options.annovar or \
	   not os.path.isfile(options.annovar):
		options.run_annovar = False
	else:
		print("Annovar found at: " + options.annovar)

	if termination:
		parser.print_help()
		sys.exit(1)

	logger = open(options.log_f_name, "wt")
	logger.write(" ".join(sys.argv) + "\n\n")
	logger.write("Reference genome file: " + options.ref_genome + "\n")
	logger.write("Output directory:      " + options.o_dir + "\n")
	logger.write("Temporary directory:   " + options.tmp_dir + "\n")
	logger.write("Filtered barcode list: " + str(options.CB_file) + "\n")
	logger.write("Path to samtools:      " + options.samtools + "\n")
	logger.write("Path to Bcftools:      " + options.bcftools + "\n")
	logger.write("Path to tabix:         " + options.tabix + "\n")
	logger.write("Path to longshot:      " + options.longshot + "\n")
	logger.write("Prevalence filter:     " + str(options.prevalence) + "\n")
	logger.write("Minimal coverage:      " + options.longshot_min_cov + "\n")
	logger.write("Minimal alt. count:    " + options.longshot_min_alt_count + "\n")
	logger.write("Minimal read quality:  " + options.min_read_quality + "\n")
	if options.run_annovar:
		logger.write("Path to Annovar:       " + options.annovar + "\n")
	logger.close()

	#=== set env variables
	os.environ["OMP_NUM_THREADS"]        = str(options.ncores)
	os.environ["OPENBLAS_NUM_THREADS"]   = str(options.ncores)
	os.environ["MKL_NUM_THREADS"]        = str(options.ncores)
	os.environ["VECLIB_MAXIMUM_THREADS"] = str(options.ncores)
	os.environ["NUMEXPR_NUM_THREADS"]    = str(options.ncores)

	import numpy as np

	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), "\n", flush = True)

	start_time = time.time()

	#===load CB list===
	print("Loading CB list...", flush = True)
	bam_list, CB_list = [], []
	if not options.CB_file:
		bam_list = glob.glob(os.path.join(options.tmp_dir, "*.curated.minimap2.bam"))
		for bam_f in bam_list:
			split_f_name = bam_f.split('/')
			BC = split_f_name[len(split_f_name) - 1].split('.')
			CB_list.append(BC[0])
	else:
		with open(os.path.join(options.o_dir, options.CB_file), "rt") as CBF:
			while True:
				CB_name = CBF.readline().rstrip()
				if not CB_name:
					break
				CB_list.append(CB_name)
		bam_list = [x + ".curated.minimap2.bam" for x in CB_list]

	print("\nMerge Longshot result...", flush = True)
	merge_longshot(CB_list, options)
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(start_time)
	misc.report_time_elapse(hours, minutes, seconds)

	print("\nCounting reads no. supporting SNVs...", flush = True)

	with poolcontext(processes = options.ncores) as pool:
		pool.map(partial(compute_read_no, options = options), CB_list)

	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(start_time)
	misc.report_time_elapse(hours, minutes, seconds)

	print("\nGenerating SNVs depth table...", flush = True)
	pileup_df, options = merge_mpileup(CB_list, options)
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(start_time)
	misc.report_time_elapse(hours, minutes, seconds)

	print("\nCorrect VCF genotype...", flush = True)
	correct_vcf(pileup_df, options)
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(start_time)
	misc.report_time_elapse(hours, minutes, seconds)

	if not options.keep_meta:
		os.system("mv " + options.o_snv_l + " " + options.tmp_dir)

	print()

	#===ANNOVAR===
	if options.annovar:
		start_time = time.time()
		print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
		print("Annotating VCF by ANNOVAR...\n", flush = True)
		proc_annovar(options)
		print("\nANNOVAR Done !\n")
		print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
		hours, minutes, seconds = misc.get_time_elapse(start_time)
		misc.report_time_elapse(hours, minutes, seconds)

		print()

	print()
