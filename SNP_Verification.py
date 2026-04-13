#!/usr/bin/env python3

#   AMRPlusPlus_SNP_Verification
#   Copyright (C) 2022  Nathalie Bonin
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see https://www.gnu.org/licenses/.

from SNP_Verification_Tools import Gene
from SNP_Verification_Processes import verify
from Bio import SeqIO
import pysam, sys, getopt, os, configparser, csv
import pandas as pd
import multiprocessing
from concurrent.futures import ProcessPoolExecutor as ppe

def parse_config():
# Define Command Line Arguments and Read Config
    configFile = "config.ini"
    config = configparser.ConfigParser()
    arguments_string = """

        -a: amrplusplus; is either 'true' or 'false'
        -c: config file; if this argument is used, must be the first listed
        -h: help
        -i: BAM input file
        -l: license disclaimer
        -o: main output folder
        -r: conditions for redistribution
        -t: threads for multiprocessing
        
        --mt_and_wt:            false by default, used in case of insertion leading to presence of both mt and wt; if true, mark as resistant; if false, mark as susceptible
        --output_reads:         false by default, output list of resistant reads per gene
        --detailed_output:      false by default, determines whether a more detailed output will be given; can be either 'false', 'true', 'all', or include a list of accessions seperated by commas
        --count_matrix:         count matrix that will be updated if amrplusplus is true
        --count_matrix_final:   the file where the updated count matrix will be found if amrplusplus is true
        
        """
    argList = []
    try:
        options, args = getopt.getopt(sys.argv[1:], "hlra:c:i:o:t:", ["mt_and_wt=", "detailed_output=", "output_reads=", "count_matrix=", "count_matrix_final="])
    except getopt.GetoptError:
        print("ERROR - this is the list of arguments recognized by the program:{}".format(arguments_string))
        sys.exit(-1)
    countMatrixFinal = False
    for i, (opt, arg) in enumerate(options):
        if opt == "-h":
            print("List of arguments:{}".format(arguments_string))
            sys.exit()
        elif opt == "-r":
            print("""
            
            AMRPlusPlus_SNP_Verification
            Copyright (C) 2022  Nathalie Bonin
            
            This program is free software: you can redistribute it and/or modify
            it under the terms of the GNU General Public License as published by
            the Free Software Foundation, either version 3 of the License, or
            (at your option) any later version.
            
            """)
            sys.exit()
        elif opt == "-l":
            print("""
            
            AMRPlusPlus_SNP_Verification
            Copyright (C) 2022  Nathalie Bonin
            
            This program is distributed in the hope that it will be useful,
            but WITHOUT ANY WARRANTY; without even the implied warranty of
            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
            GNU General Public License for more details.
            
            """)
            sys.exit()
        elif opt == "-c":
            if i > 0:
                print("ERROR: config file must be the first argument listed")
                sys.exit(-1)
            configFile = arg
            config.read(configFile)
        elif opt == "-i":
            if i == 0:
                config.read(configFile)
            config['FULL_FILE_NAMES']['SAMPLE'] = os.path.basename(arg).split('.')[0]
            config['SOURCE_FILES']['BAM_INPUT'] = arg
        elif opt == "-o":
            if i == 0:
                config.read(configFile)
            # Ensure output folder path ends with /
            output_folder = arg if arg.endswith('/') else arg + '/'
            config['FOLDERS']['MAIN_OUTPUT_FOLDER'] = output_folder
            # Store output prefix for naming files (basename without trailing slash)
            config['FULL_FILE_NAMES']['OUTPUT_PREFIX'] = os.path.basename(arg.rstrip('/'))
        elif opt == "-a":
            if i == 0:
                config.read(configFile)
            if arg not in ['true', 'false']:
                print("ERROR: '-a' argument can only be either 'true' or 'false'")
                sys.exit(-1)
            config['SETTINGS']['AMRPLUSPLUS'] = arg
        elif opt == "-t":
            if i == 0:
                config.read(configFile)
            config['SETTINGS']['THREADS'] = arg
        elif opt == "--mt_and_wt":
            if i == 0:
                config.read(configFile)
            if arg not in ['true', 'false']:
                print("ERROR: mt_and_wt argument can only be either 'true' or 'false'")
                sys.exit(-1)
            config['SETTINGS']['MT_AND_WT'] = arg
        elif opt =="--detailed_output":
            if i == 0:
                config.read(configFile)
            # FIX #1: Handle 'true', 'True', 'all' correctly
            # Only set argList if it's a comma-separated list of specific accessions
            if arg.lower() not in ['false']:
                config['SETTINGS']['DETAILED'] = "true"
                if arg.lower() not in ['true', 'all']:
                    # It's a comma-separated list of specific accessions to filter
                    argList = arg.split(',')
                # If 'true' or 'all', argList stays empty = no filtering, process all genes
        elif opt =="--output_reads":
            if i == 0:
                config.read(configFile)
            if arg not in ['true', 'false']:
                print("ERROR: output_reads argument can only be either 'true' or 'false'")
                sys.exit(-1)
            config['SETTINGS']['READS'] = arg
        elif opt == "--count_matrix":
            if i == 0:
                config.read(configFile)
            if not(countMatrixFinal):
                # Use only the basename for output, not the full path
                # This prevents issues when input is "../file.csv"
                config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'] = os.path.basename(arg)
            config['SOURCE_FILES']['COUNT_MATRIX'] = arg
        elif opt == "--count_matrix_final":
            if i == 0:
                config.read(configFile)
            config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'] = arg
            countMatrixFinal = True
        i += 1
    if len(options) == 0:
        config.read(configFile)

    # Add new information to config buffer
    config['FOLDERS']['SAMPLE_OUTPUT'] = (config['FOLDERS']['MAIN_OUTPUT_FOLDER'] + 
                                        config['FULL_FILE_NAMES']['SAMPLE']) + '/'
    config['FOLDERS']['SAMPLE_DETAILED_OUTPUT'] = (config['FOLDERS']['SAMPLE_OUTPUT'] + 
                                                config['FOLDERS']['DETAILED_FOLDER'])
    config['FOLDERS']['TEMP'] = os.path.dirname(config['TEMP_FILES']['TEMP_BAM_SORTED'])

    # Get output prefix for naming files (defaults to sample name if -o not specified)
    output_prefix = config['FULL_FILE_NAMES'].get('OUTPUT_PREFIX', config['FULL_FILE_NAMES']['SAMPLE'])
    
    # Use output prefix in count matrix filename
    config['FULL_FILE_NAMES']['COUNT_MATRIX_FINAL'] = (config['FOLDERS']['MAIN_OUTPUT_FOLDER'] + 
                                                    output_prefix + '_' + config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'])
    
    # Use output prefix in type output filenames
    config['FULL_FILE_NAMES']['NTYPE_OUTPUT'] = (config['FOLDERS']['SAMPLE_OUTPUT'] + 
                                                output_prefix + '_' + config['OUTPUT_FILES']['NORMAL_TYPE_OUTPUT'])
    config['FULL_FILE_NAMES']['FTYPE_OUTPUT'] = (config['FOLDERS']['SAMPLE_OUTPUT'] + 
                                                output_prefix + '_' + config['OUTPUT_FILES']['FRAMESHIFT_TYPE_OUTPUT'])
    config['FULL_FILE_NAMES']['HTYPE_OUTPUT'] = (config['FOLDERS']['SAMPLE_OUTPUT'] + 
                                                output_prefix + '_' + config['OUTPUT_FILES']['HYPERSUSCEPTIBLE_TYPE_OUTPUT'])
    config['FULL_FILE_NAMES']['STYPE_OUTPUT'] = (config['FOLDERS']['SAMPLE_OUTPUT'] + 
                                                output_prefix + '_' + config['OUTPUT_FILES']['SUPPRESSIBLE_TYPE_OUTPUT'])
    config['FULL_FILE_NAMES']['ITYPE_OUTPUT'] = (config['FOLDERS']['SAMPLE_OUTPUT'] + 
                                                output_prefix + '_' + config['OUTPUT_FILES']['INTRINSIC_TYPE_OUTPUT'])
    return (config, argList)

def dir_check(config):
    # Verify existance of folders 
    if not(os.path.exists(config['FOLDERS']['MAIN_OUTPUT_FOLDER'])):
        os.makedirs(config['FOLDERS']['MAIN_OUTPUT_FOLDER'])
    if not(os.path.exists(config['FOLDERS']['SAMPLE_DETAILED_OUTPUT'])):
        os.makedirs(config['FOLDERS']['SAMPLE_DETAILED_OUTPUT'])
    if not(os.path.exists(config['FOLDERS']['TEMP'])):
        os.makedirs(config['FOLDERS']['TEMP'])

def parse_snp_info(config):
    gene_dict = dict()
    # Get Variants Info
    for gene in SeqIO.parse(config['SOURCE_FILES']['SNP_INFO_FASTA'], 'fasta'):
        # Find index of last pipe ('|') before variant list
        index = -1
        for pipe in range(0, 5):
            index = gene.name[index+1:].find('|') + index + 1
        name = gene.name[:index]
        variants = gene.name[index+1:]

        # Store information required for creating Gene object
        gene_dict[name + "|RequiresSNPConfirmation"] = [name, gene.seq, variants]
    return gene_dict

def iterate(process_vars):

    gene_name = process_vars[0]
    gene_object = process_vars[1]
    config = process_vars[2]
    argList = process_vars[3]
    lock = process_vars[4]

    with pysam.AlignmentFile(config['TEMP_FILES']['TEMP_BAM_SORTED'], "r") as samfile:
        alignment_iterator = list(samfile.fetch(reference=gene_name))

    # Create Gene object
    if 'Must:' in gene_object[2]:
        if 'Nuc:' in gene_object[2]:            gene_variant = Gene.IntrinsicrRNA(gene_object[0],gene_object[1],gene_object[2])
        else:                                   gene_variant = Gene.IntrinsicProtein(gene_object[0],gene_object[1],gene_object[2])
    elif 'Hyper:' in gene_object[2]:            gene_variant = Gene.Hypersusceptible(gene_object[0],gene_object[1],gene_object[2])
    elif 'FS-' in gene_object[2]:
        if 'suppression' in gene_object[2]:     gene_variant = Gene.Suppressible(gene_object[0],gene_object[1],gene_object[2])
        elif 'MEG_6142' in gene_name:           gene_variant = Gene.NormalProtein(gene_object[0],gene_object[1],gene_object[2])
        else:                                   gene_variant = Gene.Frameshift(gene_object[0],gene_object[1],gene_object[2])
    else:
        if 'Nuc:' in gene_object[2]:            gene_variant = Gene.NormalrRNA(gene_object[0],gene_object[1],gene_object[2])
        else:                                   gene_variant = Gene.NormalProtein(gene_object[0],gene_object[1],gene_object[2])

    # Go through all alignments to gene
    DEBUGGING_MODE = config.getboolean('SETTINGS', 'DEBUGGING_MODE')
    if DEBUGGING_MODE:
        lock.acquire()
    for read in alignment_iterator:
        if (read.cigarstring == None):
            continue
        elif (len(argList) != 0) and (gene_variant.getName().split("|")[0] not in argList):
            continue
        verify(read, gene_variant, config)
        gene_variant.resetForNextRead()

    if DEBUGGING_MODE:
        lock.release()
    
    return {gene_name : gene_variant}

def process_genes(config, argList, gene_dict):
    processes = list()

    pysam.sort("-o", config['TEMP_FILES']['TEMP_BAM_SORTED'], config['SOURCE_FILES']['BAM_INPUT'])
    pysam.index(config['TEMP_FILES']['TEMP_BAM_SORTED'], config['TEMP_FILES']['TEMP_BAM_SORTED']+'.bai')

    m = multiprocessing.Manager()
    lock = m.Lock()

    for gene_name, gene_object in gene_dict.items():
        processes.append((gene_name, gene_object, config, argList, lock))

    with ppe(int(config['SETTINGS']['THREADS'])) as p:
        results = p.map(iterate, processes)

    return results


def appendGeneOutputInfo(name, output_info, csvwriter, countMatrix, config, gene=None, stats=None):
    """
    Appends gene output info to CSV and updates count matrix.
    
    PERCENTAGE-BASED CALCULATION:
    - percentage = reads_with_SNP_confirmed / reads_that_covered_any_SNP_position
    - new_count = percentage × original_AMR++_count
    
    This gives an accurate estimate of what fraction of total alignments 
    are from organisms with the resistance SNP.
    """
    new_row = [name]
    new_row.extend(output_info)
    csvwriter.writerow(new_row)

    # Update count matrix if previously found in AMR++
    if config.getboolean('SETTINGS', 'AMRPLUSPLUS') and countMatrix is not None:
        sample_col = config['FULL_FILE_NAMES']['SAMPLE']
        
        # Check if sample column exists
        if sample_col not in countMatrix.columns:
            if stats is not None and stats.get('sample_col_warning', 0) == 0:
                print(f"WARNING: Sample column '{sample_col}' not found in count matrix!")
                stats['sample_col_warning'] = 1
            return
        
        # Check if gene is in matrix
        matching_rows = countMatrix[countMatrix['gene_accession'] == name]
        
        if len(matching_rows) > 0:
            index = matching_rows.index[0]
            prevResCount = countMatrix.loc[index, sample_col]
            
            # Track that we found a matching gene
            if stats is not None:
                stats['genes_found_in_matrix'] = stats.get('genes_found_in_matrix', 0) + 1
            
            if prevResCount != 0:
                is_intrinsic = (len(output_info) == 10)
                reads_analyzed = output_info[0]  # Total reads that were analyzed for this gene
                
                # Get the actual count of reads that covered any SNP position
                reads_covering_snp = 0
                if gene is not None and hasattr(gene, 'getReadsCoveringSNP'):
                    reads_covering_snp = gene.getReadsCoveringSNP()
                
                if is_intrinsic:
                    reads_with_snp = output_info[1] + output_info[2]  # All + Some = resistant
                    reads_covering_snp_intrinsic = output_info[1] + output_info[2] + output_info[4]
                    if reads_covering_snp_intrinsic > reads_covering_snp:
                        reads_covering_snp = reads_covering_snp_intrinsic
                else:
                    reads_with_snp = output_info[1]  # Resistant count
                
                # Calculate percentage
                if reads_covering_snp > 0:
                    percentage = reads_with_snp / reads_covering_snp
                elif reads_analyzed > 0:
                    percentage = reads_with_snp / reads_analyzed
                    reads_covering_snp = reads_analyzed
                else:
                    percentage = 0.0
                
                # Apply percentage to original count
                newCount = int(round(percentage * prevResCount))
                
                # What would original code have done?
                if is_intrinsic:
                    original_would_be = output_info[1] + output_info[2] + output_info[5]
                else:
                    original_would_be = output_info[1]
                
                # Track statistics
                if stats is not None:
                    stats['total_with_counts'] = stats.get('total_with_counts', 0) + 1
                    
                    if reads_with_snp > 0:
                        stats['snp_confirmed'] = stats.get('snp_confirmed', 0) + 1
                    else:
                        stats['snp_not_confirmed'] = stats.get('snp_not_confirmed', 0) + 1
                    
                    # Collect gene-level stats for output file
                    if 'gene_stats' not in stats:
                        stats['gene_stats'] = []
                    
                    stats['gene_stats'].append({
                        'sample_name': config['FULL_FILE_NAMES']['SAMPLE'],
                        'gene_name': name,
                        'gene_type': 'Intrinsic' if is_intrinsic else 'Normal',
                        'original_amrplusplus_count': prevResCount,
                        'total_reads_analyzed': reads_analyzed,
                        'reads_covering_snp_position': reads_covering_snp,
                        'reads_with_snp_confirmed': reads_with_snp,
                        'percentage': percentage,
                        'original_code_would_set': original_would_be,
                        'new_count_percentage_based': newCount,
                        'snp_confirmed': reads_with_snp > 0
                    })
                
                countMatrix.loc[index, sample_col] = newCount
        else:
            # Gene not found in matrix
            if stats is not None:
                stats['genes_not_in_matrix'] = stats.get('genes_not_in_matrix', 0) + 1


def create_output(config, argList, gene_variant_dict):
    # Create output files and write headers
    with (open(config['FULL_FILE_NAMES']['NTYPE_OUTPUT'], "w") as outputN, 
        open(config['FULL_FILE_NAMES']['FTYPE_OUTPUT'], "w") as outputF, 
        open(config['FULL_FILE_NAMES']['HTYPE_OUTPUT'], "w") as outputH,
        open(config['FULL_FILE_NAMES']['STYPE_OUTPUT'], "w") as outputS,
        open(config['FULL_FILE_NAMES']['ITYPE_OUTPUT'], "w") as outputI):

        nf_header = ["Gene","Number of reads", "Resistant", "Missense",
                    "Insertion", "Deletion", "Previously recorded nonsense",
                    "N-tuple", "Nonstop", "12+bp indel", "12+ bp frameshift",
                    "Newly found nonsense", "Frameshift till end"]

        h_header = nf_header.copy()
        h_header.extend(['Hypersusceptible mutations + resistance-conferring mutations'])

        s_header = nf_header.copy()
        s_header.pop()
        s_header.extend(["Frameshift at end", "Suppressible frameshift at res 531", 
                        "Frameshift at res 531 that is not suppressible"])

        i_header = ["Gene", "Number of reads", "All", "Some", "None", "Mutations", "Acquired", 
                    "12+bp indel", "12+ bp frameshift", "Nonsense", "Frameshift till end"]

        nwriter = csv.writer(outputN, delimiter = ',')
        fwriter = csv.writer(outputF, delimiter= ',')
        hwriter = csv.writer(outputH, delimiter= ',')
        swriter = csv.writer(outputS, delimiter= ',')
        iwriter = csv.writer(outputI, delimiter= ',')

        nwriter.writerow(nf_header)
        fwriter.writerow(nf_header)
        hwriter.writerow(h_header)
        swriter.writerow(s_header)
        iwriter.writerow(i_header)

        # Retrieve count matrix data
        countMatrix = pd.read_csv(config['SOURCE_FILES']['COUNT_MATRIX']) if config.getboolean('SETTINGS', 'AMRPLUSPLUS') else None

        # Get output prefix for naming files
        output_prefix = config['FULL_FILE_NAMES'].get('OUTPUT_PREFIX', config['FULL_FILE_NAMES']['SAMPLE'])

        # Create new output file for resistant reads if requested
        if config.getboolean('SETTINGS', 'READS'):
            with open(config['FOLDERS']['SAMPLE_OUTPUT'] + output_prefix + '_resistant_reads.csv', "w") as readsOutput:
                readsOutput.write("Gene Header, List of Reads\n")

        # Run through results
        genes_processed = 0
        genes_with_reads = 0
        
        # Stats for tracking count changes
        stats = {}
        
        for name, gene in gene_variant_dict.items():
            if (len(argList) != 0) and (gene.getName().split("|")[0] not in argList):
                continue
            
            genes_processed += 1
            if gene.getOutputInfo()[0] > 0:
                genes_with_reads += 1
            
            # Add to correct output - pass gene object and stats for tracking
            tag = gene.getGeneTag()
            if   tag == 'N': appendGeneOutputInfo(name, gene.getOutputInfo(), nwriter, countMatrix, config, gene, stats)
            elif tag == 'F': appendGeneOutputInfo(name, gene.getOutputInfo(), fwriter, countMatrix, config, gene, stats)
            elif tag == 'H': appendGeneOutputInfo(name, gene.getOutputInfo(), hwriter, countMatrix, config, gene, stats)
            elif tag == 'S': appendGeneOutputInfo(name, gene.getOutputInfo(), swriter, countMatrix, config, gene, stats)
            else:            appendGeneOutputInfo(name, gene.getOutputInfo(), iwriter, countMatrix, config, gene, stats)

            # Print more detailed output if requested
            if config.getboolean('SETTINGS', 'DETAILED') and (gene.getOutputInfo()[0] > 0):
                with open(config['FOLDERS']['SAMPLE_DETAILED_OUTPUT'] + output_prefix + '_' + name + ".csv", "w") as detailedOutput:
                    gene.writeAdditionalInfo(detailedOutput)

            # Print resistant reads if requested
            if config.getboolean('SETTINGS', 'READS') and (gene.getOutputInfo()[0] > 0):
                with open(config['FOLDERS']['SAMPLE_OUTPUT'] + output_prefix + '_resistant_reads.csv', "a") as readsOutput:
                    gene.writeResistantReads(readsOutput)
            
            gene.clearOutputInfo()

        # FIX #2: Zero out RequiresSNPConfirmation genes that weren't in SNPinfo database
        # These genes have alignments but couldn't be verified because they lack SNP info
        if config.getboolean('SETTINGS', 'AMRPLUSPLUS'):
            sample_col = config['FULL_FILE_NAMES']['SAMPLE']
            zeroed_not_in_db = 0
            total_not_in_db = 0
            genes_zeroed_list = []
            
            for idx, row in countMatrix.iterrows():
                gene_name = row['gene_accession']
                if 'RequiresSNPConfirmation' in str(gene_name) and gene_name not in gene_variant_dict:
                    total_not_in_db += 1
                    original_count = row[sample_col]
                    if original_count != 0:
                        zeroed_not_in_db += 1
                        genes_zeroed_list.append({
                            'gene_name': gene_name,
                            'original_count': original_count,
                            'reason': 'Not in SNPinfo database'
                        })
                    countMatrix.loc[idx, sample_col] = 0
            
            # Write gene coverage stats to CSV file
            gene_stats_file = config['FOLDERS']['SAMPLE_OUTPUT'] + output_prefix + '_snp_coverage_stats.csv'
            if 'gene_stats' in stats and len(stats['gene_stats']) > 0:
                gene_stats_df = pd.DataFrame(stats['gene_stats'])
                gene_stats_df.to_csv(gene_stats_file, index=False)
            
            # Write summary stats to CSV file
            summary_file = config['FOLDERS']['SAMPLE_OUTPUT'] + output_prefix + '_snp_verification_summary.csv'
            summary_data = {
                'sample_name': [sample_col],
                'genes_in_snpinfo_database': [len(gene_variant_dict)],
                'genes_processed': [genes_processed],
                'genes_with_reads': [genes_with_reads],
                'genes_found_in_count_matrix': [stats.get('genes_found_in_matrix', 0)],
                'genes_not_found_in_matrix': [stats.get('genes_not_in_matrix', 0)],
                'genes_with_nonzero_counts': [stats.get('total_with_counts', 0)],
                'genes_with_snp_confirmed': [stats.get('snp_confirmed', 0)],
                'genes_with_snp_not_confirmed': [stats.get('snp_not_confirmed', 0)],
                'snpconfirmation_genes_not_in_db': [total_not_in_db],
                'genes_zeroed_not_in_db': [zeroed_not_in_db]
            }
            summary_df = pd.DataFrame(summary_data)
            summary_df.to_csv(summary_file, index=False)
            
            # Print brief summary to console
            print(f"\nSNP Verification Complete.")
            print(f"  Genes with SNP confirmed: {stats.get('snp_confirmed', 0)}")
            print(f"  Genes with SNP NOT confirmed: {stats.get('snp_not_confirmed', 0)}")
            print(f"  Summary written to: {summary_file}")
            print(f"  Gene stats written to: {gene_stats_file}")

    # Print count matrix
    if config.getboolean('SETTINGS', 'AMRPLUSPLUS'): countMatrix.to_csv(config['FULL_FILE_NAMES']['COUNT_MATRIX_FINAL'], index=False)
    sys.exit(0)


def main():
    print("SNP_Verification.py - Running with percentage-based count adjustment")
    
    config, argList = parse_config()
    dir_check(config)
    gene_dict = parse_snp_info(config)
    results = process_genes(config, argList, gene_dict)

    gene_variant_dict = dict()
    [gene_variant_dict.update(r) for r in results]

    create_output(config, argList, gene_variant_dict)

if __name__ == '__main__':
    main()
