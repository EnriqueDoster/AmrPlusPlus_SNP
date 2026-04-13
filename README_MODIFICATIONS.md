# SNP_Verification Modified Files - Final Version

This document describes all modifications made to the AmrPlusPlus_SNP verification software.

---

## Summary of Changes

| Category | Description |
|----------|-------------|
| Bug Fix #1 | Fixed `--detailed_output=True` parsing |
| Bug Fix #2 | Zero out unverified RequiresSNPConfirmation genes |
| Bug Fix #3 | Fixed path handling for `--count_matrix` and `-o` arguments |
| New Feature | Percentage-based count calculation with SNP position coverage tracking |
| New Feature | CSV output files for gene coverage statistics and run summary |
| New Feature | Output files now include `-o` flag name as prefix |

---

## Output File Naming

All output files now include the name from the `-o` flag as a prefix. For example, if you run:

```bash
python3 SNP_Verification.py -o my_sample_output -i reads.bam --count_matrix counts.csv
```

The output files will be named:
- `my_sample_output/my_sample_output_AMR_analytic_matrix.csv` (count matrix)
- `my_sample_output/<sample>/my_sample_output_NormalOutput.csv` (N-type output)
- `my_sample_output/<sample>/my_sample_output_snp_coverage_stats.csv`
- `my_sample_output/<sample>/my_sample_output_snp_verification_summary.csv`
- `my_sample_output/<sample>/my_sample_output_resistant_reads.csv`
- `my_sample_output/<sample>/detailed/my_sample_output_<gene_name>.csv` (detailed outputs)

This makes it easier to identify output files when processing multiple samples.

---

## Files Modified

### 1. SNP_Verification.py

**Location:** `SNP_Verification.py` (main script)

#### Fix #1: `--detailed_output` Parsing (Lines 121-131)

**Problem:** Using `--detailed_output=True` caused all genes to be skipped because "True" was treated as an accession filter list.

**Original code:**
```python
if arg != "false":
    config['SETTINGS']['DETAILED'] = "true"
    if arg != "all":
        argList = arg.split(',')  # "True" becomes ['True'] - filters out everything
```

**Fixed code (Lines 121-131):**
```python
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
```

#### Fix #2: Path Handling for `-o` Argument (Lines 99-102)

**Problem:** Output folder path without trailing `/` caused path construction issues.

**Fixed code:**
```python
elif opt == "-o":
    if i == 0:
        config.read(configFile)
    # Ensure output folder path ends with /
    output_folder = arg if arg.endswith('/') else arg + '/'
    config['FOLDERS']['MAIN_OUTPUT_FOLDER'] = output_folder
```

#### Fix #3: Path Handling for `--count_matrix` Argument (Lines 139-146)

**Problem:** Using `--count_matrix ../file.csv` created invalid output path like `output../file.csv`.

**Fixed code:**
```python
elif opt == "--count_matrix":
    if i == 0:
        config.read(configFile)
    if not(countMatrixFinal):
        # Use only the basename for output, not the full path
        # This prevents issues when input is "../file.csv"
        config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'] = os.path.basename(arg)
    config['SOURCE_FILES']['COUNT_MATRIX'] = arg
```

#### Fix #4: Directory Creation (Lines 177-184)

**Added:** Create `MAIN_OUTPUT_FOLDER` if it doesn't exist.

```python
def dir_check(config):
    # Verify existance of folders 
    if not(os.path.exists(config['FOLDERS']['MAIN_OUTPUT_FOLDER'])):
        os.makedirs(config['FOLDERS']['MAIN_OUTPUT_FOLDER'])
    if not(os.path.exists(config['FOLDERS']['SAMPLE_DETAILED_OUTPUT'])):
        os.makedirs(config['FOLDERS']['SAMPLE_DETAILED_OUTPUT'])
    if not(os.path.exists(config['FOLDERS']['TEMP'])):
        os.makedirs(config['FOLDERS']['TEMP'])
```

#### New Feature: Percentage-Based Count Calculation (Lines 260-365)

**Completely rewritten `appendGeneOutputInfo()` function.**

**Key changes:**
- Uses `gene.getReadsCoveringSNP()` to get accurate denominator
- Calculates: `percentage = reads_with_snp / reads_covering_snp_position`
- Applies: `new_count = percentage × original_amrplusplus_count`
- Collects per-gene statistics for CSV output

```python
def appendGeneOutputInfo(name, output_info, csvwriter, countMatrix, config, gene=None, stats=None):
    """
    PERCENTAGE-BASED CALCULATION:
    - percentage = reads_with_SNP_confirmed / reads_that_covered_any_SNP_position
    - new_count = percentage × original_AMR++_count
    """
    # ... (see full implementation in file)
    
    # Get the actual count of reads that covered any SNP position
    reads_covering_snp = 0
    if gene is not None and hasattr(gene, 'getReadsCoveringSNP'):
        reads_covering_snp = gene.getReadsCoveringSNP()
    
    # Calculate percentage
    if reads_covering_snp > 0:
        percentage = reads_with_snp / reads_covering_snp
    
    # Apply percentage to original count
    newCount = int(round(percentage * prevResCount))
```

#### New Feature: CSV Output Files (Lines 466-497)

**Added code to write:**
1. `snp_coverage_stats.csv` - Per-gene statistics
2. `snp_verification_summary.csv` - Run summary

```python
# Write gene coverage stats to CSV file
gene_stats_file = config['FOLDERS']['SAMPLE_OUTPUT'] + 'snp_coverage_stats.csv'
if 'gene_stats' in stats and len(stats['gene_stats']) > 0:
    gene_stats_df = pd.DataFrame(stats['gene_stats'])
    gene_stats_df.to_csv(gene_stats_file, index=False)

# Write summary stats to CSV file
summary_file = config['FOLDERS']['SAMPLE_OUTPUT'] + 'snp_verification_summary.csv'
summary_data = {
    'sample_name': [sample_col],
    'genes_in_snpinfo_database': [len(gene_variant_dict)],
    # ... (additional fields)
}
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(summary_file, index=False)
```

#### Fix #5: Zero Out Unverified Genes (Lines 444-464)

**Added:** Zero out RequiresSNPConfirmation genes that aren't in the SNPinfo database.

```python
for idx, row in countMatrix.iterrows():
    gene_name = row['gene_accession']
    if 'RequiresSNPConfirmation' in str(gene_name) and gene_name not in gene_variant_dict:
        total_not_in_db += 1
        original_count = row[sample_col]
        if original_count != 0:
            zeroed_not_in_db += 1
        countMatrix.loc[idx, sample_col] = 0
```

---

### 2. SNP_Verification_Tools/Gene.py

**Location:** `SNP_Verification_Tools/Gene.py`

#### New Properties (Lines 61-68)

Added tracking properties to the `Gene` class `__init__` method:

```python
# NEW: Track confirmed N-tuple variants across reads
# Key = variant string (e.g., "Mult:Mis:E52K;Mis:I77L")
# Value = set of confirmed mutation positions for that variant
this.confirmed_ntuple_positions = dict()

# NEW: Track reads that covered any SNP position (for percentage calculation)
this.reads_covering_snp = 0
this.current_read_covered_snp = False
```

#### Updated clearOutputInfo() (Lines 121-129)

Added clearing of new tracking properties:

```python
def clearOutputInfo(this):
    for index in range(len(this.output_info)):
        this.output_info[index] = 0
    this.additional_info.clear()
    # NEW: Also clear N-tuple tracking
    this.confirmed_ntuple_positions.clear()
    # NEW: Also clear SNP coverage tracking
    this.reads_covering_snp = 0
    this.current_read_covered_snp = False
```

#### New Methods (Lines 131-146)

Added three new methods for SNP coverage tracking:

```python
# NEW: Mark that the current read covered at least one SNP position
def markReadCoveredSNP(this):
    """Mark that current read covers at least one SNP position"""
    this.current_read_covered_snp = True

# NEW: Finalize coverage tracking for current read (call after processing each read)
def finalizeReadCoverage(this):
    """Call after processing each read to finalize SNP coverage count"""
    if this.current_read_covered_snp:
        this.reads_covering_snp += 1
    this.current_read_covered_snp = False

# NEW: Get count of reads that covered any SNP position
def getReadsCoveringSNP(this):
    """Return count of reads that covered at least one SNP position"""
    return this.reads_covering_snp
```

---

### 3. SNP_Verification_Processes/MisInDelCheck.py

**Location:** `SNP_Verification_Processes/MisInDelCheck.py`

#### Added SNP Coverage Tracking

Added `gene.markReadCoveredSNP()` calls at 4 locations:

**Line 9-10** (in `missenseCheck` function):
```python
if (mapOfInterest.get(mtInfo[1]-1,False)) != False:
    # NEW: Mark that this read covers at least one SNP position
    gene.markReadCoveredSNP()
```

**Line 76-77** (for insertions with multiple positions):
```python
if (mapOfInterest.get(mtInfo[1][0]-1,False)) == False:
    continue
# NEW: Mark that this read covers at least one SNP position
gene.markReadCoveredSNP()
```

**Line 113-114** (for single-position insertions):
```python
if (mapOfInterest.get(mtInfo[1][0]-1,False)) == False:
    continue
# NEW: Mark that this read covers at least one SNP position
gene.markReadCoveredSNP()
```

**Line 139-140** (for deletions):
```python
if (mapOfInterest.get(pos-1,False)) == False:
    continue
# NEW: Mark that this read covers at least one SNP position
gene.markReadCoveredSNP()
```

---

### 4. SNP_Verification_Processes/nTupleCheck.py

**Location:** `SNP_Verification_Processes/nTupleCheck.py`

#### Added SNP Coverage Tracking

Added `gene.markReadCoveredSNP()` calls at 5 locations:

**Line 27-28** (MEG_1731 special case):
```python
if (mapOfInterest.get(pos-1,False)) == False:
    notFound += 1
    continue
# NEW: Mark that this read covers at least one SNP position
gene.markReadCoveredSNP()
```

**Line 83-84** (N-tuple insertions, multiple positions):
```python
if (mapOfInterest.get(mtInfo[1][0]-1,False)) == False:
    continue
# NEW: Mark that this read covers at least one SNP position
gene.markReadCoveredSNP()
```

**Line 119-120** (N-tuple insertions, single position):
```python
if (mapOfInterest.get(mtInfo[1][0]-1,False)) == False:
    continue
# NEW: Mark that this read covers at least one SNP position
gene.markReadCoveredSNP()
```

**Line 144-145** (N-tuple deletions):
```python
if (mapOfInterest.get(mtInfo[1][0]-1,False)) == False:
    resBool = False
    break
# NEW: Mark that this read covers at least one SNP position
gene.markReadCoveredSNP()
```

**Line 176-177** (N-tuple missense):
```python
if (mapOfInterest.get(mtInfo[1]-1,False)) == False:
    resBool = False
    break
# NEW: Mark that this read covers at least one SNP position
gene.markReadCoveredSNP()
```

---

### 5. SNP_Verification_Processes/__init__.py

**Location:** `SNP_Verification_Processes/__init__.py`

#### Added Coverage Finalization (Lines 58-59)

Added call to finalize SNP coverage tracking after processing each read:

```python
def FinalCount(gene, read):
    gene.redefineLastTupleInfo(read)
    additionalInfo = gene.getLastTupleInfo()
    
    # NEW: Finalize SNP coverage tracking for this read
    gene.finalizeReadCoverage()
    
    def nTypeCount():
        # ... rest of function
```

---

## Output Files

### 1. snp_coverage_stats.csv

Per-gene statistics with columns:

| Column | Description |
|--------|-------------|
| sample_name | Name of the sample (for easy concatenation) |
| gene_name | Full gene name with RequiresSNPConfirmation tag |
| gene_type | "Normal" or "Intrinsic" |
| original_amrplusplus_count | Count from input AMR++ matrix |
| total_reads_analyzed | Total reads processed for this gene |
| reads_covering_snp_position | Reads that covered at least one SNP position |
| reads_with_snp_confirmed | Reads with resistance SNP confirmed |
| percentage | reads_with_snp_confirmed / reads_covering_snp_position |
| original_code_would_set | What original code would have set |
| new_count_percentage_based | New percentage-based count |
| snp_confirmed | True/False |

### 2. snp_verification_summary.csv

Run summary with columns:

| Column | Description |
|--------|-------------|
| sample_name | Name of the sample |
| genes_in_snpinfo_database | Total genes in SNPinfo database |
| genes_processed | Genes that were processed |
| genes_with_reads | Genes that had reads > 0 |
| genes_found_in_count_matrix | Genes found in input count matrix |
| genes_not_found_in_matrix | Genes not found in input count matrix |
| genes_with_nonzero_counts | Genes with non-zero counts in input |
| genes_with_snp_confirmed | Genes where SNP was confirmed |
| genes_with_snp_not_confirmed | Genes where SNP was NOT confirmed |
| snpconfirmation_genes_not_in_db | RequiresSNPConfirmation genes not in SNPinfo DB |
| genes_zeroed_not_in_db | Of those, genes with counts that were zeroed |

---

## How the Percentage-Based Calculation Works

### Example

```
Gene X in original AMR++ count matrix: 1000 alignments

SNP_Verification analyzes:
  - 500 total reads processed for this gene
  - 200 reads actually covered the SNP position(s)
  - 150 of those 200 had the resistance SNP confirmed
  - 50 had wild-type (susceptible)
  - 300 reads didn't cover any SNP position

Calculation:
  percentage = 150 / 200 = 75%
  new_count = 0.75 × 1000 = 750
```

### Why This Matters

- **Without coverage tracking**: percentage = 150/500 = 30%, giving new_count = 300 (underestimate)
- **With coverage tracking**: percentage = 150/200 = 75%, giving new_count = 750 (accurate)

Reads that don't cover the SNP region shouldn't be in the denominator because we can't determine their resistance status.

---

## Concatenating Output Files Across Samples

```bash
# Coverage stats
head -1 sample1/snp_coverage_stats.csv > all_samples_coverage.csv
tail -n +2 -q */snp_coverage_stats.csv >> all_samples_coverage.csv

# Summary stats  
head -1 sample1/snp_verification_summary.csv > all_samples_summary.csv
tail -n +2 -q */snp_verification_summary.csv >> all_samples_summary.csv
```

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| v1 | Session 1 | Fixed --detailed_output parsing bug |
| v2 | Session 2 | Added zeroing of unverified genes, intrinsic gene tracking |
| v3 | Session 3 | Added percentage-based calculation |
| v4 | Session 3 | Added SNP position coverage tracking |
| v5 | Session 3 | Added CSV output files |
| v6 | Session 3 | Changed summary to CSV format |
| v7 (Final) | Session 3 | Added sample_name column to all output rows |
