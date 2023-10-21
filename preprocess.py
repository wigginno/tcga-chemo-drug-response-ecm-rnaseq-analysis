import json
from collections import defaultdict
from typing import Set, List, Union
from glob import glob
import shutil
from pathlib import Path
import random
import string
import os.path as osp
from os import makedirs
import pickle
import subprocess

import requests
import pandas as pd

# Hard-coded paths to directories
BASE_DIR = Path("/media/bean/data/gdc")
METADATA = BASE_DIR / "meta"
CLINICAL = BASE_DIR / "clinical"
EXPRESSION = BASE_DIR / "rnaseq"
GDC_RNASEQ_MANIFEST = BASE_DIR / "meta"
GDC_TOOL_PATH = "/media/data/software/gdc-client/gdc-client"

def get_genes_from_msig_set(gene_set_name, species="human"):
    """Fetches the genes associated with a given gene set name from MSigDB.

    Args:
        gene_set_name (str): The name of the gene set (e.g, NABA_MATRISOME).
        species (str, optional): The species for which to fetch the gene set. Defaults to "human".

    Returns:
        list[str]: A list of gene symbols associated with the specified gene set name.

    Raises:
        HTTPError: If the GET request to the MSigDB results in an error.
    """
    url = f"https://www.gsea-msigdb.org/gsea/msigdb/{species}/download_geneset.jsp?geneSetName={gene_set_name}&fileType=json"
    response = requests.get(url, timeout=10)
    response.raise_for_status()
    return response.json()[gene_set_name]["geneSymbols"]

def get_matrisome(core_only=False):
    """Fetches the genes associated with the NABA_MATRISOME gene set from MSigDB."""
    gene_set_name = "NABA_CORE_MATRISOME" if core_only else "NABA_MATRISOME"
    return get_genes_from_msig_set(gene_set_name)

def delete_file(file_path: Union[str, Path]):
    if isinstance(file_path, str):
        file_path = Path(file_path)
    if file_path.exists():
        file_path.unlink()

def delete_files(file_paths: List[Union[str, Path]]):
    for file_path in file_paths:
        delete_file(file_path)

def get_clinical_files() -> set[str]:
    return set(glob(str(CLINICAL / "*.json")))

def make_gdc_rnaseq_manifest(filenames: Set[str], output_path: str):
    """Creates a manifest file for downloading RNA-seq files from the Genomic Data Commons (GDC)."""
    manifest = pd.read_csv(GDC_RNASEQ_MANIFEST / "rnaseq_manifest.txt", sep="\t")
    manifest = manifest[manifest["filename"].isin(filenames)]
    print(f"Found {len(manifest)} files in manifest.")
    manifest.to_csv(output_path, sep="\t", index=False)

def download_rnaseq(rnaseq_files: List, out_dir: str):
    """Downloads RNA-seq files by filename from GDC."""
    if isinstance(out_dir, str):
        out_dir = Path(out_dir)
    makedirs(out_dir, exist_ok=True)
    rnaseq_files = set(rnaseq_files)
    make_gdc_rnaseq_manifest(rnaseq_files, out_dir / "manifest.txt")
    subprocess.run([GDC_TOOL_PATH, "download", "-m", out_dir / "manifest.txt", "-d", out_dir],
                   check=True)

def load_json_data(fp: str) -> dict:
    with open(fp, 'r') as file:
        return json.load(file)

def save_json_data(data: dict, fp: str):
    with open(fp, "w") as f:
        json.dump(data, f, indent=4)

def remove_extra_lines(fp):
    """Remove extra lines from raw rnaseq files from GDC."""
    has_extra = False
    with open(fp, "r") as f:
        lines = f.readlines()[1:]
    has_extra = lines[1].startswith("N_unmapped")
    if has_extra:
        del lines[1:5]
        with open(fp, 'w') as f:
            f.writelines(lines)

def save_patient_barcodes_to_file(fp, barcodes):
    """Save patient barcodes to a text file."""
    if not isinstance(barcodes, pd.Series):
        barcodes = pd.Series(list(barcodes))
    barcodes.to_csv(fp, index=False, header=False)

def read_patient_barcodes_from_file(fp):
    """Read patient barcodes from a text file."""
    return pd.read_csv(fp, header=None)[0].tolist()

def get_rnaseq_files_from_barcodes(barcodes):
    barcode_to_rnaseq_files = load_json_data(osp.join(METADATA, "barcode_to_rnaseq_files.json"))
    rnaseq_files = set()
    for barcode in barcodes:
        rnaseq_files.update(barcode_to_rnaseq_files[barcode])
    return rnaseq_files

def get_barcodes_for_rnaseq_files(rnaseq_filepaths):
    rnaseq_file_to_barcode = load_json_data(osp.join(METADATA, "rnaseq_file_to_barcode.json"))
    rnaseq_filepath_to_barcode = {}
    for rnaseq_filepath in rnaseq_filepaths:
        fname = osp.basename(rnaseq_filepath)
        rnaseq_filepath_to_barcode[rnaseq_filepath] = rnaseq_file_to_barcode[fname]
    return rnaseq_filepath_to_barcode

def get_unstranded_counts(rnaseq_filepath, barcode=None):
    """Extract unstranded counts column from rna-seq file after removing extra lines."""
    df = pd.read_csv(rnaseq_filepath, sep="\t", index_col=0, usecols=["gene_id", "unstranded"])
    if barcode is not None:
        # rename unstranded column to barcode
        df.rename(columns={"unstranded": barcode}, inplace=True)
    return df

def process_rnaseq_files(expression_dir, out_dir):
    """Process raw rnaseq files from GDC into an expression matrix, then save that and gene list."""
    filepaths = glob(osp.join(expression_dir, "*/*.tsv"))
    fp_to_barcode = get_barcodes_for_rnaseq_files(filepaths)

    if isinstance(out_dir, str):
        out_dir = Path(out_dir)

    makedirs(out_dir, exist_ok=True)

    all_tmp_files = []
    barcode_to_tmp_files = defaultdict(list)
    for fp in filepaths:
        barcode = fp_to_barcode[fp]
        out_fp = out_dir / ''.join(random.choices(string.ascii_letters + string.digits, k=16))
        shutil.copyfile(fp, out_fp)
        remove_extra_lines(out_fp)
        barcode_to_tmp_files[barcode].append(out_fp)
        all_tmp_files.append(out_fp)

    barcode_to_tmp_file = {}

    count_sum = 0
    count_sums = {}

    for barcode, tmp_files in barcode_to_tmp_files.items():
        if len(tmp_files) > 1:
            filepath = max(tmp_files, key=(lambda fp: get_unstranded_counts(fp)["unstranded"].sum()))
        else:
            filepath = tmp_files[0]
        current_count_sum = get_unstranded_counts(filepath)["unstranded"].sum()
        count_sum += current_count_sum
        count_sums[barcode] = current_count_sum
        barcode_to_tmp_file[barcode] = filepath

    mean_count_sum = count_sum / len(barcode_to_tmp_file)
    threshold = mean_count_sum * 0.25
    for barcode in list(barcode_to_tmp_file.keys()):
        if count_sums[barcode] < threshold:
            del barcode_to_tmp_file[barcode]

    # load one of the files to save the gene ids and names
    fp = next(iter(barcode_to_tmp_file.values()))
    index_df = pd.read_csv(fp, sep="\t", index_col=0, usecols=["gene_id", "gene_name"])
    index_df.to_csv(osp.join(out_dir, "genes.tsv"))

    dfs = [get_unstranded_counts(fp, barcode) for barcode, fp in barcode_to_tmp_file.items()]
    expression_df = pd.concat(dfs, axis=1)
    expression_df.to_csv(osp.join(out_dir, "expression.tsv"), sep="\t")

    delete_files(all_tmp_files)
