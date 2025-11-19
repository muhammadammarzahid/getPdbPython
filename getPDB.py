import os
import time
import pandas as pd
import requests
from tqdm import tqdm
import json
import re
import subprocess

# --- CONFIGURATION ---
TTD_URL = "https://ttd.idrblab.cn/files/download/P2-02-TTD_uniprot_successful.txt"
UNIPROT_API_URL = "https://rest.uniprot.org"
SIFTS_URL = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz"
RCSB_GQL_URL = "https://data.rcsb.org/graphql"
MAX_RESOLUTION = 3.0
EXPERIMENTAL_METHODS = ["X-RAY DIFFRACTION"]
RAW_DOWNLOAD_DIR = "PDB_repository_TTD"
DOWNLOAD_LOG_FILE = "download_log_ttd.csv"
# --- CHANGE: Added separate directories for each output format ---
CONVERTED_PDB_DIR = "PDB_converted_PDB"
MAE_OUTPUT_DIR = "PDB_converted_MAE"
# --- END OF CONFIGURATION ---

# --- Using a local file for testing purposes ---
# def get_ttd_targets():
#     """
#     Reads UniProt Entry Names from a local test file.
#     """
#     print("--- Phase 1: Fetching and Cleaning Druggable Targets from a local file for testing ---")
#     try:
#         with open('small_ttd_input.txt', 'r') as f:
#             content = f.read()
#     except FileNotFoundError:
#         print("Error: small_ttd_input.txt not found.")
#         return None

#     uniprot_entry_names = set()
#     raw_id_lines = re.findall(r"^UNIPROID\s+(.*)", content, re.MULTILINE)

#     for line in raw_id_lines:
#         potential_ids = line.split(';')
#         for pid in potential_ids:
#             cleaned_id = pid.strip()
#             if cleaned_id and cleaned_id.lower() != "uniprot id":
#                 uniprot_entry_names.add(cleaned_id)

#     if not uniprot_entry_names:
#         print("Error: Could not parse any valid UniProt IDs from the local file.")
#         return None

#     print(f"Found {len(uniprot_entry_names)} unique, clean UniProt Entry Names from the local file.")
#     return uniprot_entry_names

# --- Updated to download from TTD URL and handle multiple IDs per line ---
def get_ttd_targets():
    """
    Downloads the TTD target list, cleans the data, and extracts UniProt Entry Names.
    """
    print("--- Phase 1: Fetching and Cleaning Druggable Targets from TTD ---")
    print(f"Downloading target list from {TTD_URL}...")
    try:
        response = requests.get(TTD_URL)
        response.raise_for_status()
        content = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download data from TTD. {e}")
        return None

    uniprot_entry_names = set()
    raw_id_lines = re.findall(r"^UNIPROID\s+(.*)", content, re.MULTILINE)

    for line in raw_id_lines:
        potential_ids = line.split(';')
        for pid in potential_ids:
            cleaned_id = pid.strip()
            if cleaned_id and cleaned_id.lower() != "uniprot id":
                uniprot_entry_names.add(cleaned_id)

    if not uniprot_entry_names:
        print("Error: Could not parse any valid UniProt IDs from the TTD file.")
        return None

    print(f"Found {len(uniprot_entry_names)} unique, clean UniProt Entry Names from TTD.")
    return uniprot_entry_names

    """
    Downloads the TTD target list, cleans the data, and extracts UniProt Entry Names.
    This version is robust against multiple IDs per line, invisible characters, and header artifacts.
    """
    print("--- Phase 1: Fetching and Cleaning Druggable Targets from a local file for testing ---")
    try:
        with open('small_ttd_input.txt', 'r') as f:
            content = f.read()
    except FileNotFoundError:
        print("Error: small_ttd_input.txt not found.")
        return None

    uniprot_entry_names = set()
    raw_id_lines = re.findall(r"^UNIPROID\s+(.*)", content, re.MULTILINE)

    for line in raw_id_lines:
        potential_ids = line.split(';')
        for pid in potential_ids:
            cleaned_id = pid.strip()
            if cleaned_id and cleaned_id.lower() != "uniprot id":
                uniprot_entry_names.add(cleaned_id)

    if not uniprot_entry_names:
        print("Error: Could not parse any valid UniProt IDs from the local file.")
        return None

    print(f"Found {len(uniprot_entry_names)} unique, clean UniProt Entry Names from the local file.")
    return uniprot_entry_names


def convert_uniprot_ids(entry_names):
    """
    Converts UniProt Entry Names to Accession Codes using the UniProt Search API.
    """
    print("\nConverting UniProt Entry Names to Accession Codes via UniProt Search API...")
    accession_codes = set()
    entry_list = list(entry_names)
    batch_size = 100

    for i in tqdm(range(0, len(entry_list), batch_size), desc="Converting IDs"):
        batch = [name for name in entry_list[i:i + batch_size] if name]
        query_part = " OR ".join([f"(id:{name})" for name in batch])
        params = {"query": query_part, "fields": "accession", "format": "json", "size": len(batch)}

        try:
            response = requests.get(f"{UNIPROT_API_URL}/uniprotkb/search", params=params)
            response.raise_for_status()
            results = response.json().get("results", [])
            for item in results:
                if "primaryAccession" in item:
                    accession_codes.add(item["primaryAccession"])
        except requests.exceptions.RequestException as e:
            print(f"\nWarning: API request failed for a batch: {e}")
            if e.response is not None:
                print(f"Server response: {e.response.text}")
            continue

    print(f"Successfully converted to {len(accession_codes)} unique UniProt Accession Codes.")
    return accession_codes

def get_sifts_mapping(target_uniprot_ids):
    """Downloads SIFTS data and maps the target UniProt Accession Codes to PDB IDs."""
    print("\n--- Phase 2: Mapping Targets to PDB Structures ---")
    print(f"Downloading SIFTS mapping file from {SIFTS_URL}...")
    try:
        sifts_df = pd.read_csv(
            SIFTS_URL, skiprows=1,
            names=["PDB", "CHAIN", "SP_PRIMARY", "RES_BEG", "RES_END", "PDB_BEG", "PDB_END", "SP_BEG", "SP_END"],
            usecols=["PDB", "SP_PRIMARY"], compression='gzip'
        )
    except Exception as e:
        print(f"Error: Could not download or process SIFTS mapping file. {e}")
        return None

    print("Filtering SIFTS data for our target proteins...")
    sifts_df.dropna(inplace=True)
    target_sifts_df = sifts_df[sifts_df['SP_PRIMARY'].isin(target_uniprot_ids)].drop_duplicates()
    uniprot_to_pdb_map = target_sifts_df.groupby('SP_PRIMARY')['PDB'].apply(list).to_dict()
    total_pdb_ids = sum(len(pdbs) for pdbs in uniprot_to_pdb_map.values())
    print(f"Mapped {len(uniprot_to_pdb_map)} UniProt IDs to a total of {total_pdb_ids} PDB entries.")
    return uniprot_to_pdb_map

def filter_pdbs_by_quality(pdb_ids):
    """
    Filters a list of PDB IDs based on resolution and experimental method.
    """
    print("\nQuerying RCSB PDB for structural metadata (this may take a while)...")
    high_quality_pdbs = []
    batch_size = 500
    pdb_id_batches = [pdb_ids[i:i + batch_size] for i in range(0, len(pdb_ids), batch_size)]

    for batch in tqdm(pdb_id_batches, desc="Querying PDB Metadata"):
        query = """
        query($pdb_ids: [String!]!) {
          entries(entry_ids: $pdb_ids) {
            rcsb_id
            exptl { method }
            rcsb_entry_info { resolution_combined }
          }
        }
        """
        variables = {"pdb_ids": batch}
        try:
            response = requests.post(RCSB_GQL_URL, json={"query": query, "variables": variables})
            response.raise_for_status()
            results = response.json()
            if "errors" in results and results["errors"]:
                print(f"\nWarning: GraphQL API returned an error for a batch: {results['errors'][0]['message']}")
                continue
            data = results.get("data")
            if data and "entries" in data:
                for entry in data["entries"]:
                    if entry is None: continue
                    method = (entry.get("exptl") or [{}])[0].get("method")
                    resolution = (entry.get("rcsb_entry_info", {}).get("resolution_combined") or [None])[0]
                    if method in EXPERIMENTAL_METHODS and resolution is not None and resolution <= MAX_RESOLUTION:
                        high_quality_pdbs.append({"pdb_id": entry["rcsb_id"], "method": method, "resolution": resolution})
        except (requests.exceptions.RequestException, json.JSONDecodeError) as e:
            print(f"\nWarning: API call failed for a batch. {e}")
            time.sleep(2)

    print(f"Found {len(high_quality_pdbs)} PDB structures matching the quality criteria.")
    return high_quality_pdbs

def select_best_structure(pdb_info_list):
    """Selects the best PDB structure from a list based on the highest resolution."""
    if not pdb_info_list: return None
    return min(pdb_info_list, key=lambda x: x['resolution'])

def group_and_select_best_structures(high_quality_pdb_list, uniprot_map):
    """Groups high-quality PDBs by UniProt ID and selects the one with the best resolution."""
    print("\nSelecting the best resolution structure for each UniProt ID...")
    pdb_to_uniprot_map = {pdb.upper(): uniprot for uniprot, pdbs in uniprot_map.items() for pdb in pdbs}
    uniprot_to_pdbs_info = {}
    for pdb_info in high_quality_pdb_list:
        uniprot_id = pdb_to_uniprot_map.get(pdb_info["pdb_id"].upper())
        if uniprot_id:
            if uniprot_id not in uniprot_to_pdbs_info:
                uniprot_to_pdbs_info[uniprot_id] = []
            uniprot_to_pdbs_info[uniprot_id].append(pdb_info)
    best_structures = []
    for uniprot_id, pdb_info_list in uniprot_to_pdbs_info.items():
        best_pdb_info = select_best_structure(pdb_info_list)
        if best_pdb_info:
            best_pdb_info['uniprot_id'] = uniprot_id
            best_structures.append(best_pdb_info)
    print(f"Selected {len(best_structures)} unique PDB structures to download (one per protein).")
    return best_structures

def download_pdb_files(best_pdb_list, output_dir):
    """Downloads PDBx/mmCIF files for the single best structure per protein."""
    print(f"\n--- Phase 3: Downloading Selected PDBs ---")
    print(f"Downloading PDBx/mmCIF files to '{output_dir}'...")
    os.makedirs(output_dir, exist_ok=True)
    log_data = []
    for pdb_info in tqdm(best_pdb_list, desc="Downloading PDBs"):
        pdb_id = pdb_info["pdb_id"]
        uniprot_id = pdb_info["uniprot_id"]
        download_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        new_filename = f"{uniprot_id}_{pdb_id}.cif"
        filepath = os.path.join(output_dir, new_filename)
        try:
            response = requests.get(download_url, timeout=30)
            if response.status_code == 200:
                with open(filepath, 'w', encoding='utf-8') as f: f.write(response.text)
                log_data.append({"UniProt_ID": uniprot_id, "PDB_ID": pdb_id, "Resolution": pdb_info["resolution"], "Method": pdb_info["method"], "Status": "Downloaded", "Filename": new_filename})
            else:
                log_data.append({"UniProt_ID": uniprot_id, "PDB_ID": pdb_id, "Resolution": pdb_info["resolution"], "Method": pdb_info["method"], "Status": f"Failed (HTTP {response.status_code})", "Filename": ""})
        except requests.exceptions.RequestException:
            log_data.append({"UniProt_ID": uniprot_id, "PDB_ID": pdb_id, "Resolution": pdb_info["resolution"], "Method": pdb_info["method"], "Status": "Failed (Connection Error)", "Filename": ""})
    if log_data:
        log_df = pd.DataFrame(log_data)
        log_df.to_csv(DOWNLOAD_LOG_FILE, index=False)
    print(f"\nDownload complete. Log saved to '{DOWNLOAD_LOG_FILE}'.")

# --- MODIFIED: Renamed and updated to handle both PDB and MAE conversion ---
def convert_structures(best_pdb_list, raw_pdb_dir, pdb_output_dir, mae_output_dir):
    """
    Phase 4: Converts downloaded .cif files to .pdb and .mae formats using 'structconvert'.
    """
    print(f"\n--- Phase 4: Converting Structures to PDB and MAE Formats ---")
    os.makedirs(pdb_output_dir, exist_ok=True)
    os.makedirs(mae_output_dir, exist_ok=True)
    conversion_log = []
    
    schrodinger_path = os.environ.get("SCHRODINGER")
    if not schrodinger_path:
        print("Error: SCHRODINGER environment variable not set. Cannot find 'structconvert' utility.")
        return
        
    structconvert_path = os.path.join(schrodinger_path, "utilities", "structconvert")
    if not os.path.exists(structconvert_path):
        print(f"Error: Schrödinger utility 'structconvert' not found at expected path.")
        return

    for pdb_info in tqdm(best_pdb_list, desc="Converting Structures"):
        uniprot_id = pdb_info["uniprot_id"]
        pdb_id = pdb_info["pdb_id"]
        
        raw_cif_filename = f"{uniprot_id}_{pdb_id}.cif"
        raw_cif_path = os.path.join(raw_pdb_dir, raw_cif_filename)

        if not os.path.exists(raw_cif_path):
            print(f"\nWarning: Raw input file not found for {uniprot_id} - {pdb_id}. Skipping conversion.")
            continue

        # Create UniProt-specific subdirectories for better organization
        uniprot_pdb_dir = os.path.join(pdb_output_dir, uniprot_id)
        os.makedirs(uniprot_pdb_dir, exist_ok=True)
        uniprot_mae_dir = os.path.join(mae_output_dir, uniprot_id)
        os.makedirs(uniprot_mae_dir, exist_ok=True)

        # Define output paths for both file types
        converted_pdb_path = os.path.join(uniprot_pdb_dir, f"{uniprot_id}_{pdb_id}.pdb")
        converted_mae_path = os.path.join(uniprot_mae_dir, f"{uniprot_id}_{pdb_id}.mae")

        log_entry = {
            "UniProt_ID": uniprot_id, "PDB_ID": pdb_id, 
            "PDB_Status": "Failed", "PDB_Output_File": "",
            "MAE_Status": "Failed", "MAE_Output_File": "",
            "Error": ""
        }

        # --- Step 1: Convert to PDB ---
        try:
            pdb_cmd = [structconvert_path, raw_cif_path, converted_pdb_path]
            subprocess.run(pdb_cmd, capture_output=True, text=True, check=True)
            log_entry["PDB_Status"] = "Converted"
            log_entry["PDB_Output_File"] = converted_pdb_path
        except subprocess.CalledProcessError as e:
            log_entry["PDB_Status"] = "Failed (Conversion)"
            log_entry["Error"] += f"PDB conversion failed: {e.stderr.strip()}; "
        
        # --- Step 2: Convert to MAE ---
        try:
            mae_cmd = [structconvert_path, raw_cif_path, converted_mae_path]
            subprocess.run(mae_cmd, capture_output=True, text=True, check=True)
            log_entry["MAE_Status"] = "Converted"
            log_entry["MAE_Output_File"] = converted_mae_path
        except subprocess.CalledProcessError as e:
            log_entry["MAE_Status"] = "Failed (Conversion)"
            log_entry["Error"] += f"MAE conversion failed: {e.stderr.strip()}"

        conversion_log.append(log_entry)

    if conversion_log:
        conv_log_df = pd.DataFrame(conversion_log)
        conv_log_df.to_csv("conversion_log.csv", index=False)
        print(f"\nConversion complete. Log saved to 'conversion_log.csv'.")

# --- MODIFIED main function ---
def main():
    """Main function to orchestrate the download and conversion workflow."""
    # Phase 1 & 2: Get targets, map to PDBs, and filter by quality
    target_entry_names = get_ttd_targets()
    if not target_entry_names: return
    
    accession_codes = convert_uniprot_ids(target_entry_names)
    if not accession_codes: return
    
    uniprot_to_pdb_map = get_sifts_mapping(accession_codes)
    if not uniprot_to_pdb_map: return
    
    all_pdb_ids = list(set(pdb for pdbs in uniprot_to_pdb_map.values() for pdb in pdbs))
    high_quality_pdb_list = filter_pdbs_by_quality(all_pdb_ids)
    if not high_quality_pdb_list: return
    
    best_structures_to_download = group_and_select_best_structures(high_quality_pdb_list, uniprot_to_pdb_map)
    if not best_structures_to_download: return
    
    # Phase 3: Download the selected structures as .cif files
    download_pdb_files(best_structures_to_download, RAW_DOWNLOAD_DIR)
    
    # Phase 4: Convert the downloaded .cif files to .pdb and .mae formats
    convert_structures(
        best_structures_to_download, 
        RAW_DOWNLOAD_DIR, 
        CONVERTED_PDB_DIR, 
        MAE_OUTPUT_DIR
    )

    print("\n\nWorkflow finished.")

if __name__ == "__main__":
    main()