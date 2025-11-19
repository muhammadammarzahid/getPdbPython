# TTD Druggable Target Structure Pipeline

A Python pipeline designed to automate the retrieval and processing of high-quality structural data for druggable targets. This tool fetches target lists from the Therapeutic Target Database (TTD), maps them to PDB entries via UniProt and SIFTS, filters for the best available resolutions, and converts the files into analysis-ready formats (`.pdb` and `.mae`).

## 📋 Workflow

The script executes a 4-phase automated pipeline:

1.  **Fetch & Clean:** Downloads the latest druggable target list from TTD and extracts clean UniProt Entry Names.
2.  **Map & Filter:**
      * Converts Entry Names to Accession Codes using the UniProt API.
      * Maps Accession Codes to PDB IDs using SIFTS data.
      * Queries the RCSB PDB to filter for **X-Ray Diffraction** structures with **resolution \< 3.0Å** (configurable).
      * Selects the single best-resolution structure per protein.
3.  **Download:** Retrieves raw `.cif` files from the RCSB PDB.
4.  **Convert:** Utilizes the Schrödinger Suite (`structconvert`) to generate `.pdb` and `.mae` (Maestro) files for every downloaded structure.

## 🛠 Prerequisites

### 1\. Python Dependencies

This project relies on the following Python libraries:

  * `pandas` (Data manipulation)
  * `requests` (API handling)
  * `tqdm` (Progress bars)

### 2\. External Software (Schrödinger)

**Crucial:** This script requires the `structconvert` utility from the **Schrödinger Software Suite**.

  * You must have a valid installation of Schrödinger.
  * You must set the `SCHRODINGER` environment variable to your installation path.

## 📦 Installation

1.  **Clone the repository:**

    ```bash
    git clone https://github.com/yourusername/getPdbPython.git
    cd getPdbPython
    ```

2.  **Install Python requirements:**

    ```bash
    pip install -r requirements.txt
    ```

3.  **Configure Environment Variables (Linux/macOS):**
    Point the variable to your Schrödinger installation folder (e.g., version 2024-1):

    ```bash
    export SCHRODINGER=/opt/schrodinger2024-1
    ```

## 🚀 Usage

Run the main script to start the pipeline:

```bash
python getPdbPython.py
```

### Configuration

You can modify the `--- CONFIGURATION ---` section at the top of the script to change:

  * **`MAX_RESOLUTION`**: The resolution cutoff (default: 3.0).
  * **`EXPERIMENTAL_METHODS`**: The allowed experiment types (default: X-RAY DIFFRACTION).
  * **`RAW_DOWNLOAD_DIR`**: The folder name for raw downloads.

## 📂 Output Structure

The script will create the following directories and files:

```text
├── PDB_repository_TTD/      # Raw downloaded .cif files
├── PDB_converted_PDB/       # Processed .pdb files (sorted by UniProt ID)
├── PDB_converted_MAE/       # Processed .mae files (sorted by UniProt ID)
├── download_log_ttd.csv     # CSV log of the PDB download results
└── conversion_log.csv       # CSV log of the file conversion results
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
