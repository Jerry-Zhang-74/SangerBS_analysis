# SangerBS-Analysis: Automated Analysis Tool for Sanger Bisulfite Sequencing

### Project Overview
SangerBS-Analysis is an automated software tool designed to analyze Sanger sequencing data (.ab1 format) derived from Bisulfite Sequencing (BS). It addresses the common pain points in manual BS-Seq analysis, such as low efficiency, high error rates, and difficulties in handling reverse sequencing or alignment shifts. By utilizing a physical sliding alignment algorithm and local pixel-level peak detection, the tool automatically extracts methylation percentages for target CpG sites and generates standardized analysis reports.

### Core Functions
1.  **Automated Sequence Alignment**: Implements a 3-Letter reduction algorithm to precisely align low-complexity bisulfite-converted sequences with reference genomes (either wild-type or simulated-converted).
2.  **Direction Recognition**: Automatically detects the orientation of the sequencing primers (Forward or Reverse) and maps coordinates accordingly without manual intervention.
3.  **Indel Correction**: Addresses polymer slips (Indels) common in BS sequencing via local fine-tuning alignment, ensuring peak extraction occurs at the correct physical coordinates for every CpG site.
4.  **Methylation Quantitation**: Calculates single-base resolution methylation levels by extracting calibrated fluorescence intensities of Cytosine (C) and Thymine (T) peaks using the formula `C / (C + T)`.
5.  **Experimental Quality Control (QC)**: Automatically calculates the conversion rate of non-CpG cytosines to assess the reliability of the bisulfite treatment.

### Input Requirements
To run the analysis, the following two types of files are required:
* **Reference Sequence File** (1 file):
    * Supported formats: `.dna` (SnapGene), `.fasta`, `.fa`, or `.txt`.
    * Content: Providing the original wild-type (unconverted) sequence is recommended, as the program will automatically identify internal CpG sites.
* **Sequencing Data Files** (One or more):
    * Format: `.ab1`.
    * Content: Original trace files provided by sequencing service providers.

### Output Results
Upon completion, the tool generates three types of output:
1.  **Aligned_Matrix.csv (Alignment Matrix)**:
    * A standardized table with CpG sites as rows and sample names as columns. It displays the methylation percentage for each site, ready for import into statistical software like GraphPad Prism.
2.  **Details.csv (Analysis Details)**:
    * Contains specific reference coordinates, peak heights, and raw data for every site, facilitating data traceability and verification.
3.  **Methylation_Plot.html (Interactive Chart)**:
    * A visualized methylation map supporting zoom functions and hover-data. It provides an intuitive comparison of methylation trends across different experimental groups.

### Technical Implementation
* **Language**: Python
* **Core Libraries**: Pandas (Data Processing), Numpy (Numerical Calculation), Plotly (Visualization), CustomTkinter / Streamlit (UI Interaction).
* **Algorithm Logic**: Sliding window alignment based on physical distance combined with Sanger DATA channel signal extraction.

---

### Instructions for Use
1.  Launch the application and upload your reference sequence (.dna or .fasta) in the sidebar.
2.  Batch select and upload the .ab1 sequencing files to be analyzed.
3.  Adjust parameters such as "Minimum Signal Threshold" or "Trim Start" if necessary.
4.  Click "Run Analysis." The system will automatically perform alignment and display the results matrix and methylation trends.
5.  Click the download buttons to save the generated CSV files.
