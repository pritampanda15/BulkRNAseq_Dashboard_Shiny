# RNA-seq Analysis Dashboard

This is a Shiny web application for RNA-seq data analysis. It allows users to upload count data and metadata, run **DESeq2** for differential gene expression analysis, and visualize the results using various plots.

## Features

- **Upload Count and Metadata Files**:
  - Supports `.csv`, `.tsv`, and `.txt` file formats.
  - Input count data (gene expression) and metadata (sample conditions).

- **Differential Expression Analysis**:
  - Uses the `DESeq2` package for gene-level differential expression analysis.

- **Interactive Data Visualization**:
  - **Volcano Plot**: Highlights significant genes based on log-fold change and p-values.
  - **Heatmap**: Displays the top expressed genes.
  - **PCA Plot**: Visualizes sample clustering using Principal Component Analysis.
  - **Sample Distance Heatmap**: Shows sample-to-sample distances for quality assessment.
  - **MA Plot**: Log-fold changes vs. mean counts.
  - **Dispersion Plot**: Visualizes dispersion estimates for genes.

## Installation

### Prerequisites
To run the app, ensure you have R and the following R packages installed:

1. Install R and RStudio: [Download R](https://cran.r-project.org/) | [Download RStudio](https://posit.co/downloads/)

2. Install the required R packages:

```r
install.packages(c("shiny", "shinydashboard", "ggplot2", "DT", "pheatmap"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "EnhancedVolcano"))
```

---

## How to Run the Application

1. Clone this repository to your local machine:
   ```bash
   git clone https://github.com/username/repository.git
   cd repository
   ```

2. Launch the Shiny application in R or RStudio:

   ```r
   library(shiny)
   runApp("app.R")
   ```

3. Upload your RNA-seq **counts file** and **metadata file** when prompted.

   - **Counts file**: Rows as gene IDs, columns as sample IDs.
   - **Metadata file**: Contains sample IDs (matching column names in counts file) and condition labels.

---

## Input File Format

### **Counts File Example**:
| Gene ID          | SRR1146097 | SRR1146103 | SRR1146116 |  
|------------------|------------|------------|------------|  
| ENSG00000000003  | 187        | 427        | 697        |  
| ENSG00000000005  | 3          | 17         | 119        |  

### **Metadata File Example**:
| Run          | condition    |  
|--------------|--------------|  
| SRR1146097   | psoriasis    |  
| SRR1146103   | psoriasis    |  
| SRR1146116   | normal       |  

---

## Dashboard Features

- **Data Upload Tab**:
   - Upload count data and metadata files.
   - Trigger analysis with the "Run Analysis" button.

- **Differential Expression Tab**:
   - View the results in an interactive table.

- **Plots Tab**:
   - **Volcano Plot**: Significant genes based on thresholds.
   - **Heatmap**: Top 50 highly expressed genes.
   - **PCA Plot**: Sample clustering.
   - **Sample Distance Heatmap**: Explore relationships between samples.
   - **MA Plot**: Log-fold changes vs. mean counts.
   - **Dispersion Estimates**: Visualize model dispersion.
---

## Contributing

Contributions are welcome! Feel free to:
- Submit pull requests for improvements.
- Open issues to report bugs or suggest new features.

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---


```

### **Instructions to Customize**:
1. Replace:
   - `username/repository` with your actual GitHub username and repository name.
   - `path/to/volcano_plot.png` with the path to your actual images (optional).
2. Add screenshots of the dashboard plots to your repository for better visual appeal.

### **Commit and Push**:
After saving the `README.md` file, commit and push it to GitHub:

```bash
git add README.md
git commit -m "Add README file for Shiny RNA-seq Analysis Dashboard"
git push origin main
```
