## Integration analysis of scfa and microbiome

We want to bring together the key parts of each analysis and then look at the two datasets together.  
We should not have to start all over because much of the two separate analyses in ../microbiome and ../scfa_metabolomics are done and very good

SCFA reference source: scfa_metabolomics/scfa-project2-analysis-improved-v2.Rmd
Microbiome reference source: microbiome_analysis/analysis.v4_zymo.Rmd

We do want the integrated analysis to be based on an R markdown file that when rendered will run each part of the analysis like in the microbiome and scfa analyses. It should not cherry pick the completed analyses from both of them.  

We want these sections. For each results section provide prose to introduce the section and what is attempted to achieve. For each figure and table, describe the statistical methods used.

Make all of the tables and graphs publication quality. Use color and themes from the SCFA analysis for the integrated analysis. All code should be acordianed and in the collapsed state initially opening an html report. Suppress warnings and intermediate output in the report. This can particularly be an issue for the microbiome analysis.

Try to integrate the metadata into a shared sample_data component of the phyloseq object so that the same data are being used for the integrated analysis throughout and when filtering is being done we assure consistency and rigor.

Implementation details and assumptions

- SCFA reference source: `scfa_metabolomics/scfa-project2-analysis-improved-v2.Rmd`
- Microbiome reference source: `microbiome_analysis/analysis.v4_zymo.Rmd`
- R version: 4.5.1 (packages already installed as per source analyses; required packages declared in the source Rmd files)
- Output: single `html_document` with HTML-specific features enabled (e.g., code folding/accordions, dynamic tables)
- Visuals: use themes, color palettes, and general styling from the SCFA analysis as the primary reference, applied consistently across integrated analyses
- Code visibility: all code chunks should be collapsed/accordioned by default in the HTML output; warnings and intermediate console output should be suppressed
- Metadata integration:
  - SCFA data sample identifier: `sampleid`
  - Microbiome data sample identifier: `SampleID`
  - A new, integrated master metadata table will be constructed and used to populate a shared `sample_data` component of the phyloseq object, ensuring consistent filtering and covariate usage across SCFA and microbiome analyses
- qPCR / absolute abundance:
  - qPCR data file: `microbiome_analysis/zr24558.16S_250813.zymo/00...AllSamples.Bac16Sv34/Sample_Information/absolute.abundance.csv`
  - qPCR sample identifier column: `customer_label`
  - We will use total genome copies per sample to convert relative abundances to absolute abundances for each taxon in the microbiome data
- Target taxon:
  - Focal genus: `g_Fusicantenibacter sp.` (aggregated at the genus level across all features assigned to this genus)
  - Absolute abundance plots will be generated for this genus by Control vs Case and 0 h vs 48 h within each carbohydrate condition
- Responder vs non-responder groups:
  - Defined separately for butyrate and propionate based on changes from 0 h to 48 h
  - Initial approach (to be implemented and evaluated): use a data-driven threshold (e.g., median or quantile-based split on delta SCFA) to define responder status, then examine overlap in responder classifications between butyrate and propionate
- Development strategy:
  - Implement one major section at a time (Summary Data → Metabolomics → Microbiome → Integrated analyses), reusing and adapting code from the source Rmds
  - Cache heavier modeling steps (e.g., mixed-effects models, ANCOM-BC2, MaAsLin2) so that users can delete cache directories when full recomputation is desired

1. Introduction - take the base from the scfa analysis and add information about the microbiome analysis as it fits into the project. You should use a web search and citations as necessary searching highly reliable and academic sources like pubmed.
2. Project Objectives - take the base from scfa analysis and add information about the microbiome analysis. 
3. Methods  - integrate the methods from the scfa and microbiome analysis. 
4. Summary Data - a concise table with the number of samples for control and case groups in each of the carbohydrate conditions
5. Metabolomics 
- Give a brief section introduction to what and why the analysis
- Calculate and show the analysis from scfa called 'Carbohydrate Type Effects (Combined Groups)'. This will include the Control and Case Groups separately and combined as in the scfa rmd.
- Calculate and show  the analysis from scfa callled 'SCFA Concentrations Over Time'. 
- Calculate and show from scfa analysis 'Advanced Statistical Modeling: Mixed-Effects Models for SCFA Concentrations'. Order the table from lowest to highest p_adj. The table should be a dynamic table as per the original analysis
- Include a brief summary section of the scfa results with key findings
6. Microbiome
- Give a brief section introduction to what and why the analysis
- Calculate and show 'Alpha Diversity' per the source microbiome analysis. 
  - Show 'Observed', 'Shannon', 'Simpson' at 0 and 48 hr per the source. 
  - Calculate and show the 'Combined-Group Modeling'
  - Show the 'Alpha diversity at 48h'. 
  - Calculate and show 'Beta Diversity' per the source microbiome analysis
    - Show the PCoA using Bray distance for 0H and 48H per the microbiome source
    - Show the 'PERMONOVA results' table 
    - Show the 'Beta Diversity at 48 h' from the source
    - Include the 'Baseline Clustering and Trajectories' section from the source analysis but only include graphs from the original. 
      - NEW: Calculate the PC1 distance change for each carbohydrate condition and cluster between 0H and 48H. Show the results graphically and provide statistics.
  - Calculate and show Community Composition graphs and data
    - Show 'Genus  Composition Changes from 0H to 48H by Carbohydrate Type'. Try to repair the RDC panel that didn't work properly in the source.
    - Show 'Validation: Bar plots of genus abundances by timepoint and carbohydrate'
  - Calculate and include 'Differential Abundance and Predictive Modeling' from the source
    - Show a summary table with columns taxon | Genus Species | DA Method | p value | q value (or similar)
    - Show the  'Forest Plot of ANCOM-BC2 Results'
    - Show the 'Side-by-Side Comparison: RDC vs No Carb and SDC vs No Carb' for ANCOM-BC2 analysis
    - Show 'Forest Plot of ANCOM-BC2 48h Results'
    - Show 'Forest Plot: RDC vs SDC at 48h' and 'Forest Plot: Group × Carbohydrate Interaction at 48h' from the ANCOM-BC2 analysis
    - Show 'LEfSe Visualization' but call it 'LEfSe Differential Abundance' 
    - Show the MaAsLin2 table with statistics but also add column for genus species next to taxon
    - Show 'MaAsLin2 Visualization' but this graph needs to be taller to space out the y axis text properly
    - Show the 'MaAsLin2 Volcano Plot' and use ggrepel to label the points features in 'MaAsLin2: Top Effects by Variable'
    - Show 'MaAsLin2: Top Effects by Variable'
    - Show 'MaAsLin2 Volcano Plot: 48H RDC vs  SDC'
    - Show 'MaAsLin2: Timepoint Effects (RDC vs SDC)'
  - NEW: Absolute Abundances
    - microbiome_analysis/zr24558.16S_250813.zymo/00...AllSamples.Bac16Sv34/Sample_Information/absolute.abundance.csv has the qPCR determine genome copies per sample information.  We want to use relative abundance data to convert to absolute abundances for each taxon in the original data. We then want to show g_Fusicantenibacter sp. absolute abundance plots in Control vs. Case vs. 0H and 48H in each carbohydrate condition.  
7. Integrated analysis
  - Taxon - SCFA relationships
    - DA analysis using ANCOM-BC2 and MaAsLin2 for taxon ~ butyrate + carbohydrate using relative and absolute data
    -  DA analysis using ANCOM-BC2 and MaAsLin2 for taxon ~ propionate + carbohydrate using relative and absolute data
    - Graph g_Fusicantenibacter sp correlations to butyrate and propionate concentrations at 0H and 48H using relateive and absolute abudances
  - Group analysis - looking at the source SCFA analysis 'SCFA Concentrations by Group and Carbohydrate Type' and focusing on butyrate and propionate it looks like there are non-responder and responder groups. there is some separation of samples where each SCFA is higher and almost negligable. We want to separate those groups and separately evaluate them for differences in their composition. We can focus on them at 48H.
8. Summary
- key findings 
