# Contents
Original scripts to elucidate drug modes of action by identifying transcription factors that organize the expression of chemically induced genes using large-scale ChIP-seq data.

# Author
Zhaonan Zou (Department of Drug Discovery Medicine, Kyoto University Graduate School of Medicine, JAPAN)

# 1. Performance of ChIPEA 

1) ```CTD_Chemical_Gene_list.sh``` ```CTD_Genelist_process_for_ChIPEA.sh```
    
    Extraction of the gene symbols of up- and down-regulated genes induced by a query chemical (from CTD).
2) ```CTD_ChIPEA.sh```
    
    ChIPEA with up- and down-regulated genes extracted in 1) as input.
    

    In particular, the overlaps between the transcription start site ± 5 kb regions of chemically induced genes and peak-call data of all TF-related experiments archived in ChIP-Atlas were counted, using the “intersect” command of BEDTools2 (version, v2.23.0) [73].
    
    Enrichment scores (−log10(p-values)) were calculated using the two-tailed Fisher’s exact probability test, with the null hypothesis that the two data sets (up- and down-regulated genes) overlap with the ChIP-seq peak-call data in the same proportion.
    
    Fold enrichment values were returned at the same time.

    These procedure were performed across hundreds of chemicals using "insilicoChIP" command, the core of which is CLI-based ChIPEA.
    
    * Note that GUI-based and API-based ChIPEA are also provided on the website of ChIP-Atlas (https://chip-atlas.org/enrichment_analysis).
    
      * Submission form of ChIPEA on website. 
      ![图片 1](https://user-images.githubusercontent.com/74224230/135547550-c3b0eb2c-1685-4af6-a8a0-a940e36e60bb.png)
    
        When used for identify pivotal TFs involved in drug MoAs, first, genome assembly should be set as “hg19/hg38” (hg19 was used in this paper). Then “TFs and others” need to be selected in panal “1. Antigen Class”. “2. Cell type Class” and “3. Threshold for Significance” may be changed by the user according to demand. “4. Enter dataset A” dialog box is to be filled in by the batch of up-regulated genes, and “5. Enter dataset B” is for down-regulated genes. After specify the “Distance range from TSS” in “6. Analysis description” panel, the user can press the “Submit” button to submit the parameters to the server, and ChIPEA should be initialized immediately. 
      
      * Interpretation of the results.
      ![image](https://user-images.githubusercontent.com/74224230/135548459-d5e44b6c-63d6-40d0-ad42-026b45ff513d.png)      
        The overlaps between the genomic loci (originated from panels 4 and 5 of the submission form) and reference peak-call data (specified on upper panels 1 to 3 of the submission form) are counted with bedtools intersect command (BedTools2; ver 2.23.0). The results are returned in html as well as tsv format. P-values are calculated with two-tailed Fisher’s exact probability test. The null hypothesis is that the intersection of reference peaks with submitted data in panel 4 occurs in the same proportion to those with data in panel 5 of the submission form. Q-values are calculated with the Benjamini & Hochberg method. Fold enrichment is calculated by (column 6) / (column 7) of the same row. If the ratio > 1, the rightmost column is ‘TRUE’, meaning that the proteins from column 3 binds to the data of panel 4 in a greater proportion than to those of panel 5 specified in the submission form.
    
    Please refer to Oki, S. et al, EMBO Rep, 2018 or documentation of ChIP-Atlas (https://github.com/inutano/chip-atlas/wiki) for detailed information about primary processing of ChIP-seq data and data annotation.

# 2. Identification of chemical–TF–disease triple associations.

   * ```DisGeNET_link.sh```
    Relating the chemical–TF matrix identified using ChIPEA to TF–disease associations.
    
        Data about gene/protein–disease associations derived from DisGeNET (formula 1; where Pm and Dn represent a protein and a disease) were related to the chemical–TF associations presented by ChIPEA (formula 2; where Ci, Tj and Eij represent a chemical, a TF and an enrichment score) when Tj is also included in DisGeNET as Pm (formula 3). The enrichment scores calculated by ChIPEA were also used to evaluate the probability of each chemical–disease prediction (formula 4). 
        
        ![image](https://user-images.githubusercontent.com/74224230/135550163-1f2f0997-c4a7-4614-94ce-de34206c72ff.png)
        
        If a chemical–disease pair was predicted via multiple TFs, the highest enrichment score was adopted.

# 3. DEG-connected method

   * ```Gene_pivotal.sh```
    A baseline method for predicting chemical–disease associations.
    
        Chemically induced and disease-specific expression changes in 28,268 genes were classified as up-regulated, down-regulated, or non-DEG.
        
        RefSeq genes were obtained from the UCSC genome annotation database of the human genome (download link, http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz; genome, hg19; downloaded on July 7th, 2020).
        
        All chemical and disease profiles were further arranged into two-by-two cross tabulations for each of the chemical–disease pairs in two ways (“DEGs vs. non-DEGs” and “up- vs. down-regulated genes”), as illustrated below.

        ![image](https://user-images.githubusercontent.com/74224230/135549275-09313173-fd47-4216-8b75-8ce9ca69a576.png)

        The chemical–disease associations were evaluated based on p-values calculated using the two-tailed Fisher’s exact probability test with the null hypothesis that the comparative gene expression patterns in response to a given C (chemical) and D (disease) (n1, n2, n3, and n4) were uniformly distributed. Chemical–disease pairs with smaller p-values were considered to be more firmly associated.

# 4. Validation of predicted results.

   * ```chem-tf_AUC.sh``` ```chem-tf-disease_AUC```
    
    1) Assessment of the correctness of
         a) predicted chemical–TF associations by ChIPEA using known chemical–protein interactions data obtained from KEGG DRUG as standard data.
         b) predicted chemical–disease associations by ChIPEA-based approach as well as the baseline (DEG-connected) method using chemical–disease associations obtained from CTD were used as standard data.
       
    2) Performance of ROC and PR curve
         a) Receiver operating characteristic (ROC) curve: a plot of true-positive rates as a function of false-positive rates;
         b) Precision-recall (PR) curve: a plot of precision (positive predictive value) as a function of recall (sensitivity).
   
         The evaluation was summarized using the area under the ROC curve (AUROC) score, where 1 is perfect classification and 0.5 is random classification, and the area under the PR curve (AUPR) score, where 1 is perfect inference and the ratio of positive examples in the standard data is random inference.
      
    
    3) Calculation of global AUROC and AUPR
         After 1), all predicted chemical–TF (or disease) associations were first arranged into a single m × n matrix of enrichment scores with correctness information, where m is the total number of chemicals and n is the number of TFs (or diseases). We then stored the maximum value of enrichment scores within each column (TF or disease) into a vector with n elements. We generated ROC and PR curves and summarized the results into global AUROC and AUPR scores as described in the Results section.


    
chem-tf_AUC.sh
chem-tf-disease_AUC
Summary_4_TFEA.sh
Comparison_of_databases.sh
