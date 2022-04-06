# Soil Pyruvate Experiment 
## Effect of drought on microbial metabolic activity and VOC production, assessed using multi-omics and stable isotope probing in the tropical rainforest at Biosphere 2

### **Description of Project** <br>
Due to climate change, the frequency and duration of droughts in
tropical rain forests (TRFs) are expected to increase, having a
significant impact on soil carbon dynamics . The role of microbes as
drivers of changing carbon flow, particularly in relation to metabolic
pathways and volatile organic compounds (VOCs), remains largely
unknown.  

#### **Biosphere 2 drought experiment**<br>
To examine how microbial activity, particularly in relation to carbon allocation, shifts
during drought, we utilized the controlable conditions of the glass and steelenclosed
TRF at Biosphere 2 (Fig 1). Here we created a 66-day drought in order to
study carbon cycling by plants and microbes before, during, and after drought (Water,
Atmosphere, and Life Dynamics [WALD]) . As part of this larger study, we perfomed
multi-omics and traced carbon allocation by soil microbes using position-specific (C1 or
C2) C-pyruvate labeling before and during drought. Our results will inform key
processes in tropical soil carbon cycling. For more detailed description of project, see wiki page link below<br>

**URL of \<Code\>  Section of Repository:** https://github.com/linneakh/SoilPyruvate

**URL of the Wiki Section of Repository:** https://github.com/linneakh/SoilPyruvate/wiki

![Carbon Cycle](soil-carbon-cycling.png?raw=true)

***

### Description of Experiment

**Soil C-pyruvate labeling:**
- C1 or C2 postition-specific 13C-pyruvate was added to soil
within automatic chambers located at three sites of the Biosphere 2 TRF in
order to track how carbon was allocated into 13C-CO2 and 13C-VOCs.
- 13CO2 and C-VOCs were measured using a Licor8100 coupled to Picarro G2201
and proton-transfer-reaction time-of-flight mass spectrometry (PTR-TOF-MS)
- Soil samples were collected at time = 0, 6, and 48 hr post pyruvate addition for metabolomics, metagenomics, and metatranscriptomics.

### Description of Overall Datasets

**12/13C-CO2**

**12/13C-VOCs**

**MetaT**
- All metatranscriptomic data including raw and processed files are stored on IMG
- gene expression data analyzed here was downloaded from IMG as gene copies per kegg orthology (KO)

**MetaG**
- All metagenomic data including raw and processed files are stored on IMG.
- Gene data analyzed here was downloaded from IMG as gene copies per kegg orthology (KO), unaltered downloaded gene copies are in the raw data folder.
- Raw gene copy table with control samples removed were input to deseq2 analysis.

**MetaB**
- FTICR - 
- NMR

### Description of Data Files Included here
- ClusterProfiler/modules-gene-long.csv (Table of KO ids and the module [obtained from WGCNA] it belongs to)
- Deseq2/Drought_vs_predrought_metaT_missing_1_no_feature.csv (Gene copies for metaT, obtained from raw download minus water control samples)
- Deseq2/metaT-metadata.csv (Metadata information for metaT gene copy table)
- 

***
Good practice:
Updated: 04/6/2022, L. Honeker
