#!/usr/bin/env python
# coding: utf-8

# # Download genomic information from NCBI
# 
# ### Dependencies:
# For metadata:
# * [Entrez-direct](https://anaconda.org/bioconda/entrez-direct) v.13.9. NCBI
# For genomic datasets:
# * [NCBI datasets command line tool](https://www.ncbi.nlm.nih.gov/datasets/docs/quickstarts/command-line-tools/). v. 12.12. [Available here.](https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets
# ).
# 
# This code is written in bash and was tested excecuting in unix-terminal (Ubuntu 20.04.3 LTS)

# In[ ]:


# Download all available information in Assembly database in xml format, for each class:

esearch -db assembly -query 'Cubozoa"[Organism]'| esummary >  Cubozoa_NCBI_Assemblydb.xml
esearch -db assembly -query 'Scyphozoa"[Organism]'| esummary >  Scyphozoa_NCBI_Assemblydb.xml
esearch -db assembly -query 'Hydrozoa"[Organism]'| esummary >  Hydrozoa_NCBI_Assemblydb.xml
esearch -db assembly -query 'Staurozoa"[Organism]'| esummary >  Staurozoa_NCBI_Assemblydb.xml

# Download all available information in SRA database in xml format, for each class:
esearch -db sra -query 'Cubozoa"[Organism]'| esummary >  Cubozoa_NCBI_sradb.xml
esearch -db sra -query 'Scyphozoa"[Organism]'| esummary >  Scyphozoa_NCBI_sradb.xml
esearch -db sra -query 'Hydrozoa"[Organism]'| esummary >  Hydrozoa_NCBI_sradb.xml
esearch -db sra -query 'Staurozoa"[Organism]'| esummary >  Staurozoa_NCBI_sradb.xml


# Download all available information in Bioproject database in xml format, for each class:
esearch -db bioproject -query 'Cubozoa"[Organism]'| esummary >  Cubozoa_NCBI_bioproject.xml
esearch -db bioproject -query 'Scyphozoa"[Organism]'| esummary >  Scyphozoa_NCBI_bioproject.xml
esearch -db bioproject -query 'Hydrozoa"[Organism]'| esummary >  Hydrozoa_NCBI_bioproject.xml
esearch -db bioproject -query 'Staurozoa"[Organism]'| esummary >  Staurozoa_NCBI_bioproject.xml

# Download all available information in Biosample database in xml format, for each class:
esearch -db biosample -query 'Cubozoa"[Organism]'| esummary >  Cubozoa_NCBI_biosample.xml
esearch -db biosample -query 'Scyphozoa"[Organism]'| esummary >  Scyphozoa_NCBI_biosample.xml
esearch -db biosample -query 'Hydrozoa"[Organism]'| esummary >  Hydrozoa_NCBI_biosample.xml
esearch -db biosample -query 'Staurozoa"[Organism]'| esummary >  Staurozoa_NCBI_biosample.xml

#Download all available information in Genome database in xml format, for each class:
esearch -db genome -query 'Cubozoa"[Organism]'| esummary >  Cubozoa_NCBI_genome.xml
esearch -db genome -query 'Scyphozoa"[Organism]'| esummary >  Scyphozoa_NCBI_genome.xml
esearch -db genome -query 'Hydrozoa"[Organism]'| esummary >  Hydrozoa_NCBI_genome.xml
esearch -db genome -query 'Staurozoa"[Organism]'| esummary >  Staurozoa_NCBI_genome.xml

# Download assemblies by taxon name:
datasets download genome taxon Cubozoa --filename Cubozoa
datasets download genome taxon Scyphozoa --filename Scyphozoa
datasets download genome taxon Hydrozoa --filename Hydrozoa
datasets download genome taxon Staurozoa --filename Staurozoa


# # Assembly statistics
# ### Dependencies:
# * [BBMap from BBtools Suite](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/). Bushnell B. Last modified: July 25, 2019
# 
# 
# In a folder containing all assemblies, run:

# In[ ]:


~/PATH_TO/BBMap_38.73/bbmap/statswrapper.sh in=Aurelia.Genome_v1.2_11-27-18.fasta,GCA_000004095.1.fna,GCA_000219015.1.fna,GCA_003687565.1_ASM368756v1_genomic.fna,GCA_003864495.1.fna,GCA_003991215.1_MVIv1_genomic.fna.gz,GCA_004118115.1_Hvir_v1_genomic.fna,GCA_004118135.1_Ho_v1_genomic.fna,GCA_004194395.1.fna,GCA_004194415.1.fna,GCA_008930755.1_Aala_01_genomic.fna.gz,GCA_009936425.1_ASM993642v1_genomic.fna,GCA_010016025.1_ASM1001602v1_genomic.fna,GCA_010016065.1_ASM1001606v1_genomic.fna,GCA_011634815.1_ASM1163481v1_genomic.fna,GCA_011763395.1_ASM1176339v1_genomic.fna,GCA_012295145.1_contig.fa,GCA_012295145.1_scaffold.fa,GCA_013076295.1_ASM1307629v1_genomic.fasta,GCA_013076305.1_ASM1307630v1_genomic.fna,GCA_014526335.1_ASM1452633v1_genomic.fna,GCA_014706445.1_ASM1470644v1_genomic.fna,GCA_015164055.1_ASM1516405v1_genomic.fna,GCA_018155075.1_ASM1815507v1_genomic.fna.gz,GCA_900245855.1_ASM90024585v1_genomic.fna.gz,GCA_900291935.1_Cxam_T1-A_Genome_genomic.fna,GCA_902728285.1_Clytia_hemisphaerica_genome_assembly_genomic.fna,Hm105_Dovetail_Assembly_1.0.fa,Li_et_al_genome.fastaminscaf=200 format=6 > genome_stats.txt


# # Gene statistics estimation using Agat:
# ### Dependencies:
# * [Another GFF Analysis Toolkit (AGAT)](https://github.com/NBISweden/AGAT) - Version: v0.6.0. National Bioinformatics Infrastructure Sweden (NBIS)  
# 
# 
# 
# 
# ###1. Standarization of gff files:
# To standarize gff, gff2, gff3 files for processiong, we run the following command in a folder containing all the annotation files:
# 

# In[ ]:


for f in *;  do   agat_convert_sp_gxf2gxf.pl -g $f -o $f.agat.gff;  done


# ###2. Add intron features:
# In some cases, Agat could not find introns as features in the gff files, so we ran agat_sp_add_introns.pl, when information was available
# 
# In a folder with the following annotation files:
# * genome.gff (Li et al. 2020)
# * Aurelia.Genome_Annotation_v1.2_11-28-18.gff3 (Gold et al. 2019)
# * GCA_012295145.1_contig.gff (Xia et al. 2020)
# * Chrysaora_quinquecirrha_genomic.gff (Xia et al. 2020)
# * merged_transcript_models.agat.gff3 (Leclere et al. 2019)
# 
# 

# In[ ]:


for f in *;  agat_sp_add_introns.pl --gff $f --out $f.intron.agat.gff;  done


# 
# ### 3. Calculate statistics
# Once all annotation files were standarized, we put these files in a folder and 
# ran the following command to calculate gene statistics
# 

# In[ ]:


for f in *.agat.gff; do agat_sp_statistics.pl --gff $f -o $f.statistics; done

