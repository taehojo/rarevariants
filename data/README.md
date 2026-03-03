# Data

## Included Data

### Table_S1_gene_list.csv

List of 85 AD-associated genes analyzed in this study, selected based on the ADSP Gene Verification Committee (GVC) top hits list.

| Column | Description |
|--------|-------------|
| `gene_name` | HGNC gene symbol |
| `gene_id` | Ensembl gene ID |
| `chromosome` | Chromosome number |
| `start` | Gene start position (GRCh38) |
| `end` | Gene end position (GRCh38) |
| `strand` | Strand (+/-) |
| `length_bp` | Gene length in base pairs |
| `Priority` | Evidence tier (1-4) |
| `Cell_Type` | Predominant cell type expression category |

### Table_S2_variant_data.csv

9,943 unique rare variants (AC >= 3) mapped to 85 AD gene regions with case-control statistics and AlphaGenome-predicted regulatory effect scores.

| Column | Description |
|--------|-------------|
| `variant_id` | dbSNP rsID |
| `chr_num` | Chromosome number |
| `pos` | Genomic position (GRCh38) |
| `REF` | Reference allele |
| `ALT` | Alternative allele |
| `gene_name` | Mapped gene symbol |
| `gene_id` | Ensembl gene ID |
| `population` | Population group (NHW, AA, Asian, Hispanic) |
| `case_AF` | Case allele frequency |
| `case_AC` | Case allele count |
| `case_AN` | Case allele number |
| `ctrl_AF` | Control allele frequency |
| `ctrl_AC` | Control allele count |
| `ctrl_AN` | Control allele number |
| `total_AC` | Total allele count (case_AC + ctrl_AC) |
| `cc_ratio` | Case-control allele frequency ratio |
| `log2_cc_ratio` | Log2-transformed CC ratio |
| `enrichment` | Enrichment category (case_enriched, ctrl_enriched, case_only, ctrl_only) |
| `OR` | Odds ratio |
| `fisher_p` | Fisher's exact test p-value |
| `rna_seq_effect` | AlphaGenome RNA-seq effect score |
| `cage_effect` | AlphaGenome CAGE effect score |
| `chip_histone_effect` | AlphaGenome ChIP histone effect score |
| `dnase_effect` | AlphaGenome DNase effect score |
| `splice_junctions_effect` | AlphaGenome splice junctions effect score |
| `splice_sites_effect` | AlphaGenome splice sites effect score |
| `splice_site_usage_effect` | AlphaGenome splice site usage effect score |
| `chip_tf_effect` | AlphaGenome ChIP TF effect score |
| `ag_match_type` | AlphaGenome allele match type |

## Data Not Included

### variant_cc_with_alphgenome.csv (Primary Analysis File)

The main analysis file (18,412 rows, 33 columns) contains individual-level variant data derived from ADSP WGS and cannot be publicly shared due to data use agreements.

## Accessing ADSP Data

The individual-level genetic data used in this study were obtained from the Alzheimer's Disease Sequencing Project (ADSP).

To access ADSP data:

1. Visit [NIAGADS](https://adsp.niagads.org/)
2. Submit a data access request through the NIAGADS Data Sharing Service
3. Approval is required from the ADSP Data Access Committee
4. Once approved, WGS data (Release 4) can be downloaded

**ADSP WGS R4 includes:**
- 24,595 participants (6,296 AD cases, 18,299 controls)
- Whole genome sequencing data
- Phenotype and covariate information

For questions about data access, contact NIAGADS at https://adsp.niagads.org/.
