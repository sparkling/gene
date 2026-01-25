<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  dbNSFP -> Functional Prediction Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Functional Prediction Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────────┬──────────────┐
  │ Source Field                │ Target Field                │ Transform    │
  ├─────────────────────────────┼─────────────────────────────┼──────────────┤
  │ chr                         │ chromosome                  │ Normalize    │
  │ pos                         │ position                    │ Direct       │
  │ ref                         │ reference_allele            │ Direct       │
  │ alt                         │ alternate_allele            │ Direct       │
  │ aaref                       │ dbnsfp_ref_amino_acid       │ Direct       │
  │ aaalt                       │ dbnsfp_alt_amino_acid       │ Direct       │
  │ genename                    │ gene_symbols                │ Split array  │
  │ Ensembl_geneid              │ ensembl_gene_id             │ Direct       │
  │ Ensembl_transcriptid        │ ensembl_transcript_id       │ Direct       │
  │ Ensembl_proteinid           │ ensembl_protein_id          │ Direct       │
  │ Uniprot_acc                 │ uniprot_accession           │ Direct       │
  │ HGVSc_snpEff                │ hgvs_coding                 │ Direct       │
  │ HGVSp_snpEff                │ hgvs_protein                │ Direct       │
  │ SIFT_score                  │ sift_score                  │ Direct       │
  │ SIFT_pred                   │ sift_prediction             │ Direct       │
  │ Polyphen2_HDIV_score        │ polyphen2_hdiv_score        │ Direct       │
  │ Polyphen2_HDIV_pred         │ polyphen2_hdiv_prediction   │ Direct       │
  │ Polyphen2_HVAR_score        │ polyphen2_hvar_score        │ Direct       │
  │ Polyphen2_HVAR_pred         │ polyphen2_hvar_prediction   │ Direct       │
  │ LRT_score                   │ lrt_score                   │ Direct       │
  │ MutationTaster_score        │ mutationtaster_score        │ Direct       │
  │ MutationAssessor_score      │ mutationassessor_score      │ Direct       │
  │ FATHMM_score                │ fathmm_score                │ Direct       │
  │ PROVEAN_score               │ provean_score               │ Direct       │
  │ MetaSVM_score               │ metasvm_score               │ Direct       │
  │ MetaLR_score                │ metalr_score                │ Direct       │
  │ REVEL_score                 │ revel_score                 │ Direct       │
  │ CADD_phred                  │ cadd_phred                  │ Direct       │
  │ AlphaMissense_score         │ alphamissense_pathogenicity │ Direct       │
  │ phyloP100way_vertebrate     │ phylop_100way               │ Direct       │
  │ phastCons100way_vertebrate  │ phastcons_100way            │ Direct       │
  │ GERP++_RS                   │ gerp_rs                     │ Direct       │
  │ gnomAD_exomes_AF            │ gnomad_exomes_af            │ Direct       │
  │ gnomAD_genomes_AF           │ gnomad_genomes_af           │ Direct       │
  │ 1000Gp3_AF                  │ thousand_genomes_af         │ Direct       │
  │ TOPMed_AF                   │ topmed_af                   │ Direct       │
  │ clinvar_clnsig              │ clinvar_significance        │ Direct       │
  │ clinvar_review              │ clinvar_review              │ Direct       │
  └─────────────────────────────┴─────────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, ".", "-" → null
  - Scores with multiple values (;-separated) take first value

  NOTES:
  - dbNSFP aggregates 40+ prediction scores
  - Multiple transcript annotations separated by ";"
  - First value is typically canonical/MANE transcript
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    xmlns:map="http://www.w3.org/2005/xpath-functions/map"
    xmlns:array="http://www.w3.org/2005/xpath-functions/array"
    xmlns:local="http://local.functions"
    exclude-result-prefixes="xs fn map array local">

  <xsl:output method="json" indent="yes"/>

  <!-- ============================================================
       Entry Point: Parse JSON input
       ============================================================ -->
  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <!-- ============================================================
       Main Record Transformation
       ============================================================ -->
  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== CORE VARIANT FIELDS ========== -->
      <fn:string key="chromosome">
        <xsl:value-of select="local:normalize-chromosome(fn:string[@key='chr'])"/>
      </fn:string>

      <fn:number key="position">
        <xsl:value-of select="fn:number[@key='pos']"/>
      </fn:number>

      <fn:string key="reference_allele">
        <xsl:value-of select="fn:string[@key='ref']"/>
      </fn:string>

      <fn:string key="alternate_allele">
        <xsl:value-of select="fn:string[@key='alt']"/>
      </fn:string>

      <!-- ========== GENE SYMBOLS (ARRAY) ========== -->
      <xsl:variable name="genes" select="fn:string[@key='genename']"/>
      <xsl:choose>
        <xsl:when test="local:is-null($genes)">
          <fn:null key="gene_symbols"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="gene_symbols">
            <xsl:for-each select="tokenize($genes, ';')">
              <xsl:if test="not(local:is-null(.))">
                <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
              </xsl:if>
            </xsl:for-each>
          </fn:array>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== TRANSCRIPT IDS (ARRAY) ========== -->
      <xsl:variable name="transcripts" select="fn:string[@key='Ensembl_transcriptid']"/>
      <xsl:choose>
        <xsl:when test="local:is-null($transcripts)">
          <fn:null key="transcript_ids"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="transcript_ids">
            <xsl:for-each select="tokenize($transcripts, ';')">
              <xsl:if test="not(local:is-null(.))">
                <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
              </xsl:if>
            </xsl:for-each>
          </fn:array>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DBNSFP-SPECIFIC FIELDS ========== -->

      <!-- Amino acids -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">dbnsfp_ref_amino_acid</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='aaref']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">dbnsfp_alt_amino_acid</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='aaalt']"/>
      </xsl:call-template>

      <!-- Ensembl IDs -->
      <xsl:call-template name="first-value-or-null">
        <xsl:with-param name="key">ensembl_gene_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Ensembl_geneid']"/>
      </xsl:call-template>

      <xsl:call-template name="first-value-or-null">
        <xsl:with-param name="key">ensembl_transcript_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Ensembl_transcriptid']"/>
      </xsl:call-template>

      <xsl:call-template name="first-value-or-null">
        <xsl:with-param name="key">ensembl_protein_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Ensembl_proteinid']"/>
      </xsl:call-template>

      <xsl:call-template name="first-value-or-null">
        <xsl:with-param name="key">uniprot_accession</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Uniprot_acc']"/>
      </xsl:call-template>

      <!-- HGVS -->
      <xsl:call-template name="first-value-or-null">
        <xsl:with-param name="key">hgvs_coding</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='HGVSc_snpEff']"/>
      </xsl:call-template>

      <xsl:call-template name="first-value-or-null">
        <xsl:with-param name="key">hgvs_protein</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='HGVSp_snpEff']"/>
      </xsl:call-template>

      <!-- ========== PREDICTION SCORES ========== -->

      <!-- SIFT -->
      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">sift_score</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='SIFT_score']"/>
      </xsl:call-template>

      <xsl:call-template name="first-value-or-null">
        <xsl:with-param name="key">sift_prediction</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='SIFT_pred']"/>
      </xsl:call-template>

      <!-- PolyPhen2 -->
      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">polyphen2_hdiv_score</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Polyphen2_HDIV_score']"/>
      </xsl:call-template>

      <xsl:call-template name="first-value-or-null">
        <xsl:with-param name="key">polyphen2_hdiv_prediction</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Polyphen2_HDIV_pred']"/>
      </xsl:call-template>

      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">polyphen2_hvar_score</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Polyphen2_HVAR_score']"/>
      </xsl:call-template>

      <xsl:call-template name="first-value-or-null">
        <xsl:with-param name="key">polyphen2_hvar_prediction</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Polyphen2_HVAR_pred']"/>
      </xsl:call-template>

      <!-- Other prediction scores -->
      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">lrt_score</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='LRT_score']"/>
      </xsl:call-template>

      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">mutationtaster_score</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='MutationTaster_score']"/>
      </xsl:call-template>

      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">mutationassessor_score</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='MutationAssessor_score']"/>
      </xsl:call-template>

      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">fathmm_score</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='FATHMM_score']"/>
      </xsl:call-template>

      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">provean_score</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='PROVEAN_score']"/>
      </xsl:call-template>

      <!-- Ensemble scores -->
      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">metasvm_score</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='MetaSVM_score']"/>
      </xsl:call-template>

      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">metalr_score</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='MetaLR_score']"/>
      </xsl:call-template>

      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">revel_score</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='REVEL_score']"/>
      </xsl:call-template>

      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">cadd_phred</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='CADD_phred']"/>
      </xsl:call-template>

      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">alphamissense_pathogenicity</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='AlphaMissense_score']"/>
      </xsl:call-template>

      <!-- Conservation scores -->
      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">phylop_100way</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='phyloP100way_vertebrate']"/>
      </xsl:call-template>

      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">phastcons_100way</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='phastCons100way_vertebrate']"/>
      </xsl:call-template>

      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">gerp_rs</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='GERP++_RS']"/>
      </xsl:call-template>

      <!-- Population frequencies -->
      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">gnomad_exomes_af</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='gnomAD_exomes_AF']"/>
      </xsl:call-template>

      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">gnomad_genomes_af</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='gnomAD_genomes_AF']"/>
      </xsl:call-template>

      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">thousand_genomes_af</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='1000Gp3_AF']"/>
      </xsl:call-template>

      <xsl:call-template name="first-score-or-null">
        <xsl:with-param name="key">topmed_af</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='TOPMed_AF']"/>
      </xsl:call-template>

      <!-- ClinVar -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">clinvar_significance</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='clinvar_clnsig']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">clinvar_review</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='clinvar_review']"/>
      </xsl:call-template>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">dbNSFP</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Named Templates
       ============================================================ -->

  <xsl:template name="string-or-null">
    <xsl:param name="key"/>
    <xsl:param name="value"/>
    <xsl:choose>
      <xsl:when test="local:is-null($value)">
        <fn:null key="{$key}"/>
      </xsl:when>
      <xsl:otherwise>
        <fn:string key="{$key}">
          <xsl:value-of select="$value"/>
        </fn:string>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="first-value-or-null">
    <xsl:param name="key"/>
    <xsl:param name="value"/>
    <xsl:variable name="first" select="tokenize($value, ';')[1]"/>
    <xsl:choose>
      <xsl:when test="local:is-null($first)">
        <fn:null key="{$key}"/>
      </xsl:when>
      <xsl:otherwise>
        <fn:string key="{$key}">
          <xsl:value-of select="normalize-space($first)"/>
        </fn:string>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="first-score-or-null">
    <xsl:param name="key"/>
    <xsl:param name="value"/>
    <xsl:variable name="first" select="tokenize($value, ';')[1]"/>
    <xsl:choose>
      <xsl:when test="local:is-null($first)">
        <fn:null key="{$key}"/>
      </xsl:when>
      <xsl:otherwise>
        <fn:number key="{$key}">
          <xsl:value-of select="number(normalize-space($first))"/>
        </fn:number>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <!-- ============================================================
       Reusable Functions
       ============================================================ -->

  <xsl:function name="local:normalize-chromosome" as="xs:string">
    <xsl:param name="chr" as="xs:string"/>
    <xsl:value-of select="replace($chr, '^chr', '')"/>
  </xsl:function>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or normalize-space($value) = ('', '-', '.', 'NA', 'N/A', 'null', '.')"/>
  </xsl:function>

</xsl:stylesheet>
