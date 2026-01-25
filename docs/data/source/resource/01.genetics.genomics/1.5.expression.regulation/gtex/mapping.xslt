<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  GTEx -> Expression & Regulation Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Expression & Regulation Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────────┬──────────────┐
  │ Source Field                │ Target Field                │ Transform    │
  ├─────────────────────────────┼─────────────────────────────┼──────────────┤
  │ gene_id                     │ gene_id                     │ Direct       │
  │ gene_symbol                 │ gene_symbol                 │ Direct       │
  │ variant_id                  │ variant_id                  │ Direct       │
  │ pval_nominal                │ p_value                     │ Direct       │
  │ pval_nominal                │ nominal_pvalue              │ Direct       │
  │ slope                       │ effect_size                 │ Direct       │
  │ slope                       │ slope                       │ Direct       │
  │ slope_se                    │ slope_se                    │ Direct       │
  │ pval_beta                   │ beta_adjusted_pvalue        │ Direct       │
  │ SUBJID                      │ gtex_subject_id             │ Direct       │
  │ SEX                         │ gtex_sex                    │ Direct       │
  │ AGE                         │ gtex_age                    │ Direct       │
  │ DTHHRDY                     │ gtex_death_classification   │ Direct       │
  │ SAMPID                      │ gtex_sample_id              │ Direct       │
  │ SMTS                        │ tissue_type_general         │ Direct       │
  │ SMTSD                       │ tissue_type_detailed        │ Direct       │
  │ SMRIN                       │ rna_integrity_number        │ Direct       │
  │ SMAFRZE                     │ analysis_freeze             │ Direct       │
  │ TPM                         │ tpm                         │ Direct       │
  │ read_counts                 │ read_counts                 │ Direct       │
  │ tss_distance                │ tss_distance                │ Direct       │
  │ ma_samples                  │ minor_allele_samples        │ Direct       │
  │ ma_count                    │ minor_allele_count          │ Direct       │
  │ maf                         │ minor_allele_frequency      │ Direct       │
  │ phenotype_id                │ sqtl_phenotype_id           │ Direct       │
  │ chr                         │ chromosome                  │ Direct       │
  └─────────────────────────────┴─────────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "NA" → null

  NOTES:
  - GTEx: Genotype-Tissue Expression Project
  - V8: 948 donors, 17,382 samples, 54 tissues
  - eQTL and sQTL data available
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

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== CORE FIELDS ========== -->
      <fn:string key="gene_id">
        <xsl:value-of select="fn:string[@key='gene_id']"/>
      </fn:string>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">gene_symbol</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='gene_symbol']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">variant_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='variant_id']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">chromosome</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='chr']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='pval_nominal']">
          <fn:number key="p_value">
            <xsl:value-of select="fn:number[@key='pval_nominal']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="p_value"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='slope']">
          <fn:number key="effect_size">
            <xsl:value-of select="fn:number[@key='slope']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="effect_size"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GTEX-SPECIFIC FIELDS ========== -->

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">gtex_subject_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='SUBJID']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='SEX']">
          <fn:number key="gtex_sex">
            <xsl:value-of select="fn:number[@key='SEX']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gtex_sex"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">gtex_age</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='AGE']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='DTHHRDY']">
          <fn:number key="gtex_death_classification">
            <xsl:value-of select="fn:number[@key='DTHHRDY']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gtex_death_classification"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">gtex_sample_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='SAMPID']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">tissue_type_general</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='SMTS']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">tissue_type_detailed</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='SMTSD']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='SMRIN']">
          <fn:number key="rna_integrity_number">
            <xsl:value-of select="fn:number[@key='SMRIN']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="rna_integrity_number"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">analysis_freeze</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='SMAFRZE']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='TPM']">
          <fn:number key="tpm">
            <xsl:value-of select="fn:number[@key='TPM']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="tpm"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='read_counts']">
          <fn:number key="read_counts">
            <xsl:value-of select="fn:number[@key='read_counts']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="read_counts"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='tss_distance']">
          <fn:number key="tss_distance">
            <xsl:value-of select="fn:number[@key='tss_distance']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="tss_distance"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='ma_samples']">
          <fn:number key="minor_allele_samples">
            <xsl:value-of select="fn:number[@key='ma_samples']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="minor_allele_samples"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='ma_count']">
          <fn:number key="minor_allele_count">
            <xsl:value-of select="fn:number[@key='ma_count']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="minor_allele_count"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='maf']">
          <fn:number key="minor_allele_frequency">
            <xsl:value-of select="fn:number[@key='maf']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="minor_allele_frequency"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='pval_nominal']">
          <fn:number key="nominal_pvalue">
            <xsl:value-of select="fn:number[@key='pval_nominal']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="nominal_pvalue"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='slope']">
          <fn:number key="slope">
            <xsl:value-of select="fn:number[@key='slope']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="slope"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='slope_se']">
          <fn:number key="slope_se">
            <xsl:value-of select="fn:number[@key='slope_se']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="slope_se"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='pval_beta']">
          <fn:number key="beta_adjusted_pvalue">
            <xsl:value-of select="fn:number[@key='pval_beta']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="beta_adjusted_pvalue"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">sqtl_phenotype_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='phenotype_id']"/>
      </xsl:call-template>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">GTEx</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

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

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or normalize-space($value) = ('', '-', '.', 'NA', 'N/A', 'null')"/>
  </xsl:function>

</xsl:stylesheet>
