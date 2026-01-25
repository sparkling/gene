<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  GWAS Catalog -> Expression & Regulation Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Expression & Regulation Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────────┬──────────────┐
  │ Source Field                │ Target Field                │ Transform    │
  ├─────────────────────────────┼─────────────────────────────┼──────────────┤
  │ MAPPED_GENE                 │ gene_id                     │ Direct       │
  │ rsId                        │ variant_id                  │ Direct       │
  │ rsId                        │ rsid                        │ Direct       │
  │ pvalue                      │ p_value                     │ Direct       │
  │ betaNum                     │ effect_size                 │ Direct       │
  │ betaNum                     │ beta                        │ Direct       │
  │ accessionId                 │ gwas_study_accession        │ Direct       │
  │ diseaseTrait                │ disease_trait               │ Direct       │
  │ initialSampleSize           │ initial_sample_description  │ Direct       │
  │ snpCount                    │ snp_count                   │ Direct       │
  │ imputed                     │ is_imputed                  │ Direct       │
  │ pubmedId                    │ pubmed_id                   │ Direct       │
  │ publicationDate             │ publication_date            │ Direct       │
  │ author                      │ author_info                 │ Object       │
  │ pvalueDescription           │ pvalue_description          │ Direct       │
  │ riskFrequency               │ risk_allele_frequency       │ Direct       │
  │ orPerCopyNum                │ odds_ratio                  │ Direct       │
  │ betaUnit                    │ beta_unit                   │ Direct       │
  │ range                       │ confidence_interval         │ Direct       │
  │ chromosomeName              │ chromosome                  │ Direct       │
  │ chromosomePosition          │ start_position              │ Direct       │
  │ functionalClass             │ functional_class            │ Direct       │
  │ EFO_trait                   │ efo_trait                   │ Direct       │
  │ EFO_uri                     │ efo_uri                     │ Direct       │
  │ ancestralGroups             │ ancestral_groups            │ Array        │
  └─────────────────────────────┴─────────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "NR" → null

  NOTES:
  - GWAS Catalog: NHGRI-EBI Catalog of GWAS
  - GCST accession format: GCST + digits
  - Significance threshold typically 5e-8
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
        <xsl:value-of select="(fn:string[@key='MAPPED_GENE'], fn:string[@key='gene'])[1]"/>
      </fn:string>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">variant_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='rsId']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='pvalue']">
          <fn:number key="p_value">
            <xsl:value-of select="fn:number[@key='pvalue']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="p_value"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='betaNum']">
          <fn:number key="effect_size">
            <xsl:value-of select="fn:number[@key='betaNum']"/>
          </fn:number>
        </xsl:when>
        <xsl:when test="fn:number[@key='orPerCopyNum']">
          <fn:number key="effect_size">
            <xsl:value-of select="fn:number[@key='orPerCopyNum']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="effect_size"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">chromosome</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='chromosomeName']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='chromosomePosition']">
          <fn:number key="start_position">
            <xsl:value-of select="fn:number[@key='chromosomePosition']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="start_position"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GWAS-SPECIFIC FIELDS ========== -->

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">gwas_study_accession</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='accessionId']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">disease_trait</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='diseaseTrait']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">initial_sample_description</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='initialSampleSize']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='snpCount']">
          <fn:number key="snp_count">
            <xsl:value-of select="fn:number[@key='snpCount']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="snp_count"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:boolean[@key='imputed']">
          <fn:boolean key="is_imputed">
            <xsl:value-of select="fn:boolean[@key='imputed']"/>
          </fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="is_imputed"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">pubmed_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='pubmedId']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">publication_date</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='publicationDate']"/>
      </xsl:call-template>

      <!-- Author info (preserve as object) -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='author']">
          <xsl:copy-of select="fn:map[@key='author']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="author_info"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">pvalue_description</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='pvalueDescription']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='riskFrequency']">
          <fn:number key="risk_allele_frequency">
            <xsl:value-of select="fn:number[@key='riskFrequency']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="risk_allele_frequency"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='orPerCopyNum']">
          <fn:number key="odds_ratio">
            <xsl:value-of select="fn:number[@key='orPerCopyNum']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="odds_ratio"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='betaNum']">
          <fn:number key="beta">
            <xsl:value-of select="fn:number[@key='betaNum']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="beta"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">beta_unit</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='betaUnit']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">confidence_interval</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='range']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">rsid</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='rsId']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">functional_class</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='functionalClass']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">efo_trait</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='EFO_trait']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">efo_uri</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='EFO_uri']"/>
      </xsl:call-template>

      <!-- Ancestral groups (array) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='ancestralGroups']">
          <fn:array key="ancestral_groups">
            <xsl:for-each select="fn:array[@key='ancestralGroups']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ancestral_groups"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">GWAS Catalog</fn:string>
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
    <xsl:sequence select="not($value) or normalize-space($value) = ('', '-', '.', 'NA', 'N/A', 'null', 'NR')"/>
  </xsl:function>

</xsl:stylesheet>
