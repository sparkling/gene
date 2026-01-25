<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  dbSNP -> Variant Repositories Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Variant Repositories Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ refsnp_id                   │ variant_id              │ rs prefix    │
  │ refsnp_id                   │ dbsnp_refsnp_id         │ Direct       │
  │ CHROM                       │ chromosome              │ Normalize    │
  │ POS                         │ position                │ Direct       │
  │ REF                         │ reference_allele        │ Direct       │
  │ ALT                         │ alternate_allele        │ Direct       │
  │ variant_type                │ variant_type            │ Direct       │
  │ gene_symbols                │ gene_symbols            │ Array        │
  │ hgvs_expressions            │ hgvs_expressions        │ Array        │
  │ assembly                    │ assembly                │ Direct       │
  │ create_date                 │ dbsnp_create_date       │ Direct       │
  │ last_update_date            │ dbsnp_last_update       │ Direct       │
  │ mane_select_ids             │ mane_select_ids         │ Array        │
  │ allele_count                │ allele_count            │ Direct       │
  │ total_count                 │ total_count             │ Direct       │
  │ citations                   │ citations               │ Array        │
  │ ALFA_frequencies            │ population_frequencies  │ Object       │
  │ clinical_significance       │ clinical_significance   │ Array        │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, ".", "-" → null
  - Missing numeric fields → null

  NOTES:
  - RS ID used as primary variant_id with "rs" prefix
  - Supports both JSON API format and VCF format
  - ALFA frequencies preserved as nested object
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
      <!-- ========== PRIMARY IDENTIFIER ========== -->
      <fn:string key="variant_id">
        <xsl:value-of select="concat('rs', fn:string[@key='refsnp_id'])"/>
      </fn:string>

      <!-- ========== CORE VARIANT FIELDS ========== -->
      <fn:string key="chromosome">
        <xsl:value-of select="local:normalize-chromosome(fn:string[@key='CHROM'])"/>
      </fn:string>

      <fn:number key="position">
        <xsl:value-of select="fn:number[@key='POS']"/>
      </fn:number>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='REF'])">
          <fn:null key="reference_allele"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="reference_allele">
            <xsl:value-of select="fn:string[@key='REF']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='ALT'])">
          <fn:null key="alternate_allele"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="alternate_allele">
            <xsl:value-of select="fn:string[@key='ALT']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Variant type -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='variant_type'])">
          <fn:null key="variant_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="variant_type">
            <xsl:value-of select="fn:string[@key='variant_type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GENE SYMBOLS (ARRAY) ========== -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='gene_symbols']">
          <fn:array key="gene_symbols">
            <xsl:for-each select="fn:array[@key='gene_symbols']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_symbols"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== HGVS EXPRESSIONS (ARRAY) ========== -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='hgvs_expressions']">
          <fn:array key="hgvs_expressions">
            <xsl:for-each select="fn:array[@key='hgvs_expressions']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="hgvs_expressions"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Assembly -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='assembly'])">
          <fn:null key="assembly"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="assembly">
            <xsl:value-of select="fn:string[@key='assembly']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CLINICAL SIGNIFICANCE (ARRAY) ========== -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='clinical_significance']">
          <fn:array key="clinical_significance">
            <xsl:for-each select="fn:array[@key='clinical_significance']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='clinical_significance'] and not(local:is-null(fn:string[@key='clinical_significance']))">
          <fn:array key="clinical_significance">
            <xsl:for-each select="tokenize(fn:string[@key='clinical_significance'], '[|;,]')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="clinical_significance"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DBSNP-SPECIFIC FIELDS ========== -->
      <fn:string key="dbsnp_refsnp_id">
        <xsl:value-of select="fn:string[@key='refsnp_id']"/>
      </fn:string>

      <!-- Create date -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='create_date'])">
          <fn:null key="dbsnp_create_date"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="dbsnp_create_date">
            <xsl:value-of select="fn:string[@key='create_date']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Last update date -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='last_update_date'])">
          <fn:null key="dbsnp_last_update"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="dbsnp_last_update">
            <xsl:value-of select="fn:string[@key='last_update_date']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- MANE Select IDs -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='mane_select_ids']">
          <fn:array key="mane_select_ids">
            <xsl:for-each select="fn:array[@key='mane_select_ids']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="mane_select_ids"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Allele count -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='allele_count']">
          <fn:number key="allele_count">
            <xsl:value-of select="fn:number[@key='allele_count']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="allele_count"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Total count -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='total_count']">
          <fn:number key="total_count">
            <xsl:value-of select="fn:number[@key='total_count']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="total_count"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Citations -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='citations']">
          <fn:array key="citations">
            <xsl:for-each select="fn:array[@key='citations']/fn:number">
              <fn:number><xsl:value-of select="."/></fn:number>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="citations"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Population frequencies (preserve as object) -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='ALFA_frequencies']">
          <xsl:copy-of select="fn:map[@key='ALFA_frequencies']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="population_frequencies"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">dbSNP</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="concat('rs', fn:string[@key='refsnp_id'])"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Reusable Functions
       ============================================================ -->

  <!-- Normalize chromosome: remove "chr" prefix -->
  <xsl:function name="local:normalize-chromosome" as="xs:string">
    <xsl:param name="chr" as="xs:string"/>
    <xsl:value-of select="replace($chr, '^chr', '')"/>
  </xsl:function>

  <!-- Check if value is null/empty -->
  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = ('', '-', '.', 'NA', 'N/A', 'null')"/>
  </xsl:function>

</xsl:stylesheet>
