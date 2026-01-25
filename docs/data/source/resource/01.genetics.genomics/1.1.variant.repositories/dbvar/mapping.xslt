<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  dbVar -> Variant Repositories Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Variant Repositories Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ ID (nsv/nssv)               │ variant_id              │ Direct       │
  │ CHROM                       │ chromosome              │ Normalize    │
  │ POS                         │ position                │ Direct       │
  │ REF                         │ reference_allele        │ Direct       │
  │ ALT                         │ alternate_allele        │ Direct       │
  │ SVTYPE                      │ variant_type            │ Direct       │
  │ SVTYPE                      │ sv_type                 │ Direct       │
  │ SVLEN                       │ sv_length               │ Direct       │
  │ outermost_start             │ dbvar_outermost_start   │ Direct       │
  │ outermost_stop              │ dbvar_outermost_stop    │ Direct       │
  │ inner_start                 │ dbvar_inner_start       │ Direct       │
  │ inner_stop                  │ dbvar_inner_stop        │ Direct       │
  │ variant_count               │ sv_variant_count        │ Direct       │
  │ method                      │ detection_method        │ Direct       │
  │ analysis                    │ analysis_type           │ Direct       │
  │ platform                    │ platform                │ Direct       │
  │ study                       │ study_accession         │ Direct       │
  │ IMPRECISE                   │ is_imprecise            │ Boolean      │
  │ bin_size                    │ size_category           │ Direct       │
  │ clinical_assertion          │ clinical_significance   │ Array        │
  │ gene_symbols                │ gene_symbols            │ Array        │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, ".", "-" → null
  - Missing numeric fields → null
  - SVLEN of 0 → null

  NOTES:
  - Structural variant database - coordinates may be imprecise
  - nsv = variant region, nssv = supporting variant
  - SVLEN is negative for deletions
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
        <xsl:value-of select="fn:string[@key='ID']"/>
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

      <!-- Variant type (from SVTYPE) -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='SVTYPE'])">
          <fn:null key="variant_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="variant_type">
            <xsl:value-of select="fn:string[@key='SVTYPE']"/>
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
        <xsl:when test="fn:string[@key='gene_symbols'] and not(local:is-null(fn:string[@key='gene_symbols']))">
          <fn:array key="gene_symbols">
            <xsl:for-each select="tokenize(fn:string[@key='gene_symbols'], '[|;,]')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_symbols"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CLINICAL SIGNIFICANCE (ARRAY) ========== -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='clinical_assertion']">
          <fn:array key="clinical_significance">
            <xsl:for-each select="fn:array[@key='clinical_assertion']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='clinical_assertion'] and not(local:is-null(fn:string[@key='clinical_assertion']))">
          <fn:array key="clinical_significance">
            <xsl:for-each select="tokenize(fn:string[@key='clinical_assertion'], '[|;]')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="clinical_significance"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DBVAR-SPECIFIC FIELDS ========== -->

      <!-- Outermost start -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='outermost_start']">
          <fn:number key="dbvar_outermost_start">
            <xsl:value-of select="fn:number[@key='outermost_start']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="dbvar_outermost_start"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Outermost stop -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='outermost_stop']">
          <fn:number key="dbvar_outermost_stop">
            <xsl:value-of select="fn:number[@key='outermost_stop']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="dbvar_outermost_stop"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Inner start -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='inner_start']">
          <fn:number key="dbvar_inner_start">
            <xsl:value-of select="fn:number[@key='inner_start']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="dbvar_inner_start"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Inner stop -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='inner_stop']">
          <fn:number key="dbvar_inner_stop">
            <xsl:value-of select="fn:number[@key='inner_stop']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="dbvar_inner_stop"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- SV variant count -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='variant_count']">
          <fn:number key="sv_variant_count">
            <xsl:value-of select="fn:number[@key='variant_count']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="sv_variant_count"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Detection method -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='method'])">
          <fn:null key="detection_method"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="detection_method">
            <xsl:value-of select="fn:string[@key='method']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Analysis type -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='analysis'])">
          <fn:null key="analysis_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="analysis_type">
            <xsl:value-of select="fn:string[@key='analysis']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Platform -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='platform'])">
          <fn:null key="platform"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="platform">
            <xsl:value-of select="fn:string[@key='platform']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Study accession -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='study'])">
          <fn:null key="study_accession"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="study_accession">
            <xsl:value-of select="fn:string[@key='study']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- SV type -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='SVTYPE'])">
          <fn:null key="sv_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="sv_type">
            <xsl:value-of select="fn:string[@key='SVTYPE']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- SV length -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='SVLEN'] and fn:number[@key='SVLEN'] != 0">
          <fn:number key="sv_length">
            <xsl:value-of select="fn:number[@key='SVLEN']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="sv_length"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Is imprecise -->
      <xsl:choose>
        <xsl:when test="fn:boolean[@key='IMPRECISE']">
          <fn:boolean key="is_imprecise">
            <xsl:value-of select="fn:boolean[@key='IMPRECISE']"/>
          </fn:boolean>
        </xsl:when>
        <xsl:when test="fn:string[@key='IMPRECISE']">
          <fn:boolean key="is_imprecise">
            <xsl:value-of select="fn:string[@key='IMPRECISE'] = 'true'"/>
          </fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="is_imprecise"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Size category -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='bin_size'])">
          <fn:null key="size_category"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="size_category">
            <xsl:value-of select="fn:string[@key='bin_size']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">dbVar</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='ID']"/>
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
