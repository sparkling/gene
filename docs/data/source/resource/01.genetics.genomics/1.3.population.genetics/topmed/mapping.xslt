<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  TOPMed -> Population Genetics Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Population Genetics Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ CHROM                       │ chromosome              │ Normalize    │
  │ POS                         │ position                │ Direct       │
  │ REF                         │ reference_allele        │ Direct       │
  │ ALT                         │ alternate_allele        │ Direct       │
  │ allele_freq                 │ allele_frequency        │ Direct       │
  │ allele_count                │ allele_count            │ Direct       │
  │ allele_num                  │ allele_number           │ Direct       │
  │ het_count                   │ topmed_het_count        │ Direct       │
  │ hom_count                   │ topmed_hom_count        │ Direct       │
  │ NWD_ID                      │ topmed_sample_id        │ Direct       │
  │ FILTER                      │ filter_status           │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "." → null
  - Missing counts → null

  NOTES:
  - TOPMed: Trans-Omics for Precision Medicine
  - 180,000+ deeply sequenced genomes (30X WGS)
  - Freeze 8: ~700M variants
  - Authorized access required for individual-level data
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
        <xsl:value-of select="local:normalize-chromosome(fn:string[@key='CHROM'])"/>
      </fn:string>

      <fn:number key="position">
        <xsl:value-of select="fn:number[@key='POS']"/>
      </fn:number>

      <fn:string key="reference_allele">
        <xsl:value-of select="fn:string[@key='REF']"/>
      </fn:string>

      <fn:string key="alternate_allele">
        <xsl:value-of select="fn:string[@key='ALT']"/>
      </fn:string>

      <!-- ========== GLOBAL FREQUENCY METRICS ========== -->

      <!-- Allele frequency -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='allele_freq']">
          <fn:number key="allele_frequency">
            <xsl:value-of select="fn:number[@key='allele_freq']"/>
          </fn:number>
        </xsl:when>
        <xsl:when test="fn:number[@key='AF']">
          <fn:number key="allele_frequency">
            <xsl:value-of select="fn:number[@key='AF']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="allele_frequency"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Allele count -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='allele_count']">
          <fn:number key="allele_count">
            <xsl:value-of select="fn:number[@key='allele_count']"/>
          </fn:number>
        </xsl:when>
        <xsl:when test="fn:number[@key='AC']">
          <fn:number key="allele_count">
            <xsl:value-of select="fn:number[@key='AC']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="allele_count"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Allele number -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='allele_num']">
          <fn:number key="allele_number">
            <xsl:value-of select="fn:number[@key='allele_num']"/>
          </fn:number>
        </xsl:when>
        <xsl:when test="fn:number[@key='AN']">
          <fn:number key="allele_number">
            <xsl:value-of select="fn:number[@key='AN']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="allele_number"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Filter status -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='FILTER'])">
          <fn:null key="filter_status"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="filter_status">
            <xsl:value-of select="fn:string[@key='FILTER']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== TOPMED-SPECIFIC FIELDS ========== -->

      <!-- Heterozygote count -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='het_count']">
          <fn:number key="topmed_het_count">
            <xsl:value-of select="fn:number[@key='het_count']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="topmed_het_count"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Homozygote count -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='hom_count']">
          <fn:number key="topmed_hom_count">
            <xsl:value-of select="fn:number[@key='hom_count']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="topmed_hom_count"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Sample ID (authorized access only) -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='NWD_ID'])">
          <fn:null key="topmed_sample_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="topmed_sample_id">
            <xsl:value-of select="fn:string[@key='NWD_ID']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">TOPMed</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
      </fn:map>

    </fn:map>
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
    <xsl:sequence select="not($value) or normalize-space($value) = ('', '-', '.', 'NA', 'N/A', 'null')"/>
  </xsl:function>

</xsl:stylesheet>
