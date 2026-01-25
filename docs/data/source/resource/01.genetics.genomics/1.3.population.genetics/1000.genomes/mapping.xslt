<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  1000 Genomes -> Population Genetics Unified Schema Mapping
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
  │ AF                          │ allele_frequency        │ Direct       │
  │ AC                          │ allele_count            │ Direct       │
  │ AN                          │ allele_number           │ Direct       │
  │ AFR_AF                      │ af_african              │ Direct       │
  │ AMR_AF                      │ af_american             │ Direct       │
  │ EAS_AF                      │ af_east_asian           │ Direct       │
  │ EUR_AF                      │ af_european             │ Direct       │
  │ SAS_AF                      │ af_south_asian          │ Direct       │
  │ FILTER                      │ filter_status           │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "." → null
  - Missing frequency fields → null

  NOTES:
  - 1000 Genomes Phase 3: 2,504 individuals from 26 populations
  - 5 super-populations: AFR, AMR, EAS, EUR, SAS
  - Reference genome: GRCh38
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

      <!-- ========== POPULATION-SPECIFIC FREQUENCIES ========== -->

      <!-- African -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='AFR_AF']">
          <fn:number key="af_african">
            <xsl:value-of select="fn:number[@key='AFR_AF']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="af_african"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- American/Latino -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='AMR_AF']">
          <fn:number key="af_american">
            <xsl:value-of select="fn:number[@key='AMR_AF']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="af_american"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- East Asian -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='EAS_AF']">
          <fn:number key="af_east_asian">
            <xsl:value-of select="fn:number[@key='EAS_AF']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="af_east_asian"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- European -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='EUR_AF']">
          <fn:number key="af_european">
            <xsl:value-of select="fn:number[@key='EUR_AF']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="af_european"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- South Asian -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='SAS_AF']">
          <fn:number key="af_south_asian">
            <xsl:value-of select="fn:number[@key='SAS_AF']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="af_south_asian"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">1000 Genomes</fn:string>
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
