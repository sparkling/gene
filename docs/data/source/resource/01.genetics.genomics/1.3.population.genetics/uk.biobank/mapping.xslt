<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  UK Biobank -> Population Genetics Unified Schema Mapping
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
  │ eid                         │ ukb_participant_id      │ Direct       │
  │ field_id                    │ ukb_field_id            │ Direct       │
  │ instance                    │ ukb_instance            │ Direct       │
  │ genotyping_array            │ ukb_genotyping_array    │ Direct       │
  │ batch                       │ ukb_batch               │ Direct       │
  │ sex_inferred                │ ukb_sex_inferred        │ Direct       │
  │ heterozygosity              │ ukb_heterozygosity      │ Direct       │
  │ missingness                 │ ukb_missingness         │ Direct       │
  │ genetic_ethnic_grouping     │ ukb_genetic_ancestry    │ Direct       │
  │ FILTER                      │ filter_status           │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "." → null
  - Missing fields → null

  NOTES:
  - UK Biobank: 500,000 participants with genotype and phenotype data
  - Requires registered access through UK Biobank Access Management System
  - WGS data available for ~200,000 participants
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

      <!-- ========== UK BIOBANK-SPECIFIC FIELDS ========== -->

      <!-- Participant ID (encrypted) -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='eid']">
          <fn:number key="ukb_participant_id">
            <xsl:value-of select="fn:number[@key='eid']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ukb_participant_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Field ID -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='field_id']">
          <fn:number key="ukb_field_id">
            <xsl:value-of select="fn:number[@key='field_id']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ukb_field_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Instance -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='instance']">
          <fn:number key="ukb_instance">
            <xsl:value-of select="fn:number[@key='instance']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ukb_instance"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Genotyping array -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='genotyping_array'])">
          <fn:null key="ukb_genotyping_array"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="ukb_genotyping_array">
            <xsl:value-of select="fn:string[@key='genotyping_array']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Batch -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='batch'])">
          <fn:null key="ukb_batch"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="ukb_batch">
            <xsl:value-of select="fn:string[@key='batch']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Sex inferred -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='sex_inferred']">
          <fn:number key="ukb_sex_inferred">
            <xsl:value-of select="fn:number[@key='sex_inferred']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ukb_sex_inferred"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Heterozygosity -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='heterozygosity']">
          <fn:number key="ukb_heterozygosity">
            <xsl:value-of select="fn:number[@key='heterozygosity']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ukb_heterozygosity"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Missingness -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='missingness']">
          <fn:number key="ukb_missingness">
            <xsl:value-of select="fn:number[@key='missingness']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ukb_missingness"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Genetic ancestry -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='genetic_ethnic_grouping'])">
          <fn:null key="ukb_genetic_ancestry"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="ukb_genetic_ancestry">
            <xsl:value-of select="fn:string[@key='genetic_ethnic_grouping']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">UK Biobank</fn:string>
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
