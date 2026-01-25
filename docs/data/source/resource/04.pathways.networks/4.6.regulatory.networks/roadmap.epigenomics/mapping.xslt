<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Roadmap Epigenomics -> Unified Regulatory Networks Schema Mapping
  ============================================================
  Source: ./schema.json (Roadmap Epigenomics data)
  Target: ../schema.json (Unified 4.6 Regulatory Networks Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ eid                         │ epigenome_id            │ Direct       │
  │ standardName                │ sample_name             │ Direct       │
  │ group                       │ cell_type_group         │ Direct       │
  │ anatomy                     │ anatomical_source       │ Direct       │
  │ type                        │ sample_type             │ Direct       │
  │ histone_marks               │ available_marks         │ Restructure  │
  │ chromatinState              │ roadmap_chromatin_state │ Direct       │
  │ dnaMethylation              │ roadmap_methylation     │ Direct       │
  │ peak                        │ roadmap_peak            │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings "" -> null
  - "not_available" -> null

  NOTES:
  - Epigenome IDs follow E### pattern
  - 15-state chromatin model used for state annotation
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
      <!-- ========== DIRECT FIELD MAPPINGS ========== -->

      <!-- Epigenome identifier -->
      <fn:string key="epigenome_id">
        <xsl:value-of select="fn:string[@key='eid']"/>
      </fn:string>

      <!-- Sample name -->
      <fn:string key="sample_name">
        <xsl:value-of select="fn:string[@key='standardName']"/>
      </fn:string>

      <!-- Cell type group -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='group'] and fn:string[@key='group'] != ''">
          <fn:string key="cell_type_group">
            <xsl:value-of select="fn:string[@key='group']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cell_type_group"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Anatomical source -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='anatomy'] and fn:string[@key='anatomy'] != ''">
          <fn:string key="anatomical_source">
            <xsl:value-of select="fn:string[@key='anatomy']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="anatomical_source"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Sample type -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='type'] and fn:string[@key='type'] != ''">
          <fn:string key="sample_type">
            <xsl:value-of select="fn:string[@key='type']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="sample_type"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== TRANSFORMED FIELDS ========== -->

      <!-- Available histone marks (convert to array of available marks) -->
      <xsl:if test="fn:map[@key='histone_marks']">
        <fn:array key="available_marks">
          <xsl:for-each select="fn:map[@key='histone_marks']/*[. = 'available']">
            <fn:string><xsl:value-of select="@key"/></fn:string>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- ========== ROADMAP-SPECIFIC FIELDS ========== -->

      <!-- Chromatin state -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='chromatinState']">
          <fn:map key="roadmap_chromatin_state">
            <xsl:if test="fn:map[@key='chromatinState']/fn:string[@key='chrom']">
              <fn:string key="chromosome">
                <xsl:value-of select="fn:map[@key='chromatinState']/fn:string[@key='chrom']"/>
              </fn:string>
            </xsl:if>
            <xsl:if test="fn:map[@key='chromatinState']/fn:number[@key='chromStart']">
              <fn:number key="start">
                <xsl:value-of select="fn:map[@key='chromatinState']/fn:number[@key='chromStart']"/>
              </fn:number>
            </xsl:if>
            <xsl:if test="fn:map[@key='chromatinState']/fn:number[@key='chromEnd']">
              <fn:number key="end">
                <xsl:value-of select="fn:map[@key='chromatinState']/fn:number[@key='chromEnd']"/>
              </fn:number>
            </xsl:if>
            <xsl:if test="fn:map[@key='chromatinState']/fn:string[@key='state']">
              <fn:string key="state">
                <xsl:value-of select="fn:map[@key='chromatinState']/fn:string[@key='state']"/>
              </fn:string>
            </xsl:if>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="roadmap_chromatin_state"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- DNA methylation -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='dnaMethylation']">
          <fn:map key="roadmap_methylation">
            <xsl:copy-of select="fn:map[@key='dnaMethylation']/*"/>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="roadmap_methylation"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Peak data -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='peak']">
          <fn:map key="roadmap_peak">
            <xsl:copy-of select="fn:map[@key='peak']/*"/>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="roadmap_peak"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">Roadmap Epigenomics</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='eid']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Reusable Functions
       ============================================================ -->

  <!-- Convert empty/placeholder to empty string (for null handling) -->
  <xsl:function name="local:nullify" as="xs:string">
    <xsl:param name="value" as="xs:string"/>
    <xsl:choose>
      <xsl:when test="$value = ('-', '.', '', 'NA', 'N/A', 'null', 'not_available')"></xsl:when>
      <xsl:otherwise><xsl:value-of select="$value"/></xsl:otherwise>
    </xsl:choose>
  </xsl:function>

</xsl:stylesheet>
