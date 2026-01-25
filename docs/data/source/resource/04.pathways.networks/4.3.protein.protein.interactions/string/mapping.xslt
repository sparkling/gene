<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  STRING -> Unified PPI Schema Mapping
  ============================================================
  Source: ./schema.json (STRING protein links format)
  Target: ../schema.json (Unified 4.3 PPI Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ protein1                    │ string_id_a             │ Direct       │
  │ protein2                    │ string_id_b             │ Direct       │
  │ combined_score              │ string_combined_score   │ Normalize    │
  │ nscore                      │ string_nscore           │ Normalize    │
  │ fscore                      │ string_fscore           │ Normalize    │
  │ pscore                      │ string_pscore           │ Normalize    │
  │ ascore                      │ string_ascore           │ Normalize    │
  │ escore                      │ string_escore           │ Normalize    │
  │ dscore                      │ string_dscore           │ Normalize    │
  │ tscore                      │ string_tscore           │ Normalize    │
  │ network_type                │ network_type            │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - 0 scores remain as 0 (valid value)
  - Missing optional fields -> null

  NOTES:
  - STRING scores are 0-1000 integers, normalized to 0-1 for consistency
  - Protein IDs in format "taxid.ENSP" (e.g., 9606.ENSP00000269305)
  - Combined score integrates all evidence channels
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
      <!-- ========== REQUIRED FIELDS ========== -->

      <!-- Interactor A (STRING format: taxid.ENSP) -->
      <fn:string key="interactor_a_id">
        <xsl:value-of select="fn:string[@key='protein1']"/>
      </fn:string>

      <!-- Interactor B (STRING format: taxid.ENSP) -->
      <fn:string key="interactor_b_id">
        <xsl:value-of select="fn:string[@key='protein2']"/>
      </fn:string>

      <!-- ========== STRING-SPECIFIC IDS ========== -->

      <fn:string key="string_id_a">
        <xsl:value-of select="fn:string[@key='protein1']"/>
      </fn:string>

      <fn:string key="string_id_b">
        <xsl:value-of select="fn:string[@key='protein2']"/>
      </fn:string>

      <!-- ========== COMBINED SCORE ========== -->
      <!-- Normalize 0-1000 to 0-1 -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='combined_score']">
          <fn:number key="string_combined_score">
            <xsl:value-of select="fn:number[@key='combined_score'] div 1000"/>
          </fn:number>
          <fn:number key="confidence_score">
            <xsl:value-of select="fn:number[@key='combined_score'] div 1000"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="string_combined_score"/>
          <fn:null key="confidence_score"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== EVIDENCE CHANNEL SCORES ========== -->
      <!-- All normalized from 0-1000 to 0-1 -->

      <!-- Neighborhood score (gene proximity) -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='nscore']">
          <fn:number key="string_nscore">
            <xsl:value-of select="fn:number[@key='nscore'] div 1000"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="string_nscore"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Fusion score (gene fusion events) -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='fscore']">
          <fn:number key="string_fscore">
            <xsl:value-of select="fn:number[@key='fscore'] div 1000"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="string_fscore"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Phylogenetic profile score (co-occurrence) -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='pscore']">
          <fn:number key="string_pscore">
            <xsl:value-of select="fn:number[@key='pscore'] div 1000"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="string_pscore"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Co-expression score -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='ascore']">
          <fn:number key="string_ascore">
            <xsl:value-of select="fn:number[@key='ascore'] div 1000"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="string_ascore"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Experimental score (lab evidence) -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='escore']">
          <fn:number key="string_escore">
            <xsl:value-of select="fn:number[@key='escore'] div 1000"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="string_escore"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Database score (curated sources) -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='dscore']">
          <fn:number key="string_dscore">
            <xsl:value-of select="fn:number[@key='dscore'] div 1000"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="string_dscore"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Text-mining score (literature co-mentions) -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='tscore']">
          <fn:number key="string_tscore">
            <xsl:value-of select="fn:number[@key='tscore'] div 1000"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="string_tscore"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== NETWORK TYPE ========== -->

      <xsl:choose>
        <xsl:when test="fn:string[@key='network_type'] and fn:string[@key='network_type'] != ''">
          <fn:string key="network_type">
            <xsl:value-of select="fn:string[@key='network_type']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="network_type"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ORGANISM EXTRACTION ========== -->
      <!-- Extract taxid from protein1 (format: taxid.ENSP) -->

      <xsl:choose>
        <xsl:when test="contains(fn:string[@key='protein1'], '.')">
          <fn:number key="organism_a">
            <xsl:value-of select="xs:integer(substring-before(fn:string[@key='protein1'], '.'))"/>
          </fn:number>
          <fn:number key="organism_b">
            <xsl:value-of select="xs:integer(substring-before(fn:string[@key='protein2'], '.'))"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="organism_a"/>
          <fn:null key="organism_b"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">STRING</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="concat(fn:string[@key='protein1'], '_', fn:string[@key='protein2'])"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
