<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  STITCH -> Unified Drug-Target Interactions Schema Mapping
  ============================================================
  Source: ./schema.json (STITCH chemical-protein interactions)
  Target: ../schema.json (Unified 4.4 Drug-Target Interactions Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ chemical                    │ chemical_id             │ Direct       │
  │ protein                     │ protein_id              │ Direct       │
  │ preferredName_chemical      │ chemical_name           │ Direct       │
  │ preferredName_protein       │ gene_symbol             │ Direct       │
  │ ncbiTaxonId                 │ organism_taxon          │ Direct       │
  │ combined_score              │ confidence_score        │ Scale /1000  │
  │ actions                     │ interaction_actions     │ Restructure  │
  │ experimental                │ stitch_experimental     │ Direct       │
  │ prediction                  │ stitch_prediction       │ Direct       │
  │ database                    │ stitch_database         │ Direct       │
  │ textmining                  │ stitch_textmining       │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings "" -> null
  - Score 0 -> null for channel scores

  NOTES:
  - STITCH chemical IDs have CIDm/CIDs prefixes
  - Scores are 0-1000, normalized to 0-1 for unified schema
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

      <!-- Chemical identifier -->
      <fn:string key="chemical_id">
        <xsl:value-of select="fn:string[@key='chemical']"/>
      </fn:string>

      <!-- Protein identifier -->
      <fn:string key="protein_id">
        <xsl:value-of select="fn:string[@key='protein']"/>
      </fn:string>

      <!-- Chemical name -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='preferredName_chemical'] and fn:string[@key='preferredName_chemical'] != ''">
          <fn:string key="chemical_name">
            <xsl:value-of select="fn:string[@key='preferredName_chemical']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="chemical_name"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Gene symbol -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='preferredName_protein'] and fn:string[@key='preferredName_protein'] != ''">
          <fn:string key="gene_symbol">
            <xsl:value-of select="fn:string[@key='preferredName_protein']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_symbol"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Organism taxon -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='ncbiTaxonId']">
          <fn:number key="organism_taxon">
            <xsl:value-of select="fn:number[@key='ncbiTaxonId']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="organism_taxon"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== TRANSFORMED FIELDS ========== -->

      <!-- Combined score (normalize 0-1000 to 0-1) -->
      <fn:number key="confidence_score">
        <xsl:value-of select="fn:number[@key='combined_score'] div 1000"/>
      </fn:number>

      <!-- ========== ARRAY TRANSFORMATIONS ========== -->

      <!-- Interaction actions -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='actions']">
          <fn:array key="interaction_actions">
            <xsl:for-each select="fn:array[@key='actions']/fn:map">
              <fn:map>
                <xsl:if test="fn:string[@key='mode']">
                  <fn:string key="mode">
                    <xsl:value-of select="fn:string[@key='mode']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:string[@key='action']">
                  <fn:string key="action">
                    <xsl:value-of select="fn:string[@key='action']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:boolean[@key='is_directional']">
                  <fn:boolean key="directional">
                    <xsl:value-of select="fn:boolean[@key='is_directional']"/>
                  </fn:boolean>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="interaction_actions"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== STITCH-SPECIFIC FIELDS ========== -->

      <!-- Evidence channel scores (keep as 0-1000) -->
      <xsl:if test="fn:number[@key='experimental']">
        <fn:number key="stitch_experimental">
          <xsl:value-of select="fn:number[@key='experimental']"/>
        </fn:number>
      </xsl:if>

      <xsl:if test="fn:number[@key='prediction']">
        <fn:number key="stitch_prediction">
          <xsl:value-of select="fn:number[@key='prediction']"/>
        </fn:number>
      </xsl:if>

      <xsl:if test="fn:number[@key='database']">
        <fn:number key="stitch_database">
          <xsl:value-of select="fn:number[@key='database']"/>
        </fn:number>
      </xsl:if>

      <xsl:if test="fn:number[@key='textmining']">
        <fn:number key="stitch_textmining">
          <xsl:value-of select="fn:number[@key='textmining']"/>
        </fn:number>
      </xsl:if>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">STITCH</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_chemical_id">
          <xsl:value-of select="fn:string[@key='chemical']"/>
        </fn:string>
        <fn:string key="original_protein_id">
          <xsl:value-of select="fn:string[@key='protein']"/>
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
      <xsl:when test="$value = ('-', '.', '', 'NA', 'N/A', 'null')"></xsl:when>
      <xsl:otherwise><xsl:value-of select="$value"/></xsl:otherwise>
    </xsl:choose>
  </xsl:function>

</xsl:stylesheet>
