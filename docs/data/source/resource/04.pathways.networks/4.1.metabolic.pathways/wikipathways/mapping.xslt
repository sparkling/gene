<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  WikiPathways -> Unified Metabolic Pathways Schema Mapping
  ============================================================
  Source: ./schema.json (WikiPathways GPML/API data)
  Target: ../schema.json (Unified 4.1 Metabolic Pathways Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ wpid                        │ pathway_id              │ Direct       │
  │ Name                        │ pathway_name            │ Direct       │
  │ Organism                    │ organism                │ Direct       │
  │ DataNode                    │ gene_entries            │ Restructure  │
  │ Metabolite                  │ compound_entries        │ Restructure  │
  │ Interaction                 │ reactions               │ Restructure  │
  │ Xref                        │ cross_references        │ Restructure  │
  │ GraphId                     │ wikipathways_graph_id   │ Direct       │
  │ DataNode_Type               │ wikipathways_datanode_type │ Direct    │
  │ TextLabel                   │ wikipathways_text_label │ Direct       │
  │ ArrowHead                   │ wikipathways_arrow_head │ Direct       │
  │ Group_Style                 │ wikipathways_group_style│ Direct       │
  │ State                       │ wikipathways_state      │ Array        │
  │ Shape                       │ wikipathways_shape      │ Object       │
  │ Graphics                    │ wikipathways_graphics   │ Object       │
  │ Xref_Database               │ wikipathways_xref_database│ Direct     │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings "" -> null
  - Missing optional fields -> null

  NOTES:
  - WikiPathways uses WP prefix for pathway IDs (e.g., WP554)
  - GPML format uses DataNode for all node types
  - MIM notation used for arrow heads (interaction types)
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

      <!-- Primary identifier (WP554) -->
      <fn:string key="pathway_id">
        <xsl:value-of select="fn:string[@key='wpid']"/>
      </fn:string>

      <!-- Pathway name -->
      <fn:string key="pathway_name">
        <xsl:value-of select="fn:string[@key='Name']"/>
      </fn:string>

      <!-- Organism -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='Organism'] and fn:string[@key='Organism'] != ''">
          <fn:string key="organism">
            <xsl:value-of select="fn:string[@key='Organism']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="organism"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ARRAY TRANSFORMATIONS ========== -->

      <!-- Gene entries from DataNode -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='DataNode']">
          <fn:array key="gene_entries">
            <xsl:for-each select="fn:array[@key='DataNode']/fn:map[fn:string[@key='Type'] = ('GeneProduct', 'Protein', 'Rna')]">
              <fn:map>
                <fn:string key="id">
                  <xsl:value-of select="fn:string[@key='GraphId']"/>
                </fn:string>
                <fn:string key="symbol">
                  <xsl:value-of select="fn:string[@key='TextLabel']"/>
                </fn:string>
                <xsl:if test="fn:string[@key='Type']">
                  <fn:string key="type">
                    <xsl:value-of select="fn:string[@key='Type']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_entries"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Compound entries from Metabolite nodes -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='Metabolite'] or fn:array[@key='DataNode']/fn:map[fn:string[@key='Type'] = 'Metabolite']">
          <fn:array key="compound_entries">
            <xsl:for-each select="fn:array[@key='Metabolite']/fn:map | fn:array[@key='DataNode']/fn:map[fn:string[@key='Type'] = 'Metabolite']">
              <fn:map>
                <fn:string key="id">
                  <xsl:value-of select="fn:string[@key='GraphId']"/>
                </fn:string>
                <fn:string key="name">
                  <xsl:value-of select="fn:string[@key='TextLabel']"/>
                </fn:string>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="compound_entries"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Reactions from Interaction -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='Interaction']">
          <fn:array key="reactions">
            <xsl:for-each select="fn:array[@key='Interaction']/fn:map">
              <fn:map>
                <fn:string key="id">
                  <xsl:value-of select="fn:string[@key='GraphId']"/>
                </fn:string>
                <xsl:if test="fn:string[@key='ArrowHead']">
                  <fn:string key="type">
                    <xsl:value-of select="fn:string[@key='ArrowHead']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="reactions"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Cross-references -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='Xref']">
          <fn:array key="cross_references">
            <xsl:for-each select="fn:array[@key='Xref']/fn:map">
              <fn:map>
                <fn:string key="database">
                  <xsl:value-of select="fn:string[@key='Database']"/>
                </fn:string>
                <fn:string key="id">
                  <xsl:value-of select="fn:string[@key='ID']"/>
                </fn:string>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cross_references"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== WIKIPATHWAYS-SPECIFIC FIELDS ========== -->

      <!-- Graph ID -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='GraphId'] and fn:string[@key='GraphId'] != ''">
          <fn:string key="wikipathways_graph_id">
            <xsl:value-of select="fn:string[@key='GraphId']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="wikipathways_graph_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- DataNode type -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='DataNode_Type'] and fn:string[@key='DataNode_Type'] != ''">
          <fn:string key="wikipathways_datanode_type">
            <xsl:value-of select="fn:string[@key='DataNode_Type']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="wikipathways_datanode_type"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Text label -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='TextLabel'] and fn:string[@key='TextLabel'] != ''">
          <fn:string key="wikipathways_text_label">
            <xsl:value-of select="fn:string[@key='TextLabel']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="wikipathways_text_label"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Arrow head (MIM notation) -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='ArrowHead'] and fn:string[@key='ArrowHead'] != ''">
          <fn:string key="wikipathways_arrow_head">
            <xsl:value-of select="fn:string[@key='ArrowHead']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="wikipathways_arrow_head"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Group style -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='Group_Style'] and fn:string[@key='Group_Style'] != ''">
          <fn:string key="wikipathways_group_style">
            <xsl:value-of select="fn:string[@key='Group_Style']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="wikipathways_group_style"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- State (PTM markers) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='State']">
          <fn:array key="wikipathways_state">
            <xsl:for-each select="fn:array[@key='State']/fn:map">
              <fn:map>
                <xsl:copy-of select="./*"/>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="wikipathways_state"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Shape -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='Shape']">
          <fn:map key="wikipathways_shape">
            <xsl:copy-of select="fn:map[@key='Shape']/*"/>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="wikipathways_shape"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Graphics -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='Graphics']">
          <fn:map key="wikipathways_graphics">
            <xsl:copy-of select="fn:map[@key='Graphics']/*"/>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="wikipathways_graphics"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Xref database -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='Xref_Database'] and fn:string[@key='Xref_Database'] != ''">
          <fn:string key="wikipathways_xref_database">
            <xsl:value-of select="fn:string[@key='Xref_Database']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="wikipathways_xref_database"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">WikiPathways</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='wpid']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
