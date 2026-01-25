<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  PanelApp -> Rare/Orphan Diseases Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Genomics England PanelApp)
  Target: ../schema.json (3.5 Rare/Orphan Diseases Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ panel_id                    │ panel_id                │ Direct       │
  │ panel_name                  │ panel_name              │ Direct       │
  │ panel_name                  │ disease_name            │ Direct       │
  │ panel_id                    │ disease_id              │ Prefix       │
  │ gene_symbol                 │ gene_symbols            │ Array        │
  │ hgnc_id                     │ hgnc_id                 │ Direct       │
  │ evidence_level              │ evidence_level          │ Direct       │
  │ confidence_level            │ confidence_color        │ Enum         │
  │ mode_of_inheritance         │ inheritance_patterns    │ Array        │
  │ expert_reviews              │ expert_reviews          │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NOTES:
  - PanelApp provides curated gene panels for clinical diagnostics
  - Evidence levels: 0 (Gray), 1 (Red), 2 (Amber), 3 (Green)
  - Green genes are diagnostic grade with high evidence
  - Used by NHS England 100,000 Genomes Project
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
      <!-- ========== IDENTIFIERS ========== -->

      <fn:string key="disease_id">
        <xsl:value-of select="concat('PanelApp:', fn:number[@key='panel_id'])"/>
      </fn:string>

      <xsl:choose>
        <xsl:when test="fn:number[@key='panel_id']">
          <fn:number key="panel_id">
            <xsl:value-of select="fn:number[@key='panel_id']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="panel_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PANEL/DISEASE NAME ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='panel_name'])">
          <fn:null key="panel_name"/>
          <fn:null key="disease_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="panel_name">
            <xsl:value-of select="fn:string[@key='panel_name']"/>
          </fn:string>
          <fn:string key="disease_name">
            <xsl:value-of select="fn:string[@key='panel_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GENES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='gene_symbol']/fn:string">
          <fn:array key="gene_symbols">
            <xsl:for-each select="fn:array[@key='gene_symbol']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='gene_symbol'] and not(local:is-null(fn:string[@key='gene_symbol']))">
          <fn:array key="gene_symbols">
            <fn:string><xsl:value-of select="fn:string[@key='gene_symbol']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_symbols"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== HGNC ID ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='hgnc_id'])">
          <fn:null key="hgnc_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="hgnc_id">
            <xsl:value-of select="fn:string[@key='hgnc_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== EVIDENCE LEVEL ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='evidence_level']">
          <fn:number key="evidence_level">
            <xsl:value-of select="fn:number[@key='evidence_level']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="evidence_level"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CONFIDENCE COLOR ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='confidence_level'])">
          <fn:null key="confidence_color"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="confidence_color">
            <xsl:value-of select="fn:string[@key='confidence_level']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== INHERITANCE ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='mode_of_inheritance']/fn:string">
          <fn:array key="inheritance_patterns">
            <xsl:for-each select="fn:array[@key='mode_of_inheritance']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='mode_of_inheritance'] and not(local:is-null(fn:string[@key='mode_of_inheritance']))">
          <fn:array key="inheritance_patterns">
            <fn:string><xsl:value-of select="fn:string[@key='mode_of_inheritance']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="inheritance_patterns"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== EXPERT REVIEWS ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='expert_reviews']">
          <fn:number key="expert_reviews">
            <xsl:value-of select="fn:number[@key='expert_reviews']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="expert_reviews"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">PanelApp</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:number[@key='panel_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
