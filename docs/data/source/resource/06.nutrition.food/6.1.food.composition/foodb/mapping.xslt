<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  FooDB -> Unified Food Composition Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Food Composition Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────┬─────────────────┬──────────────┐
  │ Source Field            │ Target Field    │ Transform    │
  ├─────────────────────────┼─────────────────┼──────────────┤
  │ public_id               │ food_id         │ Direct       │
  │ name                    │ food_name       │ Direct       │
  │ name_scientific         │ scientific_name │ Direct       │
  │ food_group              │ food_category[] │ Array wrap   │
  │ compounds[]             │ compound_content│ Restructure  │
  │ moldb_smiles            │ moldb_smiles    │ Direct       │
  │ moldb_inchikey          │ moldb_inchikey  │ Direct       │
  │ cas_number              │ cas_number      │ Direct       │
  │ itis_id                 │ itis_id         │ Direct       │
  │ (computed)              │ _source         │ Metadata     │
  └─────────────────────────┴─────────────────┴──────────────┘

  NULL HANDLING:
  - name_scientific: null if missing
  - itis_id: null if missing

  NOTES:
  - FooDB uses numeric food_id internally, public_id for external use
  - compound_content restructured from compounds array
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

      <!-- Primary identifier: use public_id -->
      <fn:string key="food_id">
        <xsl:value-of select="fn:string[@key='public_id']"/>
      </fn:string>

      <!-- Food name -->
      <fn:string key="food_name">
        <xsl:value-of select="fn:string[@key='name']"/>
      </fn:string>

      <!-- Scientific name (nullable) -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='name_scientific'] and fn:string[@key='name_scientific'] != ''">
          <fn:string key="scientific_name">
            <xsl:value-of select="fn:string[@key='name_scientific']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="scientific_name"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ARRAY TRANSFORMATIONS ========== -->

      <!-- Food category: wrap food_group in array -->
      <fn:array key="food_category">
        <xsl:if test="fn:string[@key='food_group']">
          <fn:string><xsl:value-of select="fn:string[@key='food_group']"/></fn:string>
        </xsl:if>
        <xsl:if test="fn:string[@key='food_subgroup']">
          <fn:string><xsl:value-of select="fn:string[@key='food_subgroup']"/></fn:string>
        </xsl:if>
      </fn:array>

      <!-- Brand: null for FooDB (not a branded product database) -->
      <fn:null key="brand"/>

      <!-- Compound content: restructure from compounds array -->
      <xsl:if test="fn:array[@key='compounds']">
        <fn:array key="compound_content">
          <xsl:for-each select="fn:array[@key='compounds']/fn:map">
            <fn:map>
              <fn:string key="compound_name">
                <xsl:value-of select="fn:string[@key='name']"/>
              </fn:string>
              <fn:string key="compound_id">
                <xsl:value-of select="fn:string[@key='compound_id']"/>
              </fn:string>
              <xsl:choose>
                <xsl:when test="fn:number[@key='content']">
                  <fn:number key="content">
                    <xsl:value-of select="fn:number[@key='content']"/>
                  </fn:number>
                </xsl:when>
                <xsl:otherwise>
                  <fn:null key="content"/>
                </xsl:otherwise>
              </xsl:choose>
              <fn:string key="unit">
                <xsl:value-of select="fn:string[@key='content_unit']"/>
              </fn:string>
            </fn:map>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- ========== CHEMICAL IDENTIFIERS ========== -->

      <!-- First compound's SMILES (if available) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='compounds']/fn:map[1]/fn:string[@key='moldb_smiles']">
          <fn:string key="moldb_smiles">
            <xsl:value-of select="fn:array[@key='compounds']/fn:map[1]/fn:string[@key='moldb_smiles']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="moldb_smiles"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ITIS ID -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='itis_id']">
          <fn:number key="itis_id">
            <xsl:value-of select="fn:number[@key='itis_id']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="itis_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SCORES (not applicable for FooDB) ========== -->
      <fn:null key="nutriscore_grade"/>
      <fn:null key="nova_group"/>
      <fn:null key="ecoscore_grade"/>
      <fn:null key="quantity"/>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">FooDB</fn:string>
        <fn:string key="version">1.0</fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='public_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
