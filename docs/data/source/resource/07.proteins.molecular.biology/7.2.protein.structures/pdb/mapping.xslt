<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  PDB -> Unified Protein Structures Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Structures Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────┬─────────────────┬──────────────┐
  │ Source Field            │ Target Field    │ Transform    │
  ├─────────────────────────┼─────────────────┼──────────────┤
  │ entry_id                │ structure_id    │ Direct       │
  │ title                   │ title           │ Direct       │
  │ experimental_method     │ method          │ Direct       │
  │ resolution              │ resolution      │ Direct       │
  │ entities                │ molecules       │ Restructure  │
  │ validation              │ quality         │ Restructure  │
  │ (computed)              │ _source         │ Metadata     │
  └─────────────────────────┴─────────────────┴──────────────┘

  NOTES:
  - PDB is the primary experimental structure database
  - Maps experimental data to unified structure schema
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    xmlns:map="http://www.w3.org/2005/xpath-functions/map"
    xmlns:array="http://www.w3.org/2005/xpath-functions/array"
    exclude-result-prefixes="xs fn map array">

  <xsl:output method="json" indent="yes"/>

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <xsl:template match="fn:map">
    <fn:map>
      <!-- Structure ID -->
      <fn:string key="structure_id">
        <xsl:value-of select="fn:string[@key='entry_id']"/>
      </fn:string>

      <!-- Title -->
      <fn:string key="title">
        <xsl:value-of select="fn:string[@key='title']"/>
      </fn:string>

      <!-- Method -->
      <fn:string key="method">
        <xsl:value-of select="fn:string[@key='experimental_method']"/>
      </fn:string>

      <!-- Structure type: experimental -->
      <fn:string key="structure_type">experimental</fn:string>

      <!-- Resolution -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='resolution']">
          <fn:number key="resolution">
            <xsl:value-of select="fn:number[@key='resolution']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="resolution"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Deposition date -->
      <xsl:if test="fn:string[@key='deposition_date']">
        <fn:string key="deposition_date">
          <xsl:value-of select="fn:string[@key='deposition_date']"/>
        </fn:string>
      </xsl:if>

      <!-- Entities/molecules -->
      <xsl:if test="fn:array[@key='entities']">
        <fn:array key="molecules">
          <xsl:for-each select="fn:array[@key='entities']/fn:map[fn:string[@key='type'] = 'polymer']">
            <fn:map>
              <fn:string key="description">
                <xsl:value-of select="fn:string[@key='description']"/>
              </fn:string>
              <xsl:if test="fn:string[@key='source_organism']">
                <fn:string key="organism">
                  <xsl:value-of select="fn:string[@key='source_organism']"/>
                </fn:string>
              </xsl:if>
              <xsl:if test="fn:string[@key='sequence']">
                <fn:string key="sequence">
                  <xsl:value-of select="fn:string[@key='sequence']"/>
                </fn:string>
              </xsl:if>
            </fn:map>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- Ligands -->
      <xsl:if test="fn:array[@key='ligands']">
        <fn:array key="ligands">
          <xsl:for-each select="fn:array[@key='ligands']/fn:map">
            <fn:map>
              <fn:string key="id">
                <xsl:value-of select="fn:string[@key='comp_id']"/>
              </fn:string>
              <fn:string key="name">
                <xsl:value-of select="fn:string[@key='name']"/>
              </fn:string>
            </fn:map>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- Quality metrics -->
      <xsl:if test="fn:map[@key='validation'] or fn:map[@key='refinement']">
        <fn:map key="quality">
          <xsl:if test="fn:map[@key='refinement']/fn:number[@key='r_work']">
            <fn:number key="r_work">
              <xsl:value-of select="fn:map[@key='refinement']/fn:number[@key='r_work']"/>
            </fn:number>
          </xsl:if>
          <xsl:if test="fn:map[@key='refinement']/fn:number[@key='r_free']">
            <fn:number key="r_free">
              <xsl:value-of select="fn:map[@key='refinement']/fn:number[@key='r_free']"/>
            </fn:number>
          </xsl:if>
          <xsl:if test="fn:map[@key='validation']/fn:number[@key='clashscore']">
            <fn:number key="clashscore">
              <xsl:value-of select="fn:map[@key='validation']/fn:number[@key='clashscore']"/>
            </fn:number>
          </xsl:if>
        </fn:map>
      </xsl:if>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">PDB</fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='entry_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
