<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  MeSH -> Disease Ontologies Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Medical Subject Headings)
  Target: ../schema.json (3.1 Disease Ontologies Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ UI                          │ mesh_ui             │ Direct       │
  │ UI                          │ id                  │ Direct       │
  │ MH                          │ name                │ Direct       │
  │ MS                          │ definition          │ Direct       │
  │ SY                          │ synonyms            │ Array        │
  │ MN                          │ tree_number         │ Direct       │
  │ ENTRY                       │ entry_term          │ Array        │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null

  NOTES:
  - MeSH uses D###### format for disease descriptors
  - Tree numbers indicate hierarchy (e.g., C01.539)
  - Used for PubMed indexing
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

      <fn:string key="id">
        <xsl:value-of select="fn:string[@key='UI']"/>
      </fn:string>

      <fn:string key="mesh_ui">
        <xsl:value-of select="fn:string[@key='UI']"/>
      </fn:string>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='MN'])">
          <fn:null key="tree_number"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="tree_number">
            <xsl:value-of select="fn:string[@key='MN']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='MH'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='MH']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='MS'])">
          <fn:null key="definition"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="definition">
            <xsl:value-of select="fn:string[@key='MS']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SYNONYMS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='SY']/fn:string">
          <fn:array key="synonyms">
            <xsl:for-each select="fn:array[@key='SY']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="synonyms"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ENTRY TERMS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='ENTRY']/fn:string">
          <fn:array key="entry_term">
            <xsl:for-each select="fn:array[@key='ENTRY']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="entry_term"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">MeSH</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='UI']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
