<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  MONDO -> Disease Ontologies Unified Schema Mapping
  ============================================================
  Source: ./schema.md (MONDO Disease Ontology)
  Target: ../schema.json (3.1 Disease Ontologies Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ id                          │ id                  │ Direct       │
  │ label                       │ name                │ Direct       │
  │ definition                  │ definition          │ Direct       │
  │ synonyms                    │ synonyms            │ Array        │
  │ is_a                        │ is_a                │ Array        │
  │ xref                        │ xref                │ Array        │
  │ exactMatch                  │ exact_match         │ Array        │
  │ closeMatch                  │ close_match         │ Array        │
  │ narrowMatch                 │ narrow_match        │ Array        │
  │ broadMatch                  │ broad_match         │ Array        │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null

  NOTES:
  - MONDO uses MONDO:####### format for IDs
  - Provides semantic mappings to other ontologies
  - Integrates OMIM, Orphanet, DOID, and others
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
        <xsl:value-of select="fn:string[@key='id']"/>
      </fn:string>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='label'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='label']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='definition'])">
          <fn:null key="definition"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="definition">
            <xsl:value-of select="fn:string[@key='definition']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SYNONYMS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='synonyms']/fn:string">
          <fn:array key="synonyms">
            <xsl:for-each select="fn:array[@key='synonyms']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="synonyms"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ONTOLOGY RELATIONSHIPS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='is_a']/fn:string">
          <fn:array key="is_a">
            <xsl:for-each select="fn:array[@key='is_a']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="is_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CROSS-REFERENCES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='xref']/fn:string">
          <fn:array key="xref">
            <xsl:for-each select="fn:array[@key='xref']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="xref"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SEMANTIC MAPPINGS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='exactMatch']/fn:string">
          <fn:array key="exact_match">
            <xsl:for-each select="fn:array[@key='exactMatch']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="exact_match"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:array[@key='closeMatch']/fn:string">
          <fn:array key="close_match">
            <xsl:for-each select="fn:array[@key='closeMatch']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="close_match"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:array[@key='narrowMatch']/fn:string">
          <fn:array key="narrow_match">
            <xsl:for-each select="fn:array[@key='narrowMatch']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="narrow_match"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:array[@key='broadMatch']/fn:string">
          <fn:array key="broad_match">
            <xsl:for-each select="fn:array[@key='broadMatch']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="broad_match"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">MONDO</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
