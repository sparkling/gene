<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Wikipedia -> Unified Knowledge Base Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Knowledge Base Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  +-------------------------+-----------------+--------------+
  | Source Field            | Target Field    | Transform    |
  +-------------------------+-----------------+--------------+
  | pageid                  | id              | Prefix       |
  | title                   | title           | Direct       |
  | extract                 | description     | Direct       |
  | wikidata_id             | wikidata_id     | Direct       |
  | categories              | categories      | Direct       |
  | infobox                 | properties      | Flatten      |
  +-------------------------+-----------------+--------------+

  NOTES:
  - Infobox data extracted as structured properties
  - References with PMIDs/DOIs preserved
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
      <!-- Primary identifier -->
      <fn:string key="id">
        <xsl:value-of select="concat('WP:', fn:number[@key='pageid'])"/>
      </fn:string>

      <!-- Title -->
      <fn:string key="title">
        <xsl:value-of select="fn:string[@key='title']"/>
      </fn:string>

      <!-- Description from extract -->
      <xsl:if test="fn:string[@key='extract']">
        <fn:string key="description">
          <xsl:value-of select="fn:string[@key='extract']"/>
        </fn:string>
      </xsl:if>

      <!-- Wikidata link -->
      <xsl:if test="fn:string[@key='wikidata_id']">
        <fn:string key="wikidata_id">
          <xsl:value-of select="fn:string[@key='wikidata_id']"/>
        </fn:string>
      </xsl:if>

      <!-- Categories -->
      <fn:array key="categories">
        <xsl:for-each select="fn:array[@key='categories']/fn:string">
          <fn:string>
            <xsl:value-of select="."/>
          </fn:string>
        </xsl:for-each>
      </fn:array>

      <!-- Properties from infobox -->
      <xsl:if test="fn:map[@key='infobox']">
        <fn:map key="properties">
          <xsl:for-each select="fn:map[@key='infobox']/fn:string">
            <fn:string key="{@key}">
              <xsl:value-of select="."/>
            </fn:string>
          </xsl:for-each>
        </fn:map>
      </xsl:if>

      <!-- Cross-references from references -->
      <fn:array key="cross_references">
        <xsl:for-each select="fn:array[@key='references']/fn:map[fn:string[@key='pmid'] or fn:string[@key='doi']]">
          <fn:map>
            <xsl:if test="fn:string[@key='pmid']">
              <fn:string key="pmid">
                <xsl:value-of select="fn:string[@key='pmid']"/>
              </fn:string>
            </xsl:if>
            <xsl:if test="fn:string[@key='doi']">
              <fn:string key="doi">
                <xsl:value-of select="fn:string[@key='doi']"/>
              </fn:string>
            </xsl:if>
          </fn:map>
        </xsl:for-each>
      </fn:array>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">Wikipedia</fn:string>
        <fn:number key="original_id">
          <xsl:value-of select="fn:number[@key='pageid']"/>
        </fn:number>
        <fn:string key="url">
          <xsl:value-of select="concat('https://en.wikipedia.org/wiki/', encode-for-uri(fn:string[@key='title']))"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
