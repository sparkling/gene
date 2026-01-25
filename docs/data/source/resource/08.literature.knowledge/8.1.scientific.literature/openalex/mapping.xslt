<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  OpenAlex -> Unified Scientific Literature Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Scientific Literature Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  +-------------------------+-----------------+--------------+
  | Source Field            | Target Field    | Transform    |
  +-------------------------+-----------------+--------------+
  | id                      | id              | Extract ID   |
  | display_name            | title           | Direct       |
  | authorships             | authors         | Normalize    |
  | doi                     | doi             | Clean        |
  | publication_year        | publication_year| Direct       |
  | cited_by_count          | citations       | Direct       |
  | open_access             | open_access     | Extract      |
  +-------------------------+-----------------+--------------+

  NOTES:
  - OpenAlex IDs are full URLs, extract W-prefixed ID
  - DOIs include https prefix, normalize to 10.xxx format
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
      <!-- Primary identifier - extract W-number -->
      <fn:string key="id">
        <xsl:value-of select="concat('OpenAlex:', replace(fn:string[@key='id'], 'https://openalex.org/', ''))"/>
      </fn:string>

      <!-- Title -->
      <fn:string key="title">
        <xsl:value-of select="fn:string[@key='display_name']"/>
      </fn:string>

      <!-- Authors -->
      <fn:array key="authors">
        <xsl:for-each select="fn:array[@key='authorships']/fn:map">
          <fn:map>
            <fn:string key="name">
              <xsl:value-of select="fn:map[@key='author']/fn:string[@key='display_name']"/>
            </fn:string>
            <xsl:if test="fn:map[@key='author']/fn:string[@key='orcid']">
              <fn:string key="orcid">
                <xsl:value-of select="replace(fn:map[@key='author']/fn:string[@key='orcid'], 'https://orcid.org/', '')"/>
              </fn:string>
            </xsl:if>
            <fn:string key="position">
              <xsl:value-of select="fn:string[@key='author_position']"/>
            </fn:string>
          </fn:map>
        </xsl:for-each>
      </fn:array>

      <!-- DOI - normalize format -->
      <xsl:if test="fn:string[@key='doi']">
        <fn:string key="doi">
          <xsl:value-of select="replace(fn:string[@key='doi'], 'https://doi.org/', '')"/>
        </fn:string>
      </xsl:if>

      <!-- Publication year -->
      <xsl:if test="fn:number[@key='publication_year']">
        <fn:number key="publication_year">
          <xsl:value-of select="fn:number[@key='publication_year']"/>
        </fn:number>
      </xsl:if>

      <!-- Work type -->
      <fn:string key="type">
        <xsl:value-of select="fn:string[@key='type']"/>
      </fn:string>

      <!-- Journal/Source -->
      <xsl:if test="fn:map[@key='primary_location']/fn:map[@key='source']/fn:string[@key='display_name']">
        <fn:string key="source">
          <xsl:value-of select="fn:map[@key='primary_location']/fn:map[@key='source']/fn:string[@key='display_name']"/>
        </fn:string>
      </xsl:if>

      <!-- Citation count -->
      <xsl:if test="fn:number[@key='cited_by_count']">
        <fn:number key="citation_count">
          <xsl:value-of select="fn:number[@key='cited_by_count']"/>
        </fn:number>
      </xsl:if>

      <!-- Open access -->
      <fn:boolean key="is_open_access">
        <xsl:value-of select="fn:map[@key='open_access']/fn:boolean[@key='is_oa'] = 'true'"/>
      </fn:boolean>

      <!-- Concepts as subjects -->
      <fn:array key="subjects">
        <xsl:for-each select="fn:array[@key='concepts']/fn:map[fn:number[@key='score'] &gt; 0.5]">
          <fn:string>
            <xsl:value-of select="fn:string[@key='display_name']"/>
          </fn:string>
        </xsl:for-each>
      </fn:array>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">OpenAlex</fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
