<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Semantic Scholar -> Unified Scientific Literature Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Scientific Literature Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  +-------------------------+-----------------+--------------+
  | Source Field            | Target Field    | Transform    |
  +-------------------------+-----------------+--------------+
  | paperId                 | id              | Prefix       |
  | title                   | title           | Direct       |
  | abstract                | abstract        | Direct       |
  | authors                 | authors         | Normalize    |
  | externalIds.DOI         | doi             | Direct       |
  | year                    | publication_year| Direct       |
  | citationCount           | citations       | Direct       |
  | tldr.text               | summary         | Direct       |
  +-------------------------+-----------------+--------------+

  NOTES:
  - Includes AI-generated TLDR summary
  - Influential citation count tracked
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
        <xsl:value-of select="concat('S2:', fn:string[@key='paperId'])"/>
      </fn:string>

      <!-- Title -->
      <fn:string key="title">
        <xsl:value-of select="fn:string[@key='title']"/>
      </fn:string>

      <!-- Abstract -->
      <xsl:if test="fn:string[@key='abstract']">
        <fn:string key="abstract">
          <xsl:value-of select="fn:string[@key='abstract']"/>
        </fn:string>
      </xsl:if>

      <!-- AI Summary (TLDR) -->
      <xsl:if test="fn:map[@key='tldr']/fn:string[@key='text']">
        <fn:string key="summary">
          <xsl:value-of select="fn:map[@key='tldr']/fn:string[@key='text']"/>
        </fn:string>
      </xsl:if>

      <!-- Authors -->
      <fn:array key="authors">
        <xsl:for-each select="fn:array[@key='authors']/fn:map">
          <fn:map>
            <fn:string key="name">
              <xsl:value-of select="fn:string[@key='name']"/>
            </fn:string>
            <fn:string key="s2_id">
              <xsl:value-of select="fn:string[@key='authorId']"/>
            </fn:string>
          </fn:map>
        </xsl:for-each>
      </fn:array>

      <!-- DOI from external IDs -->
      <xsl:if test="fn:map[@key='externalIds']/fn:string[@key='DOI']">
        <fn:string key="doi">
          <xsl:value-of select="fn:map[@key='externalIds']/fn:string[@key='DOI']"/>
        </fn:string>
      </xsl:if>

      <!-- Publication year -->
      <xsl:if test="fn:number[@key='year']">
        <fn:number key="publication_year">
          <xsl:value-of select="fn:number[@key='year']"/>
        </fn:number>
      </xsl:if>

      <!-- Venue -->
      <xsl:if test="fn:string[@key='venue']">
        <fn:string key="source">
          <xsl:value-of select="fn:string[@key='venue']"/>
        </fn:string>
      </xsl:if>

      <!-- Citation count -->
      <xsl:if test="fn:number[@key='citationCount']">
        <fn:number key="citation_count">
          <xsl:value-of select="fn:number[@key='citationCount']"/>
        </fn:number>
      </xsl:if>

      <!-- Influential citations -->
      <xsl:if test="fn:number[@key='influentialCitationCount']">
        <fn:number key="influential_citations">
          <xsl:value-of select="fn:number[@key='influentialCitationCount']"/>
        </fn:number>
      </xsl:if>

      <!-- Open access -->
      <fn:boolean key="is_open_access">
        <xsl:value-of select="fn:boolean[@key='isOpenAccess'] = 'true'"/>
      </fn:boolean>

      <!-- Fields of study -->
      <fn:array key="subjects">
        <xsl:for-each select="fn:array[@key='fieldsOfStudy']/fn:string">
          <fn:string>
            <xsl:value-of select="."/>
          </fn:string>
        </xsl:for-each>
      </fn:array>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">Semantic Scholar</fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='paperId']"/>
        </fn:string>
        <xsl:if test="fn:number[@key='corpusId']">
          <fn:number key="corpus_id">
            <xsl:value-of select="fn:number[@key='corpusId']"/>
          </fn:number>
        </xsl:if>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
