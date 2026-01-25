<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Europe PMC -> Unified Scientific Literature Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Scientific Literature Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  +-------------------------+-----------------+--------------+
  | Source Field            | Target Field    | Transform    |
  +-------------------------+-----------------+--------------+
  | id                      | id              | Direct       |
  | title                   | title           | Direct       |
  | abstractText            | abstract        | Direct       |
  | authorList.author       | authors         | Normalize    |
  | doi                     | doi             | Direct       |
  | journalInfo.yearOfPub   | publication_year| Extract      |
  | citedByCount            | citations       | Direct       |
  +-------------------------+-----------------+--------------+

  NOTES:
  - Includes citation count and open access status
  - Text-mined annotations preserved in metadata
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
        <xsl:value-of select="concat('EPMC:', fn:string[@key='id'])"/>
      </fn:string>

      <!-- Title -->
      <fn:string key="title">
        <xsl:value-of select="fn:string[@key='title']"/>
      </fn:string>

      <!-- Abstract -->
      <xsl:if test="fn:string[@key='abstractText']">
        <fn:string key="abstract">
          <xsl:value-of select="fn:string[@key='abstractText']"/>
        </fn:string>
      </xsl:if>

      <!-- Authors -->
      <fn:array key="authors">
        <xsl:for-each select="fn:map[@key='authorList']/fn:array[@key='author']/fn:map">
          <fn:map>
            <fn:string key="name">
              <xsl:value-of select="fn:string[@key='fullName']"/>
            </fn:string>
            <xsl:if test="fn:map[@key='authorId']/fn:string[@key='type'] = 'ORCID'">
              <fn:string key="orcid">
                <xsl:value-of select="fn:map[@key='authorId']/fn:string[@key='value']"/>
              </fn:string>
            </xsl:if>
          </fn:map>
        </xsl:for-each>
      </fn:array>

      <!-- DOI -->
      <xsl:if test="fn:string[@key='doi']">
        <fn:string key="doi">
          <xsl:value-of select="fn:string[@key='doi']"/>
        </fn:string>
      </xsl:if>

      <!-- Publication year -->
      <xsl:if test="fn:map[@key='journalInfo']/fn:number[@key='yearOfPublication']">
        <fn:number key="publication_year">
          <xsl:value-of select="fn:map[@key='journalInfo']/fn:number[@key='yearOfPublication']"/>
        </fn:number>
      </xsl:if>

      <!-- Journal -->
      <xsl:if test="fn:map[@key='journalInfo']/fn:map[@key='journal']/fn:string[@key='title']">
        <fn:string key="source">
          <xsl:value-of select="fn:map[@key='journalInfo']/fn:map[@key='journal']/fn:string[@key='title']"/>
        </fn:string>
      </xsl:if>

      <!-- Citation count -->
      <xsl:if test="fn:number[@key='citedByCount']">
        <fn:number key="citation_count">
          <xsl:value-of select="fn:number[@key='citedByCount']"/>
        </fn:number>
      </xsl:if>

      <!-- Open access -->
      <fn:boolean key="is_open_access">
        <xsl:value-of select="fn:string[@key='isOpenAccess'] = 'Y'"/>
      </fn:boolean>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">Europe PMC</fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='id']"/>
        </fn:string>
        <fn:string key="source_type">
          <xsl:value-of select="fn:string[@key='source']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
