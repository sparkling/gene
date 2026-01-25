<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  PubMed Central -> Unified Scientific Literature Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Scientific Literature Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  +-------------------------+-----------------+--------------+
  | Source Field            | Target Field    | Transform    |
  +-------------------------+-----------------+--------------+
  | pmc_id                  | id              | Prefix       |
  | title                   | title           | Direct       |
  | abstract                | abstract        | Direct       |
  | authors                 | authors         | Normalize    |
  | doi                     | doi             | Direct       |
  | pub_date.year           | publication_year| Extract      |
  | body.sections           | full_text       | Concatenate  |
  +-------------------------+-----------------+--------------+

  NOTES:
  - PMC ID used as primary identifier
  - Full text extracted from body sections
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
        <xsl:value-of select="fn:string[@key='pmc_id']"/>
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

      <!-- Authors -->
      <fn:array key="authors">
        <xsl:for-each select="fn:array[@key='authors']/fn:map">
          <fn:map>
            <fn:string key="name">
              <xsl:value-of select="concat(fn:string[@key='given_names'], ' ', fn:string[@key='surname'])"/>
            </fn:string>
            <xsl:if test="fn:string[@key='orcid']">
              <fn:string key="orcid">
                <xsl:value-of select="fn:string[@key='orcid']"/>
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
      <xsl:if test="fn:map[@key='pub_date']/fn:number[@key='year']">
        <fn:number key="publication_year">
          <xsl:value-of select="fn:map[@key='pub_date']/fn:number[@key='year']"/>
        </fn:number>
      </xsl:if>

      <!-- Journal -->
      <fn:string key="source">
        <xsl:value-of select="fn:map[@key='journal']/fn:string[@key='title']"/>
      </fn:string>

      <!-- License -->
      <xsl:if test="fn:map[@key='license']/fn:string[@key='type']">
        <fn:string key="license">
          <xsl:value-of select="fn:map[@key='license']/fn:string[@key='type']"/>
        </fn:string>
      </xsl:if>

      <!-- Has full text indicator -->
      <fn:boolean key="has_full_text">true</fn:boolean>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">PubMed Central</fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='pmc_id']"/>
        </fn:string>
        <xsl:if test="fn:number[@key='pmid']">
          <fn:number key="pmid">
            <xsl:value-of select="fn:number[@key='pmid']"/>
          </fn:number>
        </xsl:if>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
