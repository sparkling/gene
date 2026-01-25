<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  ClinicalTrials.gov -> Unified Regulatory Schema Mapping
  ============================================================
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    exclude-result-prefixes="xs fn">

  <xsl:output method="json" indent="yes"/>

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <xsl:template match="fn:map">
    <fn:map>
      <fn:string key="id">
        <xsl:value-of select="fn:map[@key='protocolSection']/fn:map[@key='identificationModule']/fn:string[@key='nctId']"/>
      </fn:string>

      <fn:string key="title">
        <xsl:value-of select="fn:map[@key='protocolSection']/fn:map[@key='identificationModule']/fn:string[@key='briefTitle']"/>
      </fn:string>

      <fn:string key="status">
        <xsl:value-of select="fn:map[@key='protocolSection']/fn:map[@key='statusModule']/fn:string[@key='overallStatus']"/>
      </fn:string>

      <fn:string key="study_type">
        <xsl:value-of select="fn:map[@key='protocolSection']/fn:map[@key='designModule']/fn:string[@key='studyType']"/>
      </fn:string>

      <fn:array key="conditions">
        <xsl:for-each select="fn:map[@key='protocolSection']/fn:map[@key='conditionsModule']/fn:array[@key='conditions']/fn:string">
          <fn:string><xsl:value-of select="."/></fn:string>
        </xsl:for-each>
      </fn:array>

      <fn:boolean key="has_results">
        <xsl:value-of select="fn:boolean[@key='hasResults'] = 'true'"/>
      </fn:boolean>

      <fn:map key="_source">
        <fn:string key="database">ClinicalTrials.gov</fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>

</xsl:stylesheet>
