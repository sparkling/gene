<?xml version="1.0" encoding="UTF-8"?>
<!--
  NAPRALERT -> Unified Western & Global Herbal Schema Mapping
  Source: ./schema.json
  Target: ../schema.json
  XSLT Version: 3.0
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
      <fn:string key="organism_id"><xsl:value-of select="fn:string[@key='organism_id']"/></fn:string>
      <fn:string key="scientific_name"><xsl:value-of select="fn:string[@key='scientific_name']"/></fn:string>
      <fn:string key="family"><xsl:value-of select="fn:string[@key='family']"/></fn:string>
      <fn:string key="kingdom"><xsl:value-of select="fn:string[@key='kingdom']"/></fn:string>
      <xsl:if test="fn:string[@key='compound_id']">
        <fn:string key="compound_id"><xsl:value-of select="fn:string[@key='compound_id']"/></fn:string>
        <fn:string key="compound_name"><xsl:value-of select="fn:string[@key='compound_name']"/></fn:string>
      </xsl:if>
      <xsl:if test="fn:string[@key='cas_number']">
        <fn:string key="cas_number"><xsl:value-of select="fn:string[@key='cas_number']"/></fn:string>
      </xsl:if>
      <xsl:if test="fn:map[@key='ethnobotany']">
        <fn:map key="ethnobotany">
          <xsl:copy-of select="fn:map[@key='ethnobotany']/*"/>
        </fn:map>
      </xsl:if>
      <xsl:if test="fn:map[@key='pharmacology']">
        <fn:map key="pharmacology">
          <xsl:copy-of select="fn:map[@key='pharmacology']/*"/>
        </fn:map>
      </xsl:if>
      <xsl:if test="fn:map[@key='literature']/fn:number[@key='pubmed_id']">
        <fn:number key="pubmed_id"><xsl:value-of select="fn:map[@key='literature']/fn:number[@key='pubmed_id']"/></fn:number>
      </xsl:if>
      <fn:map key="_source">
        <fn:string key="database">NAPRALERT</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:string[@key='organism_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>
</xsl:stylesheet>
