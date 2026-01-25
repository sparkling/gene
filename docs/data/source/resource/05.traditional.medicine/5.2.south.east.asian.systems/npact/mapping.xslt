<?xml version="1.0" encoding="UTF-8"?>
<!--
  NPACT -> Unified South/East Asian Systems Schema Mapping
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
      <fn:string key="npact_id"><xsl:value-of select="fn:string[@key='npact_id']"/></fn:string>
      <fn:string key="compound_name"><xsl:value-of select="fn:string[@key='compound_name']"/></fn:string>
      <xsl:if test="fn:string[@key='plant_id']">
        <fn:string key="plant_id"><xsl:value-of select="fn:string[@key='plant_id']"/></fn:string>
      </xsl:if>
      <fn:string key="botanical_name"><xsl:value-of select="fn:string[@key='scientific_name']"/></fn:string>
      <fn:string key="family"><xsl:value-of select="fn:string[@key='family']"/></fn:string>
      <xsl:if test="fn:number[@key='pubchem_cid']">
        <fn:number key="pubchem_cid"><xsl:value-of select="fn:number[@key='pubchem_cid']"/></fn:number>
      </xsl:if>
      <xsl:if test="fn:string[@key='activity_type']">
        <fn:string key="npact_activity_type"><xsl:value-of select="fn:string[@key='activity_type']"/></fn:string>
        <fn:number key="npact_activity_value"><xsl:value-of select="fn:number[@key='activity_value']"/></fn:number>
        <fn:string key="npact_activity_unit"><xsl:value-of select="fn:string[@key='activity_unit']"/></fn:string>
      </xsl:if>
      <xsl:if test="fn:string[@key='cell_line_name']">
        <fn:string key="npact_cell_line_name"><xsl:value-of select="fn:string[@key='cell_line_name']"/></fn:string>
        <fn:string key="npact_cancer_type"><xsl:value-of select="fn:string[@key='cancer_type']"/></fn:string>
      </xsl:if>
      <xsl:if test="fn:number[@key='pubmed_id']">
        <fn:number key="npact_pubmed_id"><xsl:value-of select="fn:number[@key='pubmed_id']"/></fn:number>
      </xsl:if>
      <fn:map key="_source">
        <fn:string key="database">NPACT</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:string[@key='npact_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>
</xsl:stylesheet>
