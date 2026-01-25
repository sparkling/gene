<?xml version="1.0" encoding="UTF-8"?>
<!--
  TM-MC -> Unified South/East Asian Systems Schema Mapping
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
      <fn:string key="plant_id"><xsl:value-of select="fn:string[@key='plant_id']"/></fn:string>
      <fn:string key="botanical_name"><xsl:value-of select="fn:string[@key='scientific_name']"/></fn:string>
      <fn:string key="family"><xsl:value-of select="fn:string[@key='family']"/></fn:string>
      <fn:string key="tmmc_system_code"><xsl:value-of select="fn:string[@key='system_code']"/></fn:string>
      <fn:string key="tmmc_local_name"><xsl:value-of select="fn:string[@key='local_name']"/></fn:string>
      <xsl:if test="fn:array[@key='traditional_system']">
        <fn:array key="traditional_system">
          <xsl:for-each select="fn:array[@key='traditional_system']/fn:string">
            <fn:string><xsl:value-of select="."/></fn:string>
          </xsl:for-each>
        </fn:array>
      </xsl:if>
      <xsl:if test="fn:string[@key='compound_id']">
        <fn:string key="compound_id"><xsl:value-of select="fn:string[@key='compound_id']"/></fn:string>
        <fn:string key="compound_name"><xsl:value-of select="fn:string[@key='compound_name']"/></fn:string>
      </xsl:if>
      <xsl:if test="fn:number[@key='pubchem_cid']">
        <fn:number key="pubchem_cid"><xsl:value-of select="fn:number[@key='pubchem_cid']"/></fn:number>
      </xsl:if>
      <xsl:if test="fn:array[@key='convergent_uses']">
        <fn:array key="tmmc_convergent_uses">
          <xsl:for-each select="fn:array[@key='convergent_uses']/fn:map">
            <fn:map>
              <xsl:copy-of select="./*"/>
            </fn:map>
          </xsl:for-each>
        </fn:array>
      </xsl:if>
      <xsl:if test="fn:map[@key='ayurveda_properties']">
        <fn:map key="tmmc_ayurveda_properties">
          <xsl:copy-of select="fn:map[@key='ayurveda_properties']/*"/>
        </fn:map>
      </xsl:if>
      <fn:map key="_source">
        <fn:string key="database">TM-MC</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:string[@key='plant_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>
</xsl:stylesheet>
