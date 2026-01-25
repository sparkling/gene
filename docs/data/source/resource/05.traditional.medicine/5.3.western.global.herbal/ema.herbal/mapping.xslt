<?xml version="1.0" encoding="UTF-8"?>
<!--
  EMA Herbal -> Unified Western & Global Herbal Schema Mapping
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
      <fn:string key="monograph_id"><xsl:value-of select="fn:string[@key='monograph_id']"/></fn:string>
      <fn:string key="herbal_substance"><xsl:value-of select="fn:string[@key='herbal_substance']"/></fn:string>
      <fn:string key="plant_part"><xsl:value-of select="fn:string[@key='plant_part']"/></fn:string>
      <fn:string key="common_name"><xsl:value-of select="fn:string[@key='common_name']"/></fn:string>
      <fn:string key="regulatory_status"><xsl:value-of select="fn:string[@key='status']"/></fn:string>
      <xsl:if test="fn:array[@key='therapeutic_area']">
        <fn:array key="therapeutic_areas">
          <xsl:for-each select="fn:array[@key='therapeutic_area']/fn:string">
            <fn:string><xsl:value-of select="."/></fn:string>
          </xsl:for-each>
        </fn:array>
      </xsl:if>
      <xsl:if test="fn:array[@key='indications']">
        <fn:array key="indications">
          <xsl:for-each select="fn:array[@key='indications']/fn:map">
            <fn:map>
              <xsl:copy-of select="./*"/>
            </fn:map>
          </xsl:for-each>
        </fn:array>
      </xsl:if>
      <xsl:if test="fn:map[@key='posology']">
        <fn:map key="posology">
          <xsl:copy-of select="fn:map[@key='posology']/*"/>
        </fn:map>
      </xsl:if>
      <xsl:if test="fn:array[@key='contraindications']">
        <fn:array key="contraindications">
          <xsl:for-each select="fn:array[@key='contraindications']/fn:string">
            <fn:string><xsl:value-of select="."/></fn:string>
          </xsl:for-each>
        </fn:array>
      </xsl:if>
      <fn:map key="_source">
        <fn:string key="database">EMA Herbal</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:string[@key='monograph_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>
</xsl:stylesheet>
