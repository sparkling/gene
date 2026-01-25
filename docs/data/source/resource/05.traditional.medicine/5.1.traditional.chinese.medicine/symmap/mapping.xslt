<?xml version="1.0" encoding="UTF-8"?>
<!--
  SymMap -> Unified Traditional Chinese Medicine Schema Mapping
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
      <fn:string key="symptom_id"><xsl:value-of select="fn:string[@key='symptom_id']"/></fn:string>
      <fn:string key="symptom_name"><xsl:value-of select="fn:string[@key='symptom_name']"/></fn:string>
      <fn:string key="symptom_name_cn"><xsl:value-of select="fn:string[@key='symptom_name_cn']"/></fn:string>
      <xsl:if test="fn:string[@key='herb_id']">
        <fn:string key="herb_id"><xsl:value-of select="fn:string[@key='herb_id']"/></fn:string>
        <fn:string key="herb_name"><xsl:value-of select="fn:string[@key='herb_name']"/></fn:string>
      </xsl:if>
      <xsl:if test="fn:string[@key='compound_id']">
        <fn:string key="compound_id"><xsl:value-of select="fn:string[@key='compound_id']"/></fn:string>
      </xsl:if>
      <xsl:if test="fn:string[@key='target_id']">
        <fn:string key="target_id"><xsl:value-of select="fn:string[@key='target_id']"/></fn:string>
        <fn:string key="gene_symbol"><xsl:value-of select="fn:string[@key='gene_symbol']"/></fn:string>
      </xsl:if>
      <fn:map key="_source">
        <fn:string key="database">SymMap</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:string[@key='symptom_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>
</xsl:stylesheet>
