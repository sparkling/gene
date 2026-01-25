<?xml version="1.0" encoding="UTF-8"?>
<!--
  TCMSID -> Unified Traditional Chinese Medicine Schema Mapping
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
      <fn:string key="compound_id"><xsl:value-of select="fn:string[@key='compound_id']"/></fn:string>
      <fn:string key="compound_name"><xsl:value-of select="fn:string[@key='compound_name']"/></fn:string>
      <xsl:if test="fn:string[@key='herb_id']">
        <fn:string key="herb_id"><xsl:value-of select="fn:string[@key='herb_id']"/></fn:string>
        <fn:string key="herb_name"><xsl:value-of select="fn:string[@key='herb_name']"/></fn:string>
        <fn:string key="chinese_name"><xsl:value-of select="fn:string[@key='herb_name_cn']"/></fn:string>
      </xsl:if>
      <xsl:if test="fn:string[@key='smiles']">
        <fn:string key="smiles"><xsl:value-of select="fn:string[@key='smiles']"/></fn:string>
      </xsl:if>
      <xsl:if test="fn:string[@key='inchi_key']">
        <fn:string key="inchi_key"><xsl:value-of select="fn:string[@key='inchi_key']"/></fn:string>
      </xsl:if>
      <xsl:if test="fn:number[@key='molecular_weight']">
        <fn:number key="molecular_weight"><xsl:value-of select="fn:number[@key='molecular_weight']"/></fn:number>
      </xsl:if>
      <fn:map key="_source">
        <fn:string key="database">TCMSID</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:string[@key='compound_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>
</xsl:stylesheet>
