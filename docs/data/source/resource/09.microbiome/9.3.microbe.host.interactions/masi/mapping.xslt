<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  MASI -> Unified Microbe-Host Interactions Schema Mapping
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
        <xsl:value-of select="concat('MASI:', fn:number[@key='interaction_id'])"/>
      </fn:string>

      <fn:string key="metabolite_name">
        <xsl:value-of select="fn:map[@key='metabolite']/fn:string[@key='name']"/>
      </fn:string>

      <fn:string key="metabolite_pubchem">
        <xsl:value-of select="fn:map[@key='metabolite']/fn:number[@key='pubchem_cid']"/>
      </fn:string>

      <fn:string key="metabolite_hmdb">
        <xsl:value-of select="fn:map[@key='metabolite']/fn:string[@key='hmdb_id']"/>
      </fn:string>

      <fn:string key="target_symbol">
        <xsl:value-of select="fn:map[@key='target']/fn:string[@key='symbol']"/>
      </fn:string>

      <fn:string key="target_type">
        <xsl:value-of select="fn:map[@key='target']/fn:string[@key='type']"/>
      </fn:string>

      <fn:string key="effect">
        <xsl:value-of select="fn:map[@key='interaction']/fn:string[@key='effect']"/>
      </fn:string>

      <fn:string key="pathway">
        <xsl:value-of select="fn:map[@key='pathway']/fn:string[@key='name']"/>
      </fn:string>

      <fn:number key="pmid">
        <xsl:value-of select="fn:map[@key='evidence']/fn:number[@key='pmid']"/>
      </fn:number>

      <fn:map key="_source">
        <fn:string key="database">MASI</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:number[@key='interaction_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>

</xsl:stylesheet>
