<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  gutMDisorder -> Unified Microbe-Host Interactions Schema Mapping
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
        <xsl:value-of select="concat('GUTMDISORDER:', fn:number[@key='association_id'])"/>
      </fn:string>

      <fn:number key="microbe_taxid">
        <xsl:value-of select="fn:map[@key='microbe']/fn:number[@key='taxid']"/>
      </fn:number>

      <fn:string key="microbe_name">
        <xsl:value-of select="fn:map[@key='microbe']/fn:string[@key='name']"/>
      </fn:string>

      <fn:string key="disease_name">
        <xsl:value-of select="fn:map[@key='disease']/fn:string[@key='name']"/>
      </fn:string>

      <fn:string key="disease_icd10">
        <xsl:value-of select="fn:map[@key='disease']/fn:string[@key='icd10']"/>
      </fn:string>

      <fn:string key="disease_mesh_id">
        <xsl:value-of select="fn:map[@key='disease']/fn:string[@key='mesh_id']"/>
      </fn:string>

      <fn:string key="direction">
        <xsl:value-of select="fn:map[@key='association']/fn:string[@key='direction']"/>
      </fn:string>

      <fn:number key="pmid">
        <xsl:value-of select="fn:map[@key='evidence']/fn:number[@key='pmid']"/>
      </fn:number>

      <fn:map key="_source">
        <fn:string key="database">gutMDisorder</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:number[@key='association_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>

</xsl:stylesheet>
