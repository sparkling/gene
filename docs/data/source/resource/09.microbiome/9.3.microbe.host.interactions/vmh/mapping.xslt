<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  VMH -> Unified Microbe-Host Interactions Schema Mapping
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
        <xsl:value-of select="concat('VMH:', fn:string[@key='model_id'])"/>
      </fn:string>

      <fn:string key="model_id">
        <xsl:value-of select="fn:string[@key='model_id']"/>
      </fn:string>

      <fn:string key="model_type">
        <xsl:value-of select="fn:string[@key='model_type']"/>
      </fn:string>

      <fn:string key="organism_name">
        <xsl:value-of select="fn:string[@key='organism_name']"/>
      </fn:string>

      <fn:number key="ncbi_taxid">
        <xsl:value-of select="fn:number[@key='ncbi_taxid']"/>
      </fn:number>

      <fn:number key="reactions_count">
        <xsl:value-of select="fn:number[@key='reactions_count']"/>
      </fn:number>

      <fn:number key="metabolites_count">
        <xsl:value-of select="fn:number[@key='metabolites_count']"/>
      </fn:number>

      <fn:map key="_source">
        <fn:string key="database">Virtual Metabolic Human</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:string[@key='model_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>

</xsl:stylesheet>
