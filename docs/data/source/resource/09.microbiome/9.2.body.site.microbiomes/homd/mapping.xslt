<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  HOMD -> Unified Body Site Microbiome Schema Mapping
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
        <xsl:value-of select="concat('HOMD:', fn:string[@key='homt_id'])"/>
      </fn:string>

      <fn:string key="taxon_id">
        <xsl:value-of select="fn:string[@key='homt_id']"/>
      </fn:string>

      <fn:string key="ncbi_taxid">
        <xsl:value-of select="fn:number[@key='ncbi_taxid']"/>
      </fn:string>

      <fn:string key="scientific_name">
        <xsl:value-of select="fn:string[@key='full_name']"/>
      </fn:string>

      <fn:string key="genus">
        <xsl:value-of select="fn:string[@key='genus']"/>
      </fn:string>

      <fn:string key="species">
        <xsl:value-of select="fn:string[@key='species']"/>
      </fn:string>

      <fn:string key="cultivation_status">
        <xsl:value-of select="fn:string[@key='cultivation_status']"/>
      </fn:string>

      <fn:string key="body_site">oral</fn:string>

      <fn:map key="_source">
        <fn:string key="database">Human Oral Microbiome Database</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:string[@key='homt_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>

</xsl:stylesheet>
