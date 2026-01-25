<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  MetaHIT -> Unified Gut Microbiome Schema Mapping
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
        <xsl:value-of select="concat('METAHIT:', fn:string[@key='gene_id'])"/>
      </fn:string>

      <fn:string key="gene_id">
        <xsl:value-of select="fn:string[@key='gene_id']"/>
      </fn:string>

      <fn:number key="length">
        <xsl:value-of select="fn:number[@key='length']"/>
      </fn:number>

      <fn:string key="taxonomy">
        <xsl:value-of select="fn:map[@key='taxonomy']/fn:string[@key='species']"/>
      </fn:string>

      <fn:string key="kegg_ko">
        <xsl:value-of select="fn:map[@key='functional_annotation']/fn:string[@key='kegg_ko']"/>
      </fn:string>

      <fn:string key="sample_origin">
        <xsl:value-of select="fn:string[@key='sample_origin']"/>
      </fn:string>

      <fn:map key="_source">
        <fn:string key="database">MetaHIT</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:string[@key='gene_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>

</xsl:stylesheet>
