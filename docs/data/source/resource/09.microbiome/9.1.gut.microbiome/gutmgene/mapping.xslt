<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  gutMGene -> Unified Gut Microbiome Schema Mapping
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
        <xsl:value-of select="concat('GUTMGENE:', fn:number[@key='association_id'])"/>
      </fn:string>

      <fn:string key="microbe_taxid">
        <xsl:value-of select="fn:map[@key='microbe']/fn:number[@key='taxid']"/>
      </fn:string>

      <fn:string key="microbe_name">
        <xsl:value-of select="fn:map[@key='microbe']/fn:string[@key='name']"/>
      </fn:string>

      <fn:string key="gene_id">
        <xsl:value-of select="fn:map[@key='gene']/fn:number[@key='gene_id']"/>
      </fn:string>

      <fn:string key="gene_symbol">
        <xsl:value-of select="fn:map[@key='gene']/fn:string[@key='symbol']"/>
      </fn:string>

      <fn:string key="effect_direction">
        <xsl:value-of select="fn:map[@key='effect']/fn:string[@key='direction']"/>
      </fn:string>

      <fn:string key="pmid">
        <xsl:value-of select="fn:map[@key='evidence']/fn:number[@key='pmid']"/>
      </fn:string>

      <fn:map key="_source">
        <fn:string key="database">gutMGene</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:number[@key='association_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>

</xsl:stylesheet>
