<?xml version="1.0" encoding="UTF-8"?>
<!--
  KampoDB -> Unified South/East Asian Systems Schema Mapping
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
      <fn:string key="kampo_formula_id"><xsl:value-of select="fn:string[@key='formula_id']"/></fn:string>
      <fn:string key="kampo_formula_name"><xsl:value-of select="fn:string[@key='formula_name']"/></fn:string>
      <fn:string key="kampo_formula_name_jp"><xsl:value-of select="fn:string[@key='formula_name_jp']"/></fn:string>
      <xsl:if test="fn:number[@key='crude_drug_id']">
        <fn:number key="kampo_crude_id"><xsl:value-of select="fn:number[@key='crude_drug_id']"/></fn:number>
        <fn:string key="kampo_crude_name"><xsl:value-of select="fn:string[@key='crude_drug_name']"/></fn:string>
      </xsl:if>
      <fn:string key="botanical_name"><xsl:value-of select="fn:string[@key='botanical_origin']"/></fn:string>
      <xsl:if test="fn:number[@key='compound_id']">
        <fn:number key="pubchem_cid"><xsl:value-of select="fn:number[@key='compound_id']"/></fn:number>
        <fn:string key="compound_name"><xsl:value-of select="fn:string[@key='compound_name']"/></fn:string>
      </xsl:if>
      <xsl:if test="fn:number[@key='protein_id']">
        <fn:number key="kampo_protein_id"><xsl:value-of select="fn:number[@key='protein_id']"/></fn:number>
        <fn:string key="gene_symbol"><xsl:value-of select="fn:string[@key='gene_symbol']"/></fn:string>
      </xsl:if>
      <xsl:if test="fn:map[@key='docking']">
        <fn:number key="kampo_affinity_kcal_mol"><xsl:value-of select="fn:map[@key='docking']/fn:number[@key='affinity_kcal_mol']"/></fn:number>
      </xsl:if>
      <fn:map key="_source">
        <fn:string key="database">KampoDB</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:string[@key='formula_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>
</xsl:stylesheet>
