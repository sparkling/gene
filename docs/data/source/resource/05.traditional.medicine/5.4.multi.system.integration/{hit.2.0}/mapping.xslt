<?xml version="1.0" encoding="UTF-8"?>
<!--
  HIT 2.0 -> Unified Multi-System Integration Schema Mapping
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
      <fn:string key="ingredient_id"><xsl:value-of select="fn:string[@key='ingredient_id']"/></fn:string>
      <fn:string key="ingredient_name"><xsl:value-of select="fn:string[@key='ingredient_name']"/></fn:string>
      <fn:string key="target_id"><xsl:value-of select="fn:string[@key='target_id']"/></fn:string>
      <fn:string key="gene_symbol"><xsl:value-of select="fn:string[@key='gene_symbol']"/></fn:string>
      <fn:string key="target_name"><xsl:value-of select="fn:string[@key='target_name']"/></fn:string>
      <xsl:if test="fn:number[@key='pubchem_cid']">
        <fn:number key="pubchem_cid"><xsl:value-of select="fn:number[@key='pubchem_cid']"/></fn:number>
      </xsl:if>
      <xsl:if test="fn:string[@key='smiles']">
        <fn:string key="smiles"><xsl:value-of select="fn:string[@key='smiles']"/></fn:string>
      </xsl:if>
      <xsl:if test="fn:array[@key='source_herbs']">
        <fn:array key="source_herbs">
          <xsl:for-each select="fn:array[@key='source_herbs']/fn:map">
            <fn:map>
              <xsl:copy-of select="./*"/>
            </fn:map>
          </xsl:for-each>
        </fn:array>
      </xsl:if>
      <xsl:if test="fn:map[@key='evidence']">
        <fn:map key="evidence">
          <xsl:copy-of select="fn:map[@key='evidence']/*"/>
        </fn:map>
      </xsl:if>
      <xsl:if test="fn:map[@key='binding_affinity']">
        <fn:map key="binding_affinity">
          <xsl:copy-of select="fn:map[@key='binding_affinity']/*"/>
        </fn:map>
      </xsl:if>
      <xsl:if test="fn:array[@key='traditional_system_coverage']">
        <fn:array key="traditional_systems">
          <xsl:for-each select="fn:array[@key='traditional_system_coverage']/fn:string">
            <fn:string><xsl:value-of select="."/></fn:string>
          </xsl:for-each>
        </fn:array>
      </xsl:if>
      <fn:map key="_source">
        <fn:string key="database">HIT 2.0</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:string[@key='ingredient_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>
</xsl:stylesheet>
