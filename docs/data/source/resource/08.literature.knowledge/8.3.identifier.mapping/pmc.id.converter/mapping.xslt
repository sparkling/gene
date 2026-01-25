<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  PMC ID Converter -> Unified Identifier Mapping Schema
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified ID Mapping Schema)
  XSLT Version: 3.0
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    xmlns:map="http://www.w3.org/2005/xpath-functions/map"
    xmlns:array="http://www.w3.org/2005/xpath-functions/array"
    exclude-result-prefixes="xs fn map array">

  <xsl:output method="json" indent="yes"/>

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <xsl:template match="fn:map">
    <fn:map>
      <!-- Primary identifier (prefer PMCID) -->
      <fn:string key="id">
        <xsl:choose>
          <xsl:when test="fn:string[@key='pmcid']">
            <xsl:value-of select="fn:string[@key='pmcid']"/>
          </xsl:when>
          <xsl:when test="fn:string[@key='pmid']">
            <xsl:value-of select="concat('PMID:', fn:string[@key='pmid'])"/>
          </xsl:when>
          <xsl:otherwise>unknown</xsl:otherwise>
        </xsl:choose>
      </fn:string>

      <!-- All identifiers as array -->
      <fn:array key="identifiers">
        <xsl:if test="fn:string[@key='pmcid']">
          <fn:map>
            <fn:string key="type">pmcid</fn:string>
            <fn:string key="value"><xsl:value-of select="fn:string[@key='pmcid']"/></fn:string>
          </fn:map>
        </xsl:if>
        <xsl:if test="fn:string[@key='pmid']">
          <fn:map>
            <fn:string key="type">pmid</fn:string>
            <fn:string key="value"><xsl:value-of select="fn:string[@key='pmid']"/></fn:string>
          </fn:map>
        </xsl:if>
        <xsl:if test="fn:string[@key='doi']">
          <fn:map>
            <fn:string key="type">doi</fn:string>
            <fn:string key="value"><xsl:value-of select="fn:string[@key='doi']"/></fn:string>
          </fn:map>
        </xsl:if>
      </fn:array>

      <!-- Status -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='errmsg']">
          <fn:string key="status">error</fn:string>
          <fn:string key="error"><xsl:value-of select="fn:string[@key='errmsg']"/></fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="status">success</fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">PMC ID Converter</fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
