<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  UniProt ID Mapping -> Unified Identifier Mapping Schema
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
        <xsl:value-of select="fn:map[@key='to']/fn:string[@key='primaryAccession']"/>
      </fn:string>

      <fn:array key="identifiers">
        <fn:map>
          <fn:string key="type">uniprot</fn:string>
          <fn:string key="value"><xsl:value-of select="fn:map[@key='to']/fn:string[@key='primaryAccession']"/></fn:string>
        </fn:map>
        <fn:map>
          <fn:string key="type">input</fn:string>
          <fn:string key="value"><xsl:value-of select="fn:string[@key='from']"/></fn:string>
        </fn:map>
      </fn:array>

      <fn:string key="name">
        <xsl:value-of select="fn:map[@key='to']/fn:map[@key='proteinDescription']/fn:map[@key='recommendedName']/fn:map[@key='fullName']/fn:string[@key='value']"/>
      </fn:string>

      <fn:string key="gene">
        <xsl:value-of select="fn:map[@key='to']/fn:array[@key='genes']/fn:map[1]/fn:map[@key='geneName']/fn:string[@key='value']"/>
      </fn:string>

      <fn:map key="_source">
        <fn:string key="database">UniProt ID Mapping</fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>

</xsl:stylesheet>
