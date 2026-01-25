<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  NCBI ELink -> Unified Identifier Mapping Schema
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
        <xsl:value-of select="concat(fn:string[@key='dbfrom'], ':', fn:string[@key='id'])"/>
      </fn:string>

      <fn:string key="source_db">
        <xsl:value-of select="fn:string[@key='dbfrom']"/>
      </fn:string>

      <fn:string key="source_id">
        <xsl:value-of select="fn:string[@key='id']"/>
      </fn:string>

      <fn:array key="links">
        <xsl:for-each select="fn:array[@key='linksets']/fn:map">
          <fn:map>
            <fn:string key="target_db"><xsl:value-of select="fn:string[@key='dbto']"/></fn:string>
            <fn:string key="link_name"><xsl:value-of select="fn:string[@key='linkname']"/></fn:string>
            <fn:array key="ids">
              <xsl:for-each select="fn:array[@key='ids']/fn:string">
                <fn:string><xsl:value-of select="."/></fn:string>
              </xsl:for-each>
            </fn:array>
          </fn:map>
        </xsl:for-each>
      </fn:array>

      <fn:map key="_source">
        <fn:string key="database">NCBI ELink</fn:string>
        <fn:string key="command"><xsl:value-of select="fn:string[@key='cmd']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>

</xsl:stylesheet>
