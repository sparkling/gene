<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  FDA OpenFDA -> Unified Regulatory Schema Mapping
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
        <xsl:value-of select="concat('FAERS:', fn:string[@key='safetyreportid'])"/>
      </fn:string>

      <fn:string key="report_date">
        <xsl:value-of select="fn:string[@key='receivedate']"/>
      </fn:string>

      <fn:boolean key="is_serious">
        <xsl:value-of select="fn:string[@key='serious'] = '1'"/>
      </fn:boolean>

      <fn:array key="reactions">
        <xsl:for-each select="fn:map[@key='patient']/fn:array[@key='reaction']/fn:map">
          <fn:string><xsl:value-of select="fn:string[@key='reactionmeddrapt']"/></fn:string>
        </xsl:for-each>
      </fn:array>

      <fn:array key="drugs">
        <xsl:for-each select="fn:map[@key='patient']/fn:array[@key='drug']/fn:map[fn:string[@key='drugcharacterization'] = '1']">
          <fn:map>
            <fn:string key="name"><xsl:value-of select="fn:string[@key='medicinalproduct']"/></fn:string>
            <fn:string key="indication"><xsl:value-of select="fn:string[@key='drugindication']"/></fn:string>
          </fn:map>
        </xsl:for-each>
      </fn:array>

      <fn:map key="_source">
        <fn:string key="database">FDA OpenFDA FAERS</fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>

</xsl:stylesheet>
