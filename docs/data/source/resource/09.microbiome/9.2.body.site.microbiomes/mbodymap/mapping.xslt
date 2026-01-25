<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  mBodyMap -> Unified Body Site Microbiome Schema Mapping
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
        <xsl:value-of select="concat('MBODYMAP:', fn:string[@key='sample_id'])"/>
      </fn:string>

      <fn:string key="sample_id">
        <xsl:value-of select="fn:string[@key='sample_id']"/>
      </fn:string>

      <fn:string key="body_site">
        <xsl:value-of select="fn:string[@key='site_name']"/>
      </fn:string>

      <fn:string key="site_code">
        <xsl:value-of select="fn:string[@key='site_code']"/>
      </fn:string>

      <fn:string key="site_category">
        <xsl:value-of select="fn:string[@key='site_category']"/>
      </fn:string>

      <fn:string key="body_site_ontology">
        <xsl:value-of select="fn:string[@key='anatomy_term']"/>
      </fn:string>

      <fn:string key="sequencing_type">
        <xsl:value-of select="fn:string[@key='sequencing_method']"/>
      </fn:string>

      <fn:string key="health_status">
        <xsl:value-of select="fn:string[@key='health_status']"/>
      </fn:string>

      <fn:map key="_source">
        <fn:string key="database">mBodyMap</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:string[@key='sample_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>

</xsl:stylesheet>
