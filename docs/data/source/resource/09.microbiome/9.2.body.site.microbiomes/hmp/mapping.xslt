<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  HMP Body Sites -> Unified Body Site Microbiome Schema Mapping
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
        <xsl:value-of select="concat('HMP:', fn:string[@key='sample_id'])"/>
      </fn:string>

      <fn:string key="sample_id">
        <xsl:value-of select="fn:string[@key='sample_id']"/>
      </fn:string>

      <fn:string key="subject_id">
        <xsl:value-of select="fn:string[@key='rand_subject_id']"/>
      </fn:string>

      <fn:string key="body_site">
        <xsl:value-of select="fn:string[@key='body_site']"/>
      </fn:string>

      <fn:string key="supersite">
        <xsl:value-of select="fn:string[@key='supersite']"/>
      </fn:string>

      <fn:string key="body_site_ontology">
        <xsl:value-of select="fn:string[@key='fma_body_site']"/>
      </fn:string>

      <fn:string key="sequencing_type">
        <xsl:value-of select="fn:map[@key='sequence_data']/fn:string[@key='seq_type']"/>
      </fn:string>

      <fn:string key="study">
        <xsl:value-of select="fn:string[@key='study']"/>
      </fn:string>

      <fn:map key="_source">
        <fn:string key="database">Human Microbiome Project</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:string[@key='sample_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>

</xsl:stylesheet>
