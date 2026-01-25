<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  GMrepo -> Unified Gut Microbiome Schema Mapping
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
        <xsl:value-of select="concat('GMREPO:', fn:string[@key='run_id'])"/>
      </fn:string>

      <fn:string key="sample_id">
        <xsl:value-of select="fn:string[@key='run_id']"/>
      </fn:string>

      <fn:string key="project_id">
        <xsl:value-of select="fn:string[@key='project_id']"/>
      </fn:string>

      <fn:string key="phenotype">
        <xsl:value-of select="fn:string[@key='phenotype']"/>
      </fn:string>

      <fn:string key="phenotype_mesh_id">
        <xsl:value-of select="fn:string[@key='phenotype_mesh_id']"/>
      </fn:string>

      <fn:string key="sequencing_type">
        <xsl:value-of select="fn:string[@key='sequencing_type']"/>
      </fn:string>

      <fn:string key="country">
        <xsl:value-of select="fn:string[@key='country']"/>
      </fn:string>

      <fn:map key="_source">
        <fn:string key="database">GMrepo</fn:string>
        <fn:string key="original_id"><xsl:value-of select="fn:string[@key='run_id']"/></fn:string>
      </fn:map>
    </fn:map>
  </xsl:template>

</xsl:stylesheet>
