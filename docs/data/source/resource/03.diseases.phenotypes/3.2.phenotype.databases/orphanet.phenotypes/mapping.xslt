<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Orphanet Phenotypes -> Phenotype Databases Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Orphanet Disease-Phenotype Annotations)
  Target: ../schema.json (3.2 Phenotype Databases Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ OrphaCode                   │ disease_id              │ Prefix       │
  │ Name                        │ disease_name            │ Direct       │
  │ HPO_ID                      │ hpo_ids                 │ Array        │
  │ HPO_Frequency               │ hpo_frequency           │ Object       │
  │ Diagnostic                  │ diagnostic              │ Boolean      │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NOTES:
  - Links Orphanet rare diseases to HPO phenotypes
  - HPO_Frequency provides structured frequency classification
  - Diagnostic flag indicates pathognomonic features
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    xmlns:map="http://www.w3.org/2005/xpath-functions/map"
    xmlns:array="http://www.w3.org/2005/xpath-functions/array"
    xmlns:local="http://local.functions"
    exclude-result-prefixes="xs fn map array local">

  <xsl:output method="json" indent="yes"/>

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== IDENTIFIERS ========== -->

      <fn:string key="disease_id">
        <xsl:value-of select="concat('ORPHA:', fn:number[@key='OrphaCode'])"/>
      </fn:string>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Name'])">
          <fn:null key="disease_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="disease_name">
            <xsl:value-of select="fn:string[@key='Name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== HPO PHENOTYPES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='HPO_ID']/fn:string">
          <fn:array key="hpo_ids">
            <xsl:for-each select="fn:array[@key='HPO_ID']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="hpo_ids"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== FREQUENCY ========== -->

      <xsl:choose>
        <xsl:when test="fn:map[@key='HPO_Frequency']">
          <xsl:copy-of select="fn:map[@key='HPO_Frequency']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="hpo_frequency"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DIAGNOSTIC FLAG ========== -->

      <xsl:choose>
        <xsl:when test="fn:boolean[@key='Diagnostic']">
          <fn:boolean key="diagnostic">
            <xsl:value-of select="fn:boolean[@key='Diagnostic']"/>
          </fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="diagnostic"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">Orphanet Phenotypes</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:number[@key='OrphaCode']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
