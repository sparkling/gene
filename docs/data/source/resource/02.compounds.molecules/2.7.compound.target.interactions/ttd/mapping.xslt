<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  TTD -> Compound-Target Interactions Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Therapeutic Target Database)
  Target: ../schema.json (2.7 Compound-Target Interactions Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ drug_id                     │ ttd_drug_id             │ Direct       │
  │ drug_name                   │ compound_name           │ Direct       │
  │ target_id                   │ ttd_target_id           │ Direct       │
  │ target_name                 │ target_name             │ Direct       │
  │ uniprot_id                  │ uniprot_id              │ Direct       │
  │ target_validation           │ target_validation       │ Direct       │
  │ development_status          │ development_status      │ Direct       │
  │ drug_class                  │ drug_class              │ Direct       │
  │ indication                  │ indication              │ Direct       │
  │ icd_11_code                 │ icd_11_code             │ Direct       │
  │ pathway_associations        │ pathway_associations    │ Object array │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null

  NOTES:
  - TTD uses TTDT prefix for target IDs
  - Includes therapeutic target validation status
  - Development status: Approved, Phase I-III, etc.
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

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='drug_id'])">
          <fn:null key="ttd_drug_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="ttd_drug_id">
            <xsl:value-of select="fn:string[@key='drug_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='target_id'])">
          <fn:null key="ttd_target_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="ttd_target_id">
            <xsl:value-of select="fn:string[@key='target_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DRUG/COMPOUND ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='drug_name'])">
          <fn:null key="compound_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="compound_name">
            <xsl:value-of select="fn:string[@key='drug_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='drug_class'])">
          <fn:null key="drug_class"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="drug_class">
            <xsl:value-of select="fn:string[@key='drug_class']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='development_status'])">
          <fn:null key="development_status"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="development_status">
            <xsl:value-of select="fn:string[@key='development_status']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== TARGET ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='target_name'])">
          <fn:null key="target_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="target_name">
            <xsl:value-of select="fn:string[@key='target_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='uniprot_id'])">
          <fn:null key="uniprot_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="uniprot_id">
            <xsl:value-of select="fn:string[@key='uniprot_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='target_validation'])">
          <fn:null key="target_validation"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="target_validation">
            <xsl:value-of select="fn:string[@key='target_validation']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== INDICATION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='indication'])">
          <fn:null key="indication"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="indication">
            <xsl:value-of select="fn:string[@key='indication']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='icd_11_code'])">
          <fn:null key="icd_11_code"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="icd_11_code">
            <xsl:value-of select="fn:string[@key='icd_11_code']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PATHWAY ASSOCIATIONS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='pathway_associations']/fn:map">
          <fn:array key="pathway_associations">
            <xsl:for-each select="fn:array[@key='pathway_associations']/fn:map">
              <fn:map>
                <xsl:if test="fn:string[@key='pathway_id']">
                  <fn:string key="pathway_id">
                    <xsl:value-of select="fn:string[@key='pathway_id']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:string[@key='source']">
                  <fn:string key="source">
                    <xsl:value-of select="fn:string[@key='source']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="pathway_associations"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">TTD</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='drug_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
