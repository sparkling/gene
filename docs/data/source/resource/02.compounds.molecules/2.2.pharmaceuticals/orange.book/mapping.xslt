<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Orange Book -> Pharmaceuticals Unified Schema Mapping
  ============================================================
  Source: ./schema.md (FDA Orange Book - Approved Drug Products)
  Target: ../schema.json (2.2 Pharmaceuticals Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ trade_name                  │ name                │ Direct       │
  │ appl_no                     │ appl_no             │ Direct       │
  │ appl_no                     │ application_number  │ Direct       │
  │ te_code                     │ te_code             │ Direct       │
  │ patent_no                   │ patent_no           │ Direct       │
  │ patent_expire_date          │ patent_expire_date  │ Direct       │
  │ exclusivity_code            │ exclusivity_code    │ Direct       │
  │ rld                         │ reference_drug      │ Boolean      │
  │ (derived)                   │ approval_status     │ Set approved │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - trade_name: "" → null
  - patent_no: "" → null

  NOTES:
  - Orange Book contains FDA-approved drug products with TE evaluations
  - TE codes: AA, AB, BC, etc. indicate therapeutic equivalence
  - All products in Orange Book are FDA-approved
  - RLD = Reference Listed Drug
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
        <xsl:when test="local:is-null(fn:string[@key='appl_no'])">
          <fn:null key="appl_no"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="appl_no">
            <xsl:value-of select="fn:string[@key='appl_no']"/>
          </fn:string>
          <fn:string key="application_number">
            <xsl:value-of select="fn:string[@key='appl_no']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='trade_name'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='trade_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- All Orange Book products are FDA-approved -->
      <fn:string key="approval_status">approved</fn:string>

      <!-- ========== THERAPEUTIC EQUIVALENCE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='te_code'])">
          <fn:null key="te_code"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="te_code">
            <xsl:value-of select="fn:string[@key='te_code']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PATENT INFORMATION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='patent_no'])">
          <fn:null key="patent_no"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="patent_no">
            <xsl:value-of select="fn:string[@key='patent_no']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='patent_expire_date'])">
          <fn:null key="patent_expire_date"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="patent_expire_date">
            <xsl:value-of select="fn:string[@key='patent_expire_date']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== EXCLUSIVITY ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='exclusivity_code'])">
          <fn:null key="exclusivity_code"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="exclusivity_code">
            <xsl:value-of select="fn:string[@key='exclusivity_code']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== REFERENCE LISTED DRUG ========== -->

      <xsl:choose>
        <xsl:when test="fn:string[@key='rld'] = 'Yes' or fn:string[@key='rld'] = 'Y' or fn:boolean[@key='rld'] = true()">
          <fn:boolean key="reference_drug">true</fn:boolean>
        </xsl:when>
        <xsl:when test="fn:string[@key='rld'] = 'No' or fn:string[@key='rld'] = 'N' or fn:boolean[@key='rld'] = false()">
          <fn:boolean key="reference_drug">false</fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="reference_drug"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">Orange Book</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='appl_no']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
