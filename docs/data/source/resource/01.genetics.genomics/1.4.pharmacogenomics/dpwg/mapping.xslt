<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  DPWG -> Pharmacogenomics Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Pharmacogenomics Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ gene                        │ gene                    │ Direct       │
  │ drug                        │ drugs                   │ Array        │
  │ phenotype                   │ phenotype               │ Direct       │
  │ recommendation              │ recommendation          │ Direct       │
  │ evidence_level              │ evidence_level          │ Direct       │
  │ dpwg_id                     │ dpwg_id                 │ Direct       │
  │ atc_code                    │ atc_code                │ Direct       │
  │ urgency                     │ dpwg_urgency            │ Direct       │
  │ diplotype                   │ diplotype               │ Direct       │
  │ activity_score              │ activity_score          │ Direct       │
  │ allele                      │ allele                  │ Direct       │
  │ function                    │ allele_function         │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "-" → null

  NOTES:
  - DPWG: Dutch Pharmacogenetics Working Group
  - European guidelines for pharmacogenomics implementation
  - Evidence levels: Strong, Moderate
  - Urgency: Urgent (avoid/switch immediately), Non-urgent (consider alternative)
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

  <!-- ============================================================
       Entry Point: Parse JSON input
       ============================================================ -->
  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <!-- ============================================================
       Main Record Transformation
       ============================================================ -->
  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== CORE PHARMACOGENOMICS FIELDS ========== -->

      <!-- Gene -->
      <fn:string key="gene">
        <xsl:value-of select="fn:string[@key='gene']"/>
      </fn:string>

      <!-- Drugs (array) -->
      <xsl:variable name="drug" select="fn:string[@key='drug']"/>
      <xsl:choose>
        <xsl:when test="local:is-null($drug)">
          <fn:null key="drugs"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="drugs">
            <xsl:for-each select="tokenize($drug, '[;,]')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Phenotype -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">phenotype</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='phenotype']"/>
      </xsl:call-template>

      <!-- Recommendation -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">recommendation</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='recommendation']"/>
      </xsl:call-template>

      <!-- Evidence level -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">evidence_level</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='evidence_level']"/>
      </xsl:call-template>

      <!-- ========== DPWG-SPECIFIC FIELDS ========== -->

      <!-- DPWG ID -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">dpwg_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='dpwg_id']"/>
      </xsl:call-template>

      <!-- ATC code -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">atc_code</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='atc_code']"/>
      </xsl:call-template>

      <!-- Urgency -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">dpwg_urgency</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='urgency']"/>
      </xsl:call-template>

      <!-- Diplotype -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">diplotype</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='diplotype']"/>
      </xsl:call-template>

      <!-- Activity score -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='activity_score']">
          <fn:number key="activity_score">
            <xsl:value-of select="fn:number[@key='activity_score']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="activity_score"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Allele -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">allele</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='allele']"/>
      </xsl:call-template>

      <!-- Allele function -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">allele_function</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='function']"/>
      </xsl:call-template>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">DPWG</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Named Templates
       ============================================================ -->

  <xsl:template name="string-or-null">
    <xsl:param name="key"/>
    <xsl:param name="value"/>
    <xsl:choose>
      <xsl:when test="local:is-null($value)">
        <fn:null key="{$key}"/>
      </xsl:when>
      <xsl:otherwise>
        <fn:string key="{$key}">
          <xsl:value-of select="$value"/>
        </fn:string>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <!-- ============================================================
       Reusable Functions
       ============================================================ -->

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or normalize-space($value) = ('', '-', '.', 'NA', 'N/A', 'null')"/>
  </xsl:function>

</xsl:stylesheet>
