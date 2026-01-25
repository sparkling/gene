<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  CPIC -> Pharmacogenomics Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Pharmacogenomics Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────────┬──────────────┐
  │ Source Field                │ Target Field                │ Transform    │
  ├─────────────────────────────┼─────────────────────────────┼──────────────┤
  │ gene                        │ gene                        │ Direct       │
  │ drug                        │ drugs                       │ Array        │
  │ phenotype                   │ phenotype                   │ Direct       │
  │ recommendation              │ recommendation              │ Direct       │
  │ cpic_level                  │ evidence_level              │ Map A→1A     │
  │ guideline_id                │ cpic_guideline_id           │ Direct       │
  │ publication_doi             │ cpic_publication_doi        │ Direct       │
  │ cpic_level                  │ cpic_level                  │ Direct       │
  │ allele                      │ allele                      │ Direct       │
  │ function                    │ allele_function             │ Direct       │
  │ activity_score              │ activity_score              │ Direct       │
  │ strength                    │ recommendation_strength     │ Direct       │
  │ implications                │ clinical_implications       │ Direct       │
  │ diplotype                   │ diplotype                   │ Direct       │
  └─────────────────────────────┴─────────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "-" → null

  NOTES:
  - CPIC: Clinical Pharmacogenetics Implementation Consortium
  - CPIC levels: A (required), B (recommended), C (info only), D (info only)
  - Recommendation strength: Strong, Moderate, Optional
  - Activity scores for allele function quantification
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
        <xsl:when test="fn:array[@key='drugs']">
          <fn:array key="drugs">
            <xsl:for-each select="fn:array[@key='drugs']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="not(local:is-null($drug))">
          <fn:array key="drugs">
            <xsl:for-each select="tokenize($drug, '[;,]')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="drugs"/>
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

      <!-- Evidence level (map CPIC level) -->
      <xsl:variable name="level" select="fn:string[@key='cpic_level']"/>
      <xsl:choose>
        <xsl:when test="local:is-null($level)">
          <fn:null key="evidence_level"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="evidence_level">
            <xsl:choose>
              <xsl:when test="$level = 'A'">1A</xsl:when>
              <xsl:when test="$level = 'B'">1B</xsl:when>
              <xsl:otherwise><xsl:value-of select="$level"/></xsl:otherwise>
            </xsl:choose>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CPIC-SPECIFIC FIELDS ========== -->

      <!-- CPIC guideline ID -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cpic_guideline_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='guideline_id']"/>
      </xsl:call-template>

      <!-- Publication DOI -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cpic_publication_doi</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='publication_doi']"/>
      </xsl:call-template>

      <!-- CPIC level -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cpic_level</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='cpic_level']"/>
      </xsl:call-template>

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

      <!-- Recommendation strength -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">recommendation_strength</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='strength']"/>
      </xsl:call-template>

      <!-- Clinical implications -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">clinical_implications</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='implications']"/>
      </xsl:call-template>

      <!-- Diplotype -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">diplotype</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='diplotype']"/>
      </xsl:call-template>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">CPIC</fn:string>
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
