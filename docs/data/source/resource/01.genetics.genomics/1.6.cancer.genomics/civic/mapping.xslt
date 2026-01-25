<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  CIViC -> Cancer Genomics Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Cancer Genomics Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────────┬──────────────┐
  │ Source Field                │ Target Field                │ Transform    │
  ├─────────────────────────────┼─────────────────────────────┼──────────────┤
  │ gene.name                   │ gene                        │ Direct       │
  │ name                        │ civic_variant_name          │ Direct       │
  │ name                        │ variant                     │ Direct       │
  │ id                          │ civic_id                    │ Direct       │
  │ entrez_id                   │ entrez_gene_id              │ Direct       │
  │ variant_types               │ civic_variant_types         │ Array        │
  │ coordinates                 │ civic_coordinates           │ Object       │
  │ hgvs_expressions            │ civic_hgvs_expressions      │ Array        │
  │ evidence_type               │ civic_evidence_type         │ Direct       │
  │ evidence_level              │ evidence_level              │ Direct       │
  │ evidence_direction          │ civic_evidence_direction    │ Direct       │
  │ clinical_significance       │ civic_clinical_significance │ Direct       │
  │ clinical_significance       │ clinical_significance       │ Direct       │
  │ disease                     │ civic_disease               │ Object       │
  │ disease.name                │ cancer_types                │ Array        │
  │ drugs                       │ associated_drugs            │ Array        │
  │ source                      │ civic_source                │ Object       │
  │ amp_level                   │ amp_level                   │ Direct       │
  └─────────────────────────────┴─────────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings → null

  NOTES:
  - CIViC: Clinical Interpretation of Variants in Cancer
  - Evidence types: Predictive, Diagnostic, Prognostic, Predisposing
  - Evidence levels: A, B, C, D, E
  - AMP/ASCO/CAP tier classification for somatic variants
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
      <!-- ========== CORE FIELDS ========== -->
      <fn:string key="gene">
        <xsl:value-of select="(fn:map[@key='gene']/fn:string[@key='name'], fn:string[@key='gene'])[1]"/>
      </fn:string>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">variant</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='name']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">clinical_significance</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='clinical_significance']"/>
      </xsl:call-template>

      <!-- Cancer types from disease -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='disease']/fn:string[@key='name']">
          <fn:array key="cancer_types">
            <fn:string><xsl:value-of select="fn:map[@key='disease']/fn:string[@key='name']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:array[@key='diseases']">
          <fn:array key="cancer_types">
            <xsl:for-each select="fn:array[@key='diseases']/fn:map">
              <fn:string><xsl:value-of select="fn:string[@key='name']"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="cancer_types"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">evidence_level</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='evidence_level']"/>
      </xsl:call-template>

      <!-- ========== CIVIC-SPECIFIC FIELDS ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='id']">
          <fn:number key="civic_id">
            <xsl:value-of select="fn:number[@key='id']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="civic_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">civic_variant_name</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='name']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='entrez_id']">
          <fn:number key="entrez_gene_id">
            <xsl:value-of select="fn:number[@key='entrez_id']"/>
          </fn:number>
        </xsl:when>
        <xsl:when test="fn:map[@key='gene']/fn:number[@key='entrez_id']">
          <fn:number key="entrez_gene_id">
            <xsl:value-of select="fn:map[@key='gene']/fn:number[@key='entrez_id']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="entrez_gene_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Variant types (array) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='variant_types']">
          <fn:array key="civic_variant_types">
            <xsl:for-each select="fn:array[@key='variant_types']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="civic_variant_types"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Coordinates (preserve as object) -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='coordinates']">
          <xsl:copy-of select="fn:map[@key='coordinates']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="civic_coordinates"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- HGVS expressions (array) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='hgvs_expressions']">
          <fn:array key="civic_hgvs_expressions">
            <xsl:for-each select="fn:array[@key='hgvs_expressions']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="civic_hgvs_expressions"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">civic_evidence_type</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='evidence_type']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">civic_evidence_direction</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='evidence_direction']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">civic_clinical_significance</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='clinical_significance']"/>
      </xsl:call-template>

      <!-- Disease (preserve as object) -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='disease']">
          <xsl:copy-of select="fn:map[@key='disease']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="civic_disease"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Associated drugs (array) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='drugs']">
          <fn:array key="associated_drugs">
            <xsl:for-each select="fn:array[@key='drugs']/fn:map">
              <fn:string><xsl:value-of select="fn:string[@key='name']"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="associated_drugs"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Source (preserve as object) -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='source']">
          <xsl:copy-of select="fn:map[@key='source']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="civic_source"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">amp_level</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='amp_level']"/>
      </xsl:call-template>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">CIViC</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

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

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or normalize-space($value) = ('', '-', '.', 'NA', 'N/A', 'null')"/>
  </xsl:function>

</xsl:stylesheet>
