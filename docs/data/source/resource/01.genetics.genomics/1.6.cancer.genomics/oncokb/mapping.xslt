<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  OncoKB -> Cancer Genomics Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Cancer Genomics Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────────┬──────────────┐
  │ Source Field                │ Target Field                │ Transform    │
  ├─────────────────────────────┼─────────────────────────────┼──────────────┤
  │ hugoSymbol                  │ gene                        │ Direct       │
  │ alteration                  │ variant                     │ Direct       │
  │ alteration                  │ oncokb_alteration           │ Direct       │
  │ oncogenic                   │ clinical_significance       │ Direct       │
  │ oncogenic                   │ oncokb_oncogenic            │ Direct       │
  │ oncogene                    │ oncokb_is_oncogene          │ Direct       │
  │ tsg                         │ oncokb_is_tsg               │ Direct       │
  │ geneAliases                 │ oncokb_gene_aliases         │ Array        │
  │ background                  │ oncokb_background           │ Direct       │
  │ alterationType              │ oncokb_alteration_type      │ Direct       │
  │ consequence                 │ oncokb_consequence          │ Direct       │
  │ proteinStart                │ oncokb_protein_start        │ Direct       │
  │ proteinEnd                  │ oncokb_protein_end          │ Direct       │
  │ mutationEffect              │ oncokb_mutation_effect      │ Object       │
  │ level                       │ oncokb_level                │ Direct       │
  │ level                       │ evidence_level              │ Direct       │
  │ cancerTypes                 │ cancer_types                │ Array        │
  │ drugs                       │ associated_drugs            │ Array        │
  │ pmids                       │ oncokb_pmids                │ Array        │
  │ abstracts                   │ oncokb_abstracts            │ Array        │
  │ description                 │ oncokb_description          │ Direct       │
  └─────────────────────────────┴─────────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings → null

  NOTES:
  - OncoKB: MSK Precision Oncology Knowledge Base
  - Evidence levels: 1, 2, 3A, 3B, 4, R1, R2
  - Oncogenic: Oncogenic, Likely Oncogenic, Predicted Oncogenic, etc.
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
        <xsl:value-of select="(fn:string[@key='hugoSymbol'], fn:map[@key='gene']/fn:string[@key='hugoSymbol'])[1]"/>
      </fn:string>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">variant</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='alteration']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">clinical_significance</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='oncogenic']"/>
      </xsl:call-template>

      <!-- Cancer types -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='cancerTypes']">
          <fn:array key="cancer_types">
            <xsl:for-each select="fn:array[@key='cancerTypes']/fn:map">
              <fn:string><xsl:value-of select="(fn:string[@key='mainType'], fn:string[@key='name'])[1]"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="cancer_types"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">evidence_level</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='level']"/>
      </xsl:call-template>

      <!-- ========== ONCOKB-SPECIFIC FIELDS ========== -->

      <xsl:choose>
        <xsl:when test="fn:boolean[@key='oncogene']">
          <fn:boolean key="oncokb_is_oncogene">
            <xsl:value-of select="fn:boolean[@key='oncogene']"/>
          </fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="oncokb_is_oncogene"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:boolean[@key='tsg']">
          <fn:boolean key="oncokb_is_tsg">
            <xsl:value-of select="fn:boolean[@key='tsg']"/>
          </fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="oncokb_is_tsg"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Gene aliases -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='geneAliases']">
          <fn:array key="oncokb_gene_aliases">
            <xsl:for-each select="fn:array[@key='geneAliases']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="oncokb_gene_aliases"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">oncokb_background</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='background']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">oncokb_alteration</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='alteration']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">oncokb_alteration_type</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='alterationType']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">oncokb_consequence</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='consequence']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='proteinStart']">
          <fn:number key="oncokb_protein_start">
            <xsl:value-of select="fn:number[@key='proteinStart']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="oncokb_protein_start"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='proteinEnd']">
          <fn:number key="oncokb_protein_end">
            <xsl:value-of select="fn:number[@key='proteinEnd']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="oncokb_protein_end"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">oncokb_oncogenic</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='oncogenic']"/>
      </xsl:call-template>

      <!-- Mutation effect (preserve as object) -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='mutationEffect']">
          <xsl:copy-of select="fn:map[@key='mutationEffect']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="oncokb_mutation_effect"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">oncokb_level</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='level']"/>
      </xsl:call-template>

      <!-- Associated drugs -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='drugs']">
          <fn:array key="associated_drugs">
            <xsl:for-each select="fn:array[@key='drugs']/fn:map">
              <fn:string><xsl:value-of select="fn:string[@key='drugName']"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="associated_drugs"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- PMIDs -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='pmids']">
          <fn:array key="oncokb_pmids">
            <xsl:for-each select="fn:array[@key='pmids']/fn:number">
              <fn:number><xsl:value-of select="."/></fn:number>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="oncokb_pmids"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Abstracts -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='abstracts']">
          <fn:array key="oncokb_abstracts">
            <xsl:for-each select="fn:array[@key='abstracts']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="oncokb_abstracts"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">oncokb_description</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='description']"/>
      </xsl:call-template>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">OncoKB</fn:string>
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
