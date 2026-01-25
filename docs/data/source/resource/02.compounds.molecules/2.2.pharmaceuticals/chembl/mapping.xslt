<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  ChEMBL -> Pharmaceuticals Unified Schema Mapping
  ============================================================
  Source: ./schema.md (ChEMBL Bioactivity Database)
  Target: ../schema.json (2.2 Pharmaceuticals Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ chembl_id                   │ drug_id             │ Direct       │
  │ chembl_id                   │ chembl_id           │ Direct       │
  │ pref_name                   │ name                │ Direct       │
  │ canonical_smiles            │ canonical_smiles    │ Direct       │
  │ standard_inchi_key          │ inchi_key           │ Direct       │
  │ full_mwt                    │ molecular_weight    │ Direct       │
  │ molformula                  │ molecular_formula   │ Direct       │
  │ molecule_type               │ drug_type           │ Direct       │
  │ max_phase                   │ max_phase           │ Direct       │
  │ max_phase                   │ approval_status     │ Enum map     │
  │ molregno                    │ molregno            │ Direct       │
  │ pchembl_value               │ pchembl_value       │ Direct       │
  │ natural_product             │ natural_product     │ Boolean      │
  │ np_likeness_score           │ np_likeness_score   │ Direct       │
  │ qed_weighted                │ qed_weighted        │ Direct       │
  │ confidence_score            │ confidence_score    │ Direct       │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - pref_name: "" → null
  - max_phase: -1 → null

  NOTES:
  - ChEMBL uses CHEMBL prefix for compound IDs (e.g., CHEMBL25)
  - max_phase: 0=preclinical, 1-3=clinical phases, 4=approved
  - natural_product is a 0/1 flag converted to boolean
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

      <fn:string key="drug_id">
        <xsl:value-of select="fn:string[@key='chembl_id']"/>
      </fn:string>

      <fn:string key="chembl_id">
        <xsl:value-of select="fn:string[@key='chembl_id']"/>
      </fn:string>

      <xsl:if test="fn:number[@key='molregno']">
        <fn:number key="molregno">
          <xsl:value-of select="fn:number[@key='molregno']"/>
        </fn:number>
      </xsl:if>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='pref_name'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='pref_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CHEMICAL STRUCTURE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='canonical_smiles'])">
          <fn:null key="canonical_smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="canonical_smiles">
            <xsl:value-of select="fn:string[@key='canonical_smiles']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='standard_inchi_key'])">
          <fn:null key="inchi_key"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="inchi_key">
            <xsl:value-of select="fn:string[@key='standard_inchi_key']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== MOLECULAR PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='molformula'])">
          <fn:null key="molecular_formula"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="molecular_formula">
            <xsl:value-of select="fn:string[@key='molformula']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="not(fn:number[@key='full_mwt']) or fn:number[@key='full_mwt'] &lt;= 0">
          <fn:null key="molecular_weight"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:number key="molecular_weight">
            <xsl:value-of select="fn:number[@key='full_mwt']"/>
          </fn:number>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DRUG TYPE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='molecule_type'])">
          <fn:null key="drug_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="drug_type">
            <xsl:value-of select="fn:string[@key='molecule_type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CLINICAL PHASE / APPROVAL STATUS ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='max_phase']">
          <fn:number key="max_phase">
            <xsl:value-of select="fn:number[@key='max_phase']"/>
          </fn:number>
          <fn:string key="approval_status">
            <xsl:choose>
              <xsl:when test="fn:number[@key='max_phase'] = 4">approved</xsl:when>
              <xsl:when test="fn:number[@key='max_phase'] = 3">Phase III</xsl:when>
              <xsl:when test="fn:number[@key='max_phase'] = 2">Phase II</xsl:when>
              <xsl:when test="fn:number[@key='max_phase'] = 1">Phase I</xsl:when>
              <xsl:otherwise>preclinical</xsl:otherwise>
            </xsl:choose>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="max_phase"/>
          <fn:null key="approval_status"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BIOACTIVITY METRICS ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='pchembl_value']">
          <fn:number key="pchembl_value">
            <xsl:value-of select="fn:number[@key='pchembl_value']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="pchembl_value"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='confidence_score']">
          <fn:number key="confidence_score">
            <xsl:value-of select="fn:number[@key='confidence_score']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="confidence_score"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== NATURAL PRODUCT FLAGS ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='natural_product'] = 1">
          <fn:boolean key="natural_product">true</fn:boolean>
        </xsl:when>
        <xsl:when test="fn:number[@key='natural_product'] = 0">
          <fn:boolean key="natural_product">false</fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="natural_product"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='np_likeness_score']">
          <fn:number key="np_likeness_score">
            <xsl:value-of select="fn:number[@key='np_likeness_score']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="np_likeness_score"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='qed_weighted']">
          <fn:number key="qed_weighted">
            <xsl:value-of select="fn:number[@key='qed_weighted']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="qed_weighted"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">ChEMBL</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='chembl_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
