<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  SuperCYP -> Drug Metabolism Unified Schema Mapping
  ============================================================
  Source: ./schema.md (SuperCYP CYP450 Database)
  Target: ../schema.json (2.5 Drug Metabolism and Pharmacokinetics Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ drug_id                     │ compound_id             │ Direct       │
  │ drug_name                   │ name                    │ Direct       │
  │ smiles                      │ smiles                  │ Direct       │
  │ molecular_weight            │ molecular_weight        │ Direct       │
  │ cyp_id                      │ cyp_id                  │ Direct       │
  │ interaction_type            │ interaction_type        │ Direct       │
  │ strength                    │ interaction_strength    │ Direct       │
  │ ki_value                    │ ki_value                │ Direct       │
  │ ic50_value                  │ ic50_value              │ Direct       │
  │ induction_fold              │ induction_fold          │ Direct       │
  │ evidence_level              │ evidence_level          │ Direct       │
  │ drugbank_id                 │ drugbank_id             │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null
  - Numeric values: -1 or undefined → null

  NOTES:
  - SuperCYP focuses on CYP450 drug-drug interactions
  - CYP isoforms: CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4, etc.
  - Interaction types: substrate, inhibitor, inducer
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

      <fn:string key="compound_id">
        <xsl:value-of select="fn:string[@key='drug_id']"/>
      </fn:string>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='drugbank_id'])">
          <fn:null key="drugbank_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="drugbank_id">
            <xsl:value-of select="fn:string[@key='drugbank_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='drug_name'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='drug_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='smiles'])">
          <fn:null key="smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="smiles">
            <xsl:value-of select="fn:string[@key='smiles']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="not(fn:number[@key='molecular_weight']) or fn:number[@key='molecular_weight'] &lt;= 0">
          <fn:null key="molecular_weight"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:number key="molecular_weight">
            <xsl:value-of select="fn:number[@key='molecular_weight']"/>
          </fn:number>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CYP INTERACTION DATA ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='cyp_id'])">
          <fn:null key="cyp_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="cyp_id">
            <xsl:value-of select="fn:string[@key='cyp_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='interaction_type'])">
          <fn:null key="interaction_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="interaction_type">
            <xsl:value-of select="fn:string[@key='interaction_type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='strength'])">
          <fn:null key="interaction_strength"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="interaction_strength">
            <xsl:value-of select="fn:string[@key='strength']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== KINETIC PARAMETERS ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='ki_value'] and fn:number[@key='ki_value'] > 0">
          <fn:number key="ki_value">
            <xsl:value-of select="fn:number[@key='ki_value']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ki_value"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='ic50_value'] and fn:number[@key='ic50_value'] > 0">
          <fn:number key="ic50_value">
            <xsl:value-of select="fn:number[@key='ic50_value']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ic50_value"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='induction_fold'] and fn:number[@key='induction_fold'] > 0">
          <fn:number key="induction_fold">
            <xsl:value-of select="fn:number[@key='induction_fold']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="induction_fold"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== EVIDENCE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='evidence_level'])">
          <fn:null key="evidence_level"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="evidence_level">
            <xsl:value-of select="fn:string[@key='evidence_level']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">SuperCYP</fn:string>
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
