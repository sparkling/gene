<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  SwissADME -> Drug Metabolism Unified Schema Mapping
  ============================================================
  Source: ./schema.md (SwissADME ADMET Prediction)
  Target: ../schema.json (2.5 Drug Metabolism and Pharmacokinetics Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ molecule_name               │ name                    │ Direct       │
  │ smiles                      │ smiles                  │ Direct       │
  │ inchikey                    │ inchi_key               │ Direct       │
  │ mw                          │ molecular_weight        │ Direct       │
  │ tpsa                        │ tpsa                    │ Direct       │
  │ hbd                         │ hbd                     │ Direct       │
  │ hba                         │ hba                     │ Direct       │
  │ rotatable_bonds             │ rotatable_bonds         │ Direct       │
  │ consensus_logp              │ consensus_logp          │ Direct       │
  │ lipinski_violations         │ lipinski_violations     │ Direct       │
  │ bioavailability_score       │ bioavailability_score   │ Direct       │
  │ gi_absorption               │ gi_absorption           │ Direct       │
  │ bbb_permeant                │ bbb_permeant            │ Direct       │
  │ pgp_substrate               │ pgp_substrate           │ Direct       │
  │ cyp*_inhibitor              │ cyp*_inhibitor          │ Direct       │
  │ pains_alerts                │ pains_alerts            │ Direct       │
  │ synthetic_accessibility     │ synthetic_accessibility │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null
  - Numeric values: undefined → null

  NOTES:
  - SwissADME provides ADMET predictions from structure
  - CYP inhibition predictions: Yes/No
  - Bioavailability score: 0-1 (Abbott score)
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
      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='molecule_name'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='molecule_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CHEMICAL STRUCTURE ========== -->

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
        <xsl:when test="local:is-null(fn:string[@key='inchikey'])">
          <fn:null key="inchi_key"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="inchi_key">
            <xsl:value-of select="fn:string[@key='inchikey']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PHYSICOCHEMICAL PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='mw']">
          <fn:number key="molecular_weight">
            <xsl:value-of select="fn:number[@key='mw']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="molecular_weight"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='tpsa']">
          <fn:number key="tpsa">
            <xsl:value-of select="fn:number[@key='tpsa']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="tpsa"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='hbd']">
          <fn:number key="hbd">
            <xsl:value-of select="fn:number[@key='hbd']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="hbd"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='hba']">
          <fn:number key="hba">
            <xsl:value-of select="fn:number[@key='hba']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="hba"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='rotatable_bonds']">
          <fn:number key="rotatable_bonds">
            <xsl:value-of select="fn:number[@key='rotatable_bonds']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="rotatable_bonds"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='consensus_logp']">
          <fn:number key="consensus_logp">
            <xsl:value-of select="fn:number[@key='consensus_logp']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="consensus_logp"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DRUG-LIKENESS ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='lipinski_violations']">
          <fn:number key="lipinski_violations">
            <xsl:value-of select="fn:number[@key='lipinski_violations']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="lipinski_violations"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='bioavailability_score']">
          <fn:number key="bioavailability_score">
            <xsl:value-of select="fn:number[@key='bioavailability_score']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="bioavailability_score"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PHARMACOKINETICS ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='gi_absorption'])">
          <fn:null key="gi_absorption"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gi_absorption">
            <xsl:value-of select="fn:string[@key='gi_absorption']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='bbb_permeant'])">
          <fn:null key="bbb_permeant"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="bbb_permeant">
            <xsl:value-of select="fn:string[@key='bbb_permeant']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='pgp_substrate'])">
          <fn:null key="pgp_substrate"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="pgp_substrate">
            <xsl:value-of select="fn:string[@key='pgp_substrate']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CYP INHIBITION PREDICTIONS ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='cyp1a2_inhibitor'])">
          <fn:null key="cyp1a2_inhibitor"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="cyp1a2_inhibitor">
            <xsl:value-of select="fn:string[@key='cyp1a2_inhibitor']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='cyp2c19_inhibitor'])">
          <fn:null key="cyp2c19_inhibitor"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="cyp2c19_inhibitor">
            <xsl:value-of select="fn:string[@key='cyp2c19_inhibitor']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='cyp2c9_inhibitor'])">
          <fn:null key="cyp2c9_inhibitor"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="cyp2c9_inhibitor">
            <xsl:value-of select="fn:string[@key='cyp2c9_inhibitor']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='cyp2d6_inhibitor'])">
          <fn:null key="cyp2d6_inhibitor"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="cyp2d6_inhibitor">
            <xsl:value-of select="fn:string[@key='cyp2d6_inhibitor']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='cyp3a4_inhibitor'])">
          <fn:null key="cyp3a4_inhibitor"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="cyp3a4_inhibitor">
            <xsl:value-of select="fn:string[@key='cyp3a4_inhibitor']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== MEDICINAL CHEMISTRY ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='pains_alerts']">
          <fn:number key="pains_alerts">
            <xsl:value-of select="fn:number[@key='pains_alerts']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="pains_alerts"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='synthetic_accessibility']">
          <fn:number key="synthetic_accessibility">
            <xsl:value-of select="fn:number[@key='synthetic_accessibility']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="synthetic_accessibility"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">SwissADME</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='inchikey']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
