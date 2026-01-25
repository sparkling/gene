<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Open Targets -> Disease-Gene Associations Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Open Targets Platform)
  Target: ../schema.json (3.3 Disease-Gene Associations Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ targetId                    │ target_id               │ Direct       │
  │ approvedSymbol              │ gene_symbol             │ Direct       │
  │ diseaseId                   │ ot_disease_id           │ Direct       │
  │ diseaseName                 │ disease_name            │ Direct       │
  │ score                       │ score                   │ Direct       │
  │ datatypeScores              │ datatype_scores         │ Object       │
  │ evidenceCount               │ evidence_count          │ Direct       │
  │ tractability                │ tractability            │ Object       │
  │ safety                      │ safety                  │ Object       │
  │ mechanismOfAction           │ mechanism_of_action     │ Direct       │
  │ clinicalPhase               │ clinical_phase          │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NOTES:
  - Open Targets uses EFO IDs for diseases, Ensembl IDs for targets
  - Score: 0-1 overall association score
  - Datatype scores: genetic, drug, pathway, etc.
  - Tractability: druggability assessment
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
      <!-- ========== TARGET/GENE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='targetId'])">
          <fn:null key="target_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="target_id">
            <xsl:value-of select="fn:string[@key='targetId']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='approvedSymbol'])">
          <fn:null key="gene_symbol"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gene_symbol">
            <xsl:value-of select="fn:string[@key='approvedSymbol']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DISEASE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='diseaseId'])">
          <fn:null key="ot_disease_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="ot_disease_id">
            <xsl:value-of select="fn:string[@key='diseaseId']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='diseaseName'])">
          <fn:null key="disease_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="disease_name">
            <xsl:value-of select="fn:string[@key='diseaseName']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SCORES ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='score']">
          <fn:number key="score">
            <xsl:value-of select="fn:number[@key='score']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="score"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:map[@key='datatypeScores']">
          <xsl:copy-of select="fn:map[@key='datatypeScores']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="datatype_scores"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='evidenceCount']">
          <fn:number key="evidence_count">
            <xsl:value-of select="fn:number[@key='evidenceCount']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="evidence_count"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DRUGGABILITY ========== -->

      <xsl:choose>
        <xsl:when test="fn:map[@key='tractability']">
          <xsl:copy-of select="fn:map[@key='tractability']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="tractability"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:map[@key='safety']">
          <xsl:copy-of select="fn:map[@key='safety']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="safety"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DRUG INFORMATION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='mechanismOfAction'])">
          <fn:null key="mechanism_of_action"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="mechanism_of_action">
            <xsl:value-of select="fn:string[@key='mechanismOfAction']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='clinicalPhase']">
          <fn:number key="clinical_phase">
            <xsl:value-of select="fn:number[@key='clinicalPhase']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="clinical_phase"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">Open Targets</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="concat(fn:string[@key='targetId'], '_', fn:string[@key='diseaseId'])"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
