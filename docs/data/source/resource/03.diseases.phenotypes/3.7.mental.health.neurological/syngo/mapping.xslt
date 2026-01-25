<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  SynGO -> Mental Health/Neurological Unified Schema Mapping
  ============================================================
  Source: ./schema.md (SynGO Synaptic Gene Ontology)
  Target: ../schema.json (3.7 Mental Health/Neurological Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ gene_symbol                 │ gene_symbol             │ Direct       │
  │ syngo_term                  │ syngo_term              │ Direct       │
  │ go_term                     │ go_term                 │ Direct       │
  │ uniprot_id                  │ uniprot_id              │ Direct       │
  │ synaptic_location           │ synaptic_location       │ Direct       │
  │ synaptic_function           │ synaptic_function       │ Direct       │
  │ evidence_code               │ evidence_code           │ Direct       │
  │ supporting_pmid             │ supporting_pmids        │ Array        │
  │ brain_region                │ brain_regions           │ Array        │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NOTES:
  - SynGO provides expert-curated synaptic gene annotations
  - 1,200+ annotated genes with 5,800+ annotations
  - Ontology structure: Cellular Component and Biological Process branches
  - Evidence codes follow ECO (Evidence & Conclusion Ontology)
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
      <!-- ========== GENE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='gene_symbol'])">
          <fn:null key="gene_symbol"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gene_symbol">
            <xsl:value-of select="fn:string[@key='gene_symbol']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ONTOLOGY TERMS ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='syngo_term'])">
          <fn:null key="syngo_term"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="syngo_term">
            <xsl:value-of select="fn:string[@key='syngo_term']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='go_term'])">
          <fn:null key="go_term"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="go_term">
            <xsl:value-of select="fn:string[@key='go_term']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PROTEIN IDENTIFIER ========== -->

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

      <!-- ========== SYNAPTIC LOCALIZATION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='synaptic_location'])">
          <fn:null key="synaptic_location"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="synaptic_location">
            <xsl:value-of select="fn:string[@key='synaptic_location']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SYNAPTIC FUNCTION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='synaptic_function'])">
          <fn:null key="synaptic_function"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="synaptic_function">
            <xsl:value-of select="fn:string[@key='synaptic_function']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== EVIDENCE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='evidence_code'])">
          <fn:null key="evidence_code"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="evidence_code">
            <xsl:value-of select="fn:string[@key='evidence_code']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SUPPORTING PMIDS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='supporting_pmid']/fn:number">
          <fn:array key="supporting_pmids">
            <xsl:for-each select="fn:array[@key='supporting_pmid']/fn:number">
              <fn:number><xsl:value-of select="."/></fn:number>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:number[@key='supporting_pmid']">
          <fn:array key="supporting_pmids">
            <fn:number><xsl:value-of select="fn:number[@key='supporting_pmid']"/></fn:number>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="supporting_pmids"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BRAIN REGIONS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='brain_region']/fn:string">
          <fn:array key="brain_regions">
            <xsl:for-each select="fn:array[@key='brain_region']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='brain_region'] and not(local:is-null(fn:string[@key='brain_region']))">
          <fn:array key="brain_regions">
            <fn:string><xsl:value-of select="fn:string[@key='brain_region']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="brain_regions"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">SynGO</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="concat(fn:string[@key='gene_symbol'], '_', fn:string[@key='syngo_term'])"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
