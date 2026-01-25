<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  DisGeNET -> Disease-Gene Associations Unified Schema Mapping
  ============================================================
  Source: ./schema.md (DisGeNET Disease-Gene Association Database)
  Target: ../schema.json (3.3 Disease-Gene Associations Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ geneId                      │ gene_id                 │ Direct       │
  │ geneSymbol                  │ gene_symbol             │ Direct       │
  │ diseaseId                   │ disease_id              │ Direct       │
  │ diseaseName                 │ disease_name            │ Direct       │
  │ score                       │ score                   │ Direct       │
  │ EI                          │ ei                      │ Direct       │
  │ DSI                         │ dsi                     │ Direct       │
  │ DPI                         │ dpi                     │ Direct       │
  │ associationType             │ association_type        │ Direct       │
  │ snpId                       │ snp_id                  │ Direct       │
  │ NofPmids                    │ number_of_pmids         │ Direct       │
  │ pmid                        │ pmids                   │ Array        │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NOTES:
  - DisGeNET aggregates disease-gene associations from multiple sources
  - Score: 0-1 confidence based on evidence
  - EI: Evidence Index (consistency across publications)
  - DSI: Disease Specificity Index
  - DPI: Disease Pleiotropy Index
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

      <fn:string key="gene_id">
        <xsl:value-of select="fn:string[@key='geneId']"/>
      </fn:string>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='geneSymbol'])">
          <fn:null key="gene_symbol"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gene_symbol">
            <xsl:value-of select="fn:string[@key='geneSymbol']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DISEASE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='diseaseId'])">
          <fn:null key="disease_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="disease_id">
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
        <xsl:when test="fn:number[@key='EI']">
          <fn:number key="ei">
            <xsl:value-of select="fn:number[@key='EI']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ei"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='DSI']">
          <fn:number key="dsi">
            <xsl:value-of select="fn:number[@key='DSI']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="dsi"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='DPI']">
          <fn:number key="dpi">
            <xsl:value-of select="fn:number[@key='DPI']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="dpi"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ASSOCIATION TYPE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='associationType'])">
          <fn:null key="association_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="association_type">
            <xsl:value-of select="fn:string[@key='associationType']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== VARIANT (for VDAs) ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='snpId'])">
          <fn:null key="snp_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="snp_id">
            <xsl:value-of select="fn:string[@key='snpId']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== EVIDENCE ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='NofPmids']">
          <fn:number key="number_of_pmids">
            <xsl:value-of select="fn:number[@key='NofPmids']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="number_of_pmids"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:array[@key='pmid']/fn:number">
          <fn:array key="pmids">
            <xsl:for-each select="fn:array[@key='pmid']/fn:number">
              <fn:number><xsl:value-of select="."/></fn:number>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="pmids"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">DisGeNET</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="concat(fn:string[@key='geneId'], '_', fn:string[@key='diseaseId'])"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
