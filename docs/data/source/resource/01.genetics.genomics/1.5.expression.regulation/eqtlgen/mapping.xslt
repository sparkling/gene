<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  eQTLGen -> Expression & Regulation Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Expression & Regulation Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ Gene                        │ gene_id                 │ Direct       │
  │ GeneSymbol                  │ gene_symbol             │ Direct       │
  │ SNP                         │ eqtlgen_snp             │ Direct       │
  │ SNP                         │ variant_id              │ Direct       │
  │ SNPChr                      │ eqtlgen_snp_chr         │ Direct       │
  │ SNPChr                      │ chromosome              │ String       │
  │ SNPPos                      │ eqtlgen_snp_pos         │ Direct       │
  │ AssessedAllele              │ assessed_allele         │ Direct       │
  │ OtherAllele                 │ other_allele            │ Direct       │
  │ Zscore                      │ z_score                 │ Direct       │
  │ Pvalue                      │ p_value                 │ Direct       │
  │ GeneChr                     │ gene_chr                │ Direct       │
  │ GenePos                     │ gene_tss_position       │ Direct       │
  │ NrCohorts                   │ number_cohorts          │ Direct       │
  │ NrSamples                   │ sample_size             │ Direct       │
  │ FDR                         │ fdr                     │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "NA" → null

  NOTES:
  - eQTLGen: Blood-based cis-eQTL meta-analysis
  - 31,684 samples from 37 cohorts
  - Coordinates in hg19/GRCh37
  - Z-score derived from meta-analysis
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
      <fn:string key="gene_id">
        <xsl:value-of select="fn:string[@key='Gene']"/>
      </fn:string>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">gene_symbol</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='GeneSymbol']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">variant_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='SNP']"/>
      </xsl:call-template>

      <!-- Chromosome from SNPChr -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='SNPChr']">
          <fn:string key="chromosome">
            <xsl:value-of select="fn:number[@key='SNPChr']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="chromosome"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- P-value -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='Pvalue']">
          <fn:number key="p_value">
            <xsl:value-of select="fn:number[@key='Pvalue']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="p_value"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Effect size (Z-score as proxy) -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='Zscore']">
          <fn:number key="effect_size">
            <xsl:value-of select="fn:number[@key='Zscore']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="effect_size"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== EQTLGEN-SPECIFIC FIELDS ========== -->

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">eqtlgen_snp</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='SNP']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='SNPChr']">
          <fn:number key="eqtlgen_snp_chr">
            <xsl:value-of select="fn:number[@key='SNPChr']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="eqtlgen_snp_chr"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='SNPPos']">
          <fn:number key="eqtlgen_snp_pos">
            <xsl:value-of select="fn:number[@key='SNPPos']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="eqtlgen_snp_pos"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">assessed_allele</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='AssessedAllele']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">other_allele</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='OtherAllele']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='Zscore']">
          <fn:number key="z_score">
            <xsl:value-of select="fn:number[@key='Zscore']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="z_score"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='GeneChr']">
          <fn:number key="gene_chr">
            <xsl:value-of select="fn:number[@key='GeneChr']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_chr"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='GenePos']">
          <fn:number key="gene_tss_position">
            <xsl:value-of select="fn:number[@key='GenePos']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_tss_position"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='NrCohorts']">
          <fn:number key="number_cohorts">
            <xsl:value-of select="fn:number[@key='NrCohorts']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="number_cohorts"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='NrSamples']">
          <fn:number key="sample_size">
            <xsl:value-of select="fn:number[@key='NrSamples']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="sample_size"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='FDR']">
          <fn:number key="fdr">
            <xsl:value-of select="fn:number[@key='FDR']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="fdr"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">eQTLGen</fn:string>
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
