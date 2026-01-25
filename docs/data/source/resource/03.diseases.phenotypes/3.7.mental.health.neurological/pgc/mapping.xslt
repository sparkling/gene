<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  PGC -> Mental Health/Neurological Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Psychiatric Genomics Consortium)
  Target: ../schema.json (3.7 Mental Health/Neurological Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ gene_symbol                 │ gene_symbol             │ Direct       │
  │ gwas_study_id               │ gwas_study_id           │ Direct       │
  │ disorder                    │ disorders               │ Array        │
  │ disorder_code               │ disorder_code           │ Direct       │
  │ snp_id                      │ snp_ids                 │ Array        │
  │ p_value                     │ p_value                 │ Direct       │
  │ odds_ratio                  │ odds_ratio              │ Direct       │
  │ sample_size                 │ sample_size             │ Direct       │
  │ loci_count                  │ loci_count              │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NOTES:
  - PGC is the largest psychiatric GWAS consortium
  - Disorder codes: SCZ, BIP, MDD, ADHD, ASD, etc.
  - Studies include 4M+ samples across 800+ research groups
  - P-values and odds ratios from meta-analyses
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

      <!-- ========== STUDY IDENTIFIER ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='gwas_study_id'])">
          <fn:null key="gwas_study_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gwas_study_id">
            <xsl:value-of select="fn:string[@key='gwas_study_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DISORDERS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='disorder']/fn:string">
          <fn:array key="disorders">
            <xsl:for-each select="fn:array[@key='disorder']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='disorder'] and not(local:is-null(fn:string[@key='disorder']))">
          <fn:array key="disorders">
            <fn:string><xsl:value-of select="fn:string[@key='disorder']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="disorders"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DISORDER CODE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='disorder_code'])">
          <fn:null key="disorder_code"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="disorder_code">
            <xsl:value-of select="fn:string[@key='disorder_code']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SNP IDs ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='snp_id']/fn:string">
          <fn:array key="snp_ids">
            <xsl:for-each select="fn:array[@key='snp_id']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='snp_id'] and not(local:is-null(fn:string[@key='snp_id']))">
          <fn:array key="snp_ids">
            <fn:string><xsl:value-of select="fn:string[@key='snp_id']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="snp_ids"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== STATISTICS ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='p_value']">
          <fn:number key="p_value">
            <xsl:value-of select="fn:number[@key='p_value']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="p_value"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='odds_ratio']">
          <fn:number key="odds_ratio">
            <xsl:value-of select="fn:number[@key='odds_ratio']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="odds_ratio"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='sample_size']">
          <fn:number key="sample_size">
            <xsl:value-of select="fn:number[@key='sample_size']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="sample_size"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='loci_count']">
          <fn:number key="loci_count">
            <xsl:value-of select="fn:number[@key='loci_count']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="loci_count"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">PGC</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='gwas_study_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
