<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  ImmunoBase -> Autoimmune/Inflammatory Diseases Unified Schema Mapping
  ============================================================
  Source: ./schema.md (ImmunoBase GWAS Database)
  Target: ../schema.json (3.6 Autoimmune/Inflammatory Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ gene_symbol                 │ gene_symbol             │ Direct       │
  │ disease                     │ disease_associations    │ Array        │
  │ disease_code                │ disease_code            │ Direct       │
  │ snp_id                      │ snp_ids                 │ Array        │
  │ region_id                   │ region_id               │ Direct       │
  │ credible_set                │ credible_set            │ Array        │
  │ gwas_loci                   │ gwas_loci               │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NOTES:
  - ImmunoBase covers 12 autoimmune/inflammatory diseases
  - Disease codes: T1D, RA, CEL, MS, CD, UC, etc.
  - Includes fine-mapped credible sets from ImmunoChip studies
  - SNP IDs follow dbSNP rsID format
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

      <!-- ========== DISEASE ASSOCIATIONS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='disease']/fn:string">
          <fn:array key="disease_associations">
            <xsl:for-each select="fn:array[@key='disease']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='disease'] and not(local:is-null(fn:string[@key='disease']))">
          <fn:array key="disease_associations">
            <fn:string><xsl:value-of select="fn:string[@key='disease']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="disease_associations"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DISEASE CODE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='disease_code'])">
          <fn:null key="disease_code"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="disease_code">
            <xsl:value-of select="fn:string[@key='disease_code']"/>
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

      <!-- ========== REGION ID ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='region_id'])">
          <fn:null key="region_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="region_id">
            <xsl:value-of select="fn:string[@key='region_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CREDIBLE SET ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='credible_set']/fn:string">
          <fn:array key="credible_set">
            <xsl:for-each select="fn:array[@key='credible_set']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="credible_set"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GWAS LOCI ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='gwas_loci']">
          <fn:number key="gwas_loci">
            <xsl:value-of select="fn:number[@key='gwas_loci']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gwas_loci"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">ImmunoBase</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="concat(fn:string[@key='gene_symbol'], '_', fn:string[@key='disease_code'])"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
