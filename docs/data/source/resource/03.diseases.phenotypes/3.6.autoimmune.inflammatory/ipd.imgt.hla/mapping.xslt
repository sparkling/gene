<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  IPD-IMGT/HLA -> Autoimmune/Inflammatory Diseases Unified Schema Mapping
  ============================================================
  Source: ./schema.md (IPD-IMGT/HLA Database)
  Target: ../schema.json (3.6 Autoimmune/Inflammatory Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ locus                       │ gene_symbol             │ Direct       │
  │ allele                      │ hla_allele              │ Direct       │
  │ allele_group                │ allele_group            │ Direct       │
  │ accession                   │ ipd_accession           │ Direct       │
  │ gene_locus                  │ gene_locus              │ Direct       │
  │ sequence_type               │ sequence_type           │ Direct       │
  │ disease_association         │ disease_associations    │ Array        │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NOTES:
  - IPD-IMGT/HLA is the international reference for HLA allele sequences
  - HLA allele names follow WHO nomenclature: HLA-A*02:01:01:01
  - 30,000+ alleles across HLA-A, B, C, DRB1, DQB1, DPB1
  - Critical for transplantation matching and disease association studies
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
      <!-- ========== GENE SYMBOL (LOCUS) ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='locus'])">
          <fn:null key="gene_symbol"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gene_symbol">
            <xsl:value-of select="fn:string[@key='locus']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== HLA ALLELE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='allele'])">
          <fn:null key="hla_allele"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="hla_allele">
            <xsl:value-of select="fn:string[@key='allele']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ALLELE GROUP ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='allele_group'])">
          <fn:null key="allele_group"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="allele_group">
            <xsl:value-of select="fn:string[@key='allele_group']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== IPD ACCESSION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='accession'])">
          <fn:null key="ipd_accession"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="ipd_accession">
            <xsl:value-of select="fn:string[@key='accession']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GENE LOCUS ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='gene_locus'])">
          <fn:null key="gene_locus"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gene_locus">
            <xsl:value-of select="fn:string[@key='gene_locus']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SEQUENCE TYPE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='sequence_type'])">
          <fn:null key="sequence_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="sequence_type">
            <xsl:value-of select="fn:string[@key='sequence_type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DISEASE ASSOCIATIONS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='disease_association']/fn:string">
          <fn:array key="disease_associations">
            <xsl:for-each select="fn:array[@key='disease_association']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='disease_association'] and not(local:is-null(fn:string[@key='disease_association']))">
          <fn:array key="disease_associations">
            <fn:string><xsl:value-of select="fn:string[@key='disease_association']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="disease_associations"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">IPD-IMGT/HLA</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='accession']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
