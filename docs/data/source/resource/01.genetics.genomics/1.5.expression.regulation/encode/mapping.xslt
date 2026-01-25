<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  ENCODE -> Expression & Regulation Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Expression & Regulation Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────────┬──────────────┐
  │ Source Field                │ Target Field                │ Transform    │
  ├─────────────────────────────┼─────────────────────────────┼──────────────┤
  │ gene_id                     │ gene_id                     │ Direct       │
  │ accession                   │ encode_accession            │ Direct       │
  │ assay_title                 │ encode_assay_title          │ Direct       │
  │ target                      │ encode_target               │ Object       │
  │ biosample_ontology          │ biosample_ontology          │ Object       │
  │ file_format                 │ file_format                 │ Direct       │
  │ output_type                 │ output_type                 │ Direct       │
  │ assembly                    │ assembly                    │ Direct       │
  │ cCRE_accession              │ ccre_accession              │ Direct       │
  │ cCRE_classification         │ ccre_classification         │ Direct       │
  │ chrom                       │ chromosome                  │ Direct       │
  │ chromStart                  │ start_position              │ Direct       │
  │ chromEnd                    │ end_position                │ Direct       │
  │ score                       │ encode_score                │ Direct       │
  └─────────────────────────────┴─────────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings → null

  NOTES:
  - ENCODE: Encyclopedia of DNA Elements
  - Includes ChIP-seq, ATAC-seq, RNA-seq, DNase-seq
  - cCRE = candidate cis-regulatory elements
  - Classifications: PLS, pELS, dELS, DNase-H3K4me3, CTCF-only
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
        <xsl:value-of select="(fn:string[@key='gene_id'], fn:string[@key='target.label'])[1]"/>
      </fn:string>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">chromosome</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='chrom']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='chromStart']">
          <fn:number key="start_position">
            <xsl:value-of select="fn:number[@key='chromStart']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="start_position"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='chromEnd']">
          <fn:number key="end_position">
            <xsl:value-of select="fn:number[@key='chromEnd']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="end_position"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ENCODE-SPECIFIC FIELDS ========== -->

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">encode_accession</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='accession']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">encode_assay_title</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='assay_title']"/>
      </xsl:call-template>

      <!-- Target (preserve as object) -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='target']">
          <xsl:copy-of select="fn:map[@key='target']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="encode_target"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Biosample ontology (preserve as object) -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='biosample_ontology']">
          <xsl:copy-of select="fn:map[@key='biosample_ontology']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="biosample_ontology"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">file_format</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='file_format']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">output_type</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='output_type']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">assembly</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='assembly']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">ccre_accession</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='cCRE_accession']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">ccre_classification</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='cCRE_classification']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='score']">
          <fn:number key="encode_score">
            <xsl:value-of select="fn:number[@key='score']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="encode_score"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">ENCODE</fn:string>
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
