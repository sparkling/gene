<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  SpliceAI -> Functional Prediction Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Functional Prediction Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ CHROM                       │ chromosome              │ Normalize    │
  │ POS                         │ position                │ Direct       │
  │ REF                         │ reference_allele        │ Direct       │
  │ ALT                         │ alternate_allele        │ Direct       │
  │ ALLELE                      │ spliceai_allele         │ Direct       │
  │ SYMBOL                      │ spliceai_gene_symbol    │ Direct       │
  │ SYMBOL                      │ gene_symbols            │ Array        │
  │ DS_AG                       │ spliceai_ds_ag          │ Direct       │
  │ DS_AL                       │ spliceai_ds_al          │ Direct       │
  │ DS_DG                       │ spliceai_ds_dg          │ Direct       │
  │ DS_DL                       │ spliceai_ds_dl          │ Direct       │
  │ DP_AG                       │ spliceai_dp_ag          │ Direct       │
  │ DP_AL                       │ spliceai_dp_al          │ Direct       │
  │ DP_DG                       │ spliceai_dp_dg          │ Direct       │
  │ DP_DL                       │ spliceai_dp_dl          │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "." → null
  - Delta scores of 0 retained (valid score)

  NOTES:
  - SpliceAI predicts splice-altering effects
  - DS (Delta Score) = probability of splice change
  - DP (Delta Position) = distance to splice site
  - AG=acceptor gain, AL=acceptor loss, DG=donor gain, DL=donor loss
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

  <!-- ============================================================
       Entry Point: Parse JSON input
       ============================================================ -->
  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <!-- ============================================================
       Main Record Transformation
       ============================================================ -->
  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== CORE VARIANT FIELDS ========== -->
      <fn:string key="chromosome">
        <xsl:value-of select="local:normalize-chromosome(fn:string[@key='CHROM'])"/>
      </fn:string>

      <fn:number key="position">
        <xsl:value-of select="fn:number[@key='POS']"/>
      </fn:number>

      <fn:string key="reference_allele">
        <xsl:value-of select="fn:string[@key='REF']"/>
      </fn:string>

      <fn:string key="alternate_allele">
        <xsl:value-of select="fn:string[@key='ALT']"/>
      </fn:string>

      <!-- ========== GENE SYMBOLS (ARRAY) ========== -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='SYMBOL'] and not(local:is-null(fn:string[@key='SYMBOL']))">
          <fn:array key="gene_symbols">
            <fn:string><xsl:value-of select="fn:string[@key='SYMBOL']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_symbols"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SPLICEAI-SPECIFIC FIELDS ========== -->

      <!-- Allele -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='ALLELE'])">
          <fn:null key="spliceai_allele"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="spliceai_allele">
            <xsl:value-of select="fn:string[@key='ALLELE']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Gene symbol -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='SYMBOL'])">
          <fn:null key="spliceai_gene_symbol"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="spliceai_gene_symbol">
            <xsl:value-of select="fn:string[@key='SYMBOL']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Delta Score - Acceptor Gain -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='DS_AG']">
          <fn:number key="spliceai_ds_ag">
            <xsl:value-of select="fn:number[@key='DS_AG']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="spliceai_ds_ag"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Delta Score - Acceptor Loss -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='DS_AL']">
          <fn:number key="spliceai_ds_al">
            <xsl:value-of select="fn:number[@key='DS_AL']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="spliceai_ds_al"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Delta Score - Donor Gain -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='DS_DG']">
          <fn:number key="spliceai_ds_dg">
            <xsl:value-of select="fn:number[@key='DS_DG']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="spliceai_ds_dg"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Delta Score - Donor Loss -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='DS_DL']">
          <fn:number key="spliceai_ds_dl">
            <xsl:value-of select="fn:number[@key='DS_DL']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="spliceai_ds_dl"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Delta Position - Acceptor Gain -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='DP_AG']">
          <fn:number key="spliceai_dp_ag">
            <xsl:value-of select="fn:number[@key='DP_AG']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="spliceai_dp_ag"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Delta Position - Acceptor Loss -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='DP_AL']">
          <fn:number key="spliceai_dp_al">
            <xsl:value-of select="fn:number[@key='DP_AL']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="spliceai_dp_al"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Delta Position - Donor Gain -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='DP_DG']">
          <fn:number key="spliceai_dp_dg">
            <xsl:value-of select="fn:number[@key='DP_DG']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="spliceai_dp_dg"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Delta Position - Donor Loss -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='DP_DL']">
          <fn:number key="spliceai_dp_dl">
            <xsl:value-of select="fn:number[@key='DP_DL']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="spliceai_dp_dl"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">SpliceAI</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Reusable Functions
       ============================================================ -->

  <!-- Normalize chromosome: remove "chr" prefix -->
  <xsl:function name="local:normalize-chromosome" as="xs:string">
    <xsl:param name="chr" as="xs:string"/>
    <xsl:value-of select="replace($chr, '^chr', '')"/>
  </xsl:function>

  <!-- Check if value is null/empty -->
  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = ('', '-', '.', 'NA', 'N/A', 'null')"/>
  </xsl:function>

</xsl:stylesheet>
