<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  CADD -> Functional Prediction Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Functional Prediction Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ Chrom                       │ chromosome              │ Normalize    │
  │ Pos                         │ position                │ Direct       │
  │ Ref                         │ reference_allele        │ Direct       │
  │ Alt                         │ alternate_allele        │ Direct       │
  │ RawScore                    │ cadd_raw_score          │ Direct       │
  │ PHRED                       │ cadd_phred              │ Direct       │
  │ GeneName                    │ gene_symbols            │ Array        │
  │ GeneID                      │ ensembl_gene_id         │ Direct       │
  │ FeatureID                   │ ensembl_transcript_id   │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, ".", "NA" → null

  NOTES:
  - CADD (Combined Annotation Dependent Depletion)
  - RawScore is SVM output, PHRED is -10*log10(rank/total)
  - PHRED > 20 suggests top 1% most deleterious
  - PHRED > 30 suggests top 0.1% most deleterious
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
        <xsl:value-of select="local:normalize-chromosome(fn:string[@key='Chrom'])"/>
      </fn:string>

      <fn:number key="position">
        <xsl:value-of select="fn:number[@key='Pos']"/>
      </fn:number>

      <fn:string key="reference_allele">
        <xsl:value-of select="fn:string[@key='Ref']"/>
      </fn:string>

      <fn:string key="alternate_allele">
        <xsl:value-of select="fn:string[@key='Alt']"/>
      </fn:string>

      <!-- ========== GENE SYMBOLS (ARRAY) ========== -->
      <xsl:variable name="gene" select="fn:string[@key='GeneName']"/>
      <xsl:choose>
        <xsl:when test="local:is-null($gene)">
          <fn:null key="gene_symbols"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="gene_symbols">
            <xsl:for-each select="tokenize($gene, '[;,]')">
              <xsl:if test="not(local:is-null(.))">
                <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
              </xsl:if>
            </xsl:for-each>
          </fn:array>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== TRANSCRIPT IDS (ARRAY) ========== -->
      <xsl:variable name="feature" select="fn:string[@key='FeatureID']"/>
      <xsl:choose>
        <xsl:when test="local:is-null($feature)">
          <fn:null key="transcript_ids"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="transcript_ids">
            <xsl:for-each select="tokenize($feature, '[;,]')">
              <xsl:if test="not(local:is-null(.))">
                <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
              </xsl:if>
            </xsl:for-each>
          </fn:array>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CADD-SPECIFIC FIELDS ========== -->

      <!-- Raw score -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='RawScore']">
          <fn:number key="cadd_raw_score">
            <xsl:value-of select="fn:number[@key='RawScore']"/>
          </fn:number>
        </xsl:when>
        <xsl:when test="fn:string[@key='RawScore'] and not(local:is-null(fn:string[@key='RawScore']))">
          <fn:number key="cadd_raw_score">
            <xsl:value-of select="number(fn:string[@key='RawScore'])"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cadd_raw_score"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- PHRED score -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='PHRED']">
          <fn:number key="cadd_phred">
            <xsl:value-of select="fn:number[@key='PHRED']"/>
          </fn:number>
        </xsl:when>
        <xsl:when test="fn:string[@key='PHRED'] and not(local:is-null(fn:string[@key='PHRED']))">
          <fn:number key="cadd_phred">
            <xsl:value-of select="number(fn:string[@key='PHRED'])"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cadd_phred"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Ensembl gene ID -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='GeneID'])">
          <fn:null key="ensembl_gene_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="ensembl_gene_id">
            <xsl:value-of select="fn:string[@key='GeneID']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Ensembl transcript ID (from FeatureID, take first) -->
      <xsl:variable name="transcript" select="tokenize(fn:string[@key='FeatureID'], ';')[1]"/>
      <xsl:choose>
        <xsl:when test="local:is-null($transcript)">
          <fn:null key="ensembl_transcript_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="ensembl_transcript_id">
            <xsl:value-of select="normalize-space($transcript)"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">CADD</fn:string>
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
    <xsl:sequence select="not($value) or normalize-space($value) = ('', '-', '.', 'NA', 'N/A', 'null')"/>
  </xsl:function>

</xsl:stylesheet>
