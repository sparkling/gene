<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  AlphaMissense -> Functional Prediction Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Functional Prediction Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────────┬──────────────┐
  │ Source Field                │ Target Field                │ Transform    │
  ├─────────────────────────────┼─────────────────────────────┼──────────────┤
  │ CHROM                       │ chromosome                  │ Normalize    │
  │ POS                         │ position                    │ Direct       │
  │ REF                         │ reference_allele            │ Direct       │
  │ ALT                         │ alternate_allele            │ Direct       │
  │ genome                      │ alphamissense_genome        │ Direct       │
  │ uniprot_id                  │ alphamissense_uniprot_id    │ Direct       │
  │ transcript_id               │ transcript_ids              │ Array        │
  │ protein_variant             │ alphamissense_protein_variant│ Direct      │
  │ am_pathogenicity            │ alphamissense_pathogenicity │ Direct       │
  │ am_class                    │ alphamissense_class         │ Direct       │
  │ mean_am_pathogenicity       │ alphamissense_gene_mean     │ Direct       │
  │ gene                        │ gene_symbols                │ Array        │
  └─────────────────────────────┴─────────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "." → null
  - Pathogenicity score NaN → null

  NOTES:
  - AlphaMissense provides missense variant pathogenicity predictions
  - am_class values: likely_benign, ambiguous, likely_pathogenic
  - Pathogenicity scores range 0-1
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
        <xsl:when test="fn:string[@key='gene'] and not(local:is-null(fn:string[@key='gene']))">
          <fn:array key="gene_symbols">
            <fn:string><xsl:value-of select="fn:string[@key='gene']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_symbols"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== TRANSCRIPT IDS (ARRAY) ========== -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='transcript_id'] and not(local:is-null(fn:string[@key='transcript_id']))">
          <fn:array key="transcript_ids">
            <fn:string><xsl:value-of select="fn:string[@key='transcript_id']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="transcript_ids"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ALPHAMISSENSE-SPECIFIC FIELDS ========== -->

      <!-- Genome version -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='genome'])">
          <fn:null key="alphamissense_genome"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="alphamissense_genome">
            <xsl:value-of select="fn:string[@key='genome']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- UniProt ID -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='uniprot_id'])">
          <fn:null key="alphamissense_uniprot_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="alphamissense_uniprot_id">
            <xsl:value-of select="fn:string[@key='uniprot_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Protein variant -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='protein_variant'])">
          <fn:null key="alphamissense_protein_variant"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="alphamissense_protein_variant">
            <xsl:value-of select="fn:string[@key='protein_variant']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Pathogenicity score -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='am_pathogenicity']">
          <fn:number key="alphamissense_pathogenicity">
            <xsl:value-of select="fn:number[@key='am_pathogenicity']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="alphamissense_pathogenicity"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Classification -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='am_class'])">
          <fn:null key="alphamissense_class"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="alphamissense_class">
            <xsl:value-of select="fn:string[@key='am_class']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Gene mean pathogenicity -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='mean_am_pathogenicity']">
          <fn:number key="alphamissense_gene_mean">
            <xsl:value-of select="fn:number[@key='mean_am_pathogenicity']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="alphamissense_gene_mean"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">AlphaMissense</fn:string>
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
    <xsl:sequence select="not($value) or $value = ('', '-', '.', 'NA', 'N/A', 'null', 'NaN')"/>
  </xsl:function>

</xsl:stylesheet>
