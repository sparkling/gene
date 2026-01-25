<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  MSigDB -> Unified Gene Function Ontology Schema Mapping
  ============================================================
  Source: ./schema.json (MSigDB gene sets)
  Target: ../schema.json (Unified 4.5 Gene Function Ontology Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ geneSetName                 │ set_id                  │ Direct       │
  │ systematicName              │ msigdb_id               │ Direct       │
  │ geneSetCollection           │ collection              │ Direct       │
  │ geneSetSubcollection        │ subcollection           │ Direct       │
  │ description                 │ description             │ Direct       │
  │ geneSymbols                 │ gene_symbols            │ Direct       │
  │ entrezIds                   │ entrez_ids              │ Direct       │
  │ numGenes                    │ gene_count              │ Direct       │
  │ pmid                        │ pubmed_id               │ Direct       │
  │ gseaResults                 │ msigdb_gsea_results     │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings "" -> null
  - null pmid -> null

  NOTES:
  - Gene set names use UPPERCASE_WITH_UNDERSCORES convention
  - Collections use single-letter codes (H, C1-C8)
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
      <!-- ========== DIRECT FIELD MAPPINGS ========== -->

      <!-- Gene set identifier -->
      <fn:string key="set_id">
        <xsl:value-of select="fn:string[@key='geneSetName']"/>
      </fn:string>

      <!-- MSigDB systematic name -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='systematicName'] and fn:string[@key='systematicName'] != ''">
          <fn:string key="msigdb_id">
            <xsl:value-of select="fn:string[@key='systematicName']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="msigdb_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Collection -->
      <fn:string key="collection">
        <xsl:value-of select="fn:string[@key='geneSetCollection']"/>
      </fn:string>

      <!-- Subcollection -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='geneSetSubcollection'] and fn:string[@key='geneSetSubcollection'] != ''">
          <fn:string key="subcollection">
            <xsl:value-of select="fn:string[@key='geneSetSubcollection']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="subcollection"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Description -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='description'] and fn:string[@key='description'] != ''">
          <fn:string key="description">
            <xsl:value-of select="fn:string[@key='description']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="description"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Gene count -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='numGenes']">
          <fn:number key="gene_count">
            <xsl:value-of select="fn:number[@key='numGenes']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_count"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- PubMed ID -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='pmid']">
          <fn:number key="pubmed_id">
            <xsl:value-of select="fn:number[@key='pmid']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="pubmed_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ARRAY TRANSFORMATIONS ========== -->

      <!-- Gene symbols -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='geneSymbols']">
          <fn:array key="gene_symbols">
            <xsl:for-each select="fn:array[@key='geneSymbols']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_symbols"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Entrez IDs -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='entrezIds']">
          <fn:array key="entrez_ids">
            <xsl:for-each select="fn:array[@key='entrezIds']/fn:number">
              <fn:number><xsl:value-of select="."/></fn:number>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="entrez_ids"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== MSIGDB-SPECIFIC FIELDS ========== -->

      <!-- GSEA results object -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='gseaResults']">
          <fn:map key="msigdb_gsea_results">
            <xsl:copy-of select="fn:map[@key='gseaResults']/*"/>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="msigdb_gsea_results"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">MSigDB</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='geneSetName']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Reusable Functions
       ============================================================ -->

  <!-- Convert empty/placeholder to empty string (for null handling) -->
  <xsl:function name="local:nullify" as="xs:string">
    <xsl:param name="value" as="xs:string"/>
    <xsl:choose>
      <xsl:when test="$value = ('-', '.', '', 'NA', 'N/A', 'null')"></xsl:when>
      <xsl:otherwise><xsl:value-of select="$value"/></xsl:otherwise>
    </xsl:choose>
  </xsl:function>

</xsl:stylesheet>
