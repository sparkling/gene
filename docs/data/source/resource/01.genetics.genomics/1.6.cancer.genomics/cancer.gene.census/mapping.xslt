<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Cancer Gene Census -> Cancer Genomics Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Cancer Genomics Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ Gene Symbol                 │ gene                    │ Direct       │
  │ Name                        │ cgc_gene_name           │ Direct       │
  │ Entrez GeneId               │ entrez_gene_id          │ Direct       │
  │ Genome Location             │ cgc_genome_location     │ Direct       │
  │ Tier                        │ cgc_tier                │ Direct       │
  │ Hallmark                    │ cgc_hallmark            │ Direct       │
  │ Role in Cancer              │ role_in_cancer          │ Direct       │
  │ Mutation Types              │ mutation_types          │ Direct       │
  │ Translocation Partner       │ translocation_partners  │ Direct       │
  │ Tumour Types (Somatic)      │ cancer_types            │ Array        │
  │ Synonyms                    │ cgc_synonyms            │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings → null

  NOTES:
  - Cancer Gene Census: COSMIC's curated driver gene list
  - Tier 1: Strong evidence (574 genes)
  - Tier 2: Evidence consistent with role (162 genes)
  - Mutation types: Mis (missense), N (nonsense), F (frameshift), etc.
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
      <fn:string key="gene">
        <xsl:value-of select="fn:string[@key='Gene Symbol']"/>
      </fn:string>

      <!-- Cancer types from somatic tumour types -->
      <xsl:variable name="cancers" select="fn:string[@key='Tumour Types(Somatic)']"/>
      <xsl:choose>
        <xsl:when test="local:is-null($cancers)">
          <fn:array key="cancer_types"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="cancer_types">
            <xsl:for-each select="tokenize($cancers, ',')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CGC-SPECIFIC FIELDS ========== -->

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cgc_gene_name</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Name']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='Entrez GeneId']">
          <fn:number key="entrez_gene_id">
            <xsl:value-of select="fn:number[@key='Entrez GeneId']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="entrez_gene_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cgc_genome_location</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Genome Location']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='Tier']">
          <fn:number key="cgc_tier">
            <xsl:value-of select="fn:number[@key='Tier']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cgc_tier"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cgc_hallmark</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Hallmark']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">role_in_cancer</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Role in Cancer']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">mutation_types</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Mutation Types']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">translocation_partners</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Translocation Partner']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cgc_synonyms</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Synonyms']"/>
      </xsl:call-template>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">Cancer Gene Census</fn:string>
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
