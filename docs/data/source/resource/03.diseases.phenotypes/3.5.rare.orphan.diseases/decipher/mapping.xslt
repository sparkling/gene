<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  DECIPHER -> Rare/Orphan Diseases Unified Schema Mapping
  ============================================================
  Source: ./schema.md (DECIPHER CNV/Variant Database)
  Target: ../schema.json (3.5 Rare/Orphan Diseases Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ decipher_id                 │ decipher_id             │ Direct       │
  │ decipher_id                 │ disease_id              │ Prefix       │
  │ syndrome_name               │ disease_name            │ Direct       │
  │ gene_symbol                 │ gene_symbols            │ Array        │
  │ hpo_id                      │ hpo_ids                 │ Array        │
  │ cnv_type                    │ cnv_type                │ Direct       │
  │ hi_score                    │ hi_score                │ Direct       │
  │ pLI                         │ pli                     │ Direct       │
  │ LOEUF                       │ loeuf                   │ Direct       │
  │ pathogenicity               │ pathogenicity           │ Direct       │
  │ ensembl_gene                │ ensembl_gene            │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NOTES:
  - DECIPHER provides CNV and sequence variant data for rare diseases
  - HI Score: Haploinsufficiency prediction (0-1)
  - pLI: Probability of loss-of-function intolerance (0-1)
  - LOEUF: LoF observed/expected upper bound fraction
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
      <!-- ========== IDENTIFIERS ========== -->

      <fn:string key="disease_id">
        <xsl:value-of select="concat('DECIPHER:', fn:number[@key='decipher_id'])"/>
      </fn:string>

      <xsl:choose>
        <xsl:when test="fn:number[@key='decipher_id']">
          <fn:number key="decipher_id">
            <xsl:value-of select="fn:number[@key='decipher_id']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="decipher_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DISEASE NAME ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='syndrome_name'])">
          <fn:null key="disease_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="disease_name">
            <xsl:value-of select="fn:string[@key='syndrome_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GENES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='gene_symbol']/fn:string">
          <fn:array key="gene_symbols">
            <xsl:for-each select="fn:array[@key='gene_symbol']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='gene_symbol'] and not(local:is-null(fn:string[@key='gene_symbol']))">
          <fn:array key="gene_symbols">
            <fn:string><xsl:value-of select="fn:string[@key='gene_symbol']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_symbols"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== HPO PHENOTYPES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='hpo_id']/fn:string">
          <fn:array key="hpo_ids">
            <xsl:for-each select="fn:array[@key='hpo_id']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='hpo_id'] and not(local:is-null(fn:string[@key='hpo_id']))">
          <fn:array key="hpo_ids">
            <fn:string><xsl:value-of select="fn:string[@key='hpo_id']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="hpo_ids"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CNV TYPE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='cnv_type'])">
          <fn:null key="cnv_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="cnv_type">
            <xsl:value-of select="fn:string[@key='cnv_type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CONSTRAINT SCORES ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='hi_score']">
          <fn:number key="hi_score">
            <xsl:value-of select="fn:number[@key='hi_score']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="hi_score"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='pLI']">
          <fn:number key="pli">
            <xsl:value-of select="fn:number[@key='pLI']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="pli"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='LOEUF']">
          <fn:number key="loeuf">
            <xsl:value-of select="fn:number[@key='LOEUF']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="loeuf"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PATHOGENICITY ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='pathogenicity'])">
          <fn:null key="pathogenicity"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="pathogenicity">
            <xsl:value-of select="fn:string[@key='pathogenicity']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ENSEMBL GENE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='ensembl_gene'])">
          <fn:null key="ensembl_gene"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="ensembl_gene">
            <xsl:value-of select="fn:string[@key='ensembl_gene']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">DECIPHER</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:number[@key='decipher_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
