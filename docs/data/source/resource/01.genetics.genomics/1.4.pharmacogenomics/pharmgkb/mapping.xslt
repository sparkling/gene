<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  PharmGKB -> Pharmacogenomics Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Pharmacogenomics Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────────┬──────────────┐
  │ Source Field                │ Target Field                │ Transform    │
  ├─────────────────────────────┼─────────────────────────────┼──────────────┤
  │ Symbol                      │ gene                        │ Direct       │
  │ Drug/Name                   │ drugs                       │ Array        │
  │ Phenotype                   │ phenotype                   │ Direct       │
  │ Level_of_Evidence           │ evidence_level              │ Direct       │
  │ PharmGKB_Accession          │ pharmgkb_accession          │ Direct       │
  │ Is_VIP                      │ is_vip_gene                 │ Direct       │
  │ Has_Variant_Annotation      │ has_variant_annotation      │ Direct       │
  │ Annotation_ID               │ annotation_id               │ Direct       │
  │ Phenotype_Category          │ phenotype_category          │ Direct       │
  │ Trade_Names                 │ trade_names                 │ Array        │
  │ ATC_Codes                   │ atc_codes                   │ Array        │
  │ RxNorm_Identifiers          │ rxnorm_ids                  │ Array        │
  │ DrugBank_ID                 │ drugbank_id                 │ Direct       │
  │ Has_PGx_On_Label            │ has_pgx_on_label            │ Direct       │
  │ Variant_Haplotypes          │ variant_haplotypes          │ Array        │
  │ Biogeographical_Groups      │ population_frequencies      │ Object       │
  │ Pathway_ID                  │ pathway_id                  │ Direct       │
  │ Type (PK/PD)                │ pk_pd_type                  │ Direct       │
  │ RS_ID                       │ rsid                        │ Direct       │
  └─────────────────────────────┴─────────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "-" → null

  NOTES:
  - PharmGKB: Pharmacogenomics Knowledge Base
  - Evidence levels: 1A (CPIC guideline), 1B (FDA label), 2A, 2B, 3, 4
  - PA accession format: PA + integer (e.g., PA124)
  - VIP = Very Important Pharmacogene
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
      <!-- ========== CORE PHARMACOGENOMICS FIELDS ========== -->

      <!-- Gene -->
      <fn:string key="gene">
        <xsl:value-of select="(fn:string[@key='Symbol'], fn:string[@key='Gene'])[1]"/>
      </fn:string>

      <!-- Drugs (array) -->
      <xsl:variable name="drug" select="(fn:string[@key='Drug'], fn:string[@key='Name'])[1]"/>
      <xsl:choose>
        <xsl:when test="fn:array[@key='drugs']">
          <fn:array key="drugs">
            <xsl:for-each select="fn:array[@key='drugs']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="not(local:is-null($drug))">
          <fn:array key="drugs">
            <fn:string><xsl:value-of select="$drug"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="drugs"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Phenotype -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">phenotype</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Phenotype']"/>
      </xsl:call-template>

      <!-- Evidence level -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">evidence_level</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Level_of_Evidence']"/>
      </xsl:call-template>

      <!-- ========== PHARMGKB-SPECIFIC FIELDS ========== -->

      <!-- PharmGKB accession -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">pharmgkb_accession</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='PharmGKB_Accession']"/>
      </xsl:call-template>

      <!-- Is VIP gene -->
      <xsl:choose>
        <xsl:when test="fn:boolean[@key='Is_VIP']">
          <fn:boolean key="is_vip_gene">
            <xsl:value-of select="fn:boolean[@key='Is_VIP']"/>
          </fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="is_vip_gene"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Has variant annotation -->
      <xsl:choose>
        <xsl:when test="fn:boolean[@key='Has_Variant_Annotation']">
          <fn:boolean key="has_variant_annotation">
            <xsl:value-of select="fn:boolean[@key='Has_Variant_Annotation']"/>
          </fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="has_variant_annotation"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Annotation ID -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='Annotation_ID']">
          <fn:number key="annotation_id">
            <xsl:value-of select="fn:number[@key='Annotation_ID']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="annotation_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Phenotype category -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">phenotype_category</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Phenotype_Category']"/>
      </xsl:call-template>

      <!-- Trade names (array) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='Trade_Names']">
          <fn:array key="trade_names">
            <xsl:for-each select="fn:array[@key='Trade_Names']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="trade_names"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ATC codes (array) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='ATC_Codes']">
          <fn:array key="atc_codes">
            <xsl:for-each select="fn:array[@key='ATC_Codes']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="atc_codes"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- RxNorm IDs (array) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='RxNorm_Identifiers']">
          <fn:array key="rxnorm_ids">
            <xsl:for-each select="fn:array[@key='RxNorm_Identifiers']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="rxnorm_ids"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- DrugBank ID -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">drugbank_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='DrugBank_ID']"/>
      </xsl:call-template>

      <!-- Has PGx on label -->
      <xsl:choose>
        <xsl:when test="fn:boolean[@key='Has_PGx_On_Label']">
          <fn:boolean key="has_pgx_on_label">
            <xsl:value-of select="fn:boolean[@key='Has_PGx_On_Label']"/>
          </fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="has_pgx_on_label"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Variant haplotypes (array) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='Variant_Haplotypes']">
          <fn:array key="variant_haplotypes">
            <xsl:for-each select="fn:array[@key='Variant_Haplotypes']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="variant_haplotypes"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Population frequencies (object) -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='Biogeographical_Groups']">
          <xsl:copy-of select="fn:map[@key='Biogeographical_Groups']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="population_frequencies"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Pathway ID -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">pathway_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Pathway_ID']"/>
      </xsl:call-template>

      <!-- PK/PD type -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">pk_pd_type</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Type']"/>
      </xsl:call-template>

      <!-- RS ID -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">rsid</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='RS_ID']"/>
      </xsl:call-template>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">PharmGKB</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Named Templates
       ============================================================ -->

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

  <!-- ============================================================
       Reusable Functions
       ============================================================ -->

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or normalize-space($value) = ('', '-', '.', 'NA', 'N/A', 'null')"/>
  </xsl:function>

</xsl:stylesheet>
