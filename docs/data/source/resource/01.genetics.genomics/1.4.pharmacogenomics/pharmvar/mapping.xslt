<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  PharmVar -> Pharmacogenomics Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Pharmacogenomics Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────────┬──────────────┐
  │ Source Field                │ Target Field                │ Transform    │
  ├─────────────────────────────┼─────────────────────────────┼──────────────┤
  │ gene                        │ gene                        │ Direct       │
  │ allele_name                 │ pharmvar_allele_name        │ Direct       │
  │ allele_name                 │ allele                      │ Direct       │
  │ core_allele                 │ pharmvar_core_allele        │ Direct       │
  │ activity_value              │ activity_score              │ Direct       │
  │ function                    │ allele_function             │ Direct       │
  │ defining_variants           │ defining_variants           │ Array        │
  │ ref_seq_id                  │ reference_sequence_id       │ Direct       │
  │ hgvs_g                      │ hgvs_genomic                │ Direct       │
  │ hgvs_c                      │ hgvs_coding                 │ Direct       │
  │ hgvs_p                      │ hgvs_protein                │ Direct       │
  │ rsid                        │ rsid                        │ Direct       │
  │ impact                      │ variant_impact              │ Direct       │
  │ phenotype                   │ phenotype                   │ Direct       │
  └─────────────────────────────┴─────────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "-" → null

  NOTES:
  - PharmVar: Pharmacogene Variation Consortium
  - Official repository for pharmacogene star allele nomenclature
  - Core alleles vs sub-alleles (e.g., *4 vs *4.001, *4.002)
  - Activity values: 0 (no function), 0.5 (decreased), 1 (normal), etc.
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
        <xsl:value-of select="fn:string[@key='gene']"/>
      </fn:string>

      <!-- Allele (star allele) -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">allele</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='allele_name']"/>
      </xsl:call-template>

      <!-- Allele function -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">allele_function</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='function']"/>
      </xsl:call-template>

      <!-- Activity score -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='activity_value']">
          <fn:number key="activity_score">
            <xsl:value-of select="fn:number[@key='activity_value']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="activity_score"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Phenotype -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">phenotype</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='phenotype']"/>
      </xsl:call-template>

      <!-- ========== PHARMVAR-SPECIFIC FIELDS ========== -->

      <!-- PharmVar allele name -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">pharmvar_allele_name</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='allele_name']"/>
      </xsl:call-template>

      <!-- Core allele -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">pharmvar_core_allele</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='core_allele']"/>
      </xsl:call-template>

      <!-- Defining variants (array) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='defining_variants']">
          <fn:array key="defining_variants">
            <xsl:for-each select="fn:array[@key='defining_variants']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='defining_variants'] and not(local:is-null(fn:string[@key='defining_variants']))">
          <fn:array key="defining_variants">
            <xsl:for-each select="tokenize(fn:string[@key='defining_variants'], '[;,]')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="defining_variants"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Reference sequence ID -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">reference_sequence_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='ref_seq_id']"/>
      </xsl:call-template>

      <!-- HGVS genomic -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">hgvs_genomic</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='hgvs_g']"/>
      </xsl:call-template>

      <!-- HGVS coding -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">hgvs_coding</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='hgvs_c']"/>
      </xsl:call-template>

      <!-- HGVS protein -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">hgvs_protein</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='hgvs_p']"/>
      </xsl:call-template>

      <!-- RS ID -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">rsid</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='rsid']"/>
      </xsl:call-template>

      <!-- Variant impact -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">variant_impact</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='impact']"/>
      </xsl:call-template>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">PharmVar</fn:string>
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
