<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  BRCA Exchange -> Cancer Genomics Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Cancer Genomics Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬───────────────────────────────┬──────────────┐
  │ Source Field                │ Target Field                  │ Transform    │
  ├─────────────────────────────┼───────────────────────────────┼──────────────┤
  │ Gene_Symbol                 │ gene                          │ Direct       │
  │ HGVS_Protein                │ variant                       │ Direct       │
  │ Genomic_Coordinate_hg38     │ brca_genomic_coordinate       │ Direct       │
  │ HGVS_cDNA                   │ brca_hgvs_cdna                │ Direct       │
  │ HGVS_Protein                │ brca_hgvs_protein             │ Direct       │
  │ Reference_Sequence          │ brca_reference_sequence       │ Direct       │
  │ Pathogenicity_expert        │ brca_pathogenicity_expert     │ Direct       │
  │ Pathogenicity_all           │ brca_pathogenicity_aggregated │ Direct       │
  │ Clinical_significance_ENIGMA│ clinical_significance         │ Direct       │
  │ Clinical_significance_ClinVar│ brca_clinvar_significance    │ Direct       │
  │ Date_last_evaluated_ENIGMA  │ brca_date_evaluated           │ Direct       │
  │ Source                      │ brca_source_databases         │ Direct       │
  │ SCV_ClinVar                 │ brca_clinvar_scv              │ Direct       │
  └─────────────────────────────┴───────────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "-" → null

  NOTES:
  - BRCA Exchange: BRCA1/BRCA2 variant aggregation
  - 52,000+ variants from ClinVar, LOVD, ENIGMA
  - Expert-reviewed classifications from ENIGMA consortium
  - Cancer types: Breast, Ovarian (hereditary)
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
        <xsl:value-of select="fn:string[@key='Gene_Symbol']"/>
      </fn:string>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">variant</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='HGVS_Protein']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">clinical_significance</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Clinical_significance_ENIGMA']"/>
      </xsl:call-template>

      <!-- Cancer types (BRCA = hereditary breast/ovarian) -->
      <fn:array key="cancer_types">
        <fn:string>Hereditary Breast Cancer</fn:string>
        <fn:string>Hereditary Ovarian Cancer</fn:string>
      </fn:array>

      <!-- ========== BRCA EXCHANGE-SPECIFIC FIELDS ========== -->

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">brca_genomic_coordinate</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Genomic_Coordinate_hg38']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">brca_hgvs_cdna</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='HGVS_cDNA']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">brca_hgvs_protein</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='HGVS_Protein']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">brca_reference_sequence</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Reference_Sequence']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">brca_pathogenicity_expert</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Pathogenicity_expert']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">brca_pathogenicity_aggregated</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Pathogenicity_all']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">brca_clinvar_significance</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Clinical_significance_ClinVar']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">brca_date_evaluated</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Date_last_evaluated_ENIGMA']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">brca_source_databases</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Source']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">brca_clinvar_scv</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='SCV_ClinVar']"/>
      </xsl:call-template>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">BRCA Exchange</fn:string>
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
