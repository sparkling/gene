<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  PharmGKB -> Disease-Gene Associations Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Pharmacogenomics Knowledge Base)
  Target: ../schema.json (3.3 Disease-Gene Associations Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ pgkb_id                     │ pgkb_id                 │ Direct       │
  │ Symbol                      │ gene_symbol             │ Direct       │
  │ Gene_ID                     │ gene_id                 │ Direct       │
  │ Phenotype_Category          │ disease_name            │ Direct       │
  │ Level_of_Evidence           │ evidence_level          │ Direct       │
  │ Haplotype                   │ haplotype               │ Direct       │
  │ Diplotype                   │ diplotype               │ Direct       │
  │ Phenotype                   │ phenotype               │ Direct       │
  │ Activity_Value              │ activity_value          │ Direct       │
  │ Recommendation              │ recommendation          │ Direct       │
  │ Guideline_Source            │ guideline_source        │ Enum         │
  │ Biogeographical_Frequency   │ population_frequency    │ Object       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NOTES:
  - PharmGKB uses PA prefix for identifiers
  - Evidence levels: 1A, 1B, 2A, 2B, 3, 4
  - Includes star allele haplotypes and clinical recommendations
  - Guideline sources: CPIC, DPWG, FDA
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

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='pgkb_id'])">
          <fn:null key="pgkb_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="pgkb_id">
            <xsl:value-of select="fn:string[@key='pgkb_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GENE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Symbol'])">
          <fn:null key="gene_symbol"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gene_symbol">
            <xsl:value-of select="fn:string[@key='Symbol']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Gene_ID'])">
          <fn:null key="gene_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gene_id">
            <xsl:value-of select="fn:string[@key='Gene_ID']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PHENOTYPE/DISEASE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Phenotype_Category'])">
          <fn:null key="disease_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="disease_name">
            <xsl:value-of select="fn:string[@key='Phenotype_Category']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Phenotype'])">
          <fn:null key="phenotype"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="phenotype">
            <xsl:value-of select="fn:string[@key='Phenotype']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PHARMACOGENOMICS ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Haplotype'])">
          <fn:null key="haplotype"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="haplotype">
            <xsl:value-of select="fn:string[@key='Haplotype']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Diplotype'])">
          <fn:null key="diplotype"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="diplotype">
            <xsl:value-of select="fn:string[@key='Diplotype']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Activity_Value'])">
          <fn:null key="activity_value"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="activity_value">
            <xsl:value-of select="fn:string[@key='Activity_Value']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== EVIDENCE & GUIDELINES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Level_of_Evidence'])">
          <fn:null key="evidence_level"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="evidence_level">
            <xsl:value-of select="fn:string[@key='Level_of_Evidence']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Recommendation'])">
          <fn:null key="recommendation"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="recommendation">
            <xsl:value-of select="fn:string[@key='Recommendation']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Guideline_Source'])">
          <fn:null key="guideline_source"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="guideline_source">
            <xsl:value-of select="fn:string[@key='Guideline_Source']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== POPULATION FREQUENCY ========== -->

      <xsl:choose>
        <xsl:when test="fn:map[@key='Biogeographical_Frequency']">
          <xsl:copy-of select="fn:map[@key='Biogeographical_Frequency']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="population_frequency"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">PharmGKB</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='pgkb_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
