<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  cBioPortal -> Cancer Genomics Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Cancer Genomics Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────────┬──────────────┐
  │ Source Field                │ Target Field                │ Transform    │
  ├─────────────────────────────┼─────────────────────────────┼──────────────┤
  │ Hugo_Symbol                 │ gene                        │ Direct       │
  │ Entrez_Gene_Id              │ entrez_gene_id              │ Direct       │
  │ studyId                     │ cbio_study_id               │ Direct       │
  │ cancerTypeId                │ cbio_cancer_type_id         │ Direct       │
  │ cancerTypeId                │ cancer_types                │ Array        │
  │ molecularProfileId          │ cbio_molecular_profile_id   │ Direct       │
  │ molecularAlterationType     │ cbio_alteration_type        │ Direct       │
  │ datatype                    │ cbio_datatype               │ Direct       │
  │ allSampleCount              │ cbio_sample_count           │ Direct       │
  │ referenceGenome             │ reference_genome            │ Direct       │
  │ Chromosome                  │ chromosome                  │ Direct       │
  │ Start_Position              │ start_position              │ Direct       │
  │ End_Position                │ end_position                │ Direct       │
  │ Variant_Classification      │ variant_classification      │ Direct       │
  │ Variant_Type                │ variant_type                │ Direct       │
  │ Reference_Allele            │ reference_allele            │ Direct       │
  │ Tumor_Seq_Allele2           │ tumor_allele                │ Direct       │
  │ Tumor_Sample_Barcode        │ sample_barcode              │ Direct       │
  │ HGVSp_Short                 │ variant                     │ Direct       │
  └─────────────────────────────┴─────────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings → null

  NOTES:
  - cBioPortal: Cancer Genomics Data Visualization
  - 423+ studies, 1M+ samples
  - MAF format for mutations
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
        <xsl:value-of select="(fn:string[@key='Hugo_Symbol'], fn:string[@key='hugoGeneSymbol'])[1]"/>
      </fn:string>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">variant</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='HGVSp_Short']"/>
      </xsl:call-template>

      <!-- Cancer types from cancerTypeId -->
      <xsl:variable name="cancer" select="fn:string[@key='cancerTypeId']"/>
      <xsl:choose>
        <xsl:when test="local:is-null($cancer)">
          <fn:array key="cancer_types"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="cancer_types">
            <fn:string><xsl:value-of select="$cancer"/></fn:string>
          </fn:array>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">chromosome</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Chromosome']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">reference_genome</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='referenceGenome']"/>
      </xsl:call-template>

      <!-- ========== CBIOPORTAL-SPECIFIC FIELDS ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='Entrez_Gene_Id']">
          <fn:number key="entrez_gene_id">
            <xsl:value-of select="fn:number[@key='Entrez_Gene_Id']"/>
          </fn:number>
        </xsl:when>
        <xsl:when test="fn:number[@key='entrezGeneId']">
          <fn:number key="entrez_gene_id">
            <xsl:value-of select="fn:number[@key='entrezGeneId']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="entrez_gene_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cbio_study_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='studyId']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cbio_cancer_type_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='cancerTypeId']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cbio_molecular_profile_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='molecularProfileId']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cbio_alteration_type</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='molecularAlterationType']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cbio_datatype</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='datatype']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='allSampleCount']">
          <fn:number key="cbio_sample_count">
            <xsl:value-of select="fn:number[@key='allSampleCount']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cbio_sample_count"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='Start_Position']">
          <fn:number key="start_position">
            <xsl:value-of select="fn:number[@key='Start_Position']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="start_position"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='End_Position']">
          <fn:number key="end_position">
            <xsl:value-of select="fn:number[@key='End_Position']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="end_position"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">variant_classification</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Variant_Classification']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">variant_type</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Variant_Type']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">reference_allele</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Reference_Allele']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">tumor_allele</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Tumor_Seq_Allele2']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">sample_barcode</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Tumor_Sample_Barcode']"/>
      </xsl:call-template>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">cBioPortal</fn:string>
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
