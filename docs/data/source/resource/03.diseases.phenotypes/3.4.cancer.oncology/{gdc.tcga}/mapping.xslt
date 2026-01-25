<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  GDC/TCGA -> Cancer/Oncology Unified Schema Mapping
  ============================================================
  Source: ./schema.md (GDC/TCGA Cancer Genomics)
  Target: ../schema.json (3.4 Cancer/Oncology Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ Hugo_Symbol                 │ gene_symbol             │ Direct       │
  │ case_id                     │ case_id                 │ Direct       │
  │ submitter_id                │ submitter_id            │ Direct       │
  │ project_id                  │ project_id              │ Direct       │
  │ primary_site                │ primary_site            │ Direct       │
  │ sample_type                 │ sample_type             │ Enum         │
  │ Chromosome                  │ chromosome              │ Direct       │
  │ Start_Position              │ start_position          │ Direct       │
  │ End_Position                │ end_position            │ Direct       │
  │ Reference_Allele            │ reference_allele        │ Direct       │
  │ Tumor_Seq_Allele2           │ tumor_allele            │ Direct       │
  │ Variant_Classification      │ variant_classification  │ Enum         │
  │ Variant_Type                │ variant_type            │ Enum         │
  │ HGVSc                       │ hgvsc                   │ Direct       │
  │ HGVSp                       │ hgvsp                   │ Direct       │
  │ t_depth                     │ t_depth                 │ Direct       │
  │ t_alt_count                 │ t_alt_count             │ Direct       │
  │ data_type                   │ data_type               │ Direct       │
  │ experimental_strategy       │ experimental_strategy   │ Enum         │
  │ access                      │ access                  │ Enum         │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NOTES:
  - GDC uses UUIDs for case/file IDs
  - TCGA barcodes follow pattern: TCGA-XX-XXXX
  - Sample types: Primary Tumor, Blood Derived Normal, etc.
  - Variant classifications: Missense_Mutation, Nonsense_Mutation, etc.
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
      <!-- ========== GENE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Hugo_Symbol'])">
          <fn:null key="gene_symbol"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gene_symbol">
            <xsl:value-of select="fn:string[@key='Hugo_Symbol']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CASE/SAMPLE IDENTIFIERS ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='case_id'])">
          <fn:null key="case_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="case_id">
            <xsl:value-of select="fn:string[@key='case_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='submitter_id'])">
          <fn:null key="submitter_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="submitter_id">
            <xsl:value-of select="fn:string[@key='submitter_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='project_id'])">
          <fn:null key="project_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="project_id">
            <xsl:value-of select="fn:string[@key='project_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='primary_site'])">
          <fn:null key="primary_site"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="primary_site">
            <xsl:value-of select="fn:string[@key='primary_site']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='sample_type'])">
          <fn:null key="sample_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="sample_type">
            <xsl:value-of select="fn:string[@key='sample_type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== VARIANT LOCATION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Chromosome'])">
          <fn:null key="chromosome"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="chromosome">
            <xsl:value-of select="fn:string[@key='Chromosome']"/>
          </fn:string>
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

      <!-- ========== ALLELES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Reference_Allele'])">
          <fn:null key="reference_allele"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="reference_allele">
            <xsl:value-of select="fn:string[@key='Reference_Allele']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Tumor_Seq_Allele2'])">
          <fn:null key="tumor_allele"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="tumor_allele">
            <xsl:value-of select="fn:string[@key='Tumor_Seq_Allele2']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== VARIANT CLASSIFICATION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Variant_Classification'])">
          <fn:null key="variant_classification"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="variant_classification">
            <xsl:value-of select="fn:string[@key='Variant_Classification']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Variant_Type'])">
          <fn:null key="variant_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="variant_type">
            <xsl:value-of select="fn:string[@key='Variant_Type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== HGVS NOTATION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='HGVSc'])">
          <fn:null key="hgvsc"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="hgvsc">
            <xsl:value-of select="fn:string[@key='HGVSc']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='HGVSp'])">
          <fn:null key="hgvsp"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="hgvsp">
            <xsl:value-of select="fn:string[@key='HGVSp']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SEQUENCING DEPTH ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='t_depth']">
          <fn:number key="t_depth">
            <xsl:value-of select="fn:number[@key='t_depth']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="t_depth"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='t_alt_count']">
          <fn:number key="t_alt_count">
            <xsl:value-of select="fn:number[@key='t_alt_count']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="t_alt_count"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DATA TYPE/STRATEGY ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='data_type'])">
          <fn:null key="data_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="data_type">
            <xsl:value-of select="fn:string[@key='data_type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='experimental_strategy'])">
          <fn:null key="experimental_strategy"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="experimental_strategy">
            <xsl:value-of select="fn:string[@key='experimental_strategy']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='access'])">
          <fn:null key="access"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="access">
            <xsl:value-of select="fn:string[@key='access']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">GDC/TCGA</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='case_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
