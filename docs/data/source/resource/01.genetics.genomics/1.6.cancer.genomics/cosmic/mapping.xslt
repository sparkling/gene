<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  COSMIC -> Cancer Genomics Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Cancer Genomics Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────────┬──────────────┐
  │ Source Field                │ Target Field                │ Transform    │
  ├─────────────────────────────┼─────────────────────────────┼──────────────┤
  │ Gene                        │ gene                        │ Direct       │
  │ AA_Mutation                 │ variant                     │ Direct       │
  │ COSV_ID                     │ cosmic_cosv_id              │ Direct       │
  │ COSM_ID                     │ cosmic_cosm_id              │ Direct       │
  │ Transcript                  │ cosmic_transcript           │ Direct       │
  │ CDS_Mutation                │ cosmic_cds_mutation         │ Direct       │
  │ AA_Mutation                 │ cosmic_aa_mutation          │ Direct       │
  │ Mutation_Type               │ cosmic_mutation_type        │ Direct       │
  │ Mutation_Description        │ cosmic_mutation_description │ Direct       │
  │ GRCh38_Position             │ chromosome                  │ Parse        │
  │ Sample_ID                   │ cosmic_sample_id            │ Direct       │
  │ Sample_Name                 │ cosmic_sample_name          │ Direct       │
  │ Primary_Site                │ primary_site                │ Direct       │
  │ Primary_Site                │ cancer_types                │ Array        │
  │ Primary_Histology           │ primary_histology           │ Direct       │
  │ Sample_Type                 │ sample_type                 │ Direct       │
  │ Age                         │ patient_age                 │ Direct       │
  │ census_tier                 │ cgc_tier                    │ Direct       │
  │ Role_in_Cancer              │ role_in_cancer              │ Direct       │
  │ Signature                   │ mutational_signature        │ Direct       │
  │ Type (SBS/DBS/ID)           │ signature_type              │ Direct       │
  │ Contribution                │ signature_contribution      │ Direct       │
  │ Proposed_Etiology           │ proposed_etiology           │ Direct       │
  │ Fusion_ID                   │ fusion_id                   │ Direct       │
  │ Gene_5prime                 │ fusion_gene_5prime          │ Direct       │
  │ Gene_3prime                 │ fusion_gene_3prime          │ Direct       │
  │ referenceGenome             │ reference_genome            │ Direct       │
  └─────────────────────────────┴─────────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "NS" → null

  NOTES:
  - COSMIC: Catalogue of Somatic Mutations in Cancer
  - COSV = genomic mutation ID, COSM = legacy ID
  - 17M+ coding mutations, 1.5M+ samples
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
        <xsl:value-of select="fn:string[@key='Gene']"/>
      </fn:string>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">variant</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='AA_Mutation']"/>
      </xsl:call-template>

      <!-- Cancer types from primary site -->
      <xsl:variable name="site" select="fn:string[@key='Primary_Site']"/>
      <xsl:choose>
        <xsl:when test="local:is-null($site)">
          <fn:array key="cancer_types"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="cancer_types">
            <fn:string><xsl:value-of select="$site"/></fn:string>
          </fn:array>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Chromosome from GRCh38_Position (format: chr:start-end) -->
      <xsl:variable name="pos" select="fn:string[@key='GRCh38_Position']"/>
      <xsl:choose>
        <xsl:when test="contains($pos, ':')">
          <fn:string key="chromosome">
            <xsl:value-of select="substring-before($pos, ':')"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="chromosome"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Reference genome -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='GRCh38_Position']">
          <fn:string key="reference_genome">GRCh38</fn:string>
        </xsl:when>
        <xsl:when test="fn:string[@key='GRCh37_Position']">
          <fn:string key="reference_genome">GRCh37</fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="reference_genome"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== COSMIC-SPECIFIC FIELDS ========== -->

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cosmic_cosv_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='COSV_ID']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cosmic_cosm_id</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='COSM_ID']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cosmic_transcript</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Transcript']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cosmic_cds_mutation</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='CDS_Mutation']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cosmic_aa_mutation</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='AA_Mutation']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cosmic_mutation_type</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Mutation_Type']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cosmic_mutation_description</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Mutation_Description']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='Sample_ID']">
          <fn:number key="cosmic_sample_id">
            <xsl:value-of select="fn:number[@key='Sample_ID']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cosmic_sample_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">cosmic_sample_name</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Sample_Name']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">primary_site</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Primary_Site']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">primary_histology</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Primary_Histology']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">sample_type</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Sample_Type']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='Age']">
          <fn:number key="patient_age">
            <xsl:value-of select="fn:number[@key='Age']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="patient_age"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='census_tier']">
          <fn:number key="cgc_tier">
            <xsl:value-of select="fn:number[@key='census_tier']"/>
          </fn:number>
        </xsl:when>
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
        <xsl:with-param name="key">role_in_cancer</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Role_in_Cancer']"/>
      </xsl:call-template>

      <!-- Signature fields -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">mutational_signature</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Signature']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">signature_type</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Type']"/>
      </xsl:call-template>

      <xsl:choose>
        <xsl:when test="fn:number[@key='Contribution']">
          <fn:number key="signature_contribution">
            <xsl:value-of select="fn:number[@key='Contribution']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="signature_contribution"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">proposed_etiology</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Proposed_Etiology']"/>
      </xsl:call-template>

      <!-- Fusion fields -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='Fusion_ID']">
          <fn:number key="fusion_id">
            <xsl:value-of select="fn:number[@key='Fusion_ID']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="fusion_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">fusion_gene_5prime</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Gene_5prime']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">fusion_gene_3prime</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='Gene_3prime']"/>
      </xsl:call-template>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">COSMIC</fn:string>
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
    <xsl:sequence select="not($value) or normalize-space($value) = ('', '-', '.', 'NA', 'N/A', 'null', 'NS')"/>
  </xsl:function>

</xsl:stylesheet>
