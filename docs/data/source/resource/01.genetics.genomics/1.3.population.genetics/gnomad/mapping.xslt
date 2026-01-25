<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  gnomAD -> Population Genetics Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Population Genetics Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────────┬──────────────┐
  │ Source Field                │ Target Field                │ Transform    │
  ├─────────────────────────────┼─────────────────────────────┼──────────────┤
  │ CHROM                       │ chromosome                  │ Normalize    │
  │ POS                         │ position                    │ Direct       │
  │ REF                         │ reference_allele            │ Direct       │
  │ ALT                         │ alternate_allele            │ Direct       │
  │ AF                          │ allele_frequency            │ Direct       │
  │ AC                          │ allele_count                │ Direct       │
  │ AN                          │ allele_number               │ Direct       │
  │ nhomalt                     │ gnomad_homozygote_count     │ Direct       │
  │ AC_XX                       │ gnomad_ac_xx                │ Direct       │
  │ AC_XY                       │ gnomad_ac_xy                │ Direct       │
  │ AC_afr (derived)            │ af_african                  │ Compute AF   │
  │ AC_amr (derived)            │ af_american                 │ Compute AF   │
  │ AC_eas (derived)            │ af_east_asian               │ Compute AF   │
  │ AC_nfe (derived)            │ af_european                 │ Compute AF   │
  │ AC_sas (derived)            │ af_south_asian              │ Compute AF   │
  │ AF_asj (derived)            │ gnomad_af_ashkenazi         │ Direct       │
  │ AF_fin (derived)            │ gnomad_af_finnish           │ Direct       │
  │ faf95                       │ gnomad_faf95                │ Direct       │
  │ faf99                       │ gnomad_faf99                │ Direct       │
  │ fafmax                      │ gnomad_fafmax               │ Direct       │
  │ FS                          │ qc_fisher_strand            │ Direct       │
  │ MQ                          │ qc_mapping_quality          │ Direct       │
  │ QD                          │ qc_quality_by_depth         │ Direct       │
  │ VQSLOD                      │ qc_vqslod                   │ Direct       │
  │ pLI                         │ constraint_pli              │ Direct       │
  │ LOEUF                       │ constraint_loeuf            │ Direct       │
  │ oe_lof                      │ constraint_oe_lof           │ Direct       │
  │ exp_lof                     │ constraint_expected_lof     │ Direct       │
  │ obs_lof                     │ constraint_observed_lof     │ Direct       │
  │ lof                         │ loftee_annotation           │ Direct       │
  │ lof_filter                  │ loftee_filter               │ Direct       │
  │ lof_flags                   │ loftee_flags                │ Direct       │
  │ cadd_phred                  │ annotation_cadd_phred       │ Direct       │
  │ revel                       │ annotation_revel            │ Direct       │
  │ sift                        │ annotation_sift             │ Direct       │
  │ polyphen                    │ annotation_polyphen         │ Direct       │
  │ spliceai_ds_max             │ annotation_spliceai_max     │ Direct       │
  │ FILTER                      │ filter_status               │ Direct       │
  └─────────────────────────────┴─────────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "." → null
  - Missing frequency/count fields → null

  NOTES:
  - gnomAD v4: 807K individuals, 8 ancestry groups
  - Includes gene constraint metrics (pLI, LOEUF)
  - LOFTEE annotations for loss-of-function variants
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
      <!-- ========== CORE VARIANT FIELDS ========== -->
      <fn:string key="chromosome">
        <xsl:value-of select="local:normalize-chromosome(fn:string[@key='CHROM'])"/>
      </fn:string>

      <fn:number key="position">
        <xsl:value-of select="fn:number[@key='POS']"/>
      </fn:number>

      <fn:string key="reference_allele">
        <xsl:value-of select="fn:string[@key='REF']"/>
      </fn:string>

      <fn:string key="alternate_allele">
        <xsl:value-of select="fn:string[@key='ALT']"/>
      </fn:string>

      <!-- ========== GLOBAL FREQUENCY METRICS ========== -->
      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">allele_frequency</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='AF']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">allele_count</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='AC']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">allele_number</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='AN']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">filter_status</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='FILTER']"/>
      </xsl:call-template>

      <!-- ========== POPULATION-SPECIFIC FREQUENCIES ========== -->
      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">af_african</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='AF_afr']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">af_american</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='AF_amr']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">af_east_asian</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='AF_eas']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">af_european</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='AF_nfe']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">af_south_asian</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='AF_sas']"/>
      </xsl:call-template>

      <!-- ========== GNOMAD-SPECIFIC FIELDS ========== -->
      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">gnomad_homozygote_count</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='nhomalt']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">gnomad_ac_xx</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='AC_XX']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">gnomad_ac_xy</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='AC_XY']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">gnomad_af_ashkenazi</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='AF_asj']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">gnomad_af_finnish</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='AF_fin']"/>
      </xsl:call-template>

      <!-- Filtering allele frequencies -->
      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">gnomad_faf95</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='faf95']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">gnomad_faf99</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='faf99']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">gnomad_fafmax</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='fafmax']"/>
      </xsl:call-template>

      <!-- ========== QC METRICS ========== -->
      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">qc_fisher_strand</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='FS']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">qc_mapping_quality</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='MQ']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">qc_quality_by_depth</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='QD']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">qc_vqslod</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='VQSLOD']"/>
      </xsl:call-template>

      <!-- ========== CONSTRAINT METRICS ========== -->
      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">constraint_expected_lof</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='exp_lof']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">constraint_observed_lof</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='obs_lof']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">constraint_oe_lof</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='oe_lof']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">constraint_pli</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='pLI']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">constraint_loeuf</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='LOEUF']"/>
      </xsl:call-template>

      <!-- ========== LOFTEE ANNOTATIONS ========== -->
      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">loftee_annotation</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='lof']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">loftee_filter</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='lof_filter']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">loftee_flags</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='lof_flags']"/>
      </xsl:call-template>

      <!-- ========== INTEGRATED ANNOTATIONS ========== -->
      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">annotation_cadd_phred</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='cadd_phred']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">annotation_revel</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='revel']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">annotation_sift</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='sift']"/>
      </xsl:call-template>

      <xsl:call-template name="string-or-null">
        <xsl:with-param name="key">annotation_polyphen</xsl:with-param>
        <xsl:with-param name="value" select="fn:string[@key='polyphen']"/>
      </xsl:call-template>

      <xsl:call-template name="number-or-null">
        <xsl:with-param name="key">annotation_spliceai_max</xsl:with-param>
        <xsl:with-param name="node" select="fn:number[@key='spliceai_ds_max']"/>
      </xsl:call-template>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">gnomAD</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Named Templates
       ============================================================ -->

  <xsl:template name="number-or-null">
    <xsl:param name="key"/>
    <xsl:param name="node"/>
    <xsl:choose>
      <xsl:when test="$node">
        <fn:number key="{$key}">
          <xsl:value-of select="$node"/>
        </fn:number>
      </xsl:when>
      <xsl:otherwise>
        <fn:null key="{$key}"/>
      </xsl:otherwise>
    </xsl:choose>
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

  <!-- ============================================================
       Reusable Functions
       ============================================================ -->

  <xsl:function name="local:normalize-chromosome" as="xs:string">
    <xsl:param name="chr" as="xs:string"/>
    <xsl:value-of select="replace($chr, '^chr', '')"/>
  </xsl:function>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or normalize-space($value) = ('', '-', '.', 'NA', 'N/A', 'null')"/>
  </xsl:function>

</xsl:stylesheet>
