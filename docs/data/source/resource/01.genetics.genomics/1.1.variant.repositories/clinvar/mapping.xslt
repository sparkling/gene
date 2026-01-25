<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  ClinVar -> Variant Repositories Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Variant Repositories Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ VariationID                 │ variant_id              │ VCV prefix   │
  │ AlleleID                    │ clinvar_allele_id       │ Direct       │
  │ VariationID                 │ clinvar_variation_id    │ Direct       │
  │ Chromosome                  │ chromosome              │ Direct       │
  │ Start                       │ position                │ Direct       │
  │ ReferenceAllele             │ reference_allele        │ Direct       │
  │ AlternateAllele             │ alternate_allele        │ Direct       │
  │ Type                        │ variant_type            │ Direct       │
  │ ClinicalSignificance/CLNSIG │ clinical_significance   │ Split array  │
  │ GeneSymbol                  │ gene_symbols            │ Split array  │
  │ Name                        │ hgvs_expressions        │ Direct       │
  │ Assembly                    │ assembly                │ Direct       │
  │ ReviewStatus/CLNREVSTAT     │ review_status           │ Direct       │
  │ NumberSubmitters            │ number_submitters       │ Direct       │
  │ DateLastEvaluated           │ date_last_evaluated     │ Direct       │
  │ PhenotypeIDS                │ phenotype_ids           │ Split array  │
  │ SomaticClinicalImpact       │ somatic_clinical_impact │ Direct       │
  │ Oncogenicity                │ oncogenicity            │ Direct       │
  │ Origin                      │ origin                  │ Direct       │
  │ CLNDN                       │ disease_names           │ Split array  │
  │ MC                          │ molecular_consequence   │ Direct       │
  │ RS                          │ dbsnp_refsnp_id         │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings, "-", "." → null
  - Missing numeric fields → null

  NOTES:
  - Supports both VCF INFO field format and variant_summary TSV format
  - VCV prefix added to VariationID for variant_id
  - Multiple values in ClinicalSignificance separated by "|" or ";"
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
      <!-- ========== PRIMARY IDENTIFIER ========== -->
      <fn:string key="variant_id">
        <xsl:value-of select="concat('VCV', format-number(number(fn:number[@key='VariationID']), '000000000'))"/>
      </fn:string>

      <!-- ========== CORE VARIANT FIELDS ========== -->
      <fn:string key="chromosome">
        <xsl:value-of select="local:normalize-chromosome(fn:string[@key='Chromosome'])"/>
      </fn:string>

      <fn:number key="position">
        <xsl:value-of select="fn:number[@key='Start']"/>
      </fn:number>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='ReferenceAllele'])">
          <fn:null key="reference_allele"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="reference_allele">
            <xsl:value-of select="fn:string[@key='ReferenceAllele']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='AlternateAllele'])">
          <fn:null key="alternate_allele"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="alternate_allele">
            <xsl:value-of select="fn:string[@key='AlternateAllele']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Variant type -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Type'])">
          <fn:null key="variant_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="variant_type">
            <xsl:value-of select="fn:string[@key='Type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CLINICAL SIGNIFICANCE (ARRAY) ========== -->
      <xsl:variable name="clinsig" select="(fn:string[@key='ClinicalSignificance'], fn:string[@key='CLNSIG'])[1]"/>
      <xsl:choose>
        <xsl:when test="local:is-null($clinsig)">
          <fn:null key="clinical_significance"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="clinical_significance">
            <xsl:for-each select="tokenize($clinsig, '[|;/]')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GENE SYMBOLS (ARRAY) ========== -->
      <xsl:variable name="genes" select="fn:string[@key='GeneSymbol']"/>
      <xsl:choose>
        <xsl:when test="local:is-null($genes)">
          <fn:null key="gene_symbols"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="gene_symbols">
            <xsl:for-each select="tokenize($genes, '[|;,]')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:otherwise>
      </xsl:choose>

      <!-- HGVS expressions -->
      <xsl:variable name="hgvs" select="fn:string[@key='Name']"/>
      <xsl:choose>
        <xsl:when test="local:is-null($hgvs)">
          <fn:null key="hgvs_expressions"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="hgvs_expressions">
            <fn:string><xsl:value-of select="$hgvs"/></fn:string>
          </fn:array>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Assembly -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Assembly'])">
          <fn:null key="assembly"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="assembly">
            <xsl:value-of select="fn:string[@key='Assembly']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CLINVAR-SPECIFIC FIELDS ========== -->
      <fn:number key="clinvar_allele_id">
        <xsl:value-of select="fn:number[@key='AlleleID']"/>
      </fn:number>

      <fn:number key="clinvar_variation_id">
        <xsl:value-of select="fn:number[@key='VariationID']"/>
      </fn:number>

      <!-- Review status -->
      <xsl:variable name="review" select="(fn:string[@key='ReviewStatus'], fn:string[@key='CLNREVSTAT'])[1]"/>
      <xsl:choose>
        <xsl:when test="local:is-null($review)">
          <fn:null key="review_status"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="review_status">
            <xsl:value-of select="$review"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Number of submitters -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='NumberSubmitters']">
          <fn:number key="number_submitters">
            <xsl:value-of select="fn:number[@key='NumberSubmitters']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="number_submitters"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Date last evaluated -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='DateLastEvaluated'])">
          <fn:null key="date_last_evaluated"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="date_last_evaluated">
            <xsl:value-of select="fn:string[@key='DateLastEvaluated']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Phenotype IDs -->
      <xsl:variable name="phenoIds" select="fn:string[@key='PhenotypeIDS']"/>
      <xsl:choose>
        <xsl:when test="local:is-null($phenoIds)">
          <fn:null key="phenotype_ids"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="phenotype_ids">
            <xsl:for-each select="tokenize($phenoIds, '[|;,]')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Somatic clinical impact -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='SomaticClinicalImpact'])">
          <fn:null key="somatic_clinical_impact"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="somatic_clinical_impact">
            <xsl:value-of select="fn:string[@key='SomaticClinicalImpact']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Oncogenicity -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Oncogenicity'])">
          <fn:null key="oncogenicity"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="oncogenicity">
            <xsl:value-of select="fn:string[@key='Oncogenicity']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Origin -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Origin'])">
          <fn:null key="origin"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="origin">
            <xsl:value-of select="fn:string[@key='Origin']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Disease names -->
      <xsl:variable name="diseases" select="fn:string[@key='CLNDN']"/>
      <xsl:choose>
        <xsl:when test="local:is-null($diseases)">
          <fn:null key="disease_names"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="disease_names">
            <xsl:for-each select="tokenize($diseases, '[|;]')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Molecular consequence -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='MC'])">
          <fn:null key="molecular_consequence"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="molecular_consequence">
            <xsl:value-of select="fn:string[@key='MC']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- dbSNP RS ID (from RS field if present) -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='RS'] and not(local:is-null(fn:string[@key='RS']))">
          <fn:string key="dbsnp_refsnp_id">
            <xsl:value-of select="fn:string[@key='RS']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="dbsnp_refsnp_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SOURCE METADATA ========== -->
      <fn:map key="_source">
        <fn:string key="database">ClinVar</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="concat('VCV', format-number(number(fn:number[@key='VariationID']), '000000000'))"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Reusable Functions
       ============================================================ -->

  <!-- Normalize chromosome: remove "chr" prefix -->
  <xsl:function name="local:normalize-chromosome" as="xs:string">
    <xsl:param name="chr" as="xs:string"/>
    <xsl:value-of select="replace($chr, '^chr', '')"/>
  </xsl:function>

  <!-- Check if value is null/empty -->
  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = ('', '-', '.', 'NA', 'N/A', 'null', 'na', 'not provided', 'not specified')"/>
  </xsl:function>

</xsl:stylesheet>
