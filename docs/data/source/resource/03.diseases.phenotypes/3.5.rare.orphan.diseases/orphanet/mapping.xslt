<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Orphanet -> Rare/Orphan Diseases Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Orphanet Rare Disease Database)
  Target: ../schema.json (3.5 Rare/Orphan Diseases Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ OrphaCode                   │ orpha_code              │ Direct       │
  │ OrphaCode                   │ disease_id              │ Prefix       │
  │ Name                        │ disease_name            │ Direct       │
  │ Gene_Symbol                 │ gene_symbols            │ Array        │
  │ HPO_ID                      │ hpo_ids                 │ Array        │
  │ Inheritance                 │ inheritance_patterns    │ Array        │
  │ DisorderType                │ disorder_type           │ Object       │
  │ AverageAgeOfOnset           │ average_age_of_onset    │ Array        │
  │ Prevalence                  │ prevalence              │ Direct       │
  │ PrevalenceClass             │ prevalence_class        │ Direct       │
  │ PrevalenceGeographic        │ prevalence_geographic   │ Direct       │
  │ GeneAssociationType         │ gene_association_type   │ Direct       │
  │ GeneAssociationStatus       │ gene_association_status │ Enum         │
  │ HGNC_ID                     │ hgnc_id                 │ Direct       │
  │ GeneLocus                   │ gene_locus              │ Direct       │
  │ ExternalReference           │ external_references     │ Array        │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NOTES:
  - Orphanet is the primary rare disease knowledge base
  - ORPHAcode is the unique disease identifier
  - Includes gene associations, HPO phenotypes, epidemiology
  - External references link to OMIM, ICD, UMLS, MeSH
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
        <xsl:value-of select="concat('ORPHA:', fn:number[@key='OrphaCode'])"/>
      </fn:string>

      <xsl:choose>
        <xsl:when test="fn:number[@key='OrphaCode']">
          <fn:number key="orpha_code">
            <xsl:value-of select="fn:number[@key='OrphaCode']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="orpha_code"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DISEASE NAME ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Name'])">
          <fn:null key="disease_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="disease_name">
            <xsl:value-of select="fn:string[@key='Name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GENES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='Gene_Symbol']/fn:string">
          <fn:array key="gene_symbols">
            <xsl:for-each select="fn:array[@key='Gene_Symbol']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='Gene_Symbol'] and not(local:is-null(fn:string[@key='Gene_Symbol']))">
          <fn:array key="gene_symbols">
            <fn:string><xsl:value-of select="fn:string[@key='Gene_Symbol']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_symbols"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== HPO PHENOTYPES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='HPO_ID']/fn:string">
          <fn:array key="hpo_ids">
            <xsl:for-each select="fn:array[@key='HPO_ID']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="hpo_ids"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== INHERITANCE ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='Inheritance']/fn:string">
          <fn:array key="inheritance_patterns">
            <xsl:for-each select="fn:array[@key='Inheritance']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="inheritance_patterns"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DISORDER TYPE ========== -->

      <xsl:choose>
        <xsl:when test="fn:map[@key='DisorderType']">
          <xsl:copy-of select="fn:map[@key='DisorderType']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="disorder_type"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== AGE OF ONSET ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='AverageAgeOfOnset']/fn:string">
          <fn:array key="average_age_of_onset">
            <xsl:for-each select="fn:array[@key='AverageAgeOfOnset']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="average_age_of_onset"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== EPIDEMIOLOGY ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Prevalence'])">
          <fn:null key="prevalence"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="prevalence">
            <xsl:value-of select="fn:string[@key='Prevalence']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='PrevalenceClass'])">
          <fn:null key="prevalence_class"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="prevalence_class">
            <xsl:value-of select="fn:string[@key='PrevalenceClass']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='PrevalenceGeographic'])">
          <fn:null key="prevalence_geographic"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="prevalence_geographic">
            <xsl:value-of select="fn:string[@key='PrevalenceGeographic']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GENE ASSOCIATION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='GeneAssociationType'])">
          <fn:null key="gene_association_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gene_association_type">
            <xsl:value-of select="fn:string[@key='GeneAssociationType']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='GeneAssociationStatus'])">
          <fn:null key="gene_association_status"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gene_association_status">
            <xsl:value-of select="fn:string[@key='GeneAssociationStatus']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GENE IDENTIFIERS ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='HGNC_ID'])">
          <fn:null key="hgnc_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="hgnc_id">
            <xsl:value-of select="fn:string[@key='HGNC_ID']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='GeneLocus'])">
          <fn:null key="gene_locus"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gene_locus">
            <xsl:value-of select="fn:string[@key='GeneLocus']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== EXTERNAL REFERENCES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='ExternalReference']/fn:map">
          <fn:array key="external_references">
            <xsl:for-each select="fn:array[@key='ExternalReference']/fn:map">
              <xsl:copy-of select="."/>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="external_references"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">Orphanet</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:number[@key='OrphaCode']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
