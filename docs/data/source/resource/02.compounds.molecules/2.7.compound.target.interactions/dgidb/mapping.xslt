<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  DGIdb -> Compound-Target Interactions Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Drug Gene Interaction Database)
  Target: ../schema.json (2.7 Compound-Target Interactions Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ drug_name                   │ compound_name           │ Direct       │
  │ gene_name                   │ target_name             │ Direct       │
  │ concept_id                  │ concept_id              │ Direct       │
  │ gene_categories             │ gene_categories         │ Array        │
  │ interaction_types           │ interaction_types       │ Array        │
  │ interaction_score           │ interaction_score       │ Direct       │
  │ sources                     │ dgidb_sources           │ Array        │
  │ approved                    │ approved                │ Boolean      │
  │ chembl_id                   │ chembl_id               │ Direct       │
  │ entrez_id                   │ entrez_id               │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null

  NOTES:
  - DGIdb aggregates drug-gene interactions from 40+ sources
  - Gene categories: KINASE, DRUGGABLE GENOME, etc.
  - Interaction types: inhibitor, agonist, antagonist, etc.
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
        <xsl:when test="local:is-null(fn:string[@key='concept_id'])">
          <fn:null key="concept_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="concept_id">
            <xsl:value-of select="fn:string[@key='concept_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='chembl_id'])">
          <fn:null key="chembl_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="chembl_id">
            <xsl:value-of select="fn:string[@key='chembl_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='entrez_id']">
          <fn:number key="entrez_id">
            <xsl:value-of select="fn:number[@key='entrez_id']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="entrez_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DRUG/COMPOUND ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='drug_name'])">
          <fn:null key="compound_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="compound_name">
            <xsl:value-of select="fn:string[@key='drug_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:boolean[@key='approved']">
          <fn:boolean key="approved">
            <xsl:value-of select="fn:boolean[@key='approved']"/>
          </fn:boolean>
        </xsl:when>
        <xsl:when test="fn:string[@key='approved'] = 'true' or fn:string[@key='approved'] = 'True'">
          <fn:boolean key="approved">true</fn:boolean>
        </xsl:when>
        <xsl:when test="fn:string[@key='approved'] = 'false' or fn:string[@key='approved'] = 'False'">
          <fn:boolean key="approved">false</fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="approved"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== TARGET/GENE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='gene_name'])">
          <fn:null key="target_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="target_name">
            <xsl:value-of select="fn:string[@key='gene_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:array[@key='gene_categories']/fn:string">
          <fn:array key="gene_categories">
            <xsl:for-each select="fn:array[@key='gene_categories']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_categories"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== INTERACTION DATA ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='interaction_types']/fn:string">
          <fn:array key="interaction_types">
            <xsl:for-each select="fn:array[@key='interaction_types']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="interaction_types"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='interaction_score']">
          <fn:number key="interaction_score">
            <xsl:value-of select="fn:number[@key='interaction_score']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="interaction_score"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DATA SOURCES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='sources']/fn:string">
          <fn:array key="dgidb_sources">
            <xsl:for-each select="fn:array[@key='sources']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="dgidb_sources"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">DGIdb</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='concept_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
