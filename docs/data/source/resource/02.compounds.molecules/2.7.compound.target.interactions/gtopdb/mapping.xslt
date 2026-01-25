<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  GtoPdb -> Compound-Target Interactions Unified Schema Mapping
  ============================================================
  Source: ./schema.md (IUPHAR/BPS Guide to PHARMACOLOGY)
  Target: ../schema.json (2.7 Compound-Target Interactions Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ ligand_id                   │ gtopdb_ligand_id        │ Direct       │
  │ ligand_name                 │ compound_name           │ Direct       │
  │ smiles                      │ smiles                  │ Direct       │
  │ inchi_key                   │ inchi_key               │ Direct       │
  │ target_id                   │ gtopdb_target_id        │ Direct       │
  │ target_name                 │ target_name             │ Direct       │
  │ uniprot_id                  │ uniprot_id              │ Direct       │
  │ systematic_name             │ systematic_name         │ Direct       │
  │ target_family               │ target_family           │ Direct       │
  │ action                      │ interaction_types       │ Array        │
  │ selectivity                 │ selectivity             │ Direct       │
  │ endogenous                  │ is_endogenous           │ Boolean      │
  │ inn                         │ inn                     │ Direct       │
  │ hgnc_symbol                 │ hgnc_symbol             │ Direct       │
  │ affinity_type               │ affinity_type           │ Direct       │
  │ affinity_value              │ affinity_value          │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null

  NOTES:
  - GtoPdb is curated by IUPHAR for pharmacology
  - Target families: GPCR, ion_channel, kinase, etc.
  - Includes endogenous ligand information
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
        <xsl:when test="fn:number[@key='ligand_id']">
          <fn:number key="gtopdb_ligand_id">
            <xsl:value-of select="fn:number[@key='ligand_id']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gtopdb_ligand_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='target_id']">
          <fn:number key="gtopdb_target_id">
            <xsl:value-of select="fn:number[@key='target_id']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gtopdb_target_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== LIGAND PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='ligand_name'])">
          <fn:null key="compound_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="compound_name">
            <xsl:value-of select="fn:string[@key='ligand_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='smiles'])">
          <fn:null key="smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="smiles">
            <xsl:value-of select="fn:string[@key='smiles']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='inchi_key'])">
          <fn:null key="inchi_key"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="inchi_key">
            <xsl:value-of select="fn:string[@key='inchi_key']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='systematic_name'])">
          <fn:null key="systematic_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="systematic_name">
            <xsl:value-of select="fn:string[@key='systematic_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='inn'])">
          <fn:null key="inn"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="inn">
            <xsl:value-of select="fn:string[@key='inn']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:boolean[@key='endogenous']">
          <fn:boolean key="is_endogenous">
            <xsl:value-of select="fn:boolean[@key='endogenous']"/>
          </fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="is_endogenous"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== TARGET PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='target_name'])">
          <fn:null key="target_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="target_name">
            <xsl:value-of select="fn:string[@key='target_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='uniprot_id'])">
          <fn:null key="uniprot_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="uniprot_id">
            <xsl:value-of select="fn:string[@key='uniprot_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='target_family'])">
          <fn:null key="target_family"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="target_family">
            <xsl:value-of select="fn:string[@key='target_family']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='hgnc_symbol'])">
          <fn:null key="hgnc_symbol"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="hgnc_symbol">
            <xsl:value-of select="fn:string[@key='hgnc_symbol']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== INTERACTION DATA ========== -->

      <xsl:choose>
        <xsl:when test="fn:string[@key='action'] and not(local:is-null(fn:string[@key='action']))">
          <fn:array key="interaction_types">
            <fn:string><xsl:value-of select="fn:string[@key='action']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:array[@key='action']/fn:string">
          <fn:array key="interaction_types">
            <xsl:for-each select="fn:array[@key='action']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="interaction_types"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='selectivity'])">
          <fn:null key="selectivity"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="selectivity">
            <xsl:value-of select="fn:string[@key='selectivity']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== AFFINITY DATA ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='affinity_type'])">
          <fn:null key="affinity_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="affinity_type">
            <xsl:value-of select="fn:string[@key='affinity_type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='affinity_value']">
          <fn:number key="affinity_value">
            <xsl:value-of select="fn:number[@key='affinity_value']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="affinity_value"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">GtoPdb</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:number[@key='ligand_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
