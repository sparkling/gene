<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  BindingDB -> Compound-Target Interactions Unified Schema Mapping
  ============================================================
  Source: ./schema.md (BindingDB Binding Affinity Database)
  Target: ../schema.json (2.7 Compound-Target Interactions Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ monomer_id                  │ monomer_id              │ Direct       │
  │ ligand_name                 │ compound_name           │ Direct       │
  │ smiles                      │ smiles                  │ Direct       │
  │ inchi_key                   │ inchi_key               │ Direct       │
  │ target_name                 │ target_name             │ Direct       │
  │ uniprot_id                  │ uniprot_id              │ Direct       │
  │ activity_type               │ affinity_type           │ Direct       │
  │ activity_value              │ affinity_value          │ Direct       │
  │ activity_unit               │ affinity_unit           │ Direct       │
  │ activity_relation           │ affinity_relation       │ Direct       │
  │ ph                          │ assay_ph                │ Direct       │
  │ temperature                 │ assay_temperature       │ Direct       │
  │ pdb_ids                     │ pdb_ids                 │ Array        │
  │ pmid                        │ pmid                    │ Direct       │
  │ bindingdb_id                │ bindingdb_id            │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null
  - Numeric values: undefined → null

  NOTES:
  - BindingDB contains measured binding affinities
  - Activity types: IC50, Ki, Kd, EC50
  - Includes assay conditions (pH, temperature)
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
        <xsl:when test="fn:number[@key='monomer_id']">
          <fn:number key="monomer_id">
            <xsl:value-of select="fn:number[@key='monomer_id']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="monomer_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='bindingdb_id']">
          <fn:number key="bindingdb_id">
            <xsl:value-of select="fn:number[@key='bindingdb_id']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="bindingdb_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== COMPOUND PROPERTIES ========== -->

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

      <!-- ========== BINDING AFFINITY ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='activity_type'])">
          <fn:null key="affinity_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="affinity_type">
            <xsl:value-of select="fn:string[@key='activity_type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='activity_value']">
          <fn:number key="affinity_value">
            <xsl:value-of select="fn:number[@key='activity_value']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="affinity_value"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='activity_unit'])">
          <fn:null key="affinity_unit"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="affinity_unit">
            <xsl:value-of select="fn:string[@key='activity_unit']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='activity_relation'])">
          <fn:null key="affinity_relation"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="affinity_relation">
            <xsl:value-of select="fn:string[@key='activity_relation']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ASSAY CONDITIONS ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='ph']">
          <fn:number key="assay_ph">
            <xsl:value-of select="fn:number[@key='ph']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="assay_ph"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='temperature']">
          <fn:number key="assay_temperature">
            <xsl:value-of select="fn:number[@key='temperature']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="assay_temperature"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== STRUCTURE / REFERENCES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='pdb_ids']/fn:string">
          <fn:array key="pdb_ids">
            <xsl:for-each select="fn:array[@key='pdb_ids']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="pdb_ids"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='pmid']">
          <fn:number key="pmid">
            <xsl:value-of select="fn:number[@key='pmid']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="pmid"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">BindingDB</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:number[@key='bindingdb_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
